function [p, ld_exchange] = exchange_sensors8(p, x, sig_n, sig_f, ls, f)
arguments (Input)
    p
    x
    sig_n
    sig_f (1,1) double = 1.0 % Default value if not provided
    ls (1,1) double = 1.0 % Default value if not provided
    f (1,1) double = 1.0 % Default value if not provided
end
arguments (Output)
    p
    ld_exchange
end

k = length(p);
n = length(x);

% functions for getting elements of K
K_fun = @(x) gaussKern(x,sig_f,ls);
K_fun_offdiag = @(x,x2) gaussKern(x,sig_f,ls,x2);
K_diag = @(x) gaussKern(x,sig_f,ls,'1');

alpha = @(i) 1/(sig_n^2)*K_diag(1)+1; % "on" diagonal elements
A = @(i,j) 1/(sig_n^2)*K_fun_offdiag(x(i,:),x(j,:));  % off diagonal elements
A_p = (1/(sig_n^2))*K_fun(x(p,:)) + eye(k);
L = chol(A_p, 'lower'); 

%% initialize rhs_all
p_temp = p(2:end);
test_ind = setdiff(1:n, p_temp);
rhs_all = A(p_temp, test_ind); 


converged = false;
iteration_num = 1;
no_swap_counter = 0; % Tracks consecutive iterations without a swap

while ~converged
    % Save the sensor currently being evaluated before p shifts
    p_orig_front = p(1); 
    
    % Fast Cholesky Downdate (removes the first sensor from L)
    col_to_add = L(2:end,1);
    L_fast = L(2:end,2:end);
    L_fast = (cholupdate(L_fast', col_to_add, '+'))';
    
    % Compute beta of the unswapped baseline sensor
    orig_ind = find(test_ind == p_orig_front);
    b_orig = L_fast \ rhs_all(:, orig_ind); 
    beta_orig = sqrt(alpha(p_orig_front) - b_orig'*b_orig);
    
    %% --- parallel ---
    num_candidates = length(test_ind);
    B = 1024; 
    
    num_blocks = ceil(num_candidates / B);
    rhs_chunks = cell(num_blocks, 1);
    idx_chunks = cell(num_blocks, 1);
    
    for b = 1:num_blocks
        start_idx = (b-1)*B + 1;
        end_idx = min(b*B, num_candidates);
        rhs_chunks{b} = rhs_all(:, start_idx:end_idx);
        idx_chunks{b} = start_idx:end_idx; 
    end
    
    pool = gcp(); 
    futures = parallel.FevalFuture.empty(num_blocks, 0);
    
    % Dispatch all blocks to the workers
    for b = 1:num_blocks
        futures(b) = parfeval(pool, @evaluate_block, 2, ...
            L_fast, rhs_chunks{b}, alpha(test_ind(idx_chunks{b})), beta_orig, f);
    end
    
    % Arrays to collect the best candidate from EACH block
    block_max_betas = -inf(num_blocks, 1);
    block_winner_local_idx = zeros(num_blocks, 1);
    
    %% --- Exhaustive Manager Polling Loop (Best-Improvement) ---
    % We must read EVERY future; no short-circuiting allowed.
    while any(~[futures.Read])
        [completed_idx, local_max_beta, local_i] = fetchNext(futures);
        
        if ~isempty(completed_idx)
            block_max_betas(completed_idx) = local_max_beta;
            block_winner_local_idx(completed_idx) = local_i;
        end
    end
    
    % Global Reduction: Find the absolute best candidate across all blocks
    [global_max_beta, best_block_idx] = max(block_max_betas);
    
    found_improvement = false;
    best_sensor_ind = p_orig_front; 
    b_best = b_orig;
    
    % Verify if the absolute best candidate clears the improvement threshold
    if global_max_beta > beta_orig * sqrt(f) + 1e-9
        found_improvement = true;
        
        % Map back to the global index space
        winning_local_idx = block_winner_local_idx(best_block_idx);
        global_block_idx = idx_chunks{best_block_idx}(winning_local_idx);
        best_sensor_ind = test_ind(global_block_idx);
        
        % Recompute the winning b-vector for the Cholesky update
        b_best = L_fast \ rhs_all(:, global_block_idx);
    end
    
    %% --- swap ---
    if found_improvement
        p = [p_temp best_sensor_ind];
        L = [L_fast zeros(k-1,1); b_best' global_max_beta];
        no_swap_counter = 0; % Reset counter since an update occurred
        disp("Iteration " + iteration_num + ": Best swap executed.");
    else 
        p = circshift(p, -1);
        L = [L_fast zeros(k-1,1); b_orig' beta_orig];
        best_sensor_ind = p(end); 
        no_swap_counter = no_swap_counter + 1; 
        disp("Iteration " + iteration_num + ": No swap found.");
    end
    
    %% --- Update rhs_all ---
    rhs_all(1, :) = []; 
    idx_rm = find(test_ind == best_sensor_ind);
    rhs_all(:, idx_rm) = [];
    test_ind(idx_rm) = []; 
    
    new_col = A(p_temp(2:end), p_temp(1));
    rhs_all = [rhs_all, new_col];
    test_ind = [test_ind, p_temp(1)]; 
    
    new_row = A(best_sensor_ind, test_ind);
    rhs_all = [rhs_all; new_row];
    
    p_temp = [p_temp(2:end), best_sensor_ind];
    
    iteration_num = iteration_num + 1;
    
    if no_swap_counter >= k
        converged = true;
    end
end

ld_exchange = slogdet(K_fun(x(p,:)), sig_n)
end