function [p, ld_exchange, total_swaps] = exchange_sensors6(p, x, sig_n, sig_f, ls, f)
% Function that takes initial placements and improves selection by improving
% the determinant by a factor f
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
    total_swaps
end

k = length(p);
n = length(x);

% functions for getting elements of K
K_fun = @(x) gaussKern(x, sig_f, ls);
K_fun_offdiag = @(x, x2) gaussKern(x, sig_f, ls, x2);
K_diag = @(x) gaussKern(x, sig_f, ls, '1');

alpha = @(i) 1/(sig_n^2)*K_diag(1) + 1; % "on" diagonal elements
A = @(i,j) 1/(sig_n^2)*K_fun_offdiag(x(i,:), x(j,:));  % off diagonal elements

A_p = (1/(sig_n^2))*K_fun(x(p,:)) + eye(k);
L = chol(A_p, 'lower'); % Compute the Cholesky factor of the covariance matrix

%% initialize rhs_all
p_temp = p(2:end);
test_ind = setdiff(1:n, p_temp);
rhs_all = A(p_temp, test_ind); 

%% sweep loop -- continues until no swaps happen in a sweep
converged = false;
total_swaps = 0;
sweep_num = 1;

while ~converged
    swap_count_this_sweep = 0;
    
    for l = 1:k
        % Save the sensor currently being evaluated before p shifts
        p_orig_front = p(1); 
        
        % fast downdate
        col_to_add = L(2:end,1);
        L_fast = L(2:end,2:end);
        L_fast = (cholupdate(L_fast', col_to_add, '+'))';
        
        %% --- VECTORIZED EXCHANGE SEEKING STEP ---
        b_all = L_fast \ rhs_all;
        
        A_test_diag = alpha(test_ind); 
        beta_all = sqrt(A_test_diag' - sum(b_all.^2, 1));
        
        % Compute beta of unswapped selection
        orig_ind = find(test_ind == p_orig_front);
        b_orig = b_all(:, orig_ind);
        beta_orig = sqrt(alpha(p_orig_front) - b_orig'*b_orig);
        
        % Swap decision
        [max_beta, i] = max(beta_all);
        best_sensor_ind = test_ind(i); % Map back to absolute sensor ID
        
        if max_beta > beta_orig*sqrt(f) + 1e-9
            b_best = b_all(:, i); 
            p = [p_temp best_sensor_ind];
            L = [L_fast zeros(k-1,1); b_best' max_beta];
            swap_count_this_sweep = swap_count_this_sweep + 1;
        else 
            p = circshift(p, -1);
            L = [L_fast zeros(k-1,1); b_orig' beta_orig];
            
            % Overwrite with the retained original sensor
            best_sensor_ind = p(end); 
        end
        
        %% update rhs_all
        
        rhs_all(1, :) = []; %remove top row
        
        % remove the column for the sensor pushed to the back
        idx_rm = find(test_ind == best_sensor_ind);
        rhs_all(:, idx_rm) = [];
        test_ind(idx_rm) = []; 
        
        % add the new column to the end
        new_col = A(p_temp(2:end), p_temp(1));
        rhs_all = [rhs_all, new_col];
        test_ind = [test_ind, p_temp(1)]; 
        
        % add the bottom row 
        new_row = A(best_sensor_ind, test_ind);
        rhs_all = [rhs_all; new_row];
        
        % Finalize p_temp for the next iteration (or next sweep!)
        p_temp = [p_temp(2:end), best_sensor_ind];
    end
    
    % --- End of Sweep Checking ---
    disp("Swaps in sweep " + sweep_num + ": " + swap_count_this_sweep);
    total_swaps = total_swaps + swap_count_this_sweep;
    sweep_num = sweep_num + 1;
    
    if swap_count_this_sweep == 0
        converged = true;
    end
end

ld_exchange = slogdet(K_fun(x(p,:)), sig_n)
end