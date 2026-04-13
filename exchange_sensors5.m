function [p, ld_exchange, swap_count, ld_orig] = exchange_sensors5 (p, x, sig_n, sig_f, ls, f, tol, trunc)
% function takes initial placements p and performs a neighborhood search 
% and performs a swap if the determinant improves by a factor of f
arguments (Input)
    p
    x
    sig_n
    sig_f (1,1) double = 1.0 % Default value if not provided
    ls (1,1) double = 1.0 % Default value if not provided
    f (1,1) double = 1.0 % Default value if not provided
    tol (1,1) double = 0.2 % Default value if not provided
    trunc (1,1) double = 0 % Default value if not provided
end

arguments (Output)
    p
    ld_exchange
    swap_count
    ld_orig
end

swap_count = 0;

k = length(p);
n = length(x);

% functions for getting elements of K
K_fun = @(x) gaussKern(x,sig_f,ls);
K_fun_offdiag = @(x,x2) gaussKern(x,sig_f,ls,x2);
K_diag = @(x) gaussKern(x,sig_f,ls,'1');

alpha = @(i) 1/(sig_n^2)*K_diag(x(i,:))+1; % "on" diagonal elements
A = @(i,j) 1/(sig_n^2)*K_fun_offdiag(x(i,:),x(j,:));  % off diagonal elements

A_p = (1/(sig_n^2))*K_fun(x(p,:))+eye(k);
L = chol(A_p, 'lower'); % Compute the Cholesky factor of the covariance matrix

ld_orig = slogdet(K_fun(x(p,:)), sig_n); % orig ld


for l = 1:k
    %% --- FAST DOWNDATE (O(k^2) Method) ---
    % gets the cholesky factor of downdated matrix

    idx = 1; % since swaps are added to end of selection index list, just take first element each time
        
    col_to_add = L(2:end,1);
    L_fast = L(2:end,2:end);
    L_fast = (cholupdate(L_fast', col_to_add, '+'))';

    %% --- 2. VECTORIZED EXCHANGE SEEKING STEP ---
    % note: adds new sensor to "end" of block
    p_temp = p;
    p_temp(1) = []; 

    % get the p(1)th row of K
    row_p1 = K_fun_offdiag(x(p(1),:), x);
    
    % find elements near sig_f^2 --- what is "near"? 
    % do this before filtering out already-selected sensors to avoid
    % messing with the indices
    test_ind = find(row_p1 > sig_f^2*(1-tol));    
    
    % if length(test_ind) is too long, sort to find the nearest sensors
    % avoided sorting row_p1 in case n is very very large and sorting is
    % somewhat expensive
    switch trunc
        case 0
            test_ind = find(row_p1 > sig_f^2*(1-tol));    
        case 1
            if length(test_ind) > k
                A_test_ind = A(p(1), test_ind);
                B = [A_test_ind; test_ind];
                [~,inx]=sort(B(1,:), "descend");
                B = B(:,inx);
                test_ind = B(2, 1:k);
                disp("truncated test_ind in round " + l)
            end
    end

    % compute b vectors 
    test_ind = setdiff(test_ind, p_temp); % filter out already-selected sensors
    rhs_all = A(p_temp,test_ind);
    b_all = L_fast \ rhs_all;

    % A_diag = diag(A);
    A_test_diag = alpha(test_ind); %  A_diag(test_ind); %  % Only get diagonals for candidates

    % compute all betas for this sweep
    beta_all = sqrt(A_test_diag' - sum(b_all.^2, 1));

    % compute beta of unswapped selection
    rhs_orig = A(p_temp, p(idx));
    if size(rhs_orig, 2) > 1; rhs_orig = rhs_orig'; end
    
    b_orig = L_fast \ A(p_temp, p(1));
    beta_orig = sqrt(alpha(p(idx)) - b_orig'*b_orig);

    % swap decision
    [max_beta, i] = max(beta_all);
    best_sensor_ind = test_ind(i); % Map back to absolute sensor ID

    % Add a tiny tolerance (1e-9) to prevent swapping identical values
    if max_beta > beta_orig*sqrt(f) + 1e-9
        % Grab the correct b vector using the relative index
        b_best = b_all(:, i); 
        
        p = [p_temp best_sensor_ind];
        L = [L_fast zeros(k-1,1); b_best' max_beta];
        disp("swap made in round "+l)
        swap_count = swap_count +1;
    else 
        % Push original sensor to the back of the queue
        % disp("no swap made in round "+l)
        p = circshift(p, -1);
        L = [L_fast zeros(k-1,1); b_orig' beta_orig(1)];
    end

end
disp("total swaps: "+swap_count)
ld_exchange = slogdet(K_fun(x(p,:)), sig_n)


end


