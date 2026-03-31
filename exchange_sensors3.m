function [p,ld_exchange, swap_count, ld_orig] = exchange_sensors3(p,x,sig_n,sig_f,ls,f)
% Function that takes initial placements and improves selection by improving
%  the determinant by a factor f and tests sensors further away from each
%  other
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
    % gets the cholesky factor of downdated matrix, L_fast ready to work sift through
    % n-k unselected sensors

    idx = 1; % since swaps are added to end of selection index list, just take first element each time
        
    col_to_add = L(2:end,1);
    L_fast = L(2:end,2:end);
    L_fast = (cholupdate(L_fast', col_to_add, '+'))';

    %% --- 2. VECTORIZED EXCHANGE SEEKING STEP ---
    % note: adds new sensor to "end" of block
    p_temp = p;
    p_temp(1) = []; 

    %test_ind = setdiff(1:n, p_temp);

    % compute b vectors for points on the mesh (k, 3k, 5k, \dots, 6001-k)
    test_candidates = k+1:2*k:n-k;
    test_candidates = setdiff(test_candidates, p_temp); % filter out already-selected sensors
    rhs_mesh = A(p_temp, test_candidates);
    b_mesh = L_fast \ rhs_mesh;


    A_test_diag = alpha(test_candidates); % Only get diagonals for candidates

    % compute all betas for this mesh
    beta_mesh = sqrt(A_test_diag' - sum(b_mesh.^2, 1));
    
    % find biggest beta and perform search of radius k around it
    [~, mmi_idx] = max(beta_mesh);
    center_node = test_candidates(mmi_idx); 
    
    search_nodes = (center_node - k) : (center_node + k);
    search_nodes = setdiff(search_nodes, p_temp);
    
    rhs_search = A(p_temp, search_nodes);
    b_search = L_fast \ rhs_search;
    
    A_search_diag = alpha(search_nodes);
    beta_search = sqrt(A_search_diag' - sum(b_search.^2, 1)); 


   % orig_ind = find(test_ind == p(1));
    b_orig = L_fast \ A(p_temp, p(1));
    beta_orig = sqrt(alpha(p(idx)) - b_orig'*b_orig);

    % swap decision
    [max_beta, i] = max(beta_search);
    best_sensor_ind = search_nodes(i); % Map back to absolute sensor ID

    % Add a tiny tolerance (1e-9) to prevent swapping identical values
    if max_beta > beta_orig*sqrt(f) + 1e-9
        % Grab the correct b vector using the relative index
        b_best = b_search(:, i); 
        
        p = [p_temp best_sensor_ind];
        L = [L_fast zeros(k-1,1); b_best' max_beta];
        % disp("swap made in round "+l)
        swap_count = swap_count +1;
    else 
        % Push original sensor to the back of the queue
        % disp("no swap made in round "+l)
        p = circshift(p, -1);
        L = [L_fast zeros(k-1,1); b_orig' beta_orig];
    end

end
% disp("total swaps: "+swap_count)
ld_exchange = slogdet(K_fun(x(p,:)), sig_n);
end