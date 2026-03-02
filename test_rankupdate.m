% assumes K is permutes such that selection is the (1,1) block


clear; clc;
rng(0)

% % toy example
% n = 8; k = 4;
% K = pascal(n); % test K on a small kernel instead
% % permute K
% perm = randperm(n)
% % Add a small nugget to guarantee strict positive definiteness
% K = K(perm,perm) + 1e-9*eye(n);
% p = 1:k;
% A = K(p,p);
% L = chol(A + 1e-9*eye(A), 'lower'); % Recompute Cholesky factor

% load out_pde.dat
% % n = 30;
% % 
% % x = out_pde(1:20:(20*n),1);
% % h = out_pde(1:20:(20*n),2);
% % t = out_pde(1:20:(20*n),3);
% 
% n = 6001;
% 
% x = out_pde(:,1);
% h = out_pde(:,2);
% t = out_pde(:,3);
% 
% ls = 0.5;
% sig_f = 1; % signal variance
% sig_n = 0.002*norm(h)/sqrt(n);
% 
% x_all = reshape(x, n, []);
% x = x_all(:,1); % rewrites variable x to be the domain of selection
% h_all = reshape(h, n, []);
% t_all = reshape(t, n, []);
% t = t_all(1,:); % rewrites variable t to be the time steps
% 
% y_all = (h_all + sig_n^2*randn(size(h_all)))';
% y = h;
% 
% K_fun = @(x) gaussKern(x,sig_f,ls);
% K_fun_offdiag = @(x,x2) gaussKern(x,sig_f,ls,x2);
% 
% K = K_fun(x);
% k = 30;

load droplet_setup.mat
n = N;
k = 30;

% initial gks selection
tic
pk_gks = gks(K,k); % output: Kernel of placed sensors, placed sensors
x_gks = x(pk_gks,:);
K_gks = K(pk_gks,pk_gks);
time_gks = toc;


% compute log determinant
ld_gks = slogdet(K_gks,sig_n)

p = pk_gks;
% A = (1/(sig_n^2))*K_gks+eye(k); % covariance matrix of selected sensors
% L = chol(A + 1e-9*eye(k), 'lower'); % Compute the Cholesky factor of the covariance matrix

A = (1/(sig_n^2))*K+eye(n); % covariance matrix of selected sensors
L = chol(A(p,p), 'lower'); % Compute the Cholesky factor of the covariance matrix


for l = 1:k
    %% --- FAST DOWNDATE (O(k^2) Method) ---
    % gets the cholesky factor of downdated matrix, L_fast ready to work sift through
    % n-k unselected sensors

    idx = 1; % since swaps are added to end of selection index list, just take first element each time

    % remove row and col l from A
    % B = A(p,p);
    % B(idx,:) = [];
    % B(:,idx) = [];
    % 
    % L_B = chol(B,'lower'); % for testing

    ld_test = sum(log(diag(L)));
        
        
    col_to_add = L(idx+1:end, idx);
    L_fast = L;
    L_fast(idx, :) = []; % Remove row
    L_fast(:, idx) = []; % Remove col
    
    if idx <= size(L_fast, 1)
        L_sub = L_fast(idx:end, idx:end);
        % Update the sub-block
        R_sub = cholupdate(L_sub', col_to_add, '+');
        L_fast(idx:end, idx:end) = R_sub';
    end
    
    %% --- 2. EXCHANGE SEEKING STEP ---
    % note: adds new sensor to "end" of block
    
    p_temp = p;
    p_temp(idx) = []; 
    
    % later: look into block linear solves

    % compute beta of unswapped selection
    rhs_orig = A(p_temp, p(idx));
    if size(rhs_orig, 2) > 1; rhs_orig = rhs_orig'; end
    b_orig = L_fast \ rhs_orig;
    beta_orig = sqrt(A(p(idx), p(idx)) - b_orig'*b_orig);

    % p_temp(idx) = []; 
    
    beta = [];
    for j = setdiff(1:n, p_temp)
        rhs = K(p_temp, j); 
        if size(rhs, 2) > 1; rhs = rhs'; end % ensures correct dimensions
        
        b = L_fast \ rhs;
        beta(j) = sqrt(A(j,j)-b'*b);
    end    

    if any(beta > beta_orig)
        [~,i] = max(beta); 
        b = L_fast \ A(p_temp, i);
        p_temp = [p_temp i];
        p = p_temp; % update p
    
        % update A and L
        L = [L_fast zeros(k-1,1); b' beta(i)];

        disp("swap made in round "+l)

    else 
        disp("no swap made in round "+l)
    end

    ld_swap = sum(log(diag(L)));
    
end

ld_swap