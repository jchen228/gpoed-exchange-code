clear; clc;
% rng("default")

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
% n = 50;
% 
% x = out_pde(1:20:(20*n),1);
% h = out_pde(1:20:(20*n),2);
% t = out_pde(1:20:(20*n),3);

% n = 6001;
% N = n;
% 
% % x = out_pde(:,1);
% % h = out_pde(:,2);
% % t = out_pde(:,3);
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
% k = 10;

load droplet_setup.mat

n = N;
k = 30;

%% test gks initial placement
tic
pk_gks = gks(K,k); % output: Kernel of placed sensors, placed sensors
x_gks = x(pk_gks,:);
K_gks = K(pk_gks,pk_gks);
time_gks = toc;

% compute log determinant
ld_gks = slogdet(K_gks,sig_n)

%% test against random placement
pk_rnd = datasample(1:N,30,'Replace',false); % generate k random samples of values up to N without replacement
x_rnd = x(pk_rnd,:);
K_rnd = K_fun(x_rnd);
ld_rnd = slogdet(K_rnd,sig_n)

%% test against greedy placement

[pk_g, ld_g] = greedydopt6(K_fun_offdiag,x,k,sig_n);
ld_g(end)

%% Swapping Step

% p = [1:15 5987:6001];
p = pk_g;
% p = p

K_diag = @(x) gaussKern(x,sig_f,ls,'1');

alpha = @(i) 1/(sig_n^2)*K_diag(x(i,:))+1; % "on" diagonal elements
A = @(i,j) 1/(sig_n^2)*K_fun_offdiag(x(i,:),x(j,:));  % off diagonal elements

A_p = (1/(sig_n^2))*K_fun(x(p,:))+eye(k);
L = chol(A_p, 'lower'); % Compute the Cholesky factor of the covariance matrix


for l = 1:k
    %% --- FAST DOWNDATE (O(k^2) Method) ---
    % gets the cholesky factor of downdated matrix, L_fast ready to work sift through
    % n-k unselected sensors

    idx = 1; % since swaps are added to end of selection index list, just take first element each time

    ld_test = sum(log(diag(L)));
        
        
    col_to_add = L(2:end,1);
    L_fast = L(2:end,2:end);
    L_fast = (cholupdate(L_fast', col_to_add, '+'))';

    
   
    %% --- 2. VECTORIZED EXCHANGE SEEKING STEP ---
    % note: adds new sensor to "end" of block
    p_temp = p;
    p_temp(1) = []; 

    test_ind = setdiff(1:n, p_temp);

    % compute all b vectors for this sweep to make more efficient
    rhs_all = A(p_temp,test_ind);
    b_all = L_fast \ rhs_all;

    % A_diag = diag(A);
    A_test_diag = alpha(test_ind); %  A_diag(test_ind); %  % Only get diagonals for candidates

    % compute all betas for this sweep
    beta_all = sqrt(A_test_diag' - sum(b_all.^2, 1));

    % compute beta of unswapped selection
    rhs_orig = A(p_temp, p(idx));
    if size(rhs_orig, 2) > 1; rhs_orig = rhs_orig'; end
    
    orig_ind = find(test_ind == p(1));
    b_orig = b_all(:,orig_ind);
    beta_orig = sqrt(alpha(p(idx)) - b_orig'*b_orig);

    % swap decision
    [max_beta, i] = max(beta_all);
    best_sensor_ind = test_ind(i); % Map back to absolute sensor ID

    % Add a tiny tolerance (1e-9) to prevent swapping identical values
    if max_beta > beta_orig + 1e-9 
        % Grab the correct b vector using the relative index
        b_best = b_all(:, i); 
        
        p = [p_temp best_sensor_ind];
        L = [L_fast zeros(k-1,1); b_best' max_beta];
        % disp("swap made in round "+l)
    else 
        % Push original sensor to the back of the queue
        disp("no swap made in round "+l)
        p = circshift(p, -1);
        L = [L_fast zeros(k-1,1); b_orig' beta_orig];
    end

end

slogdet(K(p,p), sig_n)
