clear; clc;

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
p = pk_gks;
% p = p
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

    ld_test = sum(log(diag(L))); % why is this here???
        
        
    col_to_add = L(2:end,1);
    L_fast = L(2:end,2:end);
    L_fast = (cholupdate(L_fast', col_to_add, '+'))';

   % col_to_add = L(idx+1:end, idx);
   % L_fast = L;
   % L_fast(idx, :) = []; % Remove row
   % L_fast(:, idx) = []; % Remove col
   % 
   % if idx <= size(L_fast, 1)
   %     L_sub = L_fast(idx:end, idx:end);
   %     % Update the sub-block
   %     R_sub = cholupdate(L_sub', col_to_add, '+');
   %     L_fast(idx:end, idx:end) = R_sub';
   % end % could be used if swaps don't take the first element each time


    
    %% --- 2. EXCHANGE SEEKING STEP ---
    % note: adds new sensor to "end" of block
    
    p_temp = p;
    p_temp(idx) = []; 

    % compute beta of unswapped selection
    rhs_orig = A(p_temp, p(idx));
    if size(rhs_orig, 2) > 1; rhs_orig = rhs_orig'; end
    b_orig = L_fast \ rhs_orig;
    beta_orig = sqrt(A(p(idx), p(idx)) - b_orig'*b_orig);

    % p_temp(idx) = []; 
    
    beta = [];
    for j = setdiff(1:n, p_temp)
        rhs = A(p_temp, j); 
        if size(rhs, 2) > 1; rhs = rhs'; end % ensures correct dimensions
        
        b = L_fast \ rhs;
        beta(j) = sqrt(A(j,j)-b'*b);
    end    

    % beta
    % beta_orig
    f = 1.025; 
    if any(beta > beta_orig*sqrt(f) + 1e-9)
        [~,i] = max(beta); 
        b = L_fast \ A(p_temp, i);
        p_temp = [p_temp i];
        p = p_temp; % update p
    
        % update A and L
        L = [L_fast zeros(k-1,1); b' beta(i)];

        % disp("swap made in round "+l)

    else 
        disp("no swap made in round "+l)
        p = circshift(p, -1);
        L = [L_fast zeros(k-1,1); b_orig' beta_orig];
    end

    % ld_swap = sum(log(diag(L)));
end

slogdet(K(p,p), sig_n)