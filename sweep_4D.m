%% 4 dimensionl example - setup
clear;
% clc;
rng(0)


D = 4; N = 10000; 
% create inout grid
% x = rand(N,D);
x = lhsdesign(N,D);
% x = [linspace(0,1,N)' linspace(0,1,N)']; % for visualizing in 2D
k = ceil(0.3*N);



% set hyperparameters - on original blob spacing
% ls = 0.11;
% sig_f = 16; % signal variance

% set hyperparameters - on adjusted blob spacing
sig_f = 6.5;
ls = 0.15;

% generate true noisy data
% y = arrayfun(@(i) grlee09(x(i,:)), 1:N)';
y = arrayfun(@(i) zhou98(x(i,:)), 1:N)';
% sig_n = 0.002*norm(y)/sqrt(N);
sig_n = 0.05*norm(y)/sqrt(N);
% y = y + sig_n^2*randn(size(y));


% plot 2D original function - use linspace
% 
% [X, Y] = meshgrid(x(:,1),x(:,2));
% Z = zeros(size(X));
% grid = [reshape(X, [], 1) reshape(Y, [], 1)];
% Z = arrayfun(@(i) zhou98(grid(i,:)), 1:length(grid))';
% Z = reshape(Z, size(X));
% s = surf(X,Y,Z);

%% Perform sweep
% ls_sweep = 0.01:0.01:0.2; %[0.08:0.01:0.15];
% sig_f_sweep = 4:0.5:8; % 14:1:18;
% 
% n_h = length(y);
% 
% m = 1;
% n = 1;
% nugget = 1e-6;
% Ainv = (sig_n^-2)*eye(n_h);
% % sweep search
% for i = sig_f_sweep
%     K_ls = @(ls) gaussKern(sea_points,i,ls); % update sig_f in K_ls
% 
%     for j = ls_sweep
%         K_fun_offdiag = @(x,x2) gaussKern(x,i,j,x2);
%         K_fun_offdiag = @(x,x2) sepKern(x,i,j*ones(1,D),x2);
%         [F,~,p] = pivotedcholesky(K_fun_offdiag, x, k, 'rpcholesky');
%         L =  chol(F(:,1:(p-2))'*F(:,1:(p-2)) + (sig_n^2)*eye(p-2),'lower'); % K + sig_n*eye = L*L'
%         Binv = L\((sig_n^-2)*F(:,1:(p-2))');
%         alpha = Ainv - Binv'*Binv;
%         fun = -(1/2)*y'*alpha*y - ((n_h - k)/2)*log(sig_n^2)+sum(log(diag(L))) - (n_h/2)*log(2*pi);
%         obj(m,n) = fun;
%         n = n+1;
%         i
%         j 
%     end
%     n = 1; % reset n
%     m = m+1;
% end
% 
% 
% maximum = max(obj,[],"all");
% [m_opt,n_opt]=find(obj==maximum)
% 
% minimum = min(min(obj));
% [m_test,n_test]=find(obj==minimum)
% 
% sig_f = sig_f_sweep(m_opt)
% ls = ls_sweep(n_opt)

%%
% kernel function: exponential kernel
K_fun = @(x) gaussKern(x,sig_f,ls);
K_fun_offdiag = @(x,x2) gaussKern(x,sig_f,ls,x2);


% separable exponential kernel - same as Gaussian Kernel when all ls vector
%   elements are same
% K_fun = @(x) sepKern(x,sig_f,ls);
% K_fun_offdiag = @(x,x2) sepKern(x,sig_f,ls,x2);

% K_funMat = funMat(K_fun,K_fun,N); % K is a function handle, use K_funMat*x to form kernel 

%% random sensor placements
rr = 10000;

ld_rnd = zeros(1,rr);
fk_rnd = zeros(rr,N);
fk_rnd_err = zeros(1,rr);
for i = 1:rr
    pk_rnd = datasample(1:N,k,'Replace',false); % generate k random samples of values up to N without replacement
    x_rnd = x(pk_rnd,:);
    K_rnd = K_fun(x_rnd);
    ld_rnd(i) = slogdet(K_rnd,sig_n);
    fk_rnd(i,:) = krr(x_rnd, y(pk_rnd), K_fun_offdiag, pk_rnd, sig_n, x);
    fk_rnd_err(i) = norm(fk_rnd(i,:)'-y)/norm(y);
end

ld_rnd_stats = datastats(ld_rnd');

fk_rnd_stats = datastats(fk_rnd_err');
% save('zhou_setup_moved_hp.mat')
