% test nystrom with H2 pack
clear; clc;

rng("default")

% load("droplet_setup.mat")

addpath(genpath('/Users/jchen228/Desktop/gitGPOED/H2Pack-Matlab-master'))

global sig_n

a = 0; b = 1; % range of data
N = 500; % size of grid / number of candidate sensors
D = 1; % dimension of data
k = ceil(0.1 * N); % number of available sensors
p =  ceil(0.2 * N); % k+20; % ceil(0.4 * N);


% % design matrix (need N x D)
% x = a + ((b-a))*rand(N,D);
x = linspace(a,b,N)';

sig_f = 1;
ls = 1;

f = @(x) exp(-x).*cos(2*pi*x); 
% f = @(x) sin(2*pi*x) + 2*cos(3*pi*x) - sin(6*pi*x);

y = f(x);

sig_n = 0.02*norm(y)/sqrt(N);

y = y + sig_n*randn(size(y));

%% Kernel formation

K_fun = @(x) gaussKern(x,sig_f,ls);
K_fun_offdiag = @(x,x2) gaussKern(x,sig_f,ls,x2);

K_true = K_fun(x);

K_funMat = funMat(K_fun,K_fun,N); % K is a function handle, use K_funMat*x to form kernel 
K = K_funMat*x;

%% NysH2


[u_k,lambda_k] = Nys(K, k, p);

err_NysH2_2 = norm(K - u_k*lambda_k*u_k')/norm(K);
err_NysH2_fro = norm(K - u_k*lambda_k*u_k','fro')/norm(K,'fro')


%%  NysGKS
% tol_setting = 1e-9; % 1e-9; % epsilon tolerance or k
% KXX_hat = zeros(length(x)); % initialize approximated matrix





% tic
% % [S,F] = rpchol_tol(N,x,ls,sig_n,sig_f,tol_setting); % 
% 
% % use F to find sensors
% 
% % [u,~,~] = svd(F,0);
% [u,~] = NysH2(x, ls, sig_f, p, p);
% u_k = u(:,1:k); % take first k columns of V
%                 % gets left singular vectors of most significant singular
%                 % values
% [~,~,pick] = qr(u_k','vector'); % p permutation matrix such that u_k' * diag(p) = q * r
% % pk_rp = p(1:k); % column subset selection
% pk_nys = pick(1:k);
% KXX_s_hat = K(pk_nys,pk_nys); % K(X_*,X_*)
% 
% 
% ld_s_hat = slogdet(KXX_s_hat);
% rp_time = toc;
% 
% fk_nys = krr(x(pk_nys,:), y(pk_nys), K_funMat, K_fun_offdiag(x(pk_nys,:), x));
% fk_nys_err = norm(fk_nys - y)/norm(y);

% % find posterior variance - rewrites v_post, var_fit
% % test_ind = find(ismember(x,x_mini(pk)));
% v_post = post_var2(pk_nys, x, ls, sig_n, eta_n2, sig_f);
% var_fit = zeros(N,1); % vector of all variances 
% var_fit = sqrt(v_post);
% % var_fit(setdiff(1:N,pk_rp)) = v_post;
% 
% figure()
% plot(x_sel, h_sel, x_sel, fk_nys, 'LineWidth',1.5)
% set(gca,'FontSize',20);
% % ylim([0, 0.5])
% hold on
% scatter(x_sel(pk_nys), h_sel(pk_nys), 'r*', 'LineWidth',1.5)
% patch([x;flipud(x)], [1*var_fit;-1*flipud(var_fit)]+ [fk_nys'; flipud(fk_nys')],'k','FaceAlpha',0.1); % Prediction intervals 
% % title('With RPCholesky','FontSize',22)
% % subtitle("ld_s = " + num2str(ld_s_hat),'FontSize',16)
% xlabel('$x$','Interpreter','Latex')
% ylabel('$h(x,t)$','Interpreter','Latex')
% lgd = legend('Original', 'RP Cholesky','Location','north');
% fontsize(lgd, 16, 'points')
% subtitle("fk opt error = " + num2str(fk_nys_err))