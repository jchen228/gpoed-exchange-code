%% 6 dimensionl example - setup
clear;
% clc;
rng(0)


D = 4; N = 10000; 
% create inout grid
% x = rand(N,D);
x = lhsdesign(N,D);
k = ceil(0.3*N);



% set hyperparameters
% ls = 1.5;
% ls = 0.05*(1:D);
% ls = 0.25*ones(1,D);
% ls = 0.025;
% ls = [0.25 0.25 0.25 0.25];
% sig_f = 15; % signal variance
sig_f = 6.5;
ls = 0.15;



% generate true noisy data
% y = arrayfun(@(i) grlee09(x(i,:)), 1:N)';
y = arrayfun(@(i) zhou98(x(i,:)), 1:N)';
% sig_n = 0.002*norm(y)/sqrt(N);
sig_n = 0.05*norm(y)/sqrt(N);
% y = y + sig_n^2*randn(size(y));


% % plot 2D original function - use linspace
% [X, Y] = meshgrid(x(:,1),x(:,2));
% Z = zeros(size(X));
% grid = [reshape(X, [], 1) reshape(Y, [], 1)];
% Z = arrayfun(@(i) zhou98(grid(i,:)), 1:length(grid))';
% Z = reshape(Z, size(X));
% s = surf(X,Y,Z);


%%
% kernel function: exponential kernel
% K_fun = @(x) gaussKern(x,sig_f,ls);
% K_fun_offdiag = @(x,x2) gaussKern(x,sig_f,ls,x2);


% separable exponential kernel
K_fun = @(x) sepKern(x,sig_f,ls);
K_fun_offdiag = @(x,x2) sepKern(x,sig_f,ls,x2);


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

rnd_stats = datastats(ld_rnd');

datastats(fk_rnd_err')
% save('zhou_setup.mat')
