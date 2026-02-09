% hyperparameter selection - Droplet

clear all; close all; clc; 
rng("default")

load out_pde.dat

N = 6001;

x = out_pde(:,1);
h = out_pde(:,2);
t = out_pde(:,3);

x_all = reshape(x, N, []);
x = x_all(:,1); % rewrites variable x to be the domain of selection
h_all = reshape(h, N, []);
t_all = reshape(t, N, []);
t = t_all(1,:); % rewrites variable t to be the time steps

selection_time = 15;
snap = [selection_time*6001-6000 selection_time*6001];
total_post_steps = length((selection_time):length(t));

x_sel = out_pde(snap(1):snap(2), 1);
h_sel = out_pde(snap(1):snap(2), 2);

k = 30;

%% hyperparameters
ls = 0.471127877613764;
sig_f = 1; % signal variance
sig_n = 0.002*norm(h_sel)/sqrt(N);

% set y variable. assumes noise is already added
y = h_sel;

%% hyperparameter sweep
ls_sweep = 0.4:0.1:0.6;
sig_f_sweep = 0.2:0.05:0.4;

%% optimize for sigma over all possible sensor positions using data from same time step

% optimizing over all sensors and one time step prior
train_time = selection_time;
snap_prior = [train_time*6001-6000 train_time*6001];
total_prior_steps = train_time;
x_train = out_pde(snap_prior(1):snap_prior(2), 1);
h_train = out_pde(snap_prior(1):snap_prior(2), 2)';

x_h = x;
n_h = N;
y_h = h_train;

% % % maximize log marginal liklihood
% K_ls = @(ls) gaussKern(x,sig_f,ls);
% L = @(q) chol(K_ls(q) + (sig_n^2)*eye(n_h),'lower'); % K + sig_n*eye = L*L'
% alpha = @(q) (L(q)\y_h');
% fun = @(q) -(-(1/2)*(alpha(q)'*alpha(q)) - sum(log(diag(L(q)))) - (n_h/2)*log(2*pi));

%%

m = 1;
n = 1;
% sweep search
for i = sig_f_sweep
    K_ls = @(ls) gaussKern(x,i,ls); % update sig_f in K_ls
    for j = ls_sweep
        L = @(q) chol(K_ls(q) + (sig_n^2)*eye(n_h),'lower'); % K + sig_n*eye = L*L'
        alpha = @(q) (L(q)\y_h');
        fun = @(q) -(1/2)*(alpha(q)'*alpha(q)) - sum(log(diag(L(q)))) - (n_h/2)*log(2*pi);
        % fun = @(q)  - sum(log(diag(L(q))));
        obj(m,n) = fun(j);
        n = n+1;
    end
    n = 1; % reset n
    m = m+1;
end

% [M,I] = max(obj)

[maximum, index] = max(obj,[],"all");
[m_opt,n_opt]=find(obj==maximum)

minimum = min(min(obj));
[m_test,n_test]=find(obj==minimum)

sig_f = sig_f_sweep(m_opt)
ls = ls_sweep(n_opt)

% sig_f = sig_f_sweep(m_test)
% ls = ls_sweep(n_test)

%% Test


 
K_fun = @(x) gaussKern(x,sig_f,ls);
K_fun_offdiag = @(x,x2) gaussKern(x,sig_f,ls,x2);

K_funMat = funMat(K_fun,K_fun,N); % K is a function handle, use K_funMat*x to form kernel 

K = K_fun(x); % fully form kernel

tic
pk_gks = gks(K,k); % output: Kernel of placed sensors, placed sensors
x_gks = x(pk_gks,:);
K_gks = K_funMat*x_gks;
time_gks = toc;

% compute log determinant
ld_gks = slogdet(K_gks,sig_n)

fk_gks = krr(x_gks, y(pk_gks), K, pk_gks, sig_n);
fk_gks_err = norm(fk_gks - y)/norm(y)

% %% built in attempt
% 
% gprMdl2 = fitrgp(x,y,'KernelFunction','squaredexponential',...
%     'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
%     struct('Optimizer','gridsearch'));

