% Droplet problem set up

clear all; close all; 
% clc; 
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

% ls = 0.679421943058468; % optimized
% ls = 0.471127877613764;
ls = 0.5;
sig_f = 1; % signal variance
sig_n = 0.002*norm(h_sel)/sqrt(N);

% eta_n2 = eta^2;
% lambda = eta_n2;

% add noise to data
y_all = (h_all + sig_n^2*randn(size(h_all)))';
y = h_sel;

%% optimize for sigma over all possible sensor positions using data from same time step

% % optimizing over all sensors and one time step prior
% train_time = selection_time;
% snap_prior = [train_time*6001-6000 train_time*6001];
% total_prior_steps = train_time;
% x_train = out_pde(snap_prior(1):snap_prior(2), 1);
% h_train = out_pde(snap_prior(1):snap_prior(2), 2)';
% 
% x_h = x;
% n_h = N;
% y_h = h_train;
% 
% 
% % % do a cholesky decomp
% % L = @(q) chol(eta_n2*K_y(q) + eye(n_h),'lower'); % K + sig_n*eye = L*L'
% % alpha = @(q) L(q)'\(L(q)\y_h');
% 
% options = optimset('PlotFcns',@optimplotfval,'Display','iter');
% 
% % fun = @(ls) -sum(log(diag(K_y(ls))));
% % fun = @(q) -(-(1/2)*y_h*alpha(q) - sum(log(diag(K_y(q)))) - (n_h/2)*log(2*pi));
% 
% 
% K_ls = @(ls) gaussKern(x,sig_f,ls);
% % fun = @(ls) -slogdet(K_ls(ls), sig_n);
% 
% L = @(q) chol((1/(sig_n^2))*K_ls(q) + eye(n_h),'lower'); % K + sig_n*eye = L*L'
% % alpha = @(q) L(q)'\(L(q)\y_h');
% % fun = @(q) -(-(1/2)*y_h*alpha(q) - sum(log(diag(K_ls(q)))) - (n_h/2)*log(2*pi));
% % fun = @(q) -(-(1/2)*y_h*alpha(q) - slogdet(K_ls(q), sig_n) - (n_h/2)*log(2*pi));
% 
% alpha = @(q) (L(q)\y_h');
% fun = @(q) -(-(1/2)*(alpha(q)'*alpha(q)) - sum(log(diag(L(q)))) - (n_h/2)*log(2*pi));
% 
% 
% 
% x0 = ls;
% lb = 1e-2*ones(size(x0)); ub = 3; 
% % ub = Inf*ones(size(x0)); 
% % q_opt = fminsearch(fun,x0,options);
% q_opt = fmincon(fun,x0,[],[],[],[],lb,ub,[],options);
% 
% ls = q_opt % length scale


%% form kernel
% KXX = kxx(x,ls,sig_n,sig_f);
K_fun = @(x) gaussKern(x,sig_f,ls);
K_fun_offdiag = @(x,x2) gaussKern(x,sig_f,ls,x2);

K_funMat = funMat(K_fun,K_fun,N); % K is a function handle, use K_funMat*x to form kernel 

K = K_fun(x); % fully form kernel
%% Greedy

tic
[pk_g, ld_g] = greedydopt6(K_fun_offdiag,x,50,sig_n);
toc

i = 1;
% for k = 20:50
for k = 30
    x_g = x(pk_g(1:k),:);
    K_g = K_fun(x_g); 
    % time_g = toc;

    fk_g = krr(x_g, y(pk_g(1:k)), K, pk_g(1:k), sig_n);
    fk_g_err(i) = norm(fk_g - y)/norm(y);

    i = i+1;
end


%% 10000 realizations of random
rr = 10000;

ld_rnd = zeros(1,rr);
fk_rnd = zeros(rr,N);
fk_rnd_err = zeros(1,rr);
for i = 1:rr
    pk_rnd = datasample(1:N,30,'Replace',false); % generate k random samples of values up to N without replacement
    x_rnd = x(pk_rnd,:);
    K_rnd = K_fun(x_rnd);
    ld_rnd(i) = slogdet(K_rnd,sig_n);
    fk_rnd(i,:) = krr(x_rnd, y(pk_rnd), K, pk_rnd, sig_n);
    fk_rnd_err(i) = norm(fk_rnd(i,:)'-y)/norm(y);
end

rnd_stats = datastats(ld_rnd')

datastats(fk_rnd_err')


%%

% save('droplet_setup.mat')