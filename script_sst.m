% SST Script

clear; clc;
rng("default")

% global sig_n
addpath(genpath('/Users/jchen228/Desktop/gitGPOED/H2Pack-Matlab-master'))

load sst_setup.mat
% load sst_setup_geo.mat
y = sst_sel_reshape_sea';

p = k + 20;



%% NysGKS Algorithm
    
    tic
    [u,~] = Nys(Kf, k, p);
    pk_nys = NysGKS(u,k);
    x_nys = sea_points(pk_nys,:);
    K_nys = rowcolkern(Kf,pk_nys,pk_nys);
    time_nys = toc;

    % compute log determinant
    ld_nys = slogdet(K_nys,sig_n);

    % load('sst_nys_selection.mat');
% 
    fk_nys = krr(x_nys, y(pk_nys), Kf, pk_nys, sig_n);
    fk_nys_err = norm(fk_nys - y)/norm(y);

% 
% combine both land mask and sea data
fk_nys_recon = zeros(size(sst_sel_reshape));
fk_nys_recon(sea_ind) = fk_nys;
fk_nys_recon = reshape(fk_nys_recon,size(sst_sel));

% Plot reconstruction
plot_sst_map(fk_nys_recon, pk_nys)

%% RPCholesky Algorithm
    
    tol = k;

    tic
    % [~,F] = rpchol_tol(K_fun_offdiag,tol,sea_points);
    F = pivotedcholesky(K_fun_offdiag, sea_points, k, 'rpcholesky');
    pk_rp = rpgks(F,k);
    x_rp = sea_points(pk_rp,:);
    K_rp = K_fun(x_rp);
    time_rp = toc;
    
    % compute log determinant
    ld_rp = slogdet(K_rp, sig_n);
    
    fk_rp = krr(x_rp, y(pk_rp), K_fun_offdiag, pk_rp, sig_n, sea_points);
    fk_rp_err = norm(fk_rp - y)/norm(y);


% combine both land mask and sea data
fk_rp_recon = zeros(size(sst_sel_reshape));
fk_rp_recon(sea_ind) = fk_rp;
fk_rp_recon = reshape(fk_rp_recon,size(sst_sel));

% Plot reconstruction
plot_sst_map(fk_rp_recon,pk_rp)

v_post = post_var2(pk_rp, sea_points, sig_n,  K_fun_offdiag);
var_fit = zeros(N,1); % vector of all variances 
var_fit(sea_ind) = v_post;
var_fit = reshape(var_fit,size(sst_sel));

plot_sst_map(var_fit, pk_rp)

%% Greedy Pivoted Cholesky Algorithm

    tic
    F_piv = pivotedcholesky(K_fun_offdiag, sea_points, k, 'greedy');
    pk_piv = rpgks(F_piv,k);
    x_piv = sea_points(pk_piv,:);
    K_piv = K_fun(x_piv);
    time_piv = toc;
    
    % compute log determinant
    ld_piv = slogdet(K_piv, sig_n);
    
    fk_piv = krr(x_piv, y(pk_piv), K_fun_offdiag, pk_piv, sig_n, sea_points);
    fk_piv_err = norm(fk_piv - y)/norm(y);


% combine both land mask and sea data
fk_piv_recon = zeros(size(sst_sel_reshape));
fk_piv_recon(sea_ind) = fk_piv;
fk_piv_recon = reshape(fk_piv_recon,size(sst_sel));

% Plot reconstruction
plot_sst_map(fk_piv_recon,pk_piv)

%% Efficient Greedy Selection
tic
[pk_g, ld_g] = greedydopt6(K_fun_offdiag,sea_points,k,sig_n);
x_g = sea_points(pk_g,:);
greedy_time = toc;



%% compare to Q-DEIM

sst_reshape = reshape(sst, [N, length(time)]);
Train = sst_reshape(sea_ind,selection_time-261:selection_time-1); % change training time first 16 years if matching Tdeim paper

meansst = mean(Train,2); % average temperature at each point over the training data
Train = bsxfun(@minus,Train,meansst); % remove mean


[Phi, S, V] = svd(Train, 'econ');

r = k;
[~,~,pivot] = qr(Phi(:,1:r)','vector');
sensors = pivot(1:r);

x_test = sst_reshape(sea_ind,selection_time)-meansst;

xls = Phi(:,1:r)*(Phi(sensors,1:r)\x_test(sensors));
% xls = Phi*Phi'*x_test;
qdeim_err = norm(xls - x_test)/norm(x_test+meansst);
qdeim_anomaly_err = norm(xls - x_test)/norm(x_test)

% semilogy(1:size(S_q,1), diag(S_q), "*")

KXX_q = K_fun(sea_points(sensors,:));
ld_q = slogdet(KXX_q,sig_n);

fk_q = zeros(N,1);
fk_q(sea_ind) = xls; % +meansst; % not adding back mean to capture anomaly
fk_q_recon = reshape(fk_q, [size(sst_sel)]);

% Q deim plot
plot_sst_map(fk_q_recon, sensors)
title("Anomalies")
