% test and compare different cholesky methods with and without svd step

rng("default")

% global sig_n
addpath(genpath('/Users/jchen228/Desktop/gitGPOED/H2Pack-Matlab-master'))

load droplet_setup.mat

%% GKS

tic
pk_gks = gks(K,k); % output: Kernel of placed sensors, placed sensors
x_gks = x(pk_gks,:);
K_gks = K_funMat*x_gks;
time_gks = toc;

% compute log determinant
ld_gks = slogdet(K_gks,sig_n);

fk_gks = krr(x_gks, y(pk_gks), K, pk_gks, sig_n);
fk_gks_err = norm(fk_gks - y)/norm(y);

%% RP Cholesky Algorithm
tic
[F_rp,pk_rp] = pivotedcholesky(K_fun_offdiag, x, k, 'rpcholesky');
x_rp = x(pk_rp,:);
K_rp = K_funMat*x_rp; 
time_rp = toc;

% compute log determinant
ld_rp = slogdet(K_rp,sig_n);

fk_rp = krr(x_rp, y(pk_rp), K, pk_rp, sig_n);
fk_rp_err = norm(fk_rp - y)/norm(y);

% test with gks
pk_rpgks = rpgks(F_rp,k);
x_rpgks = x(pk_rpgks,:);
K_rpgks = K_funMat*x_rpgks; 
time_rpgks = toc;

% compute log determinant
ld_rpgks = slogdet(K_rpgks,sig_n);

fk_rpgks = krr(x_rpgks, y(pk_rpgks), K, pk_rpgks, sig_n);
fk_rpgks_err = norm(fk_rpgks - y)/norm(y);

%% Greedy Pivoted Cholesky

tic
[F_piv, pk_piv] = pivotedcholesky(K_fun_offdiag, x, k, 'greedy');
x_piv = x(pk_piv,:);
K_piv = K_funMat*x_piv; 
time_piv = toc;

% compute log determinant
ld_piv = slogdet(K_piv,sig_n);

fk_piv = krr(x_piv, y(pk_piv), K, pk_piv, sig_n);
fk_piv_err = norm(fk_piv - y)/norm(y);

% test with gks 
pk_pivgks = rpgks(F_piv,k);
x_pivgks = x(pk_pivgks,:);
K_pivgks = K_funMat*x_pivgks; 
time_pivgks = toc;

% compute log determinant
ld_pivgks = slogdet(K_pivgks,sig_n);

fk_pivgks = krr(x_pivgks, y(pk_pivgks), K, pk_pivgks, sig_n);
fk_pivgks_err = norm(fk_pivgks - y)/norm(y);

%% Plot result
% 
% figure()
% hold on
% scatter(x(pk_rp),y(pk_rp),'filled')
% scatter(x(pk_piv),y(pk_piv),'filled')
% scatter(x(pk_rpgks),y(pk_rpgks),'filled')
% scatter(x(pk_pivgks),y(pk_pivgks),'filled')
% hold off
% % legend('rp','piv')
% legend('rp','piv','rpgks','pivgks')

%% D-Optimality comparison plot

figure()
histogram(ld_rnd)
set(gca,'FontSize',20);
hold on
hxl = xline(ld_gks, '-', LineWidth=2);
fontsize(hxl, 16, 'points')
hxl.Color = [0.6350 0.0780 0.1840];
hxl2 = xline(ld_rp, "--" , LineWidth=2);
fontsize(hxl2, 16, 'points')
hxl2.Color = [0.9290 0.6940 0.1250];
hxl3 = xline(ld_rpgks, ":" , LineWidth=2);
hxl3.Color = [0.8500 0.3250 0.0980];
fontsize(hxl3, 16, 'points')
hxl4 = xline(ld_piv, "-.", LineWidth=2);
hxl4.Color = [0.4940 0.1840 0.5560];
hxl4 = xline(ld_pivgks, "-.", LineWidth=2);
hxl4.Color = [0.4660 0.6740 0.1880];
fontsize(hxl4, 16, 'points')
xlabel('D-Optimality Criterion', 'Interpreter','latex')
ylabel('Counts', 'Interpreter','latex')
legend('','GKS = '+string(ld_gks),'RP = '+string(ld_rp),'RP+GKS = '+string(ld_rpgks),'Greedy Pivot = '+string(ld_piv),'Greedy Pivot + GKS = '+string(ld_pivgks),Location='northwest')
hold off

