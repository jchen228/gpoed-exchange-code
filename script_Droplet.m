% script for droplet unit testing
% clear; clc;
rng("default")

% global sig_n
addpath(genpath('/Users/jchen228/Desktop/gitGPOED/H2Pack-Matlab-master'))

% load droplet_setup.mat

% plot settings
lw = 3; % adjust linewidth

% j = 20:50;
j = 30;
i = 1;
for k = j
    p = k + 10;

    %% Conceptualized GKS Algorithm
    tic
    pk_gks = gks(K,k); % output: Kernel of placed sensors, placed sensors
    x_gks = x(pk_gks,:);
    K_gks = K_funMat*x_gks;
    time_gks = toc;

    % compute log determinant
    ld_gks(i) = slogdet(K_gks,sig_n);

    fk_gks = krr(x_gks, y(pk_gks), K, pk_gks, sig_n);
    fk_gks_err(i) = norm(fk_gks - y)/norm(y);

    %% find posterior variance
    % % test_ind = find(ismember(x,x_gks));
    % v_post = post_var2(pk_gks, x, sig_n, K_fun_offdiag);
    % var_fit = zeros(N,1); % vector of all variances 
    % var_fit = sqrt(v_post);
    % 
    % % plot
    % figure()
    % plot(x_sel,h_sel,x_sel, fk_gks, '--','LineWidth',lw)
    % set(gca,'FontSize',36);
    % % ylim([0, 0.5])
    % hold on
    % % scatter(x_sel, h_sel,'k','LineWidth',0.5)
    % % xline(x(pk))
    % scatter(x_gks, h_sel(pk_gks),'r*','LineWidth',lw)
    % xlabel('$x$','Interpreter','Latex')
    % ylabel('$h(x,t)$','Interpreter','Latex')
    % patch([x;flipud(x)], [var_fit;flipud(-1*var_fit)]+ [fk_gks; flipud(fk_gks)],'k','FaceAlpha',0.1); % Prediction intervals 
    % % title('Without RPCholesky','FontSize',22)
    % % subtitle("ld_s = " + num2str(ld_gks),'FontSize',16)
    % legend('Original', 'GKS','Location','north')
    % % subtitle("fk opt error = " + num2str(fk_rp_err))


    %% NysGKS Algorithm

    tic
    [u,~] = Nys(K, k, p);
    pk_nys = NysGKS(u,k);
    x_nys = x(pk_nys,:);
    K_nys = K_funMat*x_nys;
    time_nys = toc;

    % compute log determinant
    ld_nys(i) = slogdet(K_nys,sig_n);

    fk_nys = krr(x_nys, y(pk_nys), K, pk_nys, sig_n);
    fk_nys_err(i) = norm(fk_nys - y)/norm(y);

    %% RP Cholesky Algorithm

    tic
    [F, pk_rp] = pivotedcholesky(K_fun_offdiag, x, k, 'rpcholesky');
    pk_rp = rpgks(F,k);
    x_rp = x(pk_rp,:);
    K_rp = K_funMat*x_rp; 
    time_rp = toc;

    % compute log determinant
    ld_rp(i) = slogdet(K_rp,sig_n);

    fk_rp = krr(x_rp, y(pk_rp), K, pk_rp, sig_n);
    fk_rp_err(i) = norm(fk_rp - y)/norm(y);

    %% Greedy Pivoted Cholesky

    tic
    [F_piv, pk_piv] = pivotedcholesky(K_fun_offdiag, x, k, 'greedy');
    pk_piv = rpgks(F_piv,k);
    x_piv = x(pk_piv,:);
    K_piv = K_funMat*x_piv; 
    time_piv = toc;

    % compute log determinant
    ld_piv(i) = slogdet(K_piv,sig_n);

    fk_piv = krr(x_piv, y(pk_piv), K, pk_piv, sig_n);
    fk_piv_err(i) = norm(fk_piv - y)/norm(y);

    i = i+1;
end

%% plot Dopt and RelErr as k increases

figure()
plot(j,ld_g(j),'--','LineWidth',lw)
hold on
plot(j,ld_gks,'LineWidth',lw)
plot(j,ld_rp,'-.','LineWidth',lw)
plot(j, ld_nys,':','LineWidth',lw)
plot(j, ld_piv,'--o','LineWidth',lw)
set(gca,'FontSize',36);
xlabel('$k$','Interpreter','Latex')
ylabel('D-Opt','Interpreter','Latex')
% title('Idealized Algorithm as $k$ increases','FontSize',22,'Interpreter','Latex')
legend('Greedy','Conceptualized GKS', 'RP Cholesky', 'NysGKS','Pivoted Cholesky', Location='southeast')
hold off
%%

figure()
plot(j,fk_g_err,'--','LineWidth',lw)
hold on
plot(j,fk_gks_err,'LineWidth',lw)
plot(j,fk_rp_err,'-.','LineWidth',lw)
plot(j, fk_nys_err,':','LineWidth',lw)
plot(j, fk_piv_err,'--o','LineWidth',lw)
set(gca,'FontSize',36);
xlabel('$k$','Interpreter','Latex')
ylabel('Relative Reconstruction Error','Interpreter','Latex')
% title('Idealized Algorithm as $k$ increases','FontSize',22,'Interpreter','Latex')
legend('Greedy','Conceptualized GKS', 'RP Cholesky', 'NysGKS','Pivoted Cholesky', Location='northeast')
% legend
hold off

%% compare Cholesky with and without GKS
i = 11;

figure()
histogram(ld_rnd)
set(gca,'FontSize',36);
hold on

hxl2 = xline(ld_gks(i), "--" , LineWidth=lw);
fontsize(hxl2, 16, 'points')
hxl2.Color = [ 0    0.4470    0.7410];
hxl = xline(ld_rp(i), '-', LineWidth=lw);
fontsize(hxl, 16, 'points')
hxl.Color = [0.8500    0.3250    0.0980];
hxl3 = xline(ld_rp_nogks, ":" , LineWidth=lw);
hxl3.Color = [ 0.9290    0.6940    0.1250];
fontsize(hxl3, 16, 'points')
hxl4 = xline(ld_piv(i), "-.", LineWidth=lw);
hxl4.Color = [0.4940    0.1840    0.5560];
hxl4 = xline(ld_piv_nogks, "--o", LineWidth=lw);
hxl4.Color = [0.4660    0.6740    0.1880];
fontsize(hxl4, 16, 'points')
xlabel('D-Optimality Criterion', 'Interpreter','latex')
ylabel('Counts', 'Interpreter','latex')
legend('','GKS', 'RPCholesky','RPCholesky without GKS','Pivoted Cholesky','Pivoted Cholesky without GKS',Location='northwest')
hold off

%% plot figs side by side
openfig('err_k_increasing_1d.fig')
ax1=gca;
figure;
tcl=tiledlayout(1,2)
ax1.Parent=tcl;
ax1.Layout.Tile=1;


f2 = openfig('dopt_k_increasing_1d.fig');
ax2=f2.Children(2);
% ax2 = gca;
ax2.Parent=tcl;
ax2.Layout.Tile=2;

%% D optimlaity comparison plot 

i = 11;

figure()
histogram(ld_rnd)
set(gca,'FontSize',48);
hold on

hxl2 = xline(ld_g(30), "--" , LineWidth=2);
fontsize(hxl2, 16, 'points')
hxl2.Color = [ 0    0.4470    0.7410];
hxl = xline(ld_gks(i), '-', LineWidth=2);
fontsize(hxl, 16, 'points')
hxl.Color = [0.8500    0.3250    0.0980];
hxl3 = xline(ld_nys(i), ":" , LineWidth=2);
hxl3.Color = [ 0.9290    0.6940    0.1250];
fontsize(hxl3, 16, 'points')
hxl4 = xline(ld_rp(i), "-.", LineWidth=2);
hxl4.Color = [0.4940    0.1840    0.5560];
hxl4 = xline(ld_piv(i), "--o", LineWidth=2);
hxl4.Color = [0.4660    0.6740    0.1880];
fontsize(hxl4, 16, 'points')
xlabel('D-Optimality Criterion', 'Interpreter','latex')
ylabel('Counts', 'Interpreter','latex')
legend('','Greedy', 'GKS','NysGKS','RPCholesky','Pivoted Cholesky',Location='northwest')
hold off
