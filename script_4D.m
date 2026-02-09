% sensor placement and d-optimality analysis
% test against RP Cholesky and Pivoted Cholesky methods

% clear;
% clc;
% rng(0)
% 
load zhou_setup_moved_hp.mat

% ld_rnd_stats
% fk_rnd_stats


%% Greedy Pivoted Cholesky Algorithm

    tic
    F_piv = pivotedcholesky(K_fun_offdiag, x, k, 'greedy');
    pk_piv = rpgks(F_piv,k);
    x_piv = x(pk_piv,:);
    K_piv = K_fun(x_piv);
    time_piv = toc;
    
    % compute log determinant
    ld_piv = slogdet(K_piv, sig_n)
    
    fk_piv = krr(x_piv, y(pk_piv), K_fun_offdiag, pk_piv, sig_n, x);
    fk_piv_err = norm(fk_piv - y)/norm(y);

%% RPCholesky Algorithm
    
    tol = k;

    tic
    % [~,F] = rpchol_tol(K_fun_offdiag,tol,x);
    F = pivotedcholesky(K_fun_offdiag, x, k, 'rpcholesky');
    pk_rp = rpgks(F,k);
    x_rp = x(pk_rp,:);
    K_rp = K_fun(x_rp);
    time_rp = toc;

    % compute log determinant
    ld_rp = slogdet(K_rp, sig_n)

    fk_rp = krr(x_rp, y(pk_rp), K_fun_offdiag, pk_rp, sig_n, x);
    fk_rp_err = norm(fk_rp - y)/norm(y);
    % fk_rp_err = norm(fk_rp - y)/sqrt(N);

%% Plot Dopt comparison
figure()
histogram(ld_rnd)
set(gca,'FontSize',36);
hold on

hxl3 = xline(ld_rp, "--", LineWidth=2);
hxl3.Color = [0.9290    0.6940    0.1250];
fontsize(hxl3, 16, 'points')

hxl4 = xline(ld_piv, "-.", LineWidth=2);
hxl4.Color = [0.4940    0.1840    0.5560];
fontsize(hxl4, 16, 'points')


xlabel('D-Optimality Criterion', 'Interpreter','latex')
ylabel('Counts', 'Interpreter','latex')
legend('',sprintf('RPCholesky = %.2f', ld_rp),sprintf('Pivoted Cholesky = %.2f',ld_piv),Location='north')
hold off
