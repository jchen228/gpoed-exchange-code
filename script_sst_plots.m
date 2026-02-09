% script for sst over time
load("sst_selections.mat")
rng("default")

%% plot D-opt comparison

figure()
histogram(ld_rnd)
set(gca,'FontSize',36);
hold on

hxl1 = xline(ld_g(end), "-" , LineWidth=3);
hxl1.Color = [0    0.4470    0.7410];
fontsize(hxl1, 16, 'points')

hxl2 = xline(ld_nys, "--" , LineWidth=3);
hxl2.Color = [0.8500 0.3250 0.0980];
fontsize(hxl2, 16, 'points')

hxl3 = xline(ld_rp, "-.", LineWidth=3);
hxl3.Color = [0.9290    0.6940    0.1250];
fontsize(hxl3, 16, 'points')

hxl4 = xline(ld_piv, "-.", LineWidth=3);
hxl4.Color = [0.4940    0.1840    0.5560];
fontsize(hxl4, 16, 'points')

hxl5 = xline(ld_q, ":", LineWidth=3);
hxl5.Color = [0.4660    0.6740    0.1880];
fontsize(hxl5, 16, 'points')

xlabel('D-Optimality Criterion', 'Interpreter','latex')
ylabel('Counts', 'Interpreter','latex')
legend('','Efficient Greedy','NysGKS','RPCholesky','Pivoted Cholesky','POD-DEIM',Location='northwest')
hold off

%% variance plot of NysGKS placements

% v_post = post_var2(pk_nys, sea_points, sig_n, K_fun_offdiag);
% var_fit = zeros(N,1); % vector of all variances 
% var_fit(sea_ind) = sqrt(v_post);
% var_fit = reshape(var_fit, size(sst_sel));
% 
% plot_sst_map(var_fit, pk_nys)

%% RE over time

fk_err_post_save = [];
i = 1;

% t = 1207 begin training time

for step = selection_time:length(time)

    h_post = sst_reshape(sea_ind,step);
    
    % % uncomment for each post selection error
    % fk_post = krr(x_nys, h_post(pk_nys), K_fun_offdiag, pk_nys, sig_n, sea_points);
    % fk_post = krr(x_rp, h_post(pk_rp), K_fun_offdiag, pk_rp, sig_n, sea_points);
    % fk_post = krr(x_piv, h_post(pk_piv), K_fun_offdiag, pk_piv, sig_n, sea_points);
    % fk_post = krr(x_rnd, h_post(pk_rnd), K_fun_offdiag, pk_rnd, sig_n, sea_points);
    % fk_post = krr(x_g, h_post(pk_g), K_fun_offdiag, pk_g, sig_n, sea_points);

    
    % for POD-DEIM 
    x_test = h_post-meansst;
    xls = Phi*Phi'*x_test;
    % xls = Phi(:,1:r)*(Phi(sensors,1:r)\x_test(sensors));
    qdeim_err = norm(xls - x_test)/norm(x_test+meansst); % remove mean for anomaly
    fk_err_post_save = [fk_err_post_save qdeim_err];
    size(fk_err_post_save)

    % fk_err_post = norm(fk_post - h_post)/norm(h_post);
    % fk_err_post_save = [fk_err_post_save fk_err_post];

    i = i+1;
end

fk_q_err_pod = fk_err_post_save;

%% plot error in time
figure()
scatter(t(selection_time:end), fk_g_err, 48 ,'filled')
hold on
% scatter(t(selection_time:end), fk_nys_err, 48 ,'filled')
scatter(t(selection_time:end), fk_rp_err, 48 ,'filled')
scatter(t(selection_time:end), fk_piv_err, 48 ,'filled')
scatter(t(selection_time:end), fk_q_err, 48 ,'filled')
scatter(t(selection_time:end), fk_q_err_pod, 48 ,'filled')


% xline(selection_time)

set(gca,'FontSize',36)
legend('Greedy','RPCholesky','Pivoted Cholesky','POD-DEIM','POD')
xlabel('time $t$','Interpreter','latex')
ylabel('Relative Error','Interpreter','latex')

err_stats = datastats(fk_err_post_save');