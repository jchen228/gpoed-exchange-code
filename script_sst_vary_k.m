% script for just NysGKS as k varies

clear; clc;
rng("default")

% global sig_n
addpath(genpath('/Users/jchen228/Desktop/gitGPOED/H2Pack-Matlab-master'))

load sst_setup.mat
y = sst_sel_reshape_sea';



% for each k, get relative error over time data stats
%% NysGKS Algorithm

j = 1;
for k = 150:20:350
    tol = k;
    p = k + 20;
    
    tic
    [u,~] = Nys(Kf, k, p);
    pk_nys = NysGKS(u,k);
    x_nys = sea_points(pk_nys,:);
    K_nys = rowcolkern(Kf,pk_nys,pk_nys);
    time_nys = toc;

    % tic
    % [~,F] = rpchol_tol(K_fun_offdiag,tol,sea_points);
    % pk_rp = rpgks(F,k);
    % x_rp = sea_points(pk_rp,:);
    % K_rp = K_fun(x_rp);
    % time_rp = toc;

    % tic
    % F_piv = pivotedcholesky(K_fun_offdiag, sea_points, k, 'greedy');
    % pk_piv = rpgks(F_piv,k);
    % x_piv = sea_points(pk_piv,:);
    % K_piv = K_fun(x_piv);
    % time_piv = toc;

    % greedy
    % K_g = K_fun(x_g);

    % compute log determinant
    ld_nys(j) = slogdet(K_nys,sig_n);
    % ld_rp(j) = slogdet(K_rp,sig_n);
    % ld_piv(j) = slogdet(K_piv,sig_n);
    % ld_g(j) = slogdet(K_g(1:k,1:k),sig_n);

    fk_err_post_save = [];
    i = 1;
    for step = selection_time:length(time)

    h_post = sst(:,:,step);
    h_post = reshape(h_post, N, []);
    h_post = h_post(sea_ind);
    
    fk_post = krr(x_nys, h_post(pk_nys), K_fun_offdiag, pk_nys, sig_n, sea_points);
    % fk_post = krr(x_rp, h_post(pk_rp), K_fun_offdiag, pk_rp, sig_n, sea_points);
    % fk_post = krr(x_piv, h_post(pk_piv), K_fun_offdiag, pk_piv, sig_n, sea_points);
    % fk_post = krr(x_g(1:k,:), h_post(pk_g(1:k)), K_fun_offdiag, pk_g(1:k), sig_n, sea_points);



    fk_err_post = norm(fk_post - h_post)/norm(h_post);
    fk_err_post_save = [fk_err_post_save fk_err_post];

    i = i+1;
    end
    err_stats(j) = datastats(fk_err_post_save');
j = j+1
end


k_mean_nys = [ err_stats(1).mean  err_stats(2).mean  err_stats(3).mean  err_stats(4).mean  err_stats(5).mean  err_stats(6).mean ...  
            err_stats(7).mean  err_stats(8).mean err_stats(9).mean  err_stats(10).mean err_stats(11).mean];

%%
% plot d opt as k increases in increments of 20

figure()
set(gca,'FontSize',36)
subplot(1,2,1);
plot(150:20:350, ld_nys,'-o', MarkerSize=10, LineWidth=2)
% plot(150:20:350, ld_g(150:20:350),'-o', MarkerSize=10, LineWidth=2)
% hold on
% plot(150:20:350, ld_rp,'-o', MarkerSize=10, LineWidth=2)
% plot(150:20:350, ld_piv,'-o', MarkerSize=10, LineWidth=2)
xlabel('k')
ylabel('D-Optimality')
axis auto
hold off

subplot(1,2,2); 
plot(150:20:350, k_mean_nys,'-o', MarkerSize=10, LineWidth=2)
% plot(150:20:350, k_mean_g,'-o', MarkerSize=10, LineWidth=2)
% hold on
% plot(150:20:350, k_mean_rp,'-o', MarkerSize=10, LineWidth=2)
% plot(150:20:350, k_mean_piv,'-o', MarkerSize=10, LineWidth=2)
xlabel('k')
ylabel('Relative Error')
% legend('NysGKS','RPCholesky','Pivoted Cholesky')
legend('Efficient Greedy','RPCholesky','Pivoted Cholesky')
axis auto
hold off

%%
