% droplet plot over time
% plot reconstruction over time 
% frames(length(t)) = struct('cdata',[],'colormap',[]); % init struct to store frames of a movie
fk_err_post_save = [];
fk_post_save = [];
i = 1;

snap = [1 6001];
for step = 1:length(t)
    x_post = out_pde(snap(1):snap(2), 1);
    h_post = out_pde(snap(1):snap(2), 2);
    
    % fk_post = krr(x_g, h_post(pk_g(1:30)), K, pk_g(1:30), sig_n);
    % fk_post = krr(x_nys, h_post(pk_nys), K, pk_nys, sig_n);
    % fk_post = krr(x_rp, h_post(pk_rp), K, pk_rp, sig_n);
    % fk_post = krr(x_gks, h_post(pk_gks), K, pk_gks, sig_n);
    % fk_post = krr(x_piv, h_post(pk_piv), K, pk_piv, sig_n);
    
    fk_err_post = norm(fk_post - h_post)/norm(h_post);
    fk_err_post_save = [fk_err_post_save fk_err_post];

    fk_post_save = [fk_post_save fk_post];
    
    if i == length(t)
        break
    else
        i = i+1;
        snap = snap + 6001; % rewrite/update x range of snapshot for next time step
    end
end

% % plot error in time
% figure()
% scatter(t, fk_err_post_save, 48 ,'filled')
% set(gca,'FontSize',20)
% xlabel('time $t$','Interpreter','latex')
% ylabel('Relative Error','Interpreter','latex')

err_stats = datastats(fk_err_post_save')