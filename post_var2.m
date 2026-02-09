function v_post = post_var2(test_ind, X_input, sig_n, K)
%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % MATLAB function that predicts standard deviation for GP regression
% verion that does not require full kernel to be formed
% 
% % Inputs: 
% % test_ind -> indices of selected sensors 
% % K -> covariance function for all X (grid), 
% % sig_n -> noise variance (sigma_n^2)
% 
% % Outputs:
% % v_post -> posterior variance,
%
%%%%%%%%%%%%%%%%%%%%%%%%%%

sel = sort(test_ind);
unsel = setdiff(1:length(X_input), sel);

% K_train = kxx(X_input(sel,:),sigma,sig_n,sig_f); % prior covariance of selected indices
% k_s = kxx_offdiag(test_ind,X_input,sigma,sig_n,sig_f);
% K_test_diag = kxx_diag(unsel,sig_n,sig_f);

K_train = K(X_input(sel,:),X_input(sel,:));
k_s = K(X_input(sel,:),X_input(unsel,:))';
K_test_diag = K(X_input(unsel,:),'diag');

% do a cholesky decomp
L = chol(K_train + (sig_n^2)*eye(length(K_train)),'lower'); % K + sig_n*eye = L*L'
beta = L\(k_s');

v_post = zeros(length(X_input),1);
% v_post(sel) = 0;
% size(beta)
% v_post(unsel) = diag(K_test_diag - eta_2*beta'*beta)'; % diag of covariance matrix gives variance
v_post(sel) = (sig_n^2);
v_post(unsel) = diag(K_test_diag - beta'*beta)';
% size(diag(K_test_diag - beta'*beta)')
% size(K_test_diag )
% size(beta'*beta)
end

