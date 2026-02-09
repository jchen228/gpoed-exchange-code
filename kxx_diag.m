function KXX_diag = kxx_diag(N, sig_f)
% squared exponential covariance function used here (can also extend to any
% other exponential kernel)
% sig_n = 0
% INPUT:
% number of candidate sensors N
% OUTPUT:
% diag of covariance matrix KXX which is n by 1

% global sig_n
    KXX_diag = sig_f^2*ones(N,1);
    % KXX_diag = (sig_f^2+sig_n^2)*ones(N,1);
end