function pk = gks(KXX,k)
% INPUT: 
%  KXX covariance matrix
% OUTPUT:
%  KXX_s: covariance matrix of selected sensors
%  pk: selected sensors via css
    % [~,~,v] = svd(KXX,0);
    % v_k = v(:,1:k); % take first k columns of V
    %                 % gets left singular vectors of most significant singular
    %                 % values
    [~,~,v_k] = svds(KXX, k);
    [~,~,p] = qr(v_k','vector'); % p permutation matrix such that v_k' * diag(p) = q * r
    pk = p(1:k); % column subset selection
end