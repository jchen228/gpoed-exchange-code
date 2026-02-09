function pk_rp = rpgks(F,k)

[u,~,~] = svd(F,0);
u_k = u(:,1:k); % take first k columns of V
                % gets left singular vectors of most significant singular
                % values
[~,~,p] = qr(u_k','vector'); % p permutation matrix such that u_k' * diag(p) = q * r
pk_rp = p(1:k); % column subset selection


end