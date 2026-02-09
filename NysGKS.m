function pk_nys = NysGKS(u,k)
    u_k = u(:,1:k); % take first k columns of V
                    % gets left singular vectors of most significant singular
                    % values
    [~,~,piv] = qr(u_k','vector'); % p permutation vecor such that u_k' * diag(p) = q * r
    pk_nys = piv(1:k);
end