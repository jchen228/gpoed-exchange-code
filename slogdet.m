function ld = slogdet(A,sig_n)
% log determinant of sig_n^-2 * A + eye(size(A)), A is SPD
% minimization function phi_D(S)
    
    [L, flag] = chol(1/(sig_n^2) * A + eye(size(A)),'lower');
    ld = sum(log(diag(L)));
end