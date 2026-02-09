function [U_k, Lambda_k] = Nys(Kf, k, p)
% sepcific to standard normal test matrix

    rng("default")

    if isa(Kf,'funMat')
        N = Kf.Asize(1);
    else 
        N = length(Kf);
    end

    Omega = (1/sqrt(p))*randn(N,p);
   
    
    Y = Kf*Omega;
    nu = sqrt(N)*1e-6; % works for square exponential kernel
    Y_nu = Y + nu*Omega;

    C = chol(Omega'*Y_nu);
    B = Y_nu/C;

    [U, S, ~] = svd(B,0);
    Lambda = S.^2 - nu*eye(length(S));
    Lambda = max(Lambda, 0);

    U_k = U(:,1:k);
    Lambda_k = Lambda(1:k,1:k);

end