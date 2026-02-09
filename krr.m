function fk_opt = krr(xs, y, Kf, pk, sig_n, varargin)
% GP Regression step
% % INPUTS:
% xs : selected sensors (k x D)
% y : observed data 
    if isa(Kf,'funMat')
        N = Kf.Asize(1);
        L = chol(rowcolkern(Kf,pk,pk) + sig_n^2*eye(length(xs)),'lower'); % K + sig_n*eye = L*L'
        alpha = L'\(L\y);
        fk_opt = (rowcolkern(Kf,pk,1:N)' * alpha); % fit / mean reconstruction
    elseif isa(Kf, 'function_handle')
        % disp('function handle')
        x = varargin{1};
        L = chol(Kf(xs,xs) + sig_n^2*eye(length(xs)),'lower');
        alpha = L'\(L\y);
        fk_opt = (Kf(xs,x)' * alpha); % fit / mean reconstruction
    else 
        N = length(Kf);
        L = chol(Kf(pk,pk) + sig_n^2*eye(length(xs)),'lower');
        alpha = L'\(L\y);
        fk_opt = (Kf(pk,1:N)' * alpha); % fit / mean reconstruction
    end
end
