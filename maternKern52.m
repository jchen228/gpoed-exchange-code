function KXX = maternKern52(X_input,sig_f,ls,varargin)
% Matern 5/2 kernel, ARD-capable: same calling convention as gaussKern.m
% (drop-in replacement), but ls may be a scalar (isotropic) or a 1xD
% vector (ARD -- one lengthscale per input dimension), matching the
% length-scale convention MATLAB's own fitrgp uses for 'ardmatern52'.

% INPUT:
% X_input is N x D data matrix
% ls: lengthscale, scalar or 1xD vector
% optional input: X_input2 if off-diagonal block is needed, or a
% char/string flag to request the diagonal only

% OUTPUT:
% covariance matrix KXX which is n by n (or n by 1 for the diagonal)
N = max(size(X_input));

if ~isempty(varargin)
    if ischar(varargin{1}) == 1 | isstring(varargin{1}) == 1
        KXX = sig_f^2*ones(N,1); % r=0 on the diagonal regardless of kernel family
        return
    else
        X_input2 = varargin{1};
    end
else
    X_input2 = X_input;
end
ls = ls(:)'; % row vector so it broadcasts across columns whether scalar or 1xD
R = pdist2(X_input./ls, X_input2./ls);
KXX = sig_f^2 * (1 + sqrt(5)*R + (5/3)*R.^2) .* exp(-sqrt(5)*R);
end
