function KXX = sepKern(X_input,sig_f,ls, varargin) 
% separable square exponential kernel

% INPUT:
% X_input is N x D data matrix 
% length: lengthscale
% optional input: X_input 2 if off diagonal block is needed

% OUTPUT:
% covariance matrix KXX which is n by n
N = max(size(X_input));
% D = min(size(X_input));

if ~isempty(varargin) 
    if ischar(varargin{1}) == 1 | isstring(varargin{1}) == 1
        KXX = sig_f^2*ones(N,1); % returns diagonal elements as a column vector
        return
    else
        X_input2 = varargin{1};
    end
else
    X_input2 = X_input;
end


    rfcn = @(X_input) pdist2(X_input,X_input2,'seuclidean',ls);
    R = rfcn(X_input);
    KXX = sig_f^2*exp(-R.^2);
end