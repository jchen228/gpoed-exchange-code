function M = gaussian(coord, e, sig_f)
%
%   Exponential Kernel matrix, K(x,y) = exp(-e * ||x-y||^2)
%   e is 1/length scale
%
if isnumeric(coord)
    N = size(coord, 1);
    if N == 1
        M = sig_f;
        return ;
    end 
    dist = pdist(coord);
    M = squareform(exp(- e * dist.^2)) + eye(N);
elseif iscell(coord)&&length(coord)==2
    % size(coord{1})
    % size(coord{2})
    M = pdist2(coord{1}, coord{2});
    M = sig_f*exp(- e * M.^2);
end
end