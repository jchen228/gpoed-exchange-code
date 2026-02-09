function [R, P] = pivChol(A, varargin)
% NOT USED
% pivoted cholesky, swapping for largest diagonal entry
% OUT: upper triangular matrix R st P'*A*P = R'*R

n = length(A);

if ~isempty(varargin)
    stopping = varargin{1};
else
    stopping = n-1;
end

R = zeros(size(A));
piv = 1:n;

for k = 1:stopping
    % find the pivot
    B = A(k:n,k:n);
    [maximum, ~] = max(diag(B));
    l = find(diag(A) == maximum)';

    % swap rows and columns
    % if piv(k) ~= l
        % piv([k l]) = piv([l k]);
        % A = A(:, piv); % swap rows and columns
        % R = R(:, piv); % swap columns
        % A = A(piv, :);
    % end 





    % cholesky decomposition
    R(k,k) = sqrt(abs(A(k,k)));
    R(k, k+1:n) = R(k,k)^-1 * A(k, k+1:n);

    % Update A
    A(k+1:n, k+1:n) = abs(A(k+1:n, k+1:n) - R(k, k+1:n)'*R(k, k+1:n));
end
R(n,n) = sqrt(abs(A(n,n)));

P = piv;
    

end