function p = exchange_sensors_sequential_incremental(p, x, sig_n, sig_f, ls, f, kernel_fn)
% exchange logic executed inside one worker, using the incremental
% triangular solve instead of a fresh O(c*k^2) solve every iteration.
% kernel_fn: optional, defaults to @gaussKern (isotropic squared
% exponential). Pass @maternKern32/@maternKern52 for ARD Matern -- in
% that case ls should be a 1xD vector rather than a scalar.
if nargin < 7 || isempty(kernel_fn)
    kernel_fn = @gaussKern;
end

p = p(:)';
k = length(p);
n = size(x, 1);

K_fun = @(x) kernel_fn(x,sig_f,ls);
K_fun_offdiag = @(x,x2) kernel_fn(x,sig_f,ls,x2);
K_diag = @(x) kernel_fn(x,sig_f,ls,'1');

alpha = @(i) 1/(sig_n^2)*K_diag(1)+1;
A = @(i,j) 1/(sig_n^2)*K_fun_offdiag(x(i,:),x(j,:));
A_p = (1/(sig_n^2))*K_fun(x(p,:)) + eye(k);
L = chol(A_p, 'lower');

p_orig_front = p(1);
p_temp = p(2:end)';
test_ind = setdiff(1:n, p_temp)';

% Fast Cholesky Downdate (removes the first sensor from L), once
col_to_add = L(2:end,1);
L_fast = L(2:end,2:end);
L_fast = (cholupdate(L_fast', col_to_add, '+'))';

% One-time O(c k^2) solve; every subsequent iteration reuses/updates this
% incrementally instead of re-solving against all candidates from scratch
rhs_all = A(p_temp, test_ind);
Z = L_fast \ rhs_all;

S = p_temp;
cand = test_ind;

converged = false;
no_swap_counter = 0;

while ~converged
    alpha_val = alpha(cand);
    beta_all = sqrt(max(alpha_val - sum(Z.^2,1), 0));

    orig_ind = find(cand == p_orig_front, 1);
    beta_orig = beta_all(orig_ind);

    [global_max_beta, best_idx] = max(beta_all);

    if global_max_beta > beta_orig * sqrt(f) + 1e-9
        qPos = best_idx;
        no_swap_counter = 0;
    else
        qPos = orig_ind;
        no_swap_counter = no_swap_counter + 1;
    end

    q = cand(qPos);
    alpha_q = alpha(q);
    [Z, L_fast, S, cand] = incremental_triangular_solve(A, alpha_q, L_fast, Z, S, cand, qPos);

    % the sensor that just re-entered the candidate pool is always
    % appended last, and is exactly next iteration's p(1)
    p_orig_front = cand(end);

    if no_swap_counter >= k
        converged = true;
    end
end

p = [p_orig_front, S'];
end
