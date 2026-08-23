function [Znew, Lnew, Snew, candNew] = incremental_triangular_solve(A, alpha_q, L, Z, S, cand, qPos)

m = size(L,1);
k = m + 1;

q = cand(qPos);
b = Z(:,qPos);
beta = sqrt(max(alpha_q - b'*b, 0));

Ltilde = [L, zeros(m,1); b', beta];
ell21 = Ltilde(2:end,1);
L22   = Ltilde(2:end,2:end);

zetaRow = (A(q,cand) - b'*Z) / beta;
Ztilde = [Z; zetaRow];

v = L22 \ ell21;
Y = Ztilde(2:end,:) + v .* Ztilde(1,:);

kk = k - 1;
tprev = 1;
Cdiag = zeros(kk,1);
coeff = zeros(kk,1);
Cexp  = zeros(kk,kk);
for j = 1:kk
    tj = tprev + v(j)^2;
    Cdiag(j) = sqrt(tj/tprev);
    coeff(j) = v(j) / sqrt(tj*tprev);
    Cexp(j,j) = Cdiag(j);
    Cexp((j+1):kk, j) = v((j+1):kk) * v(j) / sqrt(tj*tprev);
    tprev = tj;
end

% O(k) forward substitution against C, not a full O(k^2) triangular solve
ncand = size(Y,2);
X = zeros(kk, ncand);
s = zeros(1, ncand);
for j = 1:kk
    X(j,:) = (Y(j,:) - v(j)*s) / Cdiag(j);
    s = s + coeff(j) * X(j,:);
end

Lnew = L22 * Cexp;

Snew = [S(2:end); q];
p2 = S(1);
keepMask = true(numel(cand), 1);
keepMask(qPos) = false;
newcol = Lnew \ A(Snew, p2);

Znew = [X(:,keepMask), newcol];
candNew = [cand(keepMask); p2];

end
