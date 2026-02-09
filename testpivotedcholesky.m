clear all
clc 

rng(0)

x = rand(7,2);
kernel = @(r) exp(-r.^2/0.25);
[X,Y] = meshgrid(x,x);
K = kernel(pdist2(x,x));


% klst = 10:5:60;
% 
% err = zeros(length(klst),2);
% for j = 1:length(err)
%     F = pivotedcholesky(kernel, x, klst(j), 'greedy');
%     err(j,1) = norm(K - F*F')/norm(K);
%     F = pivotedcholesky(kernel, x, klst(j), 'rpcholesky');
%     err(j,2) = norm(K - F*F')/norm(K);
% end

klst = length(x);
% err = zeros(length(klst),2);

[F,S] = pivotedcholesky(kernel, x, klst, 'greedy');
err = norm(K - F*F')/norm(K);
% F = pivotedcholesky(kernel, x, klst, 'rpcholesky');
% err = norm(K - F*F')/norm(K);



%% check conditions on R
n = length(F);
I = eye(n);
P = I(:,S);
R = F'*P;
if all(round(vecnorm(R)) == ones(size(vecnorm(R))))
    for k = 1:n-1
        conditions = R(k:n,(k+1):n);
        if any(conditions <= R(k,k) ,'all') 
            % disp(['conditions on ', num2str(k), ' diagonal element satisfied'])
            continue
        else
            disp('conditions on R not satisfied')
        end
    end
    disp('Conditions on R are satisfied')
else 
    disp('col norm not equal to one')
end

PK = P'*K*P;