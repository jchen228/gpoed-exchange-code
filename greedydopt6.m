function [sel_ind, ld_g] = greedydopt6(K_fun,x,k,sig_n)
% version 6, efficient greedy
% update: K not fully formed

% Kernel function is fed in directly 
% K_fun: kernel matrix function
% k: number of sensors to select
% sig_n: signal variance

    n = length(x);
    alpha = @(i) 1/(sig_n^2)*K_fun(x(i,:),x(i,:))+1;
    A = @(i,j) 1/(sig_n^2)*K_fun(x(i,:),x(j,:)); 

    % Initialize variables
    sel_ind = zeros(1,k);
    w = zeros(1, n);
    ld = zeros(1,k);
    ld_g = zeros(1,k);
    L = zeros(k,k);



    for j = 1:k % for each available sensor

        delta = zeros(1,n); % reset scoring vector

        for i = find(w == 0) % perform sweep

            if j == 1
                ld_temp = log(sqrt(alpha(i)));
            else
                a = A(sel_ind(1:j-1),i);
                b = L(1:j-1,1:j-1)\a;
                beta = sqrt(alpha(i) - b'*b);
                ld_temp = ld_g(j-1) + log(beta);
            end

            delta(i) = ld_temp;
            
        end
        [ld(j), sel_ind(j)] = max(delta);

        % record increasing lower triangular matrix L 
        if j == 1
            L(j,j) = sqrt(alpha(sel_ind(j)));
        else
            a = A(sel_ind(1:j-1),sel_ind(j));
            b = L(1:j-1,1:j-1)\a;
            beta = sqrt(alpha(sel_ind(j)) - b'*b);
            L(j,1:j) = [b' beta];
        end
        w(sel_ind(j)) = 1;
        ld_g(j) = ld(j);
    end
end