function [sel_ind, ld_g] = greedydopt5(K,k,sig_n)
% version 5, efficient greedy
% update: L size is preallocated

% Kernel function is fed in directly 
% K: kernel matrix
% k: number of sensors to select

    n = length(K);
    alpha = @(i) 1/(sig_n^2)*K(i,i)+1;
    A = @(i,j) 1/(sig_n^2)*K(i,j); 

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
            w(sel_ind(j)) = 1;
            ld_g(j) = ld(j);
        else
        
        a = A(sel_ind(1:j-1),sel_ind(j));
        w(sel_ind(j)) = 1;
        b = L(1:j-1,1:j-1)\a;
        beta = sqrt(alpha(i) - b'*b);
        L(j,1:j) = [b' beta];
        ld_g(j) = sum(log(diag(L(1:j,1:j))));
        end

    end

end