function [sel_ind, ld_g] = greedydopt4(K,k,sig_n)
% version 4, efficient greedy
% update: add rows to L 


% Kernel function is fed in directly 
% K: kernel matrix
% k: number of sensors to select
% 

    n = length(K);

    sel_ind = zeros(1,k);
    w = zeros(1, n);
    ld = zeros(1,k);
    ld_g = zeros(1,k);
    % L = zeros(k,k);
    L = [];
    alpha = @(i) 1/(sig_n^2)*K(i,i)+1; % assumes Each diagonal of K is the same
    A = 1/(sig_n^2)*K;

    % combos = nchoosek(grid,k);    

    for j = 1:k % for each available sensor

        delta = zeros(1,n); % reset scoring vector

        for i = find(w == 0) % perform sweep

            % block splitting part
            w_temp = w; % get current sensor placements
            w_temp(i) = 1; % place a potential sensor
            w_test = find(w_temp>0); % find indices of all placed sensors

            % record increasing lower triangular matrix L 
            if j == 1
                % ld_temp = log(sqrt(A(1,1)))
                ld_temp = log(sqrt(1/(sig_n^2)*K(1,1)+eye(size(1))));
            else
                a = A(sel_ind(1:j-1),i);
                % a = 1/(sig_n^2)*K(sel_ind(1:j-1),i);
                b = L\a;
                % beta = sqrt(alpha - b'*b);
                beta = sqrt(alpha(i) - b'*b);
                ld_temp = sum(log(diag(L))) + log(beta);
            end

            % delta(i) = slogdet(K(w_test,w_test),sig_n); % compute logdet
            delta(i) = ld_temp;
            
        end
        delta;
        [ld(j), sel_ind(j)] = max(delta);


        % record increasing lower triangular matrix L 
        if j == 1
            L = sqrt(1/(sig_n^2)*K(1,1)+eye(size(1)));
            w(sel_ind(j)) = 1;
            ld_g(j) = ld(j);
        else
        
        a = A(sel_ind(1:j-1),sel_ind(j));
        w(sel_ind(j)) = 1;
        b = L\a;
        beta = sqrt(alpha(i) - b'*b);
        L = [L zeros(j-1,1);b' beta]
        ld_g(j) = sum(log(diag(L)));
        end
    end

end