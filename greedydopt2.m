function [sel_ind, ld_g] = greedydopt2(K,k,sig_n)
% % possibly use sub2ind


% Kernel function is fed in directly 
% K: kernel matrix
% k: number of sensors to select
% 

    n = length(K);

    sel_ind = zeros(1,k);
    w = zeros(1, n);
    % y = double(~w);
    ld = zeros(1,k);
    ld_g = zeros(1,k);

    % combos = nchoosek(grid,k);    

    for j = 1:k
        delta = zeros(1,n);
        for i = find(w == 0)
            % block splitting part
            w_temp = w; % get current sensor placements
            w_temp(i) = 1; % place a potential sensor
            w_test = find(w_temp>0); % find indices of all placed sensors

            delta(i) = slogdet(K(w_test,w_test),sig_n); % compute logdet
            
        end
        [ld(j), sel_ind(j)] = max(delta); 
        w(sel_ind(j)) = 1;

        ld_g(j) = slogdet(K(sel_ind(1:j),sel_ind(1:j)),sig_n);
    end

end