function [sel_ind, ld_g] = greedydopt7(K_fun, x, k, sig_n)
% version 7, efficient greedy
% update: this is the fast greedy MAP implementation

arguments (Input)
    K_fun   % Input the same kernel function used in greedydopt6
    x
    k
    sig_n
end

arguments (Output)
    sel_ind
    ld_g
end

n = size(x, 1);
Z = 1:n; 
c = zeros(k, n); % Pre-allocate Cholesky vectors

% 1. Initialize d2 identically to alpha in greedydopt6
d2 = zeros(1, n);
for i = 1:n
    d2(i) = (1 / sig_n^2) * K_fun(x(i,:), x(i,:)) + 1;
end

% Initialize outputs
sel_ind = zeros(1, k);
ld_g = zeros(1, k);

% First selection
[~, j] = max(log(d2));
sel_ind(1) = j;
ld_g(1) = 0.5 * log(d2(j)); % Equivalent to log(sqrt(alpha))

counter = 1;

while counter < k
    j = sel_ind(counter); % The most recently selected sensor
    
    for i = setdiff(Z, sel_ind(1:counter))
        
        % 2. Use spatial coordinates x and scale properly
        cov_ji = (1 / sig_n^2) * K_fun(x(j,:), x(i,:));
        
        % Compute dot product of previously computed Cholesky elements
        if counter == 1
            dot_prod = 0;
        else
            dot_prod = c(1:counter-1, j)' * c(1:counter-1, i);
        end
        
        % Calculate e_{ji}
        e_i = (cov_ji - dot_prod) / sqrt(d2(j));
        
        % 3. Correct assignment (no shifting)
        c(counter, i) = e_i;
        
        % Update conditional variance
        d2(i) = d2(i) - e_i^2;
    end
    
    % Decision criteria: Mask out already selected sensors with -Inf
    valid_d2 = d2;
    valid_d2(sel_ind(1:counter)) = -Inf; 

    % valid_d2 = round(valid_d2, 10);
    
    % 4. Select the next best sensor
    [~, next_j] = max(log(valid_d2)); % Maximizing d2 is equivalent to max log(d2)
    
    counter = counter + 1;
    sel_ind(counter) = next_j;
    
    % 5. Accumulate log-determinant matching greedydopt6
    ld_g(counter) = ld_g(counter-1) + 0.5 * log(d2(next_j));
end

end