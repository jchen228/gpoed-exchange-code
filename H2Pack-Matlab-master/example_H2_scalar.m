clear;

% load out_pde.dat

% N = 6001;
dim = 1;

% coord = out_pde(1:N,1);
% h = out_pde(:,2);
% t = out_pde(:,3);

irange = 5;
for i = irange
    clearvars -except irange i timer timer_control max_err min_err mean_err dim

    %%  Kernel Function/Matrix Definition
    % kernel = @(coord)multiquadric(coord, 1/2, 100);
    % kernel = @(coord)reciprocal(coord, 1);
    % kernel = @(coord)logfun(coord);
    % kernel = @(coord)exponential(coord, 0.1);
    sig_f = 1; ls = 1;
    % kernel = @(coord)gaussian(coord, 0.1, 1);
    kernel = @(coord)oedkern(coord, ls, sig_f);
    % % e = 0.1
    
    %%  Point generation
    npts = i*500;
    K = ceil(.1*npts);
    % coord = npts^(1/3)*rand(npts, dim);
    coord = linspace(0,1,npts)';

    % coord = bsxfun(@rdivide, coord, sum(coord.^2, 2).^(1/2));
    
    %%  STEP1: Hierarchical partitioning 
    minSize = 200;
    htree = hierarchical_partition(coord, minSize, dim);
    %clear coord; % htree has a reorder copy of coord
    
    %%  STEP2: Proxy Point Selection
    alpha =  1;
    reltol = 1e-5;
    %   Scheme 1: proxy surface method, only works for kernel from
    %   potential theory;
    % Yp = H2__ProxyPoint_Surface(htree, dim, alpha, 600);
      % Scheme 2: numerical selection, works for general kernels, using a
      % heuristic multi-step qr approach. 
    Yp = H2__ProxyPoint_QR_nlayer(kernel, htree, alpha, reltol);
    
    
    %%  STEP3: H2 matrix construction using proxy points
    JIT_flag = true; 
    h2mat = Mat2H2_ID_Proxy(kernel, htree, Yp, 'reltol', reltol, alpha, JIT_flag);
    
    %%  STEP4: H2 matrix-vector multiplications
    % K_mat = kernel(coord);

    
    x = randn(npts, 20);
    
    for j = 1:5 % to average over 5 runs
        
        tic
        u_h2 = H2_matvec(h2mat, htree, x);
        temp_timer(j) = toc;
    
        % tic
        % u_true = K_mat*x; 
        % temp_timer_control(j) = toc;

    end
    timer(i) = mean(temp_timer);
    % timer(i) = toc/5;

    % % time full matrix formation
    % tic
    % for j = 1:5 % to average over 5 runs
    %     % tic
    %     u_true = K_mat*x; 
    %     % temp_timer_control(j) = toc;
    % 
    % end
    % % timer_control(i) = mean(temp_timer_control);
    % timer_control(i) = toc/5;
    
    %   error checking
    idx = randperm(npts, K);
    u_exact = kernel({coord(idx, :), coord}) * x;
    err = sqrt(sum((u_h2(idx, :)-u_exact).^2, 1) ) ./ sqrt(sum(u_exact.^2, 1)); 
    % relative error in 2 norm for each randomly sampled column
    
    % fprintf("min/mean/max relative errors for 20 matvecs:\n%.3e,%.3e,%.3e\n", ...
        % min(err), mean(err), max(err));

    min_err(i) = min(err); 
    mean_err(i) = mean(err);
    max_err(i) = max(err);

    % % error relative to full formed matrix for small n only
    % score = sqrt(sum((u_h2-u_true).^2, 1) ) ./ sqrt(sum(u_true.^2, 1));
    % 
    % mean_score = mean(score);

end 

% plot
% figure()
% % plot(500*irange, min_err, 500*irange, mean_err, 500*irange, max_err, LineWidth=2)
% loglog(500*irange, min_err, 500*irange, mean_err, 500*irange, max_err, LineWidth=2)
% xlabel('size of matrix n')
% ylabel('relative error')
% legend('min', 'mean', 'max')
% title('matvec computation error')


% % plot of computation time - ffor small n only
% figure()
% loglog(500*irange, timer, 500*irange, timer_control, LineWidth=2)
% legend('h2pack', 'full')
% title('time (sec) vs size (n = 500, 1000, ..., 66,500)')

% plot of computation time - fit data
% poly = polyfit(log(500*irange(2:end)),log(timer(2:end)),1);
% 
% figure()
% loglog(500*irange, timer, LineWidth=2)
% hold on
% loglog(500*irange, (500*irange).^poly(1)*exp(poly(2)), LineWidth=2)
% xlabel('size of matrix n')
% ylabel('time (s)')
% legend('h2pack','fit')
% hold off
% % title('time (sec) vs size (n = 500, 1000, ..., 66,000)')
% 
% save('h2compare_dim_3.mat','irange','min_err','mean_err','max_err', 'poly', 'timer')