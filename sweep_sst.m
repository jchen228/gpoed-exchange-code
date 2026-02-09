% sweep sst

clear; clc;
rng("default")

[ Lat, Lon, time, mask, sst ] = read_data_enso('sst.wkmean.1990-present.nc','lsmask.nc');

t0 = datetime(1800,1,1,0,0,0) + days(time(1));
tfin = datetime(1800,1,1,0,0,0) + days(time(end));
t = 1:length(sst);

temp_bounds = [-1.8 36.16];

selection_time = length(t)-260; % 5 years before last time snapshot
% selection_time = 1200+10;
tsel = datetime(1800,1,1,0,0,0) + days(time(selection_time));

total_post_steps = length((selection_time + 1):length(t));

sst_sel = sst(:,:,selection_time); 
[N1,N2] = size(sst_sel);
N = N1*N2;
x = [reshape(Lat,1,[])' reshape(Lon,1,[])'];
x = double(x);
X_all_ind = (1:N);
D = size(x,1);
sea_ind = find(mask); % get indices of sea points

sst_sel_reshape = reshape(sst_sel,1,[]);
sst_sel_reshape_sea = sst_sel_reshape(sea_ind);
sea_points = x(sea_ind,:); % get positions of sea points

k = 250;

% hyperparameters 
% ls = 15;
% sigma = ls;
% sig_n = 1e-8;
% sig_n = 0;
% sig_f = 3.5;
sig_n = 0.002*norm(sst_sel)/sqrt(N);
% eta_n2 = eta^2;
% lambda = eta_n2;



%% hyperparameter sweep
ls_sweep = 1500:50:2000; % 1000:50:2000, opt: 1550
sig_f_sweep = 7:10; %1:10, opt: 8



y_h = sst_sel_reshape_sea;
n_h = length(y_h);

m = 1;
n = 1;
nugget = 1e-6;
Ainv = (sig_n^-2)*eye(n_h);
% sweep search
for i = sig_f_sweep
    % K_ls = @(ls) gaussKern(sea_points,i,ls); % update sig_f in K_ls
  
    for j = ls_sweep
        K_fun_offdiag = @(x,x2) gaussKern(x,i,j,x2);
        [F,~] = pivotedcholesky(K_fun_offdiag, sea_points, k, 'rpcholesky');
        L =  chol(F'*F + (sig_n^2)*eye(k),'lower'); % K + sig_n*eye = L*L'
        Binv = L\((sig_n^-2)*F');
        alpha = Ainv - Binv'*Binv;
        fun = -(1/2)*y_h*alpha*y_h' - ((n_h - k)/2)*log(sig_n^2)+sum(log(diag(L))) - (n_h/2)*log(2*pi);
        obj(m,n) = fun;
        n = n+1;
        i
        j
    end
    n = 1; % reset n
    m = m+1;
end


maximum = max(obj,[],"all");
[m_opt,n_opt]=find(obj==maximum)

minimum = min(min(obj));
[m_test,n_test]=find(obj==minimum)

sig_f = sig_f_sweep(m_opt)
ls = ls_sweep(n_opt)

% sig_f = sig_f_sweep(m_test)
% ls = ls_sweep(n_test)
%% 


K_fun = @(x) gaussKern(x,sig_f,ls);
K_fun_offdiag = @(x,x2) gaussKern(x,sig_f,ls,x2);

K_funMat = funMat(K_fun,K_fun,length(sea_ind)); % K is a function handle, use K_funMat*x to form kernel 


% dim = min(size(x));
% kernel = @(coord) oedkern(coord, ls, sig_f);
% minSize = 200;
% coord = sea_points;
% htree = hierarchical_partition(coord, minSize, dim);
% alpha =  1;
% reltol = 1e-8;
% Yp = H2__ProxyPoint_QR_nlayer(kernel, htree, alpha, reltol);
% JIT_flag = true; 
% h2mat = Mat2H2_ID_Proxy(kernel, htree, Yp, 'reltol', reltol, alpha, JIT_flag);
% 
% Kf = funMat( @(x) H2_matvec(h2mat, htree, x),  @(x) H2_matvec(h2mat, htree, x), length(sea_ind));
% 


% save('sst_setup.mat')