% set up data for SST
clear all
clc
close all
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
ls = 1550;
sig_f = 8;
sig_n = 0.002*norm(sst_sel)/sqrt(N);

K_fun = @(x) gaussKern(x,sig_f,ls);
K_fun_offdiag = @(x,x2) gaussKern(x,sig_f,ls,x2);

K_funMat = funMat(K_fun,K_fun,length(sea_ind)); % K is a function handle, use K_funMat*x to form kernel 

%%
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

%% random sensor placements - for testing
y = sst_sel_reshape_sea';

p = k + 20;


rr = 1;

fk_rnd = zeros(rr,length(sea_points));
fk_rnd_err = zeros(1,rr);
pk_rnd = datasample(1:length(sea_points),k,'Replace',false); % generate k random samples of values up to N without replacement
x_rnd = sea_points(pk_rnd,:);
K_rnd = K_fun(x_rnd);
ld_rnd = slogdet(K_rnd,sig_n);
fk_rnd = krr(x_rnd, y(pk_rnd), K_fun_offdiag, pk_rnd, sig_n, sea_points);
fk_rnd_err = norm(fk_rnd-y)/norm(y)

fk_rnd_recon = zeros(size(sst_sel_reshape));
fk_rnd_recon(sea_ind) = fk_rnd;
fk_rnd_recon = reshape(fk_rnd_recon,size(sst_sel));


plot_sst_map(fk_rnd_recon, pk_rnd)

%% Random Placements
% 
% rr = 10000;
% 
% ld_rnd = zeros(1,rr);
% fk_rnd = zeros(rr,length(sea_points));
% fk_rnd_err = zeros(1,rr);
% for i = 1:rr
%     pk_rnd = datasample(1:length(sea_points),k,'Replace',false); % generate k random samples of values up to N without replacement
%     x_rnd = sea_points(pk_rnd,:);
%     K_rnd = K_fun(x_rnd);
%     ld_rnd(i) = slogdet(K_rnd,sig_n);
%     fk_rnd(i,:) = krr(x_rnd, y(pk_rnd), K_fun_offdiag, pk_rnd, sig_n, sea_points);
%     fk_rnd_err(i) = norm(fk_rnd(i,:)'-y)/norm(y);
%     i
% end

% save('sst_setup_geo.mat')