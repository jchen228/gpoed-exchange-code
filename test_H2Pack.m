% test H2Pack for accuracy
clear;
rng("default")
addpath(genpath('/Users/jchen228/Desktop/gitGPOED/H2Pack-Matlab-master'))

global sig_n

a = 0; b = 1; % range of data
N = 5000; % size of grid / number of candidate sensors
D = 1; % dimension of data
k = ceil(0.1 * N); % number of available sensors
p =  ceil(0.2 * N); % k+20; % ceil(0.4 * N);


% % design matrix (need N x D)
x = a + ((b-a))*rand(N,D);
% x = linspace(a,b,N)';

sig_f = 1;
ls = 1;

%% Kernel formation
% 
% K_fun = @(x) gaussKern(x,sig_f,ls);
% K_fun_offdiag = @(x,x2) gaussKern(x,sig_f,ls,x2);
% 
% % z = rand(size(x));
% 
% K_true = K_fun(x);
% % K_true_z = K_true * z;
% 
% K_funMat = funMat(K_fun,K_fun,N); % K is a function handle, use K_funMat*x to form kernel 
% K = K_funMat*x;
%%


Omega = (1/sqrt(p))*randn(size(x,1),p);

dim = min(size(x));
kernel = @(coord) oedkern(coord, ls, sig_f);
minSize = 200;
coord = x;
htree = hierarchical_partition(coord, minSize, dim);
alpha =  1;
reltol = 1e-8;
Yp = H2__ProxyPoint_QR_nlayer(kernel, htree, alpha, reltol);
JIT_flag = true; 
h2mat = Mat2H2_ID_Proxy(kernel, htree, Yp, 'reltol', reltol, alpha, JIT_flag);
% Y = H2_matvec(h2mat, htree, Omega); % A*Omega step
%%

% K_fun = @(x) gaussKern(x,sig_f,ls);
% K_true = K_fun(x); 
% K_true_Omega = K_true * Omega;
% 
% % Test Matrix-vector product
% err_h2 = norm(K_true_Omega-Y,'fro')/norm(K_true_Omega,'fro');


%% Test funMAt
Kf = funMat( @(x) H2_matvec(h2mat, htree, x),  @(x) H2_matvec(h2mat, htree, x), [N N]  );
% Y = Kf*Omega;

% norm(K_true_Omega-Yf)/norm(K_true_Omega)

%%

[u_k,lambda_k] = Nys(Kf, k, p);

% err_NysH2_2 = norm(Kf - u_k*lambda_k*u_k')/norm(K);
% err_NysH2_fro = norm(Kf - u_k*lambda_k*u_k','fro')/norm(K,'fro')

