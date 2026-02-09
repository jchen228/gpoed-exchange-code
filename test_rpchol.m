% test RPCholesky

clear; clc;

rng("default")

% load("droplet_setup.mat")
% y = h_sel';
% sig_n = eta;
% p = k + 10;


addpath(genpath('/Users/jchen228/Desktop/gitGPOED/H2Pack-Matlab-master'))

global sig_n


a = 0; b = 1; % range of data
N = 100; % size of grid / number of candidate sensors
D = 1; % dimension of data
k = ceil(0.1 * N); % number of available sensors
p =  ceil(0.2 * N); % k+20; % ceil(0.4 * N);


% % design matrix (need N x D)
% x = a + ((b-a))*rand(N,D);
x = linspace(a,b,N)';

sig_f = 5;
ls = 1;

f = @(x) exp(-x).*cos(2*pi*x); 
% f = @(x) sin(2*pi*x) + 2*cos(3*pi*x) - sin(6*pi*x);

y = f(x);

sig_n = 0.02*norm(y)/sqrt(N);

y = y + sig_n*randn(size(y));

%% Kernel formation

K_fun = @(x) gaussKern(x,sig_f,ls);
K_fun_offdiag = @(x,x2) gaussKern(x,sig_f,ls,x2);

K_true = K_fun(x);

K_funMat = funMat(K_fun,K_fun,N); % K is a function handle, use K_funMat*x to form kernel 
K = K_funMat*x;

%% RP Cholesky (no gks yet)
tol = k;

[S, F] = rpchol_tol(x,sig_f,K_fun_offdiag,tol); % breaks at high resolutions, dependent on diagonal of K, sig_f


%% RP with GKS

% pk_rp = rpgks(F,k);
