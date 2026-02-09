% test Randomized Nystrom approximation
clear;
rng("default")

a = 0; b = 1; % range of data
N = 1000; % size of grid / number of candidate sensors
D = 1; % dimension of data
k = ceil(0.1 * N); % number of available sensors
p =  ceil(0.2 * N); % k+20; % ceil(0.4 * N);

% % design matrix (need N x D)
x = a + ((b-a))*rand(N,D);
% x = linspace(a,b,N)';

sig_f = 1;
ls = 1;

%% 

K = gaussKern(x, sig_f, ls);
% K = diag(diag(-1*ones(1,N-1),1)+ 1*ones(1,N)) + diag(-1*ones(1,N-1),-1);

[U, Lambda] = RandNys(K,k,p); % returns rank k approximation for now

K_hat = U*Lambda*U';

norm_2 = norm(K_hat - K)/norm(K);
norm_fro = norm(K_hat - K,'fro')/norm(K,'fro')
