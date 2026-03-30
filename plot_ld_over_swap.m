% plot swaps until no improvement by a factor f
clear; clc;

load droplet_setup.mat
addpath '/Users/jchen228/Desktop/JUQ_supplement'

n = N;
k = 30;

%% test gks initial placement
tic
pk_gks = gks(K,k); % output: Kernel of placed sensors, placed sensors
x_gks = x(pk_gks,:);
K_gks = K(pk_gks,pk_gks);
time_gks = toc;

% compute log determinant
ld_gks = slogdet(K_gks,sig_n)

%% test against random placement
pk_rnd = datasample(1:N,30,'Replace',false); % generate k random samples of values up to N without replacement
x_rnd = x(pk_rnd,:);
K_rnd = K_fun(x_rnd);
ld_rnd = slogdet(K_rnd,sig_n)

%% test against greedy placement

[pk_g, ld_g] = greedydopt6(K_fun_offdiag,x,k,sig_n);
ld_g(end)

%%  count swap steps
[sweep_gks,ld_all_gks] = count_sweeps(pk_gks,ld_gks,x,sig_n,sig_f,ls,f,@exchange_sensors2)
[sweep_g,ld_all_g] = count_sweeps(pk_g,ld_g(end),x,sig_n,sig_f,ls,f,@exchange_sensors2);
[sweep_rnd,ld_all_rnd] = count_sweeps(pk_rnd,ld_rnd,x,sig_n,sig_f,ls,f,@exchange_sensors2);

%% Plot results
lw = 1.5; %linewidth
figure(2)
plot(1:sweep_gks, ld_all_gks, '-o',LineWidth=lw)
hold on
plot(1:sweep_g, ld_all_g, '-o',LineWidth=lw)
plot(1:sweep_rnd, ld_all_rnd, '-o',LineWidth=lw)
legend('gks','greedy','rnd')