%% 以下为老师给出的代码
clear,clc
close all
rng(20); 
N = 1e3;
%%%%%%%%%%%%%%%%%%%%%%%%%%
ind = randperm(N);
ind = ind(1:p);
ind = sort(ind);
%%%%%%%%%%%%%%%%%%%%%%%%%%
omg = [0.1 0.2 0.65]*2*pi;
A = [1 0.5 1];
alpha = [0.01 0.02 0.005];
% alpha = [0.005 0.005 0.005];
%%
n = 0: N-1;
omg = omg(:);
alpha = alpha(:);
A = A(:).';
n = n(:).';
xn_ideal = A*exp((1i*omg-alpha)*n);

percent = 20;
p = ceil(percent/100*N);
[mask,ind] = SinPoisson(p, N);