function [mask, v] = SinPoisson(p, N)
% 用于非均匀采样的函数
% N: total number of sampling points，样本总点数
% p: number of NUS points，采样后的数据点数
rng(0);
v = zeros(1,N);
d1 = N/p;
adj=2.0*(d1-1);
t = 0;
while 1
    i = 0;
    n = 0;
    while i<N
        v(n+1) = i;
        i = i+1;
        t = adj*(sin((i+0.5)/(N+1)*1.5707963268));
        rand_d = psn(t);
        i = i+rand_d;
        n=n+1;
    end
    if n==p; break; end
    if n>p
        adj = adj*1.02;
    else
        adj = adj/1.02;
    end
end

v = v(1:p);
v=v+1;
mask = zeros(1, N);
mask(v) = 1; % mask为掩码
end

function k = psn(t)
L = exp(-t);
P = 1;
k = 0;
while P>=L
    u = rand(1);
    P = P*u;
    k = k+1;
end
k=k-1;
end
