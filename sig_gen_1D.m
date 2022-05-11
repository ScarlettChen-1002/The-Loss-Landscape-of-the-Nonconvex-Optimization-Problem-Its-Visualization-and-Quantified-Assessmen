function y = sig_gen_1D(theta,n)
% 根据θ(α)产生估计值信号
w = theta(1);
a = theta(2);
y = exp((1i*w-a)*n); 