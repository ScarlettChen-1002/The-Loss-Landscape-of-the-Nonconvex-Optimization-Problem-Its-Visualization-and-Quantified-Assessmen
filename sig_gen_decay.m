function y = sig_gen_decay(theta, n)
% 根据θ(α)产生纯衰减无振荡的预测信号
% 此时p*=1
sgm = theta(1);
A = theta(2);
y = exp((-sgm)*n)*A; 