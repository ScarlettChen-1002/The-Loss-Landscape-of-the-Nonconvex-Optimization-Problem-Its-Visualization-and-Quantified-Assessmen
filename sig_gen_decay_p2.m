function y = sig_gen_decay_p2(theta, n)
% 根据θ(α)产生纯衰减无振荡的预测信号
% 此时p*=2
sgm1 = theta(1);
sgm2 = theta(2);
A1 = theta(3);
A2 = theta(4);
y = exp((-sgm1)*n)*A1 + exp((-sgm2)*n)*A2; 