function y = sig_gen_2D(theta,n)
% 根据θ(α,β)产生估计值信号
w1 = theta(1);
sgm1 = theta(2);
A1 = theta(3); 
% y 是p*个估计信号累加而成的信号
% 此时p=2
y = exp((1i*w1-sgm1)*n)*A1; 