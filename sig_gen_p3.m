function y = sig_gen_p3(theta,n)
% 根据θ(α,β)产生估计值信号
w1 = theta(1);
w2 = theta(2);
w3 = theta(3);
sgm1 = theta(4);
sgm2 = theta(5);
sgm3 = theta(6);
A1 = theta(7); 
A2 = theta(8);
A3 = theta(9);
% y 是p*个估计信号累加而成的信号
% 此时p*=2
y = exp((1i*w1-sgm1)*n)*A1 + exp((1i*w2-sgm2)*n)*A2 + exp((1i*w3-sgm3)*n)*A3; 