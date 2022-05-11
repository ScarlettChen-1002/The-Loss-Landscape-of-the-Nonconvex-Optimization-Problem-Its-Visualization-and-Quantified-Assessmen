function y = sig_gen_1D2(theta,n)
w1 = theta(1);
w2 = theta(2);
ww = [w1,w2];
s1 = theta(3);
s2 = theta(4);
ss = [s1,s2];
y = exp(n(:)*(ww*1i-ss));