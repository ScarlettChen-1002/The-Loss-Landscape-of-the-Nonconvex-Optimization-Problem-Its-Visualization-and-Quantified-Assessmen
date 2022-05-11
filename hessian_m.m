function HL = hessian_m(omg, sgm, A, y, N)
% 计算L(θ)的海森矩阵▽^2(L(θ))
% 损失函数L要关于ω、σ、A三个自变量求偏导
% 所以▽^2(L(θ))是一个三阶矩阵
HL11 = 0; HL12 = 0; HL13 = 0;
HL21 = 0; HL22 = 0; HL23 = 0;
HL31 = 0; HL32 = 0; HL33 = 0;

for n = 1:N-1
    HL11_n = -2*A*n^2*exp(1i*omg*n-sgm*n)*(2*A*(exp(1i*omg*n-sgm*n))-y(n));
    HL12_n = -2*1i*A*n^2*exp(1i*omg*n-sgm*n)*(2*A*exp(1i*omg*n-sgm*n)-y(n));
    HL13_n = 2*1i*n*exp(1i*omg*n-sgm*n)*(2*A*exp(1i*omg*n-sgm*n)-y(n));

    HL21_n = -2*1i*A*n^2*exp(1i*omg*n-sgm*n)*(2*A*exp(1i*omg*n-sgm*n)-y(n));
    HL22_n = 2*A*n^2*exp(1i*omg*n-sgm*n)*(2*A*exp(1i*omg*n-sgm*n)-y(n));
    HL23_n = -2*n*exp(1i*omg*n-sgm*n)*(2*A*exp(1i*omg*n-sgm*n)-y(n));

    HL31_n = 2*1i*n*exp(1i*omg*n-sgm*n)*(2*A*exp(1i*omg*n-sgm*n)-y(n));
    HL32_n = -2*n*exp(1i*omg*n-sgm*n)*(2*A*exp(1i*omg*n-sgm*n)-y(n));
    HL33_n = 2*exp(2*(1i*omg*n-sgm*n));

    HL11 =HL11+HL11_n;
    HL12 =HL12+HL12_n;
    HL13 =HL13+HL13_n;
    HL21 =HL21+HL21_n;
    HL22 =HL22+HL22_n;
    HL23 =HL23+HL23_n;
    HL31 =HL31+HL31_n;
    HL32 =HL32+HL32_n;
    HL33 =HL33+HL33_n;
end

HL = [HL11, HL12, HL13; HL21, HL22, HL23; HL31, HL32, HL33];
% 输出3*3的海森矩阵
end
