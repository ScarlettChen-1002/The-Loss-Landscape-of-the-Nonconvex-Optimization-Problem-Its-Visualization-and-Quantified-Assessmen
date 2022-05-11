%% 作出ε-锐度随N变化的图
clear, clc
close all
% w0 = 0.05*2*pi; % 当讨论纯衰减无振荡信号时，ω0=0
w0 = 0;
sgm0 = 0.02; % 衰减参数σ0
A0 = 2; % 幅度参数A0
% theta0 = [w0; sgm0; A0];
theta0 = [sgm0; A0];

S = zeros(1, 128); % 存储不同N值对应的锐度
for k = 1:128
    N = k;
    n = 0:N-1;
    xn = exp(n(:)*(w0*1i-sgm0))*A0(:); % 采样信号

% delta = [0.04; 0.002; 0.5];
% eta = [0.01; 0.01; 2]; % 带衰减的复指数

delta = [0.002; 0.5];
eta = [0.007; 0.4]; % 纯衰减无振荡

    alpha = -1:0.01:1;
    beta = alpha;
    E = zeros(length(alpha),length(beta));
    for k1 = 1:length(alpha)
        for k2 = 1:length(beta)
            theta_v = theta0+alpha(k1)*delta+beta(k2)*eta; % 加权平均得到参数θ(α,β)
%             G = sig_gen_2D(theta_v,n);  % 带衰减的复指数
            G = sig_gen_decay(theta_v,n); % 纯衰减无振荡
            E(k1,k2) = norm(xn(:)-G(:),2); % 二维损失函数，为2-范数
        end
    end

    HL = hessian_m(w0, sgm0, A0, xn, N); %求损失函数的海森矩阵▽^2(L(θ))
    % 求ε-锐度
    S(k) = norm(HL, 2); % 海森矩阵的谱范数
end

%% 作图
figure
m = 1:128;
plot(m,S,'LineWidth',3)
xlabel('总点数N')
% title('锐度，带衰减的复指数信号')
title('锐度，纯衰减无振荡信号')
grid on

%% y = e^(-n)
clear,clc
close all
n = -3:0.01:3;
yn = exp(-n);
figure
plot(n,yn,LineWidth=2)
xlabel('n')
title('y(n) = e^(^-^n^)')
grid on
hold on
plot(-2.5, yn(n==-2.5), 'ro', 'LineWidth', 1)
text(-2.5, yn(n==-2.5),  ' \leftarrow  y(-2.5)=12.18', 'Color','red','FontSize',14)
hold on
plot(2.5, yn(n==2.5), 'ro', 'LineWidth', 1)
text(1, yn(n==2.5)+0.8,  ' y(2.5)=0.0821  \rightarrow', 'Color','red','FontSize',14)




