%% 计算损失函数的ε-sharpness
%% p*=p，此时p=1，p*=1
clear, clc
% close all
% w0 = 0.05*2*pi; % 初始化参数ω0，ω0=0则为纯衰减无振荡信号
w0 = 0;
sgm0 = 0.02; % 初始化参数σ0
A0 = 2; % 初始化参数A0
N = 128;
n = 0:N-1;
xn = exp(n(:)*(w0*1i-sgm0))*A0(:); % 采样信号

%% 展示x(n)
% figure
% subplot(4,1,1)
% subplot(4,1,2)
subplot(4,1,3)
% subplot(4,1,4)
plot(n,real(xn),LineWidth=2)
% xlabel({'n'; '(a)'}, 'FontSize', 8)
% xlabel({'n'; '(b)'}, 'FontSize', 8)
xlabel({'n'; '(c)'}, 'FontSize', 8)
% xlabel({'n'; '(d)'}, 'FontSize', 8)

% title('带衰减的复指数信号的实部，N_1=16')
% title('带衰减的复指数信号的实部，N_2=32')
% title('带衰减的复指数信号的实部，N_3=64')
% title('带衰减的复指数信号的实部，N_4=128')

% title('纯衰减无振荡的信号，N_1=16')
% title('纯衰减无振荡的信号，N_2=32')
title('纯衰减无振荡的信号，N_3=64')
% title('纯衰减无振荡的信号，N_4=128')
axis([0,128,0,2])
grid on


%% 
% theta0 = [w0; sgm0; A0]; % 带衰减的复指数
theta0 = [sgm0; A0]; % 纯衰减无振荡
% ↓↓↓在以下修改方向向量δ和η的值↓↓↓

% delta = [0.04; 0.002; 0.5];
% eta = [0.01; 0.01; 2]; % 带衰减的复指数

delta = [0.002; 0.5];
eta = [0.007; 0.4]; % 纯衰减无振荡

% ↑↑↑在以上修改方向向量δ和η的值↑↑↑
alpha = -1:0.01:1;
beta = alpha;
E = zeros(length(alpha),length(beta));
for k1 = 1:length(alpha)
    for k2 = 1:length(beta)
        theta_v = theta0+alpha(k1)*delta+beta(k2)*eta; % 加权平均得到参数θ(α,β)
%         G = sig_gen_2D(theta_v,n); % 根据θ(α,β)产生的信号，即预测信号
        G = sig_gen_decay(theta_v,n); % 纯衰减无振荡
        E(k1,k2) = norm(xn(:)-G(:),2); % 二维损失函数，为2-范数
    end
end

%% 求损失函数的海森矩阵（Hessian Martix）和锐度
HL = hessian_m(w0, sgm0, A0, xn, N); %求损失函数的海森矩阵▽^2(L(θ))
S = norm(HL, 2); % 求海森矩阵的谱范数，也即锐度

% 展示损失函数图像
% figure
% subplot(1,4,1) % 绘制网格图
% subplot(1,4,2)
% subplot(1,4,3)
subplot(1,4,4)
mesh(beta,alpha,E)

% xlabel({'\beta'; '(a)'}, 'FontSize', 8 )
% xlabel({'\beta'; '(b)'}, 'FontSize', 8 )
% xlabel({'\beta'; '(c)'}, 'FontSize', 8 )
xlabel({'\beta'; '(d)'}, 'FontSize', 8 )

ylabel('\alpha')
zlim([0,2])
view(-30,5)
title({'损失函数网格图，ω_0=0'; ['N=', num2str(N), '，锐度=', num2str(S)]} )



