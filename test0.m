%% START 只有两个参数ω和σ的情况
clear,clc
w0 = 0.01*2*pi; sgm0 = 0.1; % 初始化ω和σ
N = 100;
n = 0:N-1;
yn = exp((1i*w0-sgm0)*n); % 产生理想值/采样值
ww = (0:0.001:0.1)*2*pi; 
aa = 0:0.01:1;
E = zeros(length(ww),length(aa));
for k1 = 1:length(ww)
    for k2 = 1:length(aa)
        G = exp((1i*ww(k1)-aa(k2))*n); % 产生真实值
        E(k1,k2) = norm(yn(:)-G(:),2); % 损失函数，为2-范数
        % norm(v)返回向量v的欧几里德范数。此范数也称为2-范数、向量模或欧几里德长度。
    end
end

figure
mesh(aa,ww,E)
% mesh()将矩阵E中的值绘制为由X和Y定义的x-y平面中的网格上方的高度。边颜色因Z指定的高度而异。

%% 从只有2个参数变为有3个参数，ω、σ和A
clear,clc
w0 = [0.01, 0.03, 0.05, 0.07]*2*pi; % 初始化参数ω
sgm0 = [0.1, 0.08, 0.06, 0.04]; % 初始化参数σ
A0 = [0.1, 0.4, 0.8, 1]; % 初始化参数A
N = 100;
n = 0:N-1;
yn = exp(n(:)*(w0*1i-sgm0))*A0(:); % 理想值/采样值

%% 1D %%
% 一维线性插值法
clear,clc
w0 = 0.01*2*pi; % 初始化参数ω
sgm0 = 0.1; % 初始化参数σ
N = 100;
n = 0:N-1;
yn = exp(n(:)*(w0*1i-sgm0)); % 理想值/采样值
theta1 = [0; 0]; % 参数向量θ
theta2 = [0.03*2*pi; 0.3]; % 参数向量θ'
alpha = 0:0.01:1; % 标量参数α
theta_vec = theta1*(1-alpha)+theta2*alpha; % 加权平均θ(α)
E = zeros(length(alpha), 1);
for k = 1:length(alpha)
    G = sig_gen_1D(theta_vec(:,k),n); % 根据θ(α)产生的信号，即估计值
    E(k) = norm(yn(:)-G(:),2); % 一维损失函数，为2-范数
end
figure
plot(alpha, E)
xlabel('标量参数\alpha')
ylabel('损失函数')
title('一维线性插值法得到损失函数')
grid on

%% 1D %%
% 一维线性插值法，将w0和sgm0变为向量
clear,clc
w0 = [0.01, 0.02]*2*pi; % 初始化参数ω，ω1为0.02π，ω2为0.04π
sgm0 = [0.1, 0.05]; % 初始化参数σ，σ1为0.1，σ2为0.05
N = 100;
n = 0:N-1;
yn = exp(n(:)*(w0*1i-sgm0)); % 理想值/采样值
theta1 = [0.05*2*pi; 0.01*2*pi; 0; 0.2]; % 参数向量θ
theta2 = [0; 0.03*2*pi; 0.05; 0.1]; % 参数向量θ'
alpha = 0:0.01:1; % 标量参数α
theta_vec = theta1*(1-alpha)+theta2*alpha; % 加权平均θ(α)
E = zeros(length(alpha), 1);
for k = 1:length(alpha)
    G = sig_gen_1D2(theta_vec(:,k), n); % 根据θ(α)产生的信号，即估计值
    E(k) = norm(yn(:)-G(:),2); % 一维损失函数，为2-范数
end
figure
plot(alpha, E)
xlabel('标量参数\alpha')
ylabel('损失函数')
title({'一维线性插值法得到损失函数'; '当损失函数为0，对应的即为最优解'})
grid on

%% p*=p，此时p=1，p*=1
% 等高线图和随机方向法
clear,clc
close all
w0 = 0.01*2*pi; % 初始化参数ω
sgm0 = 0.1; % 初始化参数σ
A0 = 2; % 增加一个参数A，并对其初始化
N = 100;
n = 0:N-1;
yn = exp(n(:)*(w0*1i-sgm0))*A0(:); % 理想值/采样值
theta0 = [w0; sgm0; A0]; % 选择一个中心点θ*
% ↓↓↓在以下修改方向向量δ和η的值↓↓↓
delta = [1; 0; 0];
eta = [0; 0; 0];
% ↑↑↑在以上修改方向向量δ和η的值↑↑↑
alpha = -1:0.01:1;
beta = alpha;
E = zeros(length(alpha),length(beta));
for k1 = 1:length(alpha)
    for k2 = 1:length(beta)
        theta_v = theta0+alpha(k1)*delta+beta(k2)*eta; % 加权平均得到参数θ(α,β)
        G = sig_gen_2D(theta_v,n); % 根据θ(α,β)产生的信号，即估计值
        E(k1,k2) = norm(yn(:)-G(:),2); % 二维损失函数，为2-范数
    end
end

figure
subplot(1,2,1) % 绘制网格图
mesh(beta,alpha,E)
xlabel('\beta')
ylabel('\alpha')
title('损失函数网格图')
subplot(1,2,2) % 绘制等高线图
contour(alpha,beta,E,LineWidth=1.5)
xlabel('\beta')
ylabel('\alpha')
title('损失函数等高线图')

%% p*＞p，此时p=1，p*=2
clear,clc
close all
w0 = 0.1*2*pi; % 初始化参数ω
% 此时p=1
sgm0 = 0.1; % 初始化参数σ
A0 = 2; % 初始化参数A
N = 100;
n = 0:N-1;
yn = exp(n(:)*(w0*1i-sgm0))*A0(:); % 采样信号
theta0 = [w0; 0; sgm0; 0;  A0; 0]; % 选择一个中心点θ*
% 此时p*=2；θ*的维数6*1
% ↓↓↓在以下修改方向向量δ和η的值↓↓↓
delta = [1*pi; -2*pi; -0.01; 0.03; -0.4; 1.5];
eta = [-1*pi; 4*pi; 0.02; 0.02; 0.3; -1.6];
% ↑↑↑在以上修改方向向量δ和η的值↑↑↑
alpha = -1:0.01:1;
beta = alpha;
E = zeros(length(alpha),length(beta));
for k1 = 1:length(alpha)
    for k2 = 1:length(beta)
        theta_v = theta0+alpha(k1)*delta+beta(k2)*eta; % 加权平均得到参数θ(α,β)
        G = sig_gen_p2(theta_v,n); % 根据θ(α,β)产生的信号，即估计值
        E(k1,k2) = norm(yn(:)-G(:),2); % 二维损失函数，为2-范数
    end
end

figure
subplot(1,2,1) % 绘制网格图
mesh(beta,alpha,E)
xlabel('\beta')
ylabel('\alpha')
title('损失函数网格图')
subplot(1,2,2) % 绘制等高线图
contour(alpha,beta,E)
xlabel('\beta')
ylabel('\alpha')
title('损失函数等高线图')
