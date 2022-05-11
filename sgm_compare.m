%% 比较当σ为0、σ较大、σ较小时损失函数曲面的形态
%% 令采样信号σ=0，p*=p=1
clear, clc
close all
w0 = 0.04*2*pi; % 初始化参数ω0
sgm0 = 0; % 令σ=0
A0 = 2; % 初始化参数A0
N = 100;
n = 0:N-1;
xn = exp(n(:)*(w0*1i-sgm0))*A0(:); % 采样信号

figure % 绘制无衰减信号x(n)
plot(n, real(xn),'LineWidth',2)
xlabel('n')
title('无衰减的理想信号x(n)的实部')
grid on

%% 
close all
theta0 = [w0; sgm0; A0]; % 选择一个中心点θ*
% ↓↓↓在以下修改方向向量δ和η的值↓↓↓
delta = [0.04; 0.006; 1];
eta = [0.01; 0.007; 2];
% ↑↑↑在以上修改方向向量δ和η的值↑↑↑
alpha = -1:0.01:1;
beta = alpha;
E = zeros(length(alpha),length(beta));
for k1 = 1:length(alpha)
    for k2 = 1:length(beta)
        theta_v = theta0+alpha(k1)*delta+beta(k2)*eta; % 加权平均得到参数θ(α,β)
        G = sig_gen_2D(theta_v,n); % 根据θ(α,β)产生的信号，即预测信号
        E(k1,k2) = norm(xn(:)-G(:),2); % 二维损失函数，为2-范数
    end
end
[ind1, ind2] = find(E == min(min(E)));

figure
subplot(1,2,1) % 绘制网格图
mesh(beta,alpha,E)
xlabel('\beta')
ylabel('\alpha')
title({'损失函数网格图'; 'p*=p=1, x(n)无衰减'} )
view(-30,10)
hold on
plot3(beta(ind2), alpha(ind1), E(ind1, ind2), 'ro', 'LineWidth', 2)
text(beta(ind2), alpha(ind1), E(ind1, ind2), ' \leftarrow 最优解', 'Color','red','FontSize',14)

subplot(1,2,2) % 绘制等高线图
contour(alpha,beta,E,LineWidth=2,ShowText="on")
xlabel('\beta')
ylabel('\alpha')
title({'损失函数等高线图'; 'p*=p=1, x(n)无衰减'} )
hold on
plot(beta(ind2), alpha(ind1), 'r*', 'LineWidth', 1)
text(beta(ind2), alpha(ind1),  ' \leftarrow 最优解', 'Color','red','FontSize',14)


%% 令采样信号σ=0.01，p*=p=1
clear,clc
close all
w0 = 0.04*2*pi; % 初始化参数ω0
sgm0 = 0.01; % 令σ=0.01
A0 = 2; % 初始化参数A0
N = 100;
n = 0:N-1;
xn = exp(n(:)*(w0*1i-sgm0))*A0(:); % 采样信号

figure % 绘制小衰减信号x(n)
plot(n, real(xn),LineWidth=2)
xlabel('n')
title('小衰减的理想信号x(n)的实部，σ=0.01')
grid on

%%
close all
theta0 = [w0; sgm0; A0]; % 选择一个中心点θ*
% ↓↓↓在以下修改方向向量δ和η的值↓↓↓
delta = [0.04; 0.006; 1];
eta = [0.01; 0.007; 2];
% ↑↑↑在以上修改方向向量δ和η的值↑↑↑
alpha = -1:0.01:1;
beta = alpha;
E = zeros(length(alpha),length(beta));
for k1 = 1:length(alpha)
    for k2 = 1:length(beta)
        theta_v = theta0+alpha(k1)*delta+beta(k2)*eta; % 加权平均得到参数θ(α,β)
        G = sig_gen_2D(theta_v,n); % 根据θ(α,β)产生的信号，即预测信号
        E(k1,k2) = norm(xn(:)-G(:),2); % 二维损失函数，为2-范数
    end
end
[ind1, ind2] = find(E == min(min(E)));

figure
subplot(1,2,1) % 绘制网格图
mesh(beta,alpha,E)
xlabel('\beta')
ylabel('\alpha')
title({'损失函数网格图'; 'p*=p=1, σ=0.01'} )
view(-30,10)
hold on
plot3(beta(ind2), alpha(ind1), E(ind1, ind2), 'ro', 'LineWidth', 2)
text(beta(ind2), alpha(ind1), E(ind1, ind2), ' \leftarrow 最优解', 'Color','red','FontSize',14)

subplot(1,2,2) % 绘制等高线图
contour(alpha,beta,E,LineWidth=2,ShowText="on")
xlabel('\beta')
ylabel('\alpha')
title({'损失函数等高线图'; 'p*=p=1, 衰减因子σ=0.01, 小衰减'} )
hold on
plot(beta(ind2), alpha(ind1), 'r*', 'LineWidth', 1)
text(beta(ind2), alpha(ind1),  ' \leftarrow 最优解', 'Color','red','FontSize',14)




%% 令采样信号σ=0.1，p*=p=1
clear,clc
close all
w0 = 0.01*2*pi; % 初始化参数ω0
sgm0 = 0.1; % 令σ=0.5
A0 = 2; % 初始化参数A0
N = 100;
n = 0:N-1;
xn = exp(n(:)*(w0*1i-sgm0))*A0(:); % 采样信号

figure % 绘制大衰减信号x(n)
plot(n, real(xn),LineWidth=2)
xlabel('n')
title('大衰减的理想信号x(n)的实部，σ=0.1')
grid on

%%
close all
theta0 = [w0; sgm0; A0]; % 选择一个中心点θ*
% ↓↓↓在以下修改方向向量δ和η的值↓↓↓
delta = [0.04; 0.006; 1];
eta = [0.01; 0.007; 2];
% ↑↑↑在以上修改方向向量δ和η的值↑↑↑
alpha = -1:0.01:1;
beta = alpha;
E = zeros(length(alpha),length(beta));
for k1 = 1:length(alpha)
    for k2 = 1:length(beta)
        theta_v = theta0+alpha(k1)*delta+beta(k2)*eta; % 加权平均得到参数θ(α,β)
        G = sig_gen_2D(theta_v,n); % 根据θ(α,β)产生的信号，即预测信号
        E(k1,k2) = norm(xn(:)-G(:),2); % 二维损失函数，为2-范数
    end
end
[ind1, ind2] = find(E == min(min(E)));

figure
subplot(1,2,1) % 绘制网格图
mesh(beta,alpha,E)
xlabel('\beta')
ylabel('\alpha')
title({'损失函数网格图'; 'p*=p=1, 衰减因子σ=0.1, 大衰减'} )
view(-30,10)
hold on
plot3(beta(ind2), alpha(ind1), E(ind1, ind2), 'ro', 'LineWidth', 2)
text(beta(ind2), alpha(ind1), E(ind1, ind2), ' \leftarrow 最优解', 'Color','red','FontSize',14)

subplot(1,2,2) % 绘制等高线图
contour(alpha,beta,E,LineWidth=2,ShowText="on")
xlabel('\beta')
ylabel('\alpha')
title({'损失函数等高线图'; 'p*=p=1, 衰减因子σ=0.1, 大衰减'} )
hold on
plot(beta(ind2), alpha(ind1), 'r*', 'LineWidth', 1)
text(beta(ind2), alpha(ind1),  ' \leftarrow 最优解', 'Color','red','FontSize',14)





%% 令采样信号σ=0，p*=p=2
clear, clc
close all
w0 = [0.02, 0.05]*2*pi;
sgm0 = [0, 0]; % 修改衰减系数
A0 = [2, 4];
N = 100;
n = 0:N-1;
xn = A0*exp((w0(:)*1i-sgm0(:))*n); % 采样信号
% A0维度1*2，w0和sgm0展成2*1的列向量，n的维度1*100，矩阵乘法得到x(n)维度1*100
% 实现的是x(n)为求和信号

figure % 绘制无衰减信号x(n)
plot(n, real(xn),LineWidth=2)
xlabel('n')
title('无衰减的理想信号x(n)的实部')
% title({'小衰减的理想信号x(n)的实部'; 'σ_1=0.01，σ_2=0.03'})
% title({'大衰减的理想信号x(n)的实部'; 'σ_1=0.1，σ_2=0.15'})
grid on

%%
close all
theta0 = [w0(:); sgm0(:); A0(:)]; % 选择一个中心点θ*
% 此时p*=2；θ*的维数6*1

% ↓↓↓在以下修改方向向量δ和η的值↓↓↓
delta = [0.5; 0.8; 0.001; 0.003; 1; 0.7];
eta = [0.7; 0.4; 0.002; 0.008; 0.5; 0.6];
% ↑↑↑在以上修改方向向量δ和η的值↑↑↑

alpha = -1:0.01:1;
beta = alpha;
E = zeros(length(alpha),length(beta));
THETA = zeros(length(alpha),length(beta), 6);
for k1 = 1:length(alpha)
    for k2 = 1:length(beta)
        theta_v = theta0+alpha(k1)*delta+beta(k2)*eta; % 加权平均得到参数θ(α,β)
        THETA(k1, k2, 1) = theta_v(1);
        THETA(k1, k2, 2) = theta_v(2);
        THETA(k1, k2, 3) = theta_v(3);
        THETA(k1, k2, 4) = theta_v(4);
        THETA(k1, k2, 5) = theta_v(5);
        THETA(k1, k2, 6) = theta_v(6);
        G = sig_gen_p2(theta_v,n); % 根据θ(α,β)产生的信号，即估计值
        E(k1,k2) = norm(xn(:)-G(:),2); % 二维损失函数，为2-范数
    end
end
[ind1, ind2] = find(E == min(min(E)));

figure
subplot(1,2,1) % 绘制网格图
mesh(beta,alpha,E)
xlabel('\beta')
ylabel('\alpha')
title({'损失函数网格图'; 'p*=p=2, x(n)无衰减'} )
% title({'损失函数网格图'; 'p*=p=2, x(n)小衰减, σ_1=0.01，σ_2=0.03'} )
% title({'损失函数网格图'; 'p*=p=2, x(n)大衰减, σ_1=0.1，σ_2=0.15'} )
view(-30,10)
hold on
plot3(beta(ind2), alpha(ind1), E(ind1, ind2), 'ro', 'LineWidth', 2)
text(beta(ind2), alpha(ind1), E(ind1, ind2), ' \leftarrow 最优解', 'Color','red','FontSize',14)

subplot(1,2,2) % 绘制等高线图
contour(alpha,beta,E,LineWidth=1.2,ShowText="on")
xlabel('\beta')
ylabel('\alpha')
title({'损失函数等高线图'; 'p*=p=2, x(n)无衰减'} )
% title({'损失函数等高线图'; 'p*=p=2, x(n)小衰减, σ_1=0.01，σ_2=0.03'} )
% title({'损失函数等高线图'; 'p*=p=2, x(n)大衰减, σ_1=0.1，σ_2=0.15'} )
hold on
plot(beta(ind2), alpha(ind1), 'r*', 'LineWidth', 1)
text(beta(ind2), alpha(ind1),  ' \leftarrow 最优解', 'Color','red','FontSize',14)







