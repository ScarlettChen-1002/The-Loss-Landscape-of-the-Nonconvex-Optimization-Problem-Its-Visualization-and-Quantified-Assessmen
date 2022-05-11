%% NUS(Not-Uniform Sampling)非均匀采样
%% 此时p=p*=2
clear,clc
close all
w0 = [0.02, 0.05]*2*pi; % 初始化参数ω1、ω2
sgm0 = [0.03, 0.04]; % 初始化参数σ1、σ2
A0 = [2, 5]; % 初始化参数A1、A2

N = 100;
n = 0:N-1;
xn = A0*exp((w0(:)*1i-sgm0(:))*n); % 原采样信号
p = 25; % 非均匀采样的点数
[xn_mask, v] = SinPoisson(p, N);
xn_nus = xn.*xn_mask;

% 展示非均匀采样的结果
figure
subplot(2,1,1)
plot(n, real(xn), 'LineWidth', 1.5)
xlabel({'n'; '(a)'}, 'FontSize', 12 )
title('原信号x(n)', 'LineWidth', 2)
grid on

subplot(2,1,2)
plot(n, real(xn_nus), 'LineWidth', 1.5)
xlabel({'n'; '(b)'}, 'FontSize', 12 )
title('非均匀采样信号x\_nus(n)，采样率25%')
grid on

%% 对非均匀采样信号xn_nus进行预测
close all
theta0 = [w0(:); sgm0(:); A0(:)]; % 选择一个中心点θ*
% 此时p*=2；θ*的维数6*1

% ↓↓↓在以下修改方向向量δ和η的值↓↓↓
delta = [0.5; 0.8; 0.001; 0.003; 0.4; 0.5];
eta = [0.7; 0.4; 0.002; 0.004; 0.3; 0.6];
% ↑↑↑在以上修改方向向量δ和η的值↑↑↑
alpha = -1:0.01:1;
beta = alpha;
E = zeros(length(alpha),length(beta));
for k1 = 1:length(alpha)
    for k2 = 1:length(beta)
        theta_v = theta0+alpha(k1)*delta+beta(k2)*eta; % 加权平均得到参数θ(α,β)
        G = sig_gen_p2(theta_v,n); % 根据θ(α,β)产生的信号，即估计值
        E(k1,k2) = norm(xn_nus(:)-G(:),2); % 二维损失函数，为2-范数
    end
end
[ind1, ind2] = find(E == min(min(E)));

figure
subplot(1,2,1) % 绘制网格图
mesh(beta,alpha,E)
xlabel('\beta')
ylabel('\alpha')
title({'损失函数网格图'; 'p*=p=2, 非均匀采样，采样率25%'} )
view(-30,10)
hold on
plot3(beta(ind2), alpha(ind1), E(ind1, ind2), 'ro', 'LineWidth', 2)
text(beta(ind2), alpha(ind1), E(ind1, ind2), ' \leftarrow 最优解', 'Color','red','FontSize',14)

subplot(1,2,2) % 绘制等高线图
contour(alpha,beta,E,LineWidth=1.5,ShowText="on")
xlabel('\beta')
ylabel('\alpha')
title({'损失函数等高线图'; 'p*=p=2, 非均匀采样，采样率25%'} )
hold on
plot(beta(ind2), alpha(ind1), 'r*', 'LineWidth', 1)
text(beta(ind2), alpha(ind1),  ' \leftarrow 最优解', 'Color','red','FontSize',14)

%% 画出最优解对应的信号
best_theta = theta0+alpha(ind1)*delta+beta(ind2)*eta; % 最优解对应的θ
xrn = exp(n(:)*(best_theta(1)*1i-best_theta(3)))*best_theta(5)+...
    exp(n(:)*(best_theta(2)*1i-best_theta(4)))*best_theta(6); % 最优解对应的预测信号

figure
subplot(2,1,1)
plot(n, real(xn_nus), 'LineWidth', 2)
xlabel({'n'; '(a)'}, 'FontSize', 12 )
title({'非均匀采样信号x\_nus(n)的实部，采样率70%'; ['理想参数θ*=[', num2str(theta0(1)),'; ' ...
    , num2str(theta0(2)), '; ', num2str(theta0(3)), ';', ...
    num2str(theta0(4)), '; ', num2str(theta0(5)), ';', num2str(theta0(6)), ']']})
grid on

subplot(2,1,2)
plot(n, real(xn), 'y', 'LineWidth', 5)
grid on
hold on
plot(n, real(xrn), 'b', 'LineWidth', 1.5)
xlabel({'n'; '(b)'}, 'FontSize', 12 )
title({'蓝色：最优解对应的预测信号x_r(n)的实部；黄色：原信号x(n)'; ['最优解处θ=[', num2str(best_theta(1)),'; ' ...
    , num2str(best_theta(2)), '; ', num2str(best_theta(3)), ';', ...
    num2str(best_theta(4)), '; ', num2str(best_theta(5)), ';', num2str(best_theta(6)), ']']})
