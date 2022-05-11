%% p*<p，此时p=2，p*=1，欠参
clear, clc
close all
w0 = [0.02, 0.05]*2*pi; % 初始化参数ω1、ω2
sgm0 = [0.03, 0.04]; % 初始化参数σ1、σ2
A0 = [2, 5]; % 初始化参数A1、A2
N = 100;
n = 0:N-1;
xn = A0*exp((w0(:)*1i-sgm0(:))*n); % 采样信号
% A0维度1*2，w0和sgm0展成2*1的列向量，n的维度1*100，矩阵乘法得到x(n)维度1*100
% 实现的是x(n)为求和信号
theta0 = [w0(1); sgm0(1); A0(1)]; % 选择一个中心点θ*
% 此时p*=1，欠参；θ*的维数3*1

% ↓↓↓在以下修改方向向量δ和η的值↓↓↓
delta = [0.5; 0.001; 0.4];
eta = [0.7; 0.002; 0.3];
% ↑↑↑在以上修改方向向量δ和η的值↑↑↑
alpha = -1:0.01:1;
beta = alpha;
E = zeros(length(alpha),length(beta));
THETA = zeros(length(alpha),length(beta), 3);
for k1 = 1:length(alpha)
    for k2 = 1:length(beta)
        theta_v = theta0+alpha(k1)*delta+beta(k2)*eta; % 加权平均得到参数θ(α,β)
        THETA(k1, k2, 1) = theta_v(1);
        THETA(k1, k2, 2) = theta_v(2);
        THETA(k1, k2, 3) = theta_v(3);
        G = sig_gen_2D(theta_v,n); % 根据θ(α,β)产生采样信号
        % 此时产生信号的函数变为sig_gen_p2，对应的是p*=2的情况
        E(k1,k2) = norm(xn(:)-G(:),2); % 二维损失函数，为2-范数
    end
end
[ind1, ind2] = find(E == min(min(E)));

figure
subplot(1,2,1) % 绘制网格图
mesh(beta,alpha,E)
xlabel('\beta')
ylabel('\alpha')
title({'损失函数网格图'; 'p*=1, p=2'} )
view(-30,10)
hold on
plot3(beta(ind2), alpha(ind1), E(ind1, ind2), 'ro', 'LineWidth', 2)
text(beta(ind2), alpha(ind1), E(ind1, ind2), ' \leftarrow 最优解', 'Color','red','FontSize',14)

subplot(1,2,2) % 绘制等高线图
contour(alpha,beta,E,LineWidth=1.2)
xlabel('\beta')
ylabel('\alpha')
title({'损失函数等高线图'; 'p*=1, p=2'} )
hold on
plot(beta(ind2), alpha(ind1), 'r*', 'LineWidth', 1)
text(beta(ind2), alpha(ind1),  ' \leftarrow 最优解', 'Color','red','FontSize',14)


%% 画出最优解以及四个角对应的信号
best_theta = theta0+alpha(ind1)*delta+beta(ind2)*eta; % 最优解对应的θ
xrn = exp(n(:)*(best_theta(1)*1i-best_theta(2)))*best_theta(3); % 最优解对应的预测信号
x11 = exp(n(:)*(THETA(201, 201, 1)*1i-THETA(201, 201, 2)))*THETA(201, 201, 3); % α=1，β=1对应的信号
x12 = exp(n(:)*(THETA(201, 1, 1)*1i-THETA(201, 1, 2)))*THETA(201, 1, 3); % α=1，β=-1对应的信号
x21 = exp(n(:)*(THETA(1, 201, 1)*1i-THETA(1, 201, 2)))*THETA(201, 1, 3); % α=-1，β=1对应的信号
x22 = exp(n(:)*(THETA(1, 1, 1)*1i-THETA(1, 1, 2)))*THETA(1, 1, 3); % α=-1，β=-1对应的信号

figure
subplot(3,2,1)
plot(n, real(xn), 'LineWidth', 2)
% xlabel('n')
xlabel({'n'; '(a)'}, 'FontSize', 12 )
title({'采样信号x(n)的实部'; ['理想参数θ*=[', num2str(theta0(1)),'; ' ...
    , num2str(theta0(2)), '; ', num2str(theta0(3)), ']']})
grid on

subplot(3,2,2)
plot(n, real(xrn), 'LineWidth', 2)
xlabel({'n'; '(b)'}, 'FontSize', 12 )
title({'最优解对应的预测信号x_r(n)的实部'; ['最优解处θ=[', num2str(best_theta(1)),'; ' ...
    , num2str(best_theta(2)), '; ', num2str(best_theta(3)), ']']})
grid on

subplot(3,2,3)
plot(n, real(x11), 'LineWidth', 2)
xlabel({'n'; '(c)'}, 'FontSize', 12 )
title({'α=1，β=1对应的信号的实部'; ['此时θ=[', num2str(THETA(201, 201, 1)),'; ' ...
    , num2str(THETA(201, 201, 2)), '; ', num2str(THETA(201, 201, 3)), ']']})
grid on

subplot(3,2,4)
plot(n, real(x12), 'LineWidth', 2)
xlabel({'n'; '(d)'}, 'FontSize', 12 )
title({'α=1，β=-1对应的信号的实部'; ['此时θ=[', num2str(THETA(201, 1, 1)),'; ' ...
    , num2str(THETA(201, 1, 2)), '; ', num2str(THETA(201, 1, 3)), ']']})
grid on

subplot(3,2,5)
plot(n, real(x22), 'LineWidth', 2)
xlabel({'n'; '(e)'}, 'FontSize', 12 )
title({'α=-1，β=-1对应的信号的实部'; ['此时θ=[', num2str(THETA(1, 1, 1)),'; ' ...
    , num2str(THETA(1, 1, 2)), '; ', num2str(THETA(1, 1, 3)), ']']})
grid on

subplot(3,2,6)
plot(n, real(x21), 'LineWidth', 2)
xlabel({'n'; '(f)'}, 'FontSize', 12 )
title({'α=-1，β=1对应的信号的实部'; ['此时θ=[', num2str(THETA(1, 201, 1)),'; ' ...
    , num2str(THETA(1, 201, 2)), '; ', num2str(THETA(1, 201, 3)), ']']})
grid on
