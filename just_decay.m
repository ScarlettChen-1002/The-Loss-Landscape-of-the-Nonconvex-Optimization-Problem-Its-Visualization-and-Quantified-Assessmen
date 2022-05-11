%% 纯衰减无振荡信号的研究
%% p*=p=1
clear, clc
close all
sgm0 = 0.03; % 初始化参数σ0
A0 = 2; % 初始化参数A0
N = 100;
n = 0:N-1;
xn = exp(n(:)*(-sgm0))*A0(:); % 纯衰减无振荡的采样信号

figure % 展示采样信号x(n)
plot(n, xn, 'LineWidth', 2)
xlabel('n')
title({'纯衰减无振荡信号x(n)'; 'p*=p=1'})
grid on

%% 
close all
theta0 = [sgm0; A0]; % 选择一个中心参数向量θ*
% ↓↓↓在以下修改方向向量δ和η的值↓↓↓
% delta = [0.01; 0.5];
% eta = [0.03; 0.6];
delta = [0.02; 1];
eta = [0.003; 2];
% ↑↑↑在以上修改方向向量δ和η的值↑↑↑

alpha = -1:0.01:1;
beta = alpha;
E = zeros(length(alpha),length(beta));
for k1 = 1:length(alpha)
    for k2 = 1:length(beta)
        theta_v = theta0+alpha(k1)*delta+beta(k2)*eta; % 加权平均得到参数θ(α,β)
        G = sig_gen_decay(theta_v,n); % 根据θ(α,β)产生的信号，即预测信号
        E(k1,k2) = norm(xn(:)-G(:),2); % 二维损失函数，为2-范数
    end
end
[ind1, ind2] = find(E == min(min(E)));

figure
subplot(1,2,1) % 绘制网格图
mesh(beta,alpha,E)
xlabel('\beta')
ylabel('\alpha')
title({'损失函数网格图'; 'p*=p=1，纯衰减信号'} )
view(-30,10)
hold on
plot3(beta(ind2), alpha(ind1), E(ind1, ind2), 'ro', 'LineWidth', 2)
text(beta(ind2), alpha(ind1), E(ind1, ind2), ' \leftarrow 最优解', 'Color','red','FontSize',14)

subplot(1,2,2) % 绘制等高线图
contour(alpha,beta,E)
xlabel('\beta')
ylabel('\alpha')
title({'损失函数等高线图'; 'p*=p=1，纯衰减信号'} )
hold on
plot(beta(ind2), alpha(ind1), 'r*', 'LineWidth', 1)
text(beta(ind2), alpha(ind1),  ' \leftarrow 最优解', 'Color','red','FontSize',14)

%% 不同的总点数N对损失函数的影响
%% 展示三个N下不同的x(n)
clear, clc
close all
sgm0 = 0.03; % 初始化参数σ0
A0 = 2; % 初始化参数A0
N1 = 16; N2 = 32; N3 = 64;
n1 = 0:N1-1; n2 = 0:N2-1; n3 = 0:N3-1;
xn1 = exp(n1(:)*(-sgm0))*A0(:); % N1=32的x1(n)
xn2 = exp(n2(:)*(-sgm0))*A0(:); % N2=64的x2(n)
xn3 = exp(n3(:)*(-sgm0))*A0(:); % N3=128的x3(n)

figure 
subplot(311)
plot(n1,xn1, 'LineWidth', 2)
xlabel({'n'; '(a)'}, 'FontSize', 12 )
title('x_1(n)，N_1 = 16')
axis([0, 70, 0, 2])
grid on

subplot(312)
plot(n2,xn2, 'LineWidth', 2)
xlabel({'n'; '(b)'}, 'FontSize', 12 )
title('x_2(n)，N_2 = 32')
axis([0, 70, 0, 2])
grid on

subplot(313)
plot(n3,xn3, 'LineWidth', 2)
xlabel({'n'; '(c)'}, 'FontSize', 12 )
title('x_3(n)，N_3 = 64')
axis([0, 70, 0, 2])
grid on

%% N1 = 16
clear, clc
close all
sgm0 = 0.03; % 初始化参数σ0
A0 = 2; % 初始化参数A0
N1 = 16;
n = 0:N1-1;
xn = exp(n(:)*(-sgm0))*A0(:); % N1=16的x1(n)

theta0 = [sgm0; A0]; % 选择一个中心参数向量θ*
delta = [0.02; 1];
eta = [0.003; 2];

alpha = -1:0.01:1;
beta = alpha;
E = zeros(length(alpha),length(beta));
for k1 = 1:length(alpha)
    for k2 = 1:length(beta)
        theta_v = theta0+alpha(k1)*delta+beta(k2)*eta; % 加权平均得到参数θ(α,β)
        G = sig_gen_decay(theta_v,n); % 根据θ(α,β)产生的信号，即预测信号
        E(k1,k2) = norm(xn(:)-G(:),2); % 二维损失函数，为2-范数
    end
end
[ind1, ind2] = find(E == min(min(E)));

figure
subplot(3,2,1) % 绘制网格图
mesh(beta,alpha,E)
xlabel({'\beta'; '(a)'})
ylabel('\alpha')
title({'损失函数网格图'; 'p*=p=1，纯衰减信号，N_1=16'} )
view(-30,10)
hold on
plot3(beta(ind2), alpha(ind1), E(ind1, ind2), 'ro', 'LineWidth', 2)
text(beta(ind2), alpha(ind1), E(ind1, ind2), ' \leftarrow 最优解', 'Color','red','FontSize',14)

subplot(3,2,2) % 绘制等高线图
contour(alpha,beta,E,'ShowText','on')
xlabel({'\beta'; '(b)'})
ylabel('\alpha')
title({'损失函数等高线图'; 'p*=p=1，纯衰减信号，N_1=16'} )
hold on
plot(beta(ind2), alpha(ind1), 'r*', 'LineWidth', 1)
text(beta(ind2), alpha(ind1),  ' \leftarrow 最优解', 'Color','red','FontSize',14)


%% N2 = 32
clear, clc
% close all
sgm0 = 0.03; % 初始化参数σ0
A0 = 2; % 初始化参数A0
N2 = 32;
n = 0:N2-1;
xn = exp(n(:)*(-sgm0))*A0(:); % 纯衰减无振荡的采样信号

theta0 = [sgm0; A0]; % 选择一个中心参数向量θ*
delta = [0.02; 1];
eta = [0.003; 2];

alpha = -1:0.01:1;
beta = alpha;
E = zeros(length(alpha),length(beta));
for k1 = 1:length(alpha)
    for k2 = 1:length(beta)
        theta_v = theta0+alpha(k1)*delta+beta(k2)*eta; % 加权平均得到参数θ(α,β)
        G = sig_gen_decay(theta_v,n); % 根据θ(α,β)产生的信号，即预测信号
        E(k1,k2) = norm(xn(:)-G(:),2); % 二维损失函数，为2-范数
    end
end
[ind1, ind2] = find(E == min(min(E)));

subplot(3,2,3) % 绘制网格图
mesh(beta,alpha,E)
xlabel({'\beta'; '(c)'})
ylabel('\alpha')
title({'损失函数网格图'; 'p*=p=1，纯衰减信号，N_2=32'} )
view(-30,10)
hold on
plot3(beta(ind2), alpha(ind1), E(ind1, ind2), 'ro', 'LineWidth', 2)
text(beta(ind2), alpha(ind1), E(ind1, ind2), ' \leftarrow 最优解', 'Color','red','FontSize',14)

subplot(3,2,4) % 绘制等高线图
contour(alpha,beta,E,'ShowText','on')
xlabel({'\beta'; '(d)'})
ylabel('\alpha')
title({'损失函数等高线图'; 'p*=p=1，纯衰减信号，N_2=32'} )
hold on
plot(beta(ind2), alpha(ind1), 'r*', 'LineWidth', 1)
text(beta(ind2), alpha(ind1),  ' \leftarrow 最优解', 'Color','red','FontSize',14)

%% N3 = 64
clear, clc
% close all
sgm0 = 0.03; % 初始化参数σ0
A0 = 2; % 初始化参数A0
N3 = 64;
n = 0:N3-1;
xn = exp(n(:)*(-sgm0))*A0(:); % 纯衰减无振荡的采样信号

theta0 = [sgm0; A0]; % 选择一个中心参数向量θ*
delta = [0.02; 1];
eta = [0.003; 2];

alpha = -1:0.01:1;
beta = alpha;
E = zeros(length(alpha),length(beta));
for k1 = 1:length(alpha)
    for k2 = 1:length(beta)
        theta_v = theta0+alpha(k1)*delta+beta(k2)*eta; % 加权平均得到参数θ(α,β)
        G = sig_gen_decay(theta_v,n); % 根据θ(α,β)产生的信号，即预测信号
        E(k1,k2) = norm(xn(:)-G(:),2); % 二维损失函数，为2-范数
    end
end
[ind1, ind2] = find(E == min(min(E)));

subplot(3,2,5) % 绘制网格图
mesh(beta,alpha,E)
xlabel({'\beta'; '(e)'})
ylabel('\alpha')
title({'损失函数网格图'; 'p*=p=1，纯衰减信号，N_3=64'} )
view(-30,10)
hold on
plot3(beta(ind2), alpha(ind1), E(ind1, ind2), 'ro', 'LineWidth', 2)
text(beta(ind2), alpha(ind1), E(ind1, ind2), ' \leftarrow 最优解', 'Color','red','FontSize',14)

subplot(3,2,6) % 绘制等高线图
contour(alpha,beta,E,'ShowText','on')
xlabel({'\beta'; '(f)'})
ylabel('\alpha')
title({'损失函数等高线图'; 'p*=p=1，纯衰减信号，N_2=64'} )
hold on
plot(beta(ind2), alpha(ind1), 'r*', 'LineWidth', 1)
text(beta(ind2), alpha(ind1),  ' \leftarrow 最优解', 'Color','red','FontSize',14)





%% σ的分辨率
clear,clc
close all
sgm0 = [0.1, 0.9]; % 初始化参数σ1、σ2
A0 = [2, 3]; % 初始化参数A1、A2
N = 32;
n = 0:N-1;
xn = A0*exp((-sgm0(:))*n); % 纯衰减无振荡的采样信号

theta0 = [sgm0(:); A0(:)]; % 选择一个中心点θ*
% ↓↓↓在以下修改方向向量δ和η的值↓↓↓
delta = [0.01; 0.03; 0.4; 2];
eta = [0.02; 0.02; 3; 0.6];
% ↑↑↑在以上修改方向向量δ和η的值↑↑↑

alpha = -1:0.01:1;
beta = alpha;
E = zeros(length(alpha),length(beta));
for k1 = 1:length(alpha)
    for k2 = 1:length(beta)
        theta_v = theta0+alpha(k1)*delta+beta(k2)*eta; % 加权平均得到参数θ(α,β)
        G = sig_gen_decay_p2(theta_v,n); % 根据θ(α,β)产生采样信号
        % 此时产生信号的函数变为sig_gen_p2，对应的是p*=2的情况
        E(k1,k2) = norm(xn(:)-G(:),2); % 二维损失函数，为2-范数
    end
end

figure
subplot(1,2,1) % 绘制网格图
mesh(beta,alpha,E)
view(-50,5)
zlim([0,10])
xlabel('\beta')
ylabel('\alpha')
title({'损失函数网格图'; ['p*=p=2，N=', num2str(N)]; ['σ_1=', num2str(sgm0(1)), ...
    '，σ_2=', num2str(sgm0(2)), '，Δσ=', num2str(abs(sgm0(2)-sgm0(1)))]} )
subplot(1,2,2) % 绘制等高线图
contour(alpha,beta,log(E),[-3,-2.5,-2,-1.5,-1,0,1],'ShowText','on',LineWidth=1.5)
xlabel('\beta')
ylabel('\alpha')
xlim([-0.1,0.1])
ylim([-0.1,0.1])
title({'损失函数等高线图，值=lg(E)'; ['p*=p=2，N=', num2str(N)]; ['σ_1=', num2str(sgm0(1)), ...
    '，σ_2=', num2str(sgm0(2)), '，Δσ=', num2str(abs(sgm0(2)-sgm0(1)))]} )


%% 加噪声后的影响
%% p*=p=1
clear,clc
close all
sgm0 = 0.2; % 初始化参数σ0
A0 = 6; % 初始化参数A0，为了不让噪声过度干扰，取较大的幅度
N = 32;
n = 0:N-1;
xn = A0*exp((-sgm0)*n); % 纯衰减无振荡的采样信号
noise1 = rand(1, N); % 较小的噪声信号，幅度为0~1
noise2 = randn(1, N); % 较大的噪声信号，幅度从标准正态分布得到
xn_noise1 = xn+noise1; %加小噪声后的信号
xn_noise2 = xn+noise2; %加大噪声后的信号

figure
subplot(3,1,1)
plot(n, xn,LineWidth=2)
xlabel({'n'; '(a)'}, 'FontSize', 10 )
title('原信号x(n)')
grid on

subplot(3,1,2)
plot(n, xn_noise1,LineWidth=2)
xlabel({'n'; '(b)'}, 'FontSize', 10 )
title('加小噪声后的信号x\_noise1(n)')
grid on

subplot(3,1,3)
plot(n, xn_noise2,LineWidth=2)
xlabel({'n'; '(c)'}, 'FontSize', 10 )
title('加大噪声后的信号x\_noise2(n)')
grid on

%% 
close all
theta0 = [sgm0(:); A0(:)]; % 选择一个中心参数向量θ*
% delta = [0.15; 0.3];
% eta = [0.04; 3]; % δ1，η1
% delta = [0.06; 0.05];
% eta = [0.01; 2]; % δ2，η2
delta = [0.1; 0.4];
eta = [0.05; 3]; % δ3，η3

alpha = -1:0.01:1;
beta = alpha;
E = zeros(length(alpha),length(beta));
E_noise1 = E;
E_noise2 = E;
for k1 = 1:length(alpha)
    for k2 = 1:length(beta)
        theta_v = theta0+alpha(k1)*delta+beta(k2)*eta; % 加权平均得到参数θ(α,β)
        G = sig_gen_decay(theta_v,n); % 根据θ(α,β)产生的信号，即预测信号
        E(k1,k2) = norm(xn(:)-G(:),2); % 未加噪声的损失函数
        E_noise1(k1,k2) = norm(xn_noise1(:)-G(:),2); % 加小噪声后的损失函数
        E_noise2(k1,k2) = norm(xn_noise2(:)-G(:),2); % 加大噪声后的损失函数
    end
end
[i1, i2] = find(E == min(min(E)));
[i11, i12] = find(E_noise1 == min(min(E_noise1)));
[i21, i22] = find(E_noise2 == min(min(E_noise2)));

figure
subplot(2,3,1) % 绘制网格图
mesh(beta,alpha,E)
xlabel({'\beta'; '(a)'})
ylabel('\alpha')
view(-30,5)
zlim([0,10])
title({'损失函数网格图'; '纯衰减无振荡，p*=p=1'; 'δ_3，η_3，未加噪声'} )
hold on
plot3(beta(i2), alpha(i1), E(i1, i2), 'ro', 'LineWidth', 2)
text(beta(i2), alpha(i1), E(i1, i2), ' \leftarrow 最优解', 'Color','red','FontSize',14)

subplot(2,3,4) % 绘制等高线图
contour(alpha, beta, E, LineWidth=1.5,ShowText="on")
xlabel({'\beta'; '(d)'})
ylabel('\alpha')
title({'损失函数等高线图'; '纯衰减无振荡，p*=p=1'; 'δ_3，η_3，未加噪声'} )
hold on
plot(beta(i2), alpha(i1), 'r*', 'LineWidth', 1)
text(beta(i2), alpha(i1),  ' \leftarrow 最优解', 'Color','red','FontSize',14)

subplot(2,3,2) % 绘制网格图
mesh(beta,alpha,E_noise1)
xlabel({'\beta'; '(b)'})
ylabel('\alpha')
view(-30,5)
zlim([0,10])
title({'损失函数网格图'; '纯衰减无振荡，p*=p=1'; 'δ_3，η_3，加小噪声'} )
hold on
plot3(beta(i12), alpha(i11), E_noise1(i11, i12), 'ro', 'LineWidth', 2)
text(beta(i12), alpha(i11), E_noise1(i11, i12), ' \leftarrow 最优解', 'Color','red','FontSize',14)

subplot(2,3,5) % 绘制等高线图
contour(alpha, beta, E_noise1, LineWidth=1.5,ShowText="on")
xlabel({'\beta'; '(e)'})
ylabel('\alpha')
title({'损失函数等高线图'; '纯衰减无振荡，p*=p=1'; 'δ_3，η_3，加小噪声'} )
hold on
plot(beta(i12), alpha(i11), 'r*', 'LineWidth', 1)
text(beta(i12), alpha(i11),  ' \leftarrow 最优解', 'Color','red','FontSize',14)

subplot(2,3,3) % 绘制网格图
mesh(beta,alpha,E_noise2)
xlabel({'\beta'; '(c)'})
ylabel('\alpha')
view(-30,5)
zlim([0,10])
title({'损失函数网格图'; '纯衰减无振荡，p*=p=1'; 'δ_3，η_3，加大噪声'} )
hold on
plot3(beta(i22), alpha(i21), E_noise2(i21, i22), 'ro', 'LineWidth', 2)
text(beta(i22), alpha(i21), E_noise2(i21, i22), ' \leftarrow 最优解', 'Color','red','FontSize',14)

subplot(2,3,6) % 绘制等高线图
contour(alpha, beta, E_noise2, LineWidth=1.5,ShowText="on")
xlabel({'\beta'; '(f)'})
ylabel('\alpha')
title({'损失函数等高线图'; '纯衰减无振荡，p*=p=1'; 'δ_3，η_3，加大噪声'} )
hold on
plot(beta(i22), alpha(i21), 'r*', 'LineWidth', 1)
text(beta(i22), alpha(i21),  ' \leftarrow 最优解', 'Color','red','FontSize',14)

%% 画出最优解对应的信号
best_theta0 = theta0+alpha(i1)*delta+beta(i2)*eta; % 无噪声时最优解对应的θ
best_theta1 = theta0+alpha(i11)*delta+beta(i12)*eta; % 小噪声时最优解对应的θ
best_theta2 = theta0+alpha(i21)*delta+beta(i22)*eta; % 大噪声时最优解对应的θ

x0 = exp(n(:)*(-best_theta0(1)))*best_theta0(2); % 无噪声时最优解对应的预测信号
x1 = exp(n(:)*(-best_theta1(1)))*best_theta1(2); % 小噪声时最优解对应的预测信号
x2 = exp(n(:)*(-best_theta2(1)))*best_theta2(2); % 小噪声时最优解对应的预测信号

subplot(3,1,1)
plot(n, xn, 'y', 'LineWidth', 5)
grid on
hold on
plot(n, x0, 'b', 'LineWidth', 1.5)
xlabel({'n'; '(a)'}, 'FontSize', 10)
title({'蓝色：最优解对应的预测信号；黄色：原信号x(n)'; ['最优解处θ=[', num2str(best_theta0(1)),'; ' ...
    , num2str(best_theta0(2)), ']']})

subplot(3,1,2)
plot(n, xn_noise1, 'y', 'LineWidth', 4)
grid on
hold on
plot(n, x1, 'b', 'LineWidth', 1.5)
xlabel({'n'; '(b)'}, 'FontSize', 10)
title({'蓝色：加小噪声时最优解对应的预测信号；黄色：加小噪声的信号x\_noise1(n)'; ['最优解处θ=[', num2str(best_theta1(1)),'; ' ...
    , num2str(best_theta1(2)), ']']})

subplot(3,1,3)
plot(n, xn_noise2, 'y', 'LineWidth', 4)
grid on
hold on
plot(n, x2, 'b', 'LineWidth', 1.5)
xlabel({'n'; '(c)'}, 'FontSize', 10)
title({'蓝色：加大噪声时最优解对应的预测信号；黄色：加大噪声的信号x\_noise2(n)'; ['最优解处θ=[', num2str(best_theta2(1)),'; ' ...
    , num2str(best_theta2(2)), ']']})






%% y = x^4
clear,clc
close all
x = -5:0.1:5;
fx = x.^4;
figure
plot(x, fx, 'k', 'LineWidth', 2)
xlabel('x')
title('y = x^4')
grid on



