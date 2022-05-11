%% 尝试老师在会议中运行的代码
clear, clc
close all
rng(20);

w0 = 2*pi*[0.11, 0.13, 0.2, 0.4];
sgm0 = 0.01*ones(size(w0));
A0 = rand(size(w0))+0.1;

N = 128*2;
n = 0:N-1;
wn = 0;
xn = A0*exp((1i*w0(:)-sgm0(:))*n);
xn = xn + wn*randn(size(xn)); % 带衰减的采样信号x(n)

Nf = 2^(nextpow2(N)+2);
Xw = fftshift(fft(xn, Nf));
ff = -0.5:1/Nf:0.5-1/Nf;

figure
subplot(2,1,1)
plot(n, real(xn))
xlabel('n')
title('Re[x(n)]，x(n)的实部')

subplot(2,1,2)
plot(ff, abs(Xw))
xlabel('\omega')
title('|X(e^j^\omega)|，x(n)的频谱')

%% 带上估计信号x_r(n)，算损失函数
clear, clc
close all

w0 = 2*pi*[0.11, 0.14, 0.2, 0.4]; % w0为1*4的向量，4=采样信号的求和个数p
sgm0 = 0.01*ones(size(w0));
A0 = rand(size(w0))+0.1; % A0为随机生成的幅度，但要保证≥0.1

N = 128*2;
n = 0:N-1;
wn = 0;
xn = A0*exp((1i*w0(:)-sgm0(:))*n);
xn = xn + wn*randn(size(xn)); % 采样信号x(n)，是一个带衰减的复指数信号
% A0维度1*4，w0和sgm0展成4*1的列向量，n的维度1*256，矩阵乘法得到x(n)维度1*256
% 实现的是x(n)=Σ(k=1~4)x(A_k, sgm_k, w_k)

theta0 = [w0; sgm0; A0]; % 选择一个中心点θ*
% θ*是维度为3*4的参数向量，行数3=参数的个数；列数4=估计信号的求和数量p*
% ↓↓↓在以下修改方向向量δ和η的值↓↓↓
% delta = randn(size(theta0)); % 随机的方向向量δ
delta = [0.2; 0.1; 0.3]; % δ和η是3*1的向量
% eta = randn(size(theta0)); % 随机的方向向量η
eta = [0.5; -0.2; 0.4];
% ↑↑↑在以上修改方向向量δ和η的值↑↑↑
alpha = -1:0.01:1;
beta = alpha;
E = zeros(length(alpha),length(beta));
for k1 = 1:length(alpha)
    for k2 = 1:length(beta)
        theta_v = theta0+alpha(k1)*delta+beta(k2)*eta; % 加权平均得到参数θ(α,β)
        xrn = sig_gen_2D(theta_v,n); % 根据θ(α,β)产生的信号，即估计信号
        E(k1,k2) = norm(xn(:)-xrn(:),2); % 二维损失函数，为2-范数
    end
end

figure % 绘制网格图
% subplot(1,2,1)
mesh(beta,alpha,E)
xlabel('\beta')
ylabel('\alpha')
title('损失函数网格图')
figure % 绘制等高线图
% subplot(1,2,2)
contour(alpha,beta,E)
xlabel('\beta')
ylabel('\alpha')
title('损失函数等高线图')

%% 用循环来算求和Σ
clear, clc
close all

w0 = 2*pi*[0.1, 0.2, 0.3, 0.4]; % 此时p = 4
sgm0 = 0.01*ones(size(w0));
A0 = rand(size(w0))+0.1; % A0为随机生成的幅度，但要保证≥0.1

N = 128*2;
n = 0:N-1;
xn = 0;
for i = 1:4
    xx = A0(i)*exp((1i*w0(i)-sgm0(i))*n);
    xn = xn+xx;
end
% 采样信号x(n)，是一个带衰减的复指数信号

theta0 = [w0; sgm0; A0]; % 选择一个中心点θ*，这个θ*也是参数向量
% ↓↓↓在以下修改方向向量δ和η的值↓↓↓
delta = [1; -0.01; 0.5];
eta = [0; -0.05; 0.02];
% ↑↑↑在以上修改方向向量δ和η的值↑↑↑
alpha = -1:0.01:1;
beta = alpha;
E = zeros(length(alpha),length(beta));
for k1 = 1:length(alpha)
    for k2 = 1:length(beta)
        theta_v = theta0+alpha(k1)*delta+beta(k2)*eta; % 加权平均得到参数向量θ(α,β)
        xrn = sig_gen_2D(theta_v,n); % 根据θ(α,β)产生的信号，即估计信号
        E(k1,k2) = norm(xn(:)-xrn(:),2); % 二维损失函数，为2-范数
    end
end

figure % 绘制网格图
% subplot(1,2,1)
mesh(beta,alpha,E)
xlabel('\beta')
ylabel('\alpha')
title('损失函数网格图')
figure % 绘制等高线图
% subplot(1,2,2)
contour(alpha,beta,E)
xlabel('\beta')
ylabel('\alpha')
title('损失函数等高线图')


