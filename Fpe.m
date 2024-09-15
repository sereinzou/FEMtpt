function FokkerPlanckEquation
% 定义方程参数
alpha = 1.0;
beta = 1.0;
gamma = 1.0;
L = 10.0;

% 定义方程
m = 0;
pdefun = @(x,t,u,DuDx) [gamma*DuDx(1) + alpha*DuDx(2); -DuDx(1) - beta*DuDx(2)];
icfun = @(x) [exp(-x^2); 0];
bcfun = @(xl,ul,xr,ur,t) [ul(2); ur(1)];

% 定义求解区间和离散化
x = linspace(-L/2, L/2, 201);
t = linspace(0, 1, 101);

% 调用pdepe函数求解方程
sol = pdepe(m, pdefun, icfun, bcfun, x, t);

% 绘制结果
u = sol(:,:,1);
surf(x, t, u)
xlabel('x')
ylabel('t')
zlabel('u')