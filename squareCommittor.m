close all;
clear variables;
[node,elem] = squaremesh1([-2 2 -3 3], 0.1,0.1);
%用于生成一个二维正方形网格，并将其节点坐标和单元信息分别存储在node和elem数组中
 for k = 1:4
     [node,elem] = uniformrefine(node,elem);
     %[node,elem] = uniformbisect(node,elem);
 end
 %对网格进行四次均匀细化，每次细化后的新网格将覆盖原来的网格。

%用于求解定义在该网格上的一个偏微分方程，其中cdata是一个包含该方程的参数值的结构体，Committor是一个自己定义的函数，
% 用于求解该方程的数值解，((x+1).^2+y.^2<0.1) | ((x-1).^2+y.^2<0.1)是该方程的边界条件。
pde=cdata;
[soln,eqn,info] = Committor(node,elem,pde,'((x+1).^2+y.^2<0.1) | ((x-1).^2+y.^2<0.1)');
[node1,ind]=sortrows(node); %按坐标大小排序
node_x=unique(node1(:,1));  %抽取x坐标中唯一的量
node_y=unique(node1(:,2));  %抽取y坐标中唯一的量
u=soln.u(ind); 
[x,y]=meshgrid(node_x,node_y);
%用于重新排序节点坐标，并将其按照x坐标和y坐标分别抽取出唯一的值，最终得到一个坐标矩阵(x,y)，以及对应的数值解u。


%% Compute the intensity of S1/2
% nx=size(node_x,1);
% ny=size(node_y,1);
% hx=3/(nx-1); %x轴上划分的大小
% hy=2/(ny-1); %y轴上划分的大小
% V=@(x,y) exp(-2.5*(x.^2-1).^2-5*(y.^2));
% Z=integral2(V,-inf,inf,-inf,inf);
% idx=(node1(:,1)==0);
% x0=find(idx);
% x1=x0+ny;
% a1=u(x0);a2=u(x1);
% p_x=(a1-a2)./hx;
% Jab_x=Z^(-1).*exp(-2.5*(node1(x0,1).^2-1).^2-5*(node1(x0,2).^2)).*p_x;
% subplot(1,2,1)
% plot(node1(x0,2),Jab_x,'-');

% %% Compute the streamlines 
% [A,sl,sl_colour]=streamlines(node1,u);
% subplot(1,2,1)
% plot(sl(:,1),sl(:,2),'o');
% axis([-2 2 -3 3])

%% Graph the committor function
[m,n]=size(x); 
uh=reshape(u,[m,n]); %对u进行重组
%用于将数值解u重组为一个矩阵，其大小为(m,n)，其中m和n分别是x和y坐标矩阵的大小。
%a=flipud(gray);
% subplot(1,2,2)
figure(2)

 contourf(x,y,uh);
 title('q+(x)');
% colormap(a);
%contourf(x,y,-0.1*(log(uh)));

 hold on;
syms x y ;
%k = @(x,y) 2.5*(x^2-1)^2+5*y^2-0.4;
k = @(x,y) (x+1).^2+y.^2-0.1;
fimplicit(k,[-2 2 -3 3],'k')  % 指定绘图范围
colorbar;
% figure;
% showsolution(node,elem,soln.u);