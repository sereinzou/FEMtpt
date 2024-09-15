close all;
clear variables;
[node,elem] = squaremesh1([-2 2 -3 3], 0.1,0.1);
%��������һ����ά���������񣬲�����ڵ�����͵�Ԫ��Ϣ�ֱ�洢��node��elem������
 for k = 1:4
     [node,elem] = uniformrefine(node,elem);
     %[node,elem] = uniformbisect(node,elem);
 end
 %����������Ĵξ���ϸ����ÿ��ϸ����������񽫸���ԭ��������

%������ⶨ���ڸ������ϵ�һ��ƫ΢�ַ��̣�����cdata��һ�������÷��̵Ĳ���ֵ�Ľṹ�壬Committor��һ���Լ�����ĺ�����
% �������÷��̵���ֵ�⣬((x+1).^2+y.^2<0.1) | ((x-1).^2+y.^2<0.1)�Ǹ÷��̵ı߽�������
pde=cdata;
[soln,eqn,info] = Committor(node,elem,pde,'((x+1).^2+y.^2<0.1) | ((x-1).^2+y.^2<0.1)');
[node1,ind]=sortrows(node); %�������С����
node_x=unique(node1(:,1));  %��ȡx������Ψһ����
node_y=unique(node1(:,2));  %��ȡy������Ψһ����
u=soln.u(ind); 
[x,y]=meshgrid(node_x,node_y);
%������������ڵ����꣬�����䰴��x�����y����ֱ��ȡ��Ψһ��ֵ�����յõ�һ���������(x,y)���Լ���Ӧ����ֵ��u��


%% Compute the intensity of S1/2
% nx=size(node_x,1);
% ny=size(node_y,1);
% hx=3/(nx-1); %x���ϻ��ֵĴ�С
% hy=2/(ny-1); %y���ϻ��ֵĴ�С
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
uh=reshape(u,[m,n]); %��u��������
%���ڽ���ֵ��u����Ϊһ���������СΪ(m,n)������m��n�ֱ���x��y�������Ĵ�С��
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
fimplicit(k,[-2 2 -3 3],'k')  % ָ����ͼ��Χ
colorbar;
% figure;
% showsolution(node,elem,soln.u);