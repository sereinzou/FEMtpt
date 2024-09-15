function [A,sl,s_colour]=streamlines(node,u)
%% �����node����������֮��Ľڵ㣬u��ʾ��Ӧ�ڵ��ֵ

%% Initial data
a=1;
t=0;
V=@(x,y) exp(-2.5*(x.^2-1).^2-5*(y.^2));
Z=integral2(V,-inf,inf,-inf,inf);
N=size(node,1);
node_x=unique(node(:,1));  %��ȡx������Ψһ����
node_y=unique(node(:,2));  %��ȡy������Ψһ����
nx=size(node_x,1);
ny=size(node_y,1);
hx=3/(nx-1); %x���ϻ��ֵĴ�С
hy=2/(ny-1); %y���ϻ��ֵĴ�С

%% Find the boundary of A
expr='(x+1).^2+y.^2<0.1';
x=node(:,1);y=node(:,2);
idx=eval(expr);
a_idx=find(idx);
a_idx1=a_idx-ny;
x=node(a_idx1,1);y=node(a_idx1,2);
nbda_idx=eval(expr);
bda_idx=a_idx(~nbda_idx);
sl=node(bda_idx,:); %record the point of streamlines
sz=size(sl,1);
colour=zeros(sz,1);
clear x y a_idx a_idx1 x1 y1 nbda_idx bda_idx;

%% Eluer iterative method
nn=sz; %��¼ÿ�νڵ������λ��
A=sparse(sz,N); %�����洢streamline·���Ľڵ㣬����ÿһ�б�ʾһ��·��
A=A+sparse(1:sz,1,1:sz,sz,N);
c=1; %��¼����
idx1=zeros(sz,1);
fixedsz=sz;
x0=sl(:,1);y0=sl(:,2);
na=find(~idx1);
while a==1
    t=t+0.05;  %��������
    %% Compute the gradient of committor funciton and Jab
    n_x=round((1.5-x0)./hx); n_y=round((1-y0)./hy);
    idx_x=nx-n_x; idx_y=ny-n_y;
    sln=(idx_x-1).*ny+idx_y;
    slnx=sln-ny; slny=sln-1;
    q_x=(u(sln)-u(slnx))./hx;
    q_y=(u(sln)-u(slny))./hy;
    if all(q_x==0) && all(q_y==0)
        break
    end
    Jab_x=Z^(-1).*exp(-2.5*(x0.^2-1).^2-5*y0.^2).*q_x;
    Jab_y=Z^(-1).*exp(-2.5*(x0.^2-1).^2-5*y0.^2).*q_y;
    
    %% �ж��Ƿ��Ѿ���������S1/2����������Ӧ��streamlines��ɫ���
    idx0=(-hx<x0<hx);
    if ~all(idx0==0)
        cidx=na(idx0);
        colour(cidx)=-Jab_x(idx0);
    end
    nn=nn+sz;
    sl(nn-sz+1:nn,1)=x0+Jab_x*t;
    sl(nn-sz+1:nn,2)=y0+Jab_y*t;
    x2=sl(nn-sz+1:nn,1);y2=sl(nn-sz+1:nn,2);
    c=c+1;
    na=find(~idx1);
    A=A+sparse(na,c,nn-sz+1:nn,fixedsz,N);
    idx2=(((x2+1).^2+y2.^2)<=0.1);
    idx1(na)=idx1(na)+idx2;
    sz=length(idx1)-sum(idx1);
    if sz==0
        break
    end
    x0=x2(~idx2);y0=y2(~idx2);
end
clear idx0 idx1 x2 y2 x3 y3 c na 

%% ����streamline�Ľڵ��Լ���Ӧ��colour
axis([-2 2 -3 3]);
[i,~,s]=find(A);
cs=length(s);
s_colour=zeros(cs,1);
s_colour(s)=colour(i);
for j=1:fixedsz
    hold on
    k=(i==j);
    plot(sl(k,1),sl(k,2));
end