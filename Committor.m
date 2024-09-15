function [soln,eqn,info] = Committor(node,elem,pde,varargin,option)
%% Preprocess
if ~exist('option','var'), option = []; end
% important constants
N = size(node,1); 
NT = size(elem,1);
Ndof = N;

%% Compute geometric quantities and gradient of local basis
[Dphi,area] = gradbasis(node,elem);

%% Assemble stiffness matrix
%quaduature points
time=cputime;
if ~isfield(option,'quadorder')
    option.quadorder = 2;        % default order
    if ~isempty(pde.b) && isnumeric(pde.b) % numerical diffusion
        option.quadorder = 3;    % exact for linear diffusion coefficient
    end
end
[lambda, w] = quadpts(option.quadorder);
nQuad = size(lambda,1);
%compute non-zeros
A = sparse(Ndof,Ndof);
a=pde.a;ad=pde.ad;
for i=1:3
    for j=1:3
        bij=zeros(NT,1);
        adij=zeros(NT,1);
        for p=1:nQuad
            pxy = lambda(p,1)*node(elem(:,1),:) ...
                + lambda(p,2)*node(elem(:,2),:) ...
                + lambda(p,3)*node(elem(:,3),:);
            b=pde.b(pxy);
            bij=bij...
                +lambda(p,i)*(b(:,1).*Dphi(:,1,j)+b(:,2).*Dphi(:,2,j)).*w(p);
            adij=adij...
                +(ad(1)+ad(2)).*lambda(p,i).*Dphi(:,1,j).*w(p)...
                +(ad(3)+ad(4)).*lambda(p,i).*Dphi(:,2,j).*w(p);
        end
        aij=a(1).*Dphi(:,1,i).*Dphi(:,1,j)...
            +a(2).*Dphi(:,1,i).*Dphi(:,2,j)...
            +a(3).*Dphi(:,2,i).*Dphi(:,1,j)...
            +a(4).*Dphi(:,2,i).*Dphi(:,2,j);
        Aij=(bij-aij-adij).*area;
        A = A + sparse(elem(:,i),elem(:,j),Aij,Ndof,Ndof);
    end
end
clear  b a ad bij aij adij Aij

%% Set up boundary conditions
b=zeros(N,1);
u=zeros(N,1);
expr = varargin;
x=node(:,1);
y=node(:,2);
idx = eval(expr);
fixedNode = find(idx);
freeNode=find(~idx);
bdidx = zeros(Ndof,1);
%size(fixednode);
bdidx(fixedNode) = 1;
%bdidx(fixedNode) = zeros(size(fixedNode,1),1);
Tbd = spdiags(bdidx,0,Ndof,Ndof);
T = spdiags(1-bdidx,0,Ndof,Ndof);
AD = T*A*T+ Tbd;
u(fixedNode)=pde.g_D(node(fixedNode,:));
b(fixedNode) = pde.g_D(node(fixedNode,:));
b=b-A*u;

%% Record assembling time
assembleTime = cputime - time;
if ~isfield(option,'printlevel'), option.printlevel = 1; end
if option.printlevel >= 2
    fprintf('Time to assemble matrix equation %4.2g s\n',assembleTime);
end

%% Solve the system of linear equations
u(freeNode) = AD(freeNode,freeNode)\b(freeNode);

%% Compute Du
dudx =  u(elem(:,1)).*Dphi(:,1,1) + u(elem(:,2)).*Dphi(:,1,2) ...
      + u(elem(:,3)).*Dphi(:,1,3);
dudy =  u(elem(:,1)).*Dphi(:,2,1) + u(elem(:,2)).*Dphi(:,2,2) ...
      + u(elem(:,3)).*Dphi(:,2,3);         
Du = [dudx, dudy];

%% Output
if nargout == 1
    soln = u;
else
    soln = struct('u',u,'Du',Du);
    eqn = struct('A',AD,'b',b,'freeNode',freeNode,'Lap',A);
    info.assembleTime = assembleTime;
end

end % end of Poisson
