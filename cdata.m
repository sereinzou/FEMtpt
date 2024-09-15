function pde=cdata
pde=struct('b',@b,'a','ad','g_D',@g_D,'g_N',@g_N);
    function rhs=b(p)
        x = p(:,1); y = p(:,2);
        rhs(:,1)=0.2*x.*(1-x.^2)+y.*(1+sin(x));
        rhs(:,2)=-y+2*x.*(1-x.^2).*(1+sin(x));
    end

pde.a=[0.01 0 0 0.1];

pde.ad=[0 0 0 0];

    function u=g_D(p)
        x=p(:,1); y=p(:,2);
        u=((x-1).^2+y.^2<0.1);
    end

pde.g_N=0;
end