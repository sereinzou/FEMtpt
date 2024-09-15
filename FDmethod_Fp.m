

def  get_prob_flux_sparse(f,eps,D,xrange,yrange,Nx,Ny,px=[0,0]):
    Lx = xrange[1]-xrange[0]
    Ly = yrange[1]-yrange[0]
    N  = Nx*Ny
    hx = Lx/Nx
    hy = Ly/Ny
    D11,D22 = D[0][0],D[1][1]
    def idx(i,j,Nx=Nx,Ny=Ny): return (i-1)*Ny+j-1
    def pos(i,j,hx=hx,hy=hy): return np.array([xrange[0]+i*hx,yrange[0]+j*hy])
    row = []
    col = []
    data = []
    
    % construct the system
    for i = 1:Nx+1
        for j =1:Ny+1
            % J1(i,j-1/2)/hx
            ff = f(np.vstack([pos(i,j-1/2),pos(i-1,j-1/2),pos(i-1/2,j),pos(i-1/2,j-1)]))
            if i<nx: row="row" idx="idx" col="col" idx="idx" j="j" data="data" ff="ff" ff="ff" j1="j1" if="if" i="i"></nx:>1:
                row += [idx(i,j),idx(i,j)]
                col += [idx(i-1,j),idx(i  ,j)]
                data += [-(-.5*ff[1][0]-eps*D11/hx)/hx,
                         -(-.5*ff[1][0]+eps*D11/hx)/hx]

           % J2(i-1/2,j)/hy
            if j<ny: row="row" idx="idx" col="col" idx="idx" data="data" ff="ff" ff="ff" j2="j2" if="if" j="j"></ny:>1:
                row += [idx(i,j),idx(i,j)]
                col += [idx(i,j-1),idx(i,j)]
                data += [-(-.5*ff[3][1]-eps*D22/hy)/hy,
                         -(-.5*ff[3][1]+eps*D22/hy)/hy]
    A          = csc_matrix( (data,(row,col)) )
    print(A.astype)

    % solve the system
    idx_       = idx(np.int_((px[0]-xrange[0])/hx+.5),np.int_((px[1]-yrange[0])/hy+.5))
    ei         = np.zeros(dtype=np.float64,shape=(N,1))
    ei[idx_,0] = 1
    mask       = np.reshape(ei==1,-1)
    ei         = csc_matrix(ei, dtype=np.float64)
    b          = -A@ei
    mask       = ~mask
    A_         = A[:,mask]
    prob_      = spsolve(A_[:-1],b[:-1])
    prob       = np.insert(prob_,idx_,np.array([1]),0)

    % normalization
    prob       = np.maximum(prob,0)
    Z          = prob.mean()*Lx*Ly
    prob       = prob/Z
    return (np.linspace(xrange[0],xrange[1],Nx+1)[:-1]+np.linspace(xrange[0],xrange[1],Nx+1)[1:])/2,\
           (np.linspace(yrange[0],yrange[1],Ny+1)[:-1]+np.linspace(yrange[0],yrange[1],Ny+1)[1:])/2,\
           np.transpose(prob.reshape(Nx,Ny))

class ATwoDSystem(object):
    def __init__(self,xrange=[-2,2],yrange=[-3,3]):
        self.dim    = 2
        self.sigma  = np.diag([np.sqrt(1./50),np.sqrt(1./5)])
        self.D      = self.sigma@np.transpose(self.sigma)
        self.eps    = np.linalg.norm(self.D,2)/2
        self.Dbar   = self.D/self.eps/2
        self.D_list = [self.D/2,self.D,self.D*2]
        self.xrange = xrange
        self.yrange = yrange
    def get_f(self,X):
        if np.size(X.shape)==2: x,y = X[:,0][:,None],X[:,1][:,None]
        else: x,y = X[0],X[1]
        return np.hstack([.1*2*x*(1-x**2) + (1+np.sin(x))*y,
                          -y + 2*x*(1+np.sin(x))*(1-x**2)])
    def get_FD_sol(self,D_list,xrange,yrange,Nx=200,Ny=200,px=[-1,0]):
        self.FD_Y   = []
        for i in range(len(D_list)):
            eps        = np.linalg.norm(D_list[i],2)/2
            Dbar       = D_list[i]/eps/2
            xx,yy,prob = utils.get_prob_flux_sparse(self.get_f,eps,Dbar,xrange,yrange,Nx,Ny,px)
            self.FD_Y.append(prob.reshape(-1))
            utils.plot_epslog_prob(prob,eps,xx,yy,0,0,Vmax=eps*15)
        XX,YY     = np.meshgrid(xx,yy)
        self.FD_X = np.concatenate([XX[:,:,None],YY[:,:,None]],axis=-1).reshape(-1,2)
    def get_P(self,X,D): 
        for i in range(len(self.D_list)):
            if abs(D-self.D_list[i]).max()&lt;1e-5:
                return np.maximum(griddata(self.FD_X,self.FD_Y[i],X,method='cubic',fill_value=0),0)
    def get_V(self,X,D): return -np.linalg.norm(D,2)/2*np.log(self.get_P(X,D)+1e-40)

AS = ATwoDSystem()
AS.get_FD_sol(D_list=AS.D_list,xrange=AS.xrange,yrange=AS.yrange)