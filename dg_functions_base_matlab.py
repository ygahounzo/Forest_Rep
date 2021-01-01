
from numpy import *

# Legendre polynomial

def legendre_poly(p,x):

    L1=0;L1_1=0;L1_2=0;
    L0=1;L0_1=0;L0_2=0;

    for i in range(1,p+1):
        L2 = L1
        L2_1 = L1_1
        L2_2 = L1_2
        
        L1 = L0
        L1_1 = L0_1
        L1_2 = L0_2
        
        a = (2*i-1)/i
        b = (i-1)/i
        L0 = a*x*L1 - b*L2
        L0_1 = a*(L1 + x*L1_1) - b*L2_1
        L0_2 = a*(2*L1_1 + x*L1_2) - b*L2_2
    
    return L0, L0_1, L0_2

# The Legendre-Gauss-Lobatto points and weights
#which are the roots of the Lobatto Polynomials.

def legendre_gauss_lobatto(P):

    p=P-1;
    ph = int((p+1)/2);
    xgl = zeros(P)
    wgl = zeros(P)
    for i in range(ph+1):
        x=cos((2*i+1)*pi/(2*p+1))
        for k in range(1,30):
            L0,L0_1,L0_2 =legendre_poly(p,x)

            #Get new Newton Iteration
            dx=-(1-x**2)*L0_1/(-2*x*L0_1 + (1-x**2)*L0_2)
            x=x+dx
            if (abs(dx) < 1.0e-20):
                break

        xgl[p-i]=x
        wgl[p-i]=2/(p*(p+1)*L0**2)
        
    #Check for Zero Root
    if (p == 2*ph):
        x = 0
        L0,L0_1,L0_2 = legendre_poly(p,x)
        xgl[ph] = x
        wgl[ph] = 2/(p*(p+1)*L0**2)

    #Find remainder of roots via symmetry
    for i in range(ph):
        xgl[i] =- xgl[p-i]
        wgl[i] =+ wgl[p-i]
    
    return xgl, wgl

# creation of grid points
def create_grid_dg(ngl,nelem,xgl):
    # constants
    xmin = -1
    xmax = +1
    dx = (xmax-xmin)/nelem
    # Grid points
    ip = 0
    coord = zeros((ngl, nelem))
    intma = zeros((ngl, nelem))
    coord[0,0] = xmin
    for i in range(1,nelem+1):
        x0 = xmin + (i-1)*dx
        coord[0,i-1] = x0
        intma[0,i-1] = ip
        for j in range(1, ngl):
            ip = ip + 1
            coord[j,i-1] = (xgl[j]+1)*(dx/2) + x0
            intma[j,i-1] = ip
            
    return coord, intma

# Flux matrix
def create_Fmatrix_dg(inode,npoin,nelem,ngl,u,diss):

    Fmatrix=zeros((npoin,npoin))

    for e in range(1,nelem+1):

        #Left Side
        left=e-1
        if (e == 1):
            left=nelem
        
        i0 = 0
        iN = ngl
        n_left = -1

        IE = int(inode[i0,e-1])
        IL = int(inode[iN-1,left-1])

        Fmatrix[IE,IE]=n_left*(u/2)*(1+n_left*diss)
        Fmatrix[IE,IL]=n_left*(u/2)*(1-n_left*diss)

        #Right Side
        right=e+1
        if (e == nelem):
            right=1
            
        i0=0
        iN=ngl
        n_right=1

        IE=int(inode[iN-1,e-1])
        IR=int(inode[i0,right-1])

        Fmatrix[IE,IE]=n_right*(u/2)*(1+n_right*diss)
        Fmatrix[IE,IR]=n_right*(u/2)*(1-n_right*diss)
    
    return Fmatrix

# Lagragian basis function
def lagrange_basis3(P,Q,xgl,xnq):

    #Get Quadrature roots
    #xnq,wnq = legendre_gauss_lobatto(Q)
    psi = zeros((P,Q))
    dpsi = zeros((P,Q))
    #Perform Quadrature
    for l in range(Q):
        xl = xnq[l]

        #Construct Basis
        for i in range(P):
            xi = xgl[i]
            psi[i,l] = 1
            dpsi[i,l] = 0
            for j in range(P):
                xj = xgl[j]
                #Basis
                if (j != i):
                    psi[i,l] = psi[i,l]*(xl-xj)/(xi-xj)
                    
                ddpsi=1
                if (j != i):
                    for k in range(P):
                        xk=xgl[k]
                        #Derivative of Basis
                        if (k != i and k != j):
                            ddpsi = ddpsi*(xl-xk)/(xi-xk)


                    dpsi[i,l]=dpsi[i,l] + ddpsi/(xi-xj)
                    
    return psi, dpsi

# mass matrix
def create_mass_dg(coord,nelem,ngl,nq,wnq,psi):

    #Initialize
    mass = {}
    
    for ie in range(1,nelem+1):
        mass_e = zeros((ngl,ngl))
        x = zeros(ngl)
        #Store Coordinates
        for i in range(ngl):
            x[i]=coord[i,ie-1]

        dx=x[ngl-1]-x[0]
        jac=dx/2
        #print(jac)
        #Do LGL Integration
        for l in range(nq):
            wq=wnq[l]#*jac
            
            for i in range(ngl):
                h_i=psi[i,l]
                #print(h_i)
                for j in range(ngl):
                    h_j=psi[j,l]
                    #print(h_j)
                    mass_e[i,j] = mass_e[i,j] + wq*h_i*h_j
                    
        mass[ie] = jac*mass_e
    
    return mass

# differentiation matrix
def create_diff_matrix_dg(ngl,nq,wnq,psi,dpsi):

    #Initialize
    Dmatrix = zeros((ngl,ngl))

    #Integrate Divergence of Flux
    #for ie=1:nelem

    #LGL Integration
    for i in range(0,ngl):
        for j in range(0,ngl):
            for k in range(nq):
                wk=wnq[k]
                dhdx_ik=dpsi[i,k]
                h_jk=psi[j,k]
                Dmatrix[i,j]=Dmatrix[i,j] + wk*dhdx_ik*h_jk
    
    return Dmatrix

# Laplacian matrix
def create_Lmatrix_dg(coord,nelem,ngl,nq,wnq,dpsi):

    #Initialize
    lpla = zeros((ngl,ngl))
    laplacian_matrix = {}
    

    for e in range(1,nelem+1):
        x = zeros(ngl)
        #Store Coordinates
        for i in range(ngl):
            x[i]=coord[i,e-1]

        dx=x[ngl-1]-x[0]
        jac=dx/2
        dksi_dx=2/dx

        #LGL Integration
        for i in range(ngl):
            for j in range(ngl):
                for k in range(nq):
                    wq=wnq[k]*jac
                    dhdx_ik=dpsi[i,k]*dksi_dx
                    dhdx_jk=dpsi[j,k]*dksi_dx
                    lpla[i,j] = lpla[i,j] + wq*dhdx_ik*dhdx_jk
        
        laplacian_matrix[e] = lpla
        
    return laplacian_matrix

# filter
def apply_filter_dg(qp,f,nelem,ngl):

    #Initialize
    rhs = zeros((ngl,nelem))
    q_e = zeros(ngl)
    qf = zeros(ngl)

    #Integrate Divergence of Flux
    for ie in range(1,nelem+1):

        print("new iteration")
        print(ie)
        print(nelem)
        print(len(qp))
        print(len(q_e))

        #Store Local Solution
        for i in range(ngl):
            q_e[i]=qp[i,ie-1]

        #Form Filtered Variable
        qf=f*q_e

        #Store It
        for i in range(ngl):
            rhs[i,ie-1] = qf[i]
    
    return rhs

# Exact solution
def exact_solution_dg(coord,nelem,ngl,time,icase):

    #Set some constants
    w = 1
    h = 0
    xc = 0
    xmin=-1
    xmax=+1
    xl=xmax-xmin
    sigma0=0.125
    rc=0.125
    sigma = sigma0
    u = w*xl
    #icase = 0

    #Initialize
    qe = zeros((ngl,nelem))

    timec = time - floor(time)

    #Generate Grid Points
    for ie in range(1,nelem+1):
        for i in range(ngl):
            x = coord[i,ie-1]
            xbar = xc + u*timec
            if (xbar > xmax): 
                xbar = xmin + (xbar-xmax)
          
            r = x-xbar
            if (icase == 1):
                #qe(i,ie)=sigma0/sigma*exp( -(x-xbar)^2/(4*sigma^2) );
                qe[i,ie-1] = exp(-128*(x-xbar)**2)
                #qe(i,ie)=sigma0/sigma*exp( -(x-xbar)^2/(sigma^2) );
            elif (icase == 2):
                 if (abs(r) <= rc):
                    qe[i,ie-1]=1
                
            elif (icase == 3):
                qe[i,ie]=sigma0/sigma*exp(-(x-xbar)**2/(2*sigma**2));
            elif (icase == 4):
                if (abs(r) <= rc):
                    qe[i,ie-1]=1
                    
            elif (icase == 5):
                if (x <= xc):
                    qe[i,ie-1]=1

            elif (icase == 6):
                qe[i,ie-1] = sin(2*pi*x)

    return qe

# source function
def source_function_dg(coord,nelem,ngl,time,icase):

    #Set some constants
    w=1;
    h=0;
    xmin=-1;
    xmax=+1;
    xl=xmax-xmin;
    sigma0=0.125;
    xc=-0.7;
    #xc=-0.75;
    rc=1.e-6;
    #rc=0.125;
    fc=10;
    sigma = sigma0
    u=w*xl;
    #icase;

    #Initialize
    fe=zeros((ngl,nelem))

    timec=time - floor(time)

    #Generate Grid Points
    for ie in range(1,nelem+1):
        for i in range(ngl):
            x = coord[i,ie-1]
            r = abs(x-xc)
            if (icase <= 2): 
                fe[i,ie-1]=0
            elif (icase == 3):
                #fe(i,ie)=sigma0/sigma*exp( -(x-xc)^2/(2*sigma^2) );
                if (r <= rc): 
                    fe[i,ie-1] =fc

            elif (icase == 4):
                #fe(i,ie)=sigma0/sigma*exp( -(x-xc)^2/(2*sigma^2) );
                if (r <= rc): 
                    fe[i,ie-1] =fc
                    
            elif (icase == 5): 
                fe[i,ie-1] = 0
            elif (icase == 6): 
                fe[i,ie-1]=0
                
    return fe
