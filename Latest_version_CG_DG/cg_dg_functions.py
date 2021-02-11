
from numpy import *

def Legendre_deriv(Q, x):
    
    '''
    This function compute the Legendre polynomial and its derivative
    
    Inputs:
    -------
            Q  : Integration order(N+1: for exact, N: for inexact integration)
            x  : value of x at which the polynomial is evaluated
            
    Outputs:
    -------
           L1  : Value of the polynomial at x
           dL1 : First derivative
           ddLi: Second derivative
    '''
    
    L0 = 1; dL0 = 0
    L1 = x; dL1 = 1
    ddL0 = 0; ddL1 = 0
    
    for i in range(2, Q+1):
        
        Li = ((2*i-1)/i)*x*L1 - ((i-1)/i)*L0  # iteration of the polynomials
        dLi = i*L1 + x*dL1
        ddLi = (i+1.0)*dL1 + x*ddL1
        
        L0,L1 = L1,Li
        
        dL0,dL1 = dL1,dLi
       
        ddL0,ddL1 = ddL1,ddLi
        
    return L1, dL1, ddL1

def Lobatto_deriv(Q, x):
    
    '''
    This function compute the Lobatto polynomial and its derivative
    
    Inputs:
    -------
            Q  : Integration order(N+1: for exact, N: for inexact integration)
            x  : value of x at which the polynomial is evaluated
            
    Outputs:
    -------
           B  : Value of the polynomial at x
           dB : First derivative
    '''
    
    L,dL, ddL = Legendre_deriv(Q-1, x)
    B = (1.0-x**2)*dL                      # lobatto polynomial
    dB = -2.0*x*dL + (1.0-x**2)*ddL        # derivative of lobatto polynomial   
    
    return B, dB

def Lobatto_p(Q):
    
    '''
    This function compute the Lobatto points
    
    Input:
    -------
            Q  : Integration order(N+1: for exact, N: for inexact integration)
            
    Output:
    -------
           X: array containing the Lobatto points
    '''
    
    X = []                                      # Array that contains legendre points
    K = 100                                     # Order of approximation of Newton method
    e = 1e-20                                   # tolerance
    for i in range(Q):
        xik = cos(((2*i+1)/(2*Q-1))*pi)         # Chebchev points

        for k in range(K):
            out1, out2 = Lobatto_deriv(Q, xik)
            xikk = xik - out1/out2              # approximation of the solution using Newton

            if abs(xikk-xik) < e:

                break

            xik = xikk

        X.append(xikk)
        
    return array(X[::-1])

# Lagrange basis for single value x
def LagrangeBasis(N, i, xl, Xr):
    
    '''
    This function compute the Lagrange polynomial(basis function) and its derivative
    
    Inputs:
    -------
            N  : polynomial order
            i  : ith polynomial 
            xl : values at which the polynomial is evaluated
            Xr : Lobatto points or the roots of the generating polynomial used to construct the basis function
            
    Outputs:
    -------
           L   : Value of the polynomial
           dL  : Derivative
    '''
    L = 1; dL = 0
        
    for j in range(N+1):
            
        prod = 1
        
        if (j != i):
            L = L*(xl-Xr[j])/(Xr[i]-Xr[j])
                
            for k in range(N+1):
                if (k!=i  and k!=j):
                    
                    prod = prod*(xl-Xr[k])/(Xr[i]-Xr[k])
        
            dL = dL+prod/(Xr[i]-Xr[j])
            
    return L, dL

# Lagrange basis for an array that contains value of x
def LagrangeBasis_deriv(N,Q,Xn, Xq):

    l_basis = zeros((N+1,Q+1))
    dl_basis = zeros((N+1,Q+1))
    
    for k in range(Q+1):
        xl = Xq[k]
        
        for i in range(N+1):
            # Call of LagrangeBasis function
            l_basis[i,k], dl_basis[i,k] = LagrangeBasis(N, i, xl, Xn)
            
    return l_basis, dl_basis

# intma function
def intma_cdg(N, Ne, method_type):
    
    '''
    This function compute the intma array for the CG or DG
    
    Inputs:
    -------
            N          : polynomial order
            Ne         : number of element
            method_type: CG or DG
            
    Output:
    -------
           intma: (matrix) that contains intma values
    '''
    
    intma = zeros((N+1,Ne))
    
    # intma for CG
    if (method_type == 'cg'):
        for e in range(1,Ne+1):
        
        
            t = (e-1)*N
            r = N*e
            intmm = []
            for s in range(t, r+1):
                intmm.append(s)
            intma[:,e-1] = array(intmm)
    
    # intma for DG
    if (method_type == 'dg'):
        for e in range(1,Ne+1):
        
            t = int((e-1)*N)
            r = int(e*N)

            intmm = []
            for s in range(t, r+1):
                it = e-1+s
                intmm.append(it)
            intma[:,e-1] = array(intmm)
        
    return intma


#funtion that compute weight values based on quadrature rule
def weight(Q):
    
    '''
    This function compute the weight for the integration
    
    Inputs:
    -------
            Q : Integration order(N+1: for exact, N: for inexact integration)
            
    Output:
    -------
           w : (array) that contains the weight values
    '''
    
    xi = Lobatto_p(Q+1)
    w = zeros(Q+1)
    for i in range(Q+1):
        
        out1, out2, out3 = Legendre_deriv(Q, xi[i])
        w[i] = 2/(Q*(Q+1)*(out1)**2)
        
    return w 

# grid points
def grid_dg(Q,Ne, xe, ax, bx):
    
    '''
    This function compute the weight for the integration
    
    Inputs:
    -------
            Q     : Integration order(N+1: for exact, N: for inexact integration)
            Ne    : Number of elements in the domain
            xe    : grid points within the element
            ax, bx: boundaries
            
    Output:
    -------
           grid: (matrix) that contains all the grid points
    '''
    
    grid = zeros((Q,Ne))

    xel = linspace(ax,bx,Ne+1)

    for e in range(1,Ne+1):
        ae=xel[e-1] ; be=xel[e]

        xsi=((be-ae)/2)*(xe-1) + be

        for i in range(1,Q+1):

            grid[i-1,e-1]=xsi[i-1]
    return grid

# Element mass matrix

def Element_matrix(N,Q, wght,l_basis):
    
    '''
    This function compute the element mass matrix
    
    Inputs:
    -------
            Q      : Integration order(N+1: for exact, N: for inexact integration)
            N      : Polynomial order
            wght   : weights
            l_basis: basis function values
            
    Output:
    -------
           Me: Element mass matrix
    '''
    
    Me = zeros((Q, Q))       # initialisation of the matrix
    
    for k in range(Q+1):
        wk = wght[k]
        
        for i in range(N+1):
            xi = l_basis[i,k]
            
            for j in range(N+1):
                xj = l_basis[j,k]
                
                Me[i,j] = Me[i,j]+ wk*xi*xj

    Me = (1/2)*Me

    return Me

#Differentiation element matrix
def Element_Diff_matrix(N,Q,wght,l_basis,dl_basis):
    
    '''
    This function compute the element differentiation matrix
    
    Inputs:
    -------
            Q       : Integration order(N+1: for exact, N: for inexact integration)
            N       : Polynomial order
            wght    : weights
            l_basis : basis function values
            dl_basis: derivative values of the basis function
            
    Output:
    -------
           De: Element differentiation matrix
    '''
    
    De = zeros((Q, Q))
    for k in range(Q+1):
        wk = wght[k]
        
        for i in range(N+1):
            dxi = dl_basis[i,k]
            
            for j in range(N+1):
                xj = l_basis[j,k]
                
                De[i,j] = De[i,j]+ wk*dxi*xj

    return De

# DSS operator

def DSS_operator(A, Ne, N, Np, intma, periodic, coord, mtype = None):
    
    '''
    This function is the Direct Stiffness Summation operator
    
    Inputs:
    -------
            A       : Matrix ( Me or De, ...)
            N       : Polynomial order
            Ne      : Number of elements
            Np      : Number of global grid points
            intma   : intma array
            periodic: Array that helps to deal with the boundaries and periodicity
            coord   : all the grid points
            mtype   : method used (CG or DG)
            
    Output:
    -------
           M: Global matrix
    '''
    
    M = zeros((Np, Np))
    
    for e in range(1,Ne+1):
        x = coord[:,e-1]
        dx=x[N]-x[0]
        
        for j in range(N+1):
            
            J = int(intma[j,e-1])
            J = int(periodic[J])
            #J = int(J)
            
            for i in range(N+1):
                
                I = int(intma[i, e-1])
                I = int(periodic[I])
                #I = int(I)
                
                # diff = differentiation matrix
                if (mtype == 'diff'):
                    M[I,J] = M[I,J] + A[i,j]
                
                else:
                    M[I,J] = M[I,J] + dx*A[i,j]
                    
    return M


# Flux matrix
def flux_matrix(Ne, intma, N, Np, diss,u):
    
    '''
    This function compute the flux matrix for the DG method
    
    Inputs:
    -------
            N       : Polynomial order
            Ne      : Number of elements
            Np      : Number of global grid points
            intma   : intma array
            u       : velocity
            diss    : 0(centered flux), 1(Rusanov flux)
            
    Output:
    -------
           FM: Flux matrix
    '''
    
    FM = zeros((Np,Np))
    n_left = -1 
    n_right = 1
    for e in range(1, Ne + 1):
        
        # left side
        L = e - 1
        if e == 1:
            L = Ne
        
        I = int(intma[0,e-1])
        J = int(intma[N,L-1])

        FM[I,I] = (u/2)*n_left*(1 + n_left*diss)
        FM[I,J] = (u/2)*n_left*(1 - n_left*diss)
    
        # right side
        R = e+1

        if e == Ne:
            R = 1

        I = int(intma[N,e-1])
        J = int(intma[0,R-1])

        FM[I,I] = (u/2)*n_right*(1 + n_right*diss)
        FM[I,J] = (u/2)*n_right*(1 - n_right*diss)
    
    return FM

# Exact solution
def exact_solution_dg(coord,Ne,Q,time, ax, bx, u, case):
    
    '''
    This function compute the exact solution
    
    Inputs:
    -------
            Q       : Integration order(N+1: for exact, N: for inexact integration)
            N       : Polynomial order
            Ne      : Number of elements
            coord   : all the grid points
            u       : velocity
            ax, bx  : Left and right boundaries of the physical domain
            case    : For the initial condition type(icase = 1, for gaussian or 2, for sine)
            time    
            
    Output:
    -------
           qe : exact values
    '''
    
    #Initialize
    qe = zeros((Q,Ne))

    timec = time - floor(time)

    #Generate Grid Points
    for e in range(1, Ne + 1):
        for i in range(Q):
            x = coord[i,e - 1]
            xbar = u*timec
            if (xbar > bx): 
                xbar = ax + (xbar - bx)
            
            if(case == 1):
                qe[i,e - 1] = exp(-128*(x - xbar)**2)
            
            elif(case == 2):
                qe[i, e-1] = sin(2*pi*(x-xbar))
                
    return qe


def cg_dgSolver(N,Q,nel,Np, ax, bx, integration_type, method_type, icase,diss,u,Courant_max, time_final, kstages):
    
    '''
    This function is CG/DG solver for 1D wave equqtion
    
    Inputs:
    -------
            N              : Polynomial order
            Q              : Integration order(N+1: for exact, N: for inexact integration)
            nel            : Number of element
            Np             : Global grid points(nel*N+1 for CG, nel*(N+1) for DG)
            ax, bx         : Left and right boundaries of the physical domain
            intergatio_type: Exact or Inexact integration
            method_type    : CG or DG
            icase          : For the initial condition type(icase = 1, for gaussian or 2, for sine)
            diss           : 0(centered flux), 1(Rusanov flux)
            u              : velocity
            Courant_max    : CFL
            time_final     : Ending time for the computation
            kstages        : Type for the time integration(2 = RK2, 3 = RK3, 4 = RK4)
    
    Outputs:
    --------
    
            qexact         : Exact solution
            q              : Computed solution
            coord          : All grid points
            intma          : Intma(CG/DG)
    '''

    #Compute Interpolation and Integration Points
    xgl = Lobatto_p(Q)
    wgl = weight(N)
    
    if (integration_type == 1):
        noq = Q
    elif (integration_type == 2):
        noq = Q+1

    nq = noq
    xnq = Lobatto_p(nq)
    wnq = weight(Q)
    
    #Create Grid
    coord = grid_dg(Q, nel,xgl, ax, bx)
    # create intma
    intma = intma_cdg(N, nel, method_type)
    #space stuff
    dx = coord[1,0] - coord[0,0]
    #time stuff
    dt = Courant_max*dx/u
    ntime = int(time_final/dt)
    dt = time_final/ntime
    Courant = u*dt/dx

    # Compute Exact Solution
    time = 0
    qe = exact_solution_dg(coord,nel,Q,time, ax, bx, u, icase)

    # Periodic BC Pointers
    periodic = empty(Np)

    if method_type == 'cg':
        for i in range(Np):
            periodic[i] = i
        periodic[Np-1] = periodic[0]

    elif method_type == 'dg':
        for i in range(Np):
            periodic[i] = i
    
    # Lagrange basis and its derivatives
    l_basis, dl_basis = LagrangeBasis_deriv(N,Q,xgl, xnq)
    
    # Form Element Mass and Differentiation Matrices
    Me = Element_matrix(N,Q,wnq,l_basis)                             
    De = Element_Diff_matrix(N,Q,wnq,l_basis,dl_basis)               
    # Form Global Mass and Differentiation Matrices
    GMmatrix = DSS_operator(Me, nel, N, Np, intma, periodic, coord)
    GDmatrix = DSS_operator(u*De, nel, N, Np, intma, periodic, coord, 'diff')

    # Flux matrix for DG method
    Fmatrix = zeros((Np,Np))
    if method_type == 'cg':
        GMmatrix[Np-1,Np-1] = 1

    elif method_type == 'dg':
        Fmatrix = flux_matrix(nel, intma, N, Np, diss,u)
    
    # RHS matrix
    Rmatrix = linalg.solve(GMmatrix,(GDmatrix - Fmatrix))

    #Form Long Exact Solution Vector
    qexact = zeros(Np)
    for e in range(1,nel+1):
        for i in range(Q):
            ip = int(intma[i,e-1])
            qexact[ip] = qe[i,e-1]

    #Initialize State Vector
    qexact = transpose(qexact)
    q1 = qexact
    q = qexact
    qp = qexact

    #Time Integration
    for itime in range(ntime):
        time = time + dt

        for ik in range(1,kstages+1):
            if (kstages == 2):
                if(ik == 1):
                    a0 = 1
                    a1 = 0
                    beta = 1
                elif(ik == 2):
                    a0 = 0.5
                    a1 = 0.5
                    beta = 0.5
            elif(kstages == 3):
                if(ik == 1):
                    a0=1
                    a1=0
                    beta=1

                if(ik == 2):
                    a0=3/4
                    a1=1/4
                    beta=1/4 
                if(ik == 3):
                    a0 = 1/3
                    a1 = 2/3
                    beta = 2/3

            elif(kstages == 4):
                if(ik == 1):
                    a0=1
                    a1=0
                    beta=1/2

                if(ik == 2):
                    a0=0
                    a1=1
                    beta=1/2

                if(ik == 3):
                    a0 = 2/3
                    a1 = 1/3
                    beta = 1/6

                if (ik == 4):
                    a0 = 0
                    a1 = 1
                    beta = 1/2

            dtt = dt*beta

            qp = a0*q + a1*q1 + dtt*Rmatrix@qp

            #apply Periodic bc
            if (method_type == 'cg'):
                ip = periodic[Np-1]
                qp[Np-1] = qp[int(ip)] 

            #Update
            q1 = qp

        #Update Q
        q = qp

    #Compute Exact Solution
    qe = exact_solution_dg(coord,nel,Q,time, ax, bx, u, icase)

    #Form Long Exact Solution Vector
    for e in range(1,nel+1):
        for i in range(Q):
            ip = int(intma[i,e-1])
            qexact[ip] = qe[i,e-1]

    return qexact, q, coord, intma
