#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 10 18:57:31 2021

@author: mathuser
"""

from numpy import *
from time import perf_counter
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import inv
from scipy.sparse.linalg import spsolve


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
    for i in range(Q+1):
        xik = cos(((2*i+1)/(2*(Q+1)-1))*pi)         # Chebchev points

        for k in range(K):
            out1, out2 = Lobatto_deriv(Q+1, xik)
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
    
    xi = Lobatto_p(Q)
    w = zeros(Q+1)
    for i in range(Q+1):
        
        out1, out2, out3 = Legendre_deriv(Q, xi[i])
        w[i] = 2/(Q*(Q+1)*(out1)**2)
        
    return w 

# grid points
def grid_dg(N,Ne, xe, ax, bx):
    
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
    
    grid = zeros((N+1,Ne))

    xel = linspace(ax,bx,Ne+1)

    for e in range(1,Ne+1):
        
        ae = xel[e-1] ; be = xel[e]

        xsi = ((be-ae)/2)*(xe-1) + be
        
        for i in range(N+1):

            grid[i,e-1] = xsi[i]
            
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
    
    Me = zeros((N+1, N+1))       # initialisation of the matrix
    
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
def Element_Diff_matrix(N, Q, wght, l_basis, dl_basis):
    
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
    
    De = zeros((N+1, N+1))
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
            
            for i in range(N+1):
                
                I = int(intma[i, e-1])
                I = int(periodic[I])
                
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
def exact_solution(coord,Ne,N,time, ax, bx, u, case, initial_condition):
    
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
    qe = zeros((N+1,Ne))

    timec = time - floor(time)

    #Generate Grid Points
    for e in range(1, Ne + 1):
        for i in range(N+1):
            x = coord[i,e - 1]
            xbar = u*timec
            if (xbar > bx): 
                xbar = ax + (xbar - bx)
            
            qe[i,e - 1] = initial_condition(x - xbar, case)
                
    return qe

# Explicit RK methods
def RKstages(ik, kstages):
    
    if(kstages == 1):
        a0 = 1; a1 = 0; beta = 1
    
    elif(kstages == 2):
        if(ik == 1):
            a0 = 1; a1 = 0;  beta = 1
        elif(ik == 2):
            a0 = 0.5; a1 = 0.5; beta = 0.5
    elif(kstages == 3):
        if(ik == 1):
            a0 = 1; a1 = 0; beta = 1

        if(ik == 2):
            a0 = 3/4; a1 = 1/4; beta = 1/4 
        if(ik == 3):
            
            a0 = 1/3; a1 = 2/3; beta = 2/3

    elif(kstages == 4):
        if(ik == 1):
            a0=1; a1=0; beta=1/2

        if(ik == 2):
            a0=0; a1 = 1; beta = 1/2

        if(ik == 3):
            a0 = 2/3; a1 = 1/3; beta = 1/6

        if (ik == 4):
            a0 = 0; a1 = 1; beta = 1/2
            
    return a0, a1, beta

# Implicit RK methods

def IRK_coefficients(ti_method):
    
  
    if(ti_method == 1):
        stages = 2
        alpha = zeros((stages,stages))
        beta = zeros(stages)
        alpha[1,0] = 0
        alpha[1,1] = 1
        beta[:] = alpha[stages-1,:]
    elif(ti_method == 2):
        stages = 3
        alpha = zeros((stages,stages))
        beta = zeros(stages)
        alpha[1,0] = 1 - 1/sqrt(2)
        alpha[1,1] = alpha[1,0]
        alpha[2,0] = 1/(2*sqrt(2))
        alpha[2,1] = alpha[2,0]
        alpha[2,2] = alpha[1,1]
        beta[:]=alpha[stages-1,:]
    elif(ti_method == 3):
        stages = 4
        alpha = zeros((stages,stages))
        beta = zeros(stages)
        alpha[1,0] = 1767732205903.0/4055673282236.0
        alpha[1,1] = 1767732205903.0/4055673282236.0
        alpha[2,0] = 2746238789719.0/10658868560708.0
        alpha[2,1] = -640167445237.0/6845629431997.0
        alpha[2,2] = alpha[1,1]
        alpha[3,0] = 1471266399579.0/7840856788654.0
        alpha[3,1] = -4482444167858.0/7529755066697.0
        alpha[3,2] = 11266239266428.0/11593286722821.0
        alpha[3,3] = alpha[1,1]
        beta[:] = alpha[stages-1,:]
    elif(ti_method == 4):
        stages = 6
        alpha = zeros((stages,stages))
        beta = zeros(stages)
        alpha[1,0] = 1.0/4.0
        alpha[1,1] = 1.0/4.0
        alpha[2,0] = 8611.0/62500.0
        alpha[2,1] = -1743.0/31250.0
        alpha[2,2] = alpha[1,1]
        alpha[3,0] = 5012029.0/34652500.0
        alpha[3,1] = -654441.0/2922500.0
        alpha[3,2] = 174375.0/388108.0
        alpha[3,3] = alpha[1,1]
        alpha[4,0] = 15267082809.0/155376265600.0
        alpha[4,1] = -71443401.0/120774400.0
        alpha[4,2] = 730878875.0/902184768.0
        alpha[4,3] = 2285395.0/8070912.0
        alpha[4,4] = alpha[1,1]
        alpha[5,0] = 82889.0/524892.0
        alpha[5,1] = 0.0
        alpha[5,2] = 15625.0/83664.0
        alpha[5,3] = 69875.0/102672.0
        alpha[5,4] = -2260.0/8211.0
        alpha[5,5] = alpha[1,1]
        beta[:] = alpha[stages-1,:]
    elif(ti_method == 5):
        stages = 8
        alpha = zeros((stages,stages))
        beta = zeros(stages)
        alpha[1,0] = 41.0/200.0
        alpha[1,1] = 41.0/200.0
        alpha[2,0] = 41.0/400.0
        alpha[2,1] = -567603406766.0/11931857230679.0
        alpha[2,2] = alpha[1,1]
        alpha[3,0] = 683785636431.0/9252920307686.0
        alpha[3,1] = 0.0
        alpha[3,2] = -110385047103.0/1367015193373.0
        alpha[3,3] = alpha[1,1]
        alpha[4,0] =  3016520224154.0/10081342136671.0
        alpha[4,1] = 0.0
        alpha[4,2] = 30586259806659.0/12414158314087.0
        alpha[4,3] = -22760509404356.0/11113319521817.0
        alpha[4,4] = alpha[1,1]
        alpha[5,0] = 218866479029.0/1489978393911.0
        alpha[5,1] = 0.0
        alpha[5,2] = 638256894668.0/5436446318841.0
        alpha[5,3] = -1179710474555.0/5321154724896.0
        alpha[5,4] = -60928119172.0/8023461067671.0
        alpha[5,5] = alpha[1,1]
        alpha[6,0] = 1020004230633.0/5715676835656.0
        alpha[6,1] = 0.0
        alpha[6,2] = 25762820946817.0/25263940353407.0
        alpha[6,3] = -2161375909145.0/9755907335909.0
        alpha[6,4] = -211217309593.0/5846859502534.0
        alpha[6,5] = -4269925059573.0/7827059040749.0
        alpha[6,6] = alpha[2,2]
        alpha[7,0] = -872700587467.0/9133579230613.0
        alpha[7,1] = 0.0
        alpha[7,2] = 0.0
        alpha[7,3] = 22348218063261.0/9555858737531.0
        alpha[7,4] = -1143369518992.0/8141816002931.0
        alpha[7,5] = - 39379526789629.0/19018526304540.0
        alpha[7,6] = 32727382324388.0/42900044865799.0
        alpha[7,7] = alpha[1,1]
        beta[:] = alpha[stages-1,:]
    elif(ti_method == 6):
        stages = 10
        alpha = zeros((stages,stages))
        beta = zeros(stages)
        alpha[1,1] = 0.2928932188134525
        alpha[2,1] = 0.3812815664617709
        alpha[3,1] = 0.4696699141100893
        alpha[4,1] = 0.5580582617584078
        alpha[5,1] = 0.6464466094067263
        alpha[6,1] = 0.7348349570550446
        alpha[7,1] = 0.8232233047033631
        alpha[8,1] = 0.9116116523516815
        alpha[9,1] = 0.7071067811865476
        alpha[9,9] = 0.2928932188134525
        beta[:] = alpha[stages-1,:];

    return alpha, beta, stages
                    
            


def cg_dgSolver(N,Q,nel, Np, ax, bx, integration_type, method_type, icase,diss,u,CFL, time_final,\
                kstages, initial_condition,ti_method,time_integration_type):
    
    '''
    This function is CG/DG solver for 1D wave equqtion
    
    Inputs:
    -------
            N              : Polynomial order
            Q              : Integration order(N+1: for exact, N: for inexact integration)
            nel            : Number of element
            nel0           : The first number of element 
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
            time_step      : function that compute the time step and number of time( double time per element)
    Outputs:
    --------
    
            qexact         : Exact solution
            q              : Computed solution
            coord          : All grid points
            intma          : Intma(CG/DG)
    '''

    #Compute Interpolation and Integration Points
    time_list = []
    
    xgl = Lobatto_p(N)
    wgl = weight(N)
    
    xnq = Lobatto_p(Q)
    wnq = weight(Q)
    
    # preparation time 
    t0 = perf_counter()
    
    # create intma
    intma = intma_cdg(N, nel, method_type)
    
    #Create Grid and space stuff
    
    coord = grid_dg(N, nel,xgl, ax, bx)
    
    dx = coord[1,0] - coord[0,0]
    #time stuff
    
    #dt = CFL*dx/u
    #ntime = int(floor(time_final/dt)) + 1
    #dt = time_final/ntime
    
    #dx = (bx-ax)/Np
    dt_est = CFL*dx/u;
    ntime = int(floor(time_final/dt_est))
    dt = time_final/ntime
    
    # Compute Exact Solution
    t_time = 0
    qe = exact_solution(coord,nel,N,t_time, ax, bx, u, icase, initial_condition)

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
        
    if(time_integration_type == 'explicit'):
        
        # RHS matrix
        Rmatrix = linalg.solve(GMmatrix,(GDmatrix - Fmatrix))
    
        Rmatrix = csr_matrix(Rmatrix)
        
    elif(time_integration_type == 'implicit'):
        
        DFmatrix = Fmatrix - GDmatrix
    
        GMmatrix = csr_matrix(GMmatrix)
        DFmatrix = csr_matrix(DFmatrix)
        
        # RK coefficients for implicit integration
        alpha,beta, stages = IRK_coefficients(ti_method)
        
    t1 = perf_counter()
    
        
    time_list.append(t1-t0)
    
    
    #Form Long Exact Solution Vector
    qexact = zeros(Np)
    for e in range(1,nel+1):
        for i in range(N+1):
            ip = int(intma[i,e-1])
            qexact[ip] = qe[i,e-1]

    #Initialize State Vector
    #qexact = transpose(qexact)
    # integration time
    t2 = perf_counter()
        
    # Time integration
    if(time_integration_type == 'explicit'):
        q1 = qexact
        q = qexact
        qp = qexact
        
        
        
        
        #Time Integration
        for itime in range(ntime):
            t_time = t_time + dt
    
            for ik in range(1,kstages+1):
                
                a0, a1, beta = RKstages(ik, kstages)
    
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
    elif(time_integration_type == 'implicit'):
        
        q = qexact
        
        Q = zeros((Np,stages))
        QQ = zeros((Np,stages))
        R = zeros((Np,stages))
        
        #ntime = 5
        for itime in range(ntime):
            
            t_time = t_time + dt
                
            Q[:,0] = q[:]
            R[:,0] = DFmatrix@Q[:,0]
            
            for i in range(1,stages):
                R_sum = zeros(Np)
                for j in range(i):
                    R_sum = R_sum + alpha[i,j]*R[:,j]
                    #print('aij = ',alpha[i,j])
                RR = GMmatrix@q + dt*R_sum
                
                Amatrix = GMmatrix - dt*alpha[i,i]*DFmatrix
                #print('aii = ',alpha[i,i])
                
                Q[:,i] = spsolve(Amatrix, RR)
                R[:,i] = DFmatrix@Q[:,i]
                
            R_sum = zeros(Np)
            for i in range(stages):
                R_sum = R_sum + beta[i]*R[:,i]
                #print('bi = ',beta[i])
            
            qp = q + dt*spsolve(GMmatrix,R_sum)
            
            if (method_type == 'cg'):
                ip = periodic[Np-1]
                qp[Np-1] = qp[int(ip)]
            
            q = qp
    
    t3 = perf_counter()
    
    time_list.append(t3-t2)
    
    #Compute Exact Solution
    qe = exact_solution(coord,nel,N,t_time, ax, bx, u, icase, initial_condition)
    
    #Form Long Exact Solution Vector
    for e in range(1,nel+1):
        for i in range(N+1):
            ip = int(intma[i,e-1])
            qexact[ip] = qe[i,e-1]

    return qexact, q, coord, intma, array(time_list), dt, ntime