
from numpy import *

def Lobatto_deriv(Q, x):
    F = [] # array containing the polynomials
    dF = []
    ddF = []


    f0 = 1; F.append(f0)  # first polynomial
    df0 = 0; dF.append(df0)
    f1 = x; F.append(f1)  # second polynomial
    df1 = 1; dF.append(df1)
    ddF = [0,0]

    B = []
    dB = []


    for i in range(2, Q+1):
        fi = ((2*i-1)/i)*x*F[i-1] - ((i-1)/i)*F[i-2]  # iteration of the polynomials
        dfi = i*F[i-1] + x*dF[i-1]                    # first derivative
        F.append(fi)
        dF.append(dfi)
        ddfi = (i+1.0)*dF[i-1] + x*ddF[i-1]           # second derivative
        ddF.append(ddfi)

        bi = (1.0-x**2)*dF[i-1]                       # lobatto polynomial
        dbi = -2.0*x*dF[i-1] + (1.0-x**2)*ddF[i-1]    # derivative of lobatto polynomial
        B.append(bi)
        dB.append(dbi)
        
    return B[-1], dB[-1]


def Legendre_deriv(Q, x):
    L = [] # array containing the polynomials
    dL = []


    f0 = 1; L.append(f0)  # first polynomial
    df0 = 0; dL.append(df0)
    f1 = x; L.append(f1)  # second polynomial
    df1 = 1; dL.append(df1)
    

    for i in range(2, Q+1):
        fi = ((2*i-1)/i)*x*L[i-1] - ((i-1)/i)*L[i-2]  # iteration of the polynomials
        dfi = i*L[i-1] + x*dL[i-1]
        L.append(fi)
        dL.append(dfi)
        
    return L[-1], dL[-1]

def Lobatto_p(Q):
    X = []  # Array that contains legendre points
    K = 100  # Order of approximation of Newton method
    e = 10**(-20) #error
    for i in range(0, Q+1):
        xi0 = cos(((2*i+1)/(2*Q+2))*pi)   # Chebchev points
        

        xik = xi0

        for k in range(0, K):
            out1, out2 = Lobatto_deriv(Q+1, xik)
            xikk = xik - out1/out2   # approximation of the solution using Newton

            if abs(xikk-xik) < e:

                break

            xik = xikk

        X.append(xikk)
    return X


def LagrangeBasis(N, i, xl, Xr):
    
    L = 1
    dL = 0
        
    for j in range(N+1):
            
        prod = 1
        
        if (j != i):
            L = L*(xl-Xr[j])/(Xr[i]-Xr[j])
                
            for k in range(N+1):
                if (k!=i  and k!=j):
                    prod = prod*(xl-Xr[k])/(Xr[i]-Xr[k])
        
            dL = dL+prod/(Xr[i]-Xr[j])
    return L, dL

def intmaDG(e, N):
    
    t = (e-1)*N
    r = e*N
    
    intmm = []
    for s in range(t, r+1):
        it = e-1+s
        intmm.append(it)
        
    return array(intmm)


#funtion that compute weight values based on quadrature rule
def weight(Q):
    xi = Lobatto_p(Q)
    w = zeros(Q+1)
    for i in range(Q+1):
        out1, out2 = Legendre_deriv(Q, xi[i])
        w[i] = 2/(Q*(Q+1)*(out1)**2)
        
    return w 


# Element mass matrix

def Element_matrix(N,Q):
    Me = zeros((N+1, N+1))       # initialisation of the matrix
    Xr = Lobatto_p(N)               # roots
    Xi = Lobatto_p(Q)               # Lobatto points
    w = weight(Q)                   # weight values

    for i in range(N+1):

        for j in range(N+1):
            for k in range(Q+1):
                xi, a1 = LagrangeBasis(N, i, Xi[k], Xr)
                xj, a2 = LagrangeBasis(N, j, Xi[k], Xr)
                Me[i,j] = Me[i,j]+ w[k]*xi*xj



    Me = (1/2)*Me

    print('Elemental Mass Matrix')
    print(Me) #print Mass Matrix
    
    return Me

#Differentiation element matrix

def Element_Diff_matrix(N,Q):
    De = zeros((N+1, N+1))
    Xi = Lobatto_p(N)
    w = weight(N)

    for i in range(N+1):

        for j in range(N+1):
            for k in range(Q):
                xi, a1 = LagrangeBasis(N, i, Xi[k], Xi)
                a2, xj = LagrangeBasis(N, j, Xi[k], Xi)
                De[i,j] = De[i,j]- w[k]*a1*a2

    print('Elemental Differentiation Matrix')
    print(De) #print Diff Matrix
                
    return De

#Element flux matrix

def Element_Flux_matrix(N, Q):
    Fe = zeros((N+1, N+1))
    Xi = Lobatto_p(N)[::-1]

    for i in range(N+1):

        for j in range(Q):
                x1, a1 = LagrangeBasis(N, i, Xi[-1], Xi)
                x2, a2 = LagrangeBasis(N, j, Xi[-1], Xi)
                x3, a3 = LagrangeBasis(N, i, Xi[0], Xi)
                x4, a4 = LagrangeBasis(N, j, Xi[0], Xi)
                
                Fe[i,j] = x1*x2 - x3*x4

    print('Elemental Flux Matrix')
    print(Fe) #print Flux Matrix
                
    return Fe


# DSS operator

def DSS(A, Ne, intmaDG, N, Q):
    
    Np = Ne*(N+1)
    
    M = zeros((Np, Np))
    
    for e in range(1,Ne+1):
        for j in range(N+1):
            
            J = intmaDG(e, N)[j]
            
            for i in range(N+1):
                
                I = intmaDG(e, N)[i]
                
                M[I,J] = M[I,J] + A[i,j]
    return M


# Global flux vector

def GNFV(Ne, intmaDG, N, diss,u,f,q):
    
    Np = Ne*(N+1)

    fstar = zeros(Np)
    
    for e in range(1,Ne+1):
        L = e
        R = e+1
        if e == Ne:
            R = 1
            
        I = intmaDG(L,N)[N]
        J = intmaDG(R,N)[0]
        
        fs = (1/2)*(f[I]+f[J]-u*diss*(q[J]-q[I]))
        fstar[J] = fs
        fstar[I] = fs
        
    return fstar


# function that compute global residual vector

def Resi(N,Ne, De, Fe, f,q,Miv,u,diss,fstar):
    
    Np = Ne*(N+1)
    
    #fe = lambda q: u*q
    
    R = zeros(Np)                 # global residual vector
    
    for e in range(1, Ne+1):         # element loop

        # discretizing into element 

        Ie  = intmaDG(e, N)
        
        # discretizing into element 
        
        qe = q[Ie]
               
        # residual for each element

        Re = zeros(N+1)
        fe = f(qe)
        
        
        #fstar = ENFV(e, intmaDG, N, diss,u,fe,qe)
        fse = fstar[Ie]
        
        for i in range(N+1):
        
            Re[i] = dot(De[i],f(qe))
            Re[i] -= dot(Fe[i],fse)
            
            # compuataion of global residual vector

            I = Ie[i] #intma(e,N)[i]            

            R[I] = R[I] + Re[i]
                
                
            
    # reinitialisation of the global residual vector using inverse mass matrix
    
    GR = zeros(Np) 
    
    for I in range(Np):
        
        GR[I] = dot(Miv[I],R)
        
    return GR



def SolverDG(Np,N,Ne,u,diss,f,Nt, x,t,dt,bc,qinit,ax,bx):
    
    Q = N+1
    # Boundary conditions
    q0 = qinit(x)
    q0[0] = qinit(ax) 
    q0[-1] = qinit(bx)     
    
    # Global mass matrix

    Me = ((bx-ax)/Ne)*Element_matrix(N,Q)

    MG = DSS(Me, Ne, intmaDG, N, Q)
    Miv = linalg.inv(MG)
    
    # Element differentiation matrix

    De = Element_Diff_matrix(N,Q)

    # Element flux matrix

    Fe = Element_Flux_matrix(N,Q)
    
    # computation of the solution of 1D wave equation
    qn = q0 
        
    for n in range(Nt):                   # time loop
        ff = f(qn)
        fstar = GNFV(Ne, intmaDG, N, diss,u,ff,qn)

        K1 = Resi(N,Ne, De, Fe, f,qn,Miv,u,diss,fstar)

        # soultion for the wave equation at time n+1

        qh = qn + (dt/2)*K1
        
        ff1 = f(qh)
        fstar1 = GNFV(Ne, intmaDG, N, diss,u,ff1,qh)
        
        K2 = Resi(N,Ne, De, Fe, f,qh,Miv,u,diss,fstar1)
        
        Ph = qn + (dt/2)*K2
        
        ff2 = f(Ph)
        fstar2 = GNFV(Ne, intmaDG, N, diss,u,ff2,Ph)
        
        K3 = Resi(N,Ne, De, Fe, f,Ph,Miv,u,diss,fstar2)
        
        P = qn + dt*K3
        
        ff3 = f(P)
        fstar3 = GNFV(Ne, intmaDG, N, diss,u,ff3,P)
        
        K4 = Resi(N,Ne, De, Fe, f,P,Miv,u,diss,fstar3)
        

        qn1 = qn + (dt/6)*(K1+2*K2+2*K3+K4) 

        qn1[0], qn1[-1] = bc(t[n+1])
        
        qn = qn1
    
    return qn1
