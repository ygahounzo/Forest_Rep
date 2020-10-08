
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

def intma(e, N):
    t = (e-1)*N
    r = N*e
    intmm = []
    for s in range(t, r+1):
        intmm.append(s)
    
    return array(intmm)


#funtion that compute weight values based on quadrature rule
def weight(Q):
    xi = Lobatto_p(Q)
    w = zeros(Q+1)
    for i in range(Q+1):
        out1, out2 = Legendre_deriv(Q, xi[i])
        w[i] = 2/(Q*(Q+1)*(out1)**2)
        
    return w 

#weight(Q)

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

    return Me

#Element_matrix(N,Q)

def GMM(Ne, intma, N, Q):
    
    Np = Ne*N+1
    
    w = weight(Q)
    Me = Element_matrix(N,Q)
    M = zeros((Np, Np))
    
    for e in range(1,Ne+1):
        for j in range(N+1):
            
            J = intma(e, N)[j]
            
            for i in range(N+1):
                
                I = intma(e, N)[i]
                
                M[I,J] = M[I,J] + ((bx-ax)/Ne)*Me[i,j]
    return M

#GMM(Ne, intma, N, Q)

def Element_Diff_matrix(N):
    De = zeros((N+1, N+1))
    Xi = Lobatto_p(N)
    w = weight(N)

    for i in range(N+1):

        for j in range(N+1):
            for k in range(N+1):
                xi, a1 = LagrangeBasis(N, i, Xi[k], Xi)
                a2, xj = LagrangeBasis(N, j, Xi[k], Xi)
                De[i,j] = De[i,j]- w[k]*xi*xj

    return De


# function that compute global residual vector

def Resi(Ne, q, N,Miv):
    
    De = Element_Diff_matrix(N)
    
    Np = Ne*N+1
    
    fe = lambda q: u*q
    
    R = zeros(Np)                 # global residual vector
    
    for e in range(1, Ne+1):         # element loop

        # discretizing into element 

        Ie  = intma(e, N)
        
        # discretizing into element 
        
        qe = q[Ie]
               
        # residual for each element

        Re = zeros(N+1)

        for i in range(N+1):
        
            Re[i] = - dot(De[i],fe(qe))
            
            # compuataion of global residual vector

            I = Ie[i] #intma(e,N)[i]            

            R[I] = R[I] + Re[i]
    
    # reinitialisation of the global residual vector using inverse mass matrix
    
    GR = zeros(Np) 
    
    for I in range(Np):
        
        GR[I] = dot(Miv[I],R)
        
    return GR

def Solver_1DW(N,Ne,M, x,t,dt):
    
    
    #q0 = qinit(x)         # defined initial condition
    q0 = array([qinit(l) for l in x])
    # Boundary conditions

    q0[0] = qinit(ax) 
    q0[-1] = qinit(bx) 

    
    #inverse of global mass matrix
    
    MM = GMM(Ne, intma, N, Q)
    
    Miv = linalg.inv(MM)     

    # computation of the solution of 1D wave equation
    qn = q0     

    for n in range(M):                   # time loop

        K1 = Resi(Ne, qn,N, Miv)

        # soultion for the wave equation at time n+1

        qh = qn + (dt/2)*K1

        K2 = Resi(Ne,qh, N,Miv)
        
        Ph = qn + (dt/2)*K2
        
        K3 = Resi(Ne, Ph, N, Miv)
        
        P = qn + dt*K3
        
        K4 = Resi(Ne, P, N, Miv)
        

        qn1 = qn + (dt/6)*(K1+2*K2+2*K3+K4) 

        qn1[0] = qinit(ax-u*t[n+1])
        qn1[-1] = qinit(bx-u*t[n+1])
        
        qn = qn1
    
    return qn1
