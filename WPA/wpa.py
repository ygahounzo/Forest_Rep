
from numpy import *

def limiter(r,lim_choice='MC'):
    wlimiter = empty(r.shape)
    if lim_choice == None:
        wlimiter = ones(r.shape)
    elif lim_choice == 'minmod':        
        # wlimitr = dmax1(0.d0, dmin1(1.d0, r))
        wlimiter = maximum(0,minimum(1,r))
    elif lim_choice == 'superbee':
        # wlimitr = dmax1(0.d0, dmin1(1.d0, 2.d0*r), dmin1(2.d0, r))
        a1 = minimum(1,2*r)
        a2 = minimum(2,r)        
        wlimiter = maximum(0,maximum(a1,a2))
    if lim_choice == 'MC':
        # c = (1.d0 + r)/2.d0
        # wlimitr = dmax1(0.d0, dmin1(c, 2.d0, 2.d0*r))
        c = (1 + r)/2
        wlimiter = maximum(0,minimum(c,minimum(2,2*r)))
    elif lim_choice == 'vanleer':
        # wlimitr = (r + dabs(r)) / (1.d0 + dabs(r))
        wlimiter = (r + abs(r))/(1 + abs(r))
            
    return wlimiter


# ----------------------------------------------------
# Wave propagation algorithm
# 
# Adapted from the Clawpack package (www.clawpack.org)
# -----------------------------------------------------

def claw1(ax, bx, mx, Tfinal, nout, 
          meqn=1, \
          rp=None, \
          qinit=None, \
          bc=None, \
          limiter_choice='MC', \
          second_order=True):

    dx = (bx-ax)/mx
    xe = linspace(ax,bx,mx+1)  # Edge locations
    xc = xe[:-1] + dx/2       # Cell-center locations

    # For many problems we can assume mwaves=meqn.
    mwaves = meqn
        
    # Temporal mesh
    t0 = 0
    tvec = linspace(t0,Tfinal,nout+1)
    dt = Tfinal/nout
    
    assert rp is not None,    'No user supplied Riemann solver'
    assert qinit is not None, 'No user supplied initialization routine'
    assert bc is not None,    'No user supplied boundary conditions'

    # Initial the solution
    q0 = qinit(xc,meqn)    # Should be [size(xc), meqn]
    
    # Store time stolutions
    Q = empty((mx,meqn,nout+1))    # Doesn't include ghost cells
    Q[:,:,0] = q0

    q = q0
    dtdx = dt/dx
    t = t0
    for n in range(0,nout):
        t = tvec[n]
        
        # Add 2 ghost cells at each end of the domain;  
        q_ext = bc(q)

        # Get waves, speeds and fluctuations
        waves, speeds, amdq, apdq = rp(q_ext)
                
        # First order update
        q = q - dtdx*(apdq[1:-2,:] + amdq[2:-1,:])
        
        # Second order corrections (with possible limiter)
        if second_order:    
            cxx = zeros((q.shape[0]+1,meqn))  # Second order corrections defined at interfaces
            for p in range(mwaves):
                sp = speeds[p][1:-1]   # Remove unneeded ghost cell values added by Riemann solver.
                wavep = waves[p]
                if limiter_choice is not None:
                    wl = sum(wavep[:-2,:] *wavep[1:-1,:],axis=1)  # Sum along dim=1
                    wr = sum(wavep[2:,:]  *wavep[1:-1,:],axis=1)
                    w2 = sum(wavep[1:-1,:]*wavep[1:-1,:],axis=1)
                
                    # Create mask to avoid dividing by 0. 
                    m = w2 > 0
                    r = ones(w2.shape)
                    
                    r[m] = where(sp[m,0] > 0,  wl[m]/w2[m],wr[m]/w2[m])

                    wlimiter = limiter(r,lim_choice=limiter_choice)
                    wlimiter.shape = sp.shape
                    
                else:
                    wlimiter = 1
                    
                cxx += 0.5*abs(sp)*(1 - abs(sp)*dtdx)*wavep[1:-1,:]*wlimiter
                
            Fp = cxx  # Second order derivatives
            Fm = cxx
        
            # update with second order corrections
            q -= dtdx*(Fp[1:,:] - Fm[:-1,:])
        
        Q[:,:,n+1] = q
        
    return Q, xc, tvec

