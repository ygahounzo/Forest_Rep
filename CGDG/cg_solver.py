#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 10 18:34:23 2021

@author: mathuser
"""


import matplotlib.pyplot as plt
from numpy import *
import pandas as pd
from time import perf_counter


# import CG/DG functions
from cg_dg_functions import *

# Exact solution
def initial_condition(x, case):
    
    '''
    This function compute the initial condition
    
    Inputs:
    -------
            x: variable
            case: 1 for Gaussain, 2 for sinusoidal function
            
    Output:
    -------
           return the initial condition according to the case
    '''
    
    if (case == 1):
        return exp(-128*x**2)
    elif(case == 2):
        return sin(2*pi*x)

    
#time stuff
def time_step(dx, ax, bx, CFL, nel0, N, u, time_final, fixed_dt):
    
    if (fixed_dt == True):
        dx0 = (bx-ax)/(nel0*N)
        dtest = CFL*dx0/abs(u)
        ntime0 = int(time_final//dtest) + 1   # Number of time steps
        dt0 = time_final/ntime0

        f = log2(dx0/dx)
        dt = dt0/2**f
        ntime = int(ntime0*2**f)
        #print('ntime = {:d}'.format(ntime))
    else:        
        dtest = CFL*dx/abs(u)
        ntime = int(time_final/dtest) + 1    # Number of time steps
        dt = time_final/ntime
        #print('ntime = {:d}'.format(ntime))
        
    return dt, ntime


## Call of CG/DG solver

#order = array([1,2])         # polynomial order
#N_element = array([32,64])  # number of elements
#N_element = array([32,64,128,256,512,1024])

kstages = 4               # RK2, RK3, RK4   explicit methods
ti_method = 4             # IRK1, IRK2, IRK3, IRK4, IRK5: implicit methods
dt = 1e-2                 # time-step, fraction of one revolution
CFL = 0.25                # CFL number
time_final = 1.0          # final time in revolutions
integration_type = 2      # % = 1 is inexact and = 2 is exact
time_integration_type = 'implicit'       # explicit or implicit 
method_type = 'cg'        # CG or DG
iplot = False             # plot the solution
icase = 2                 # case number: 1 is a Gaussian and 2 is a sinusoidal
diss = 0 
u = 1.0
ax = 0
bx = 1.0
#ntime = time_final/dt;

#len_el = len(N_element)
#len_pol = len(order)
#l2e_norm = zeros((len_pol, len_el))
#Np_array = zeros((len_pol, len_el))

#nel0 = N_element[0]
#Nv = N_element

# Create a data type for storing results
dt_data = dtype([('N',int),('wall','d'), ('setup','d'),('integ','d'),('M',int),('dt','d'),('cfl','d'),
                  ('1-norm','d'),('2-norm','d'),('inf-norm','d')])  

def dg_simulation(file,iN, N,integration_type,Nv):
    
    #for iN,N in enumerate(order):
    #N = order[iN]
    if (integration_type == 1):
        Q = N
    elif (integration_type == 2):
        Q = N+1

    wall = 0

    for e, nel in enumerate(Nv):
        #nel = N_element[e]

        if (method_type == 'cg'):
            Np = nel*N + 1
        elif (method_type == 'dg'):
            Np = nel*(N+1)

        # Call of 1D wave solver
        '''
        outputs:
        --------
        qexact         : Exact solution
        q              : Computed solution
        coord          : All grid points
        intma          : Intma(CG/DG)
        '''
        # wall time 
        tic = perf_counter()

        qexact, q, coord, intma, t_list, dt, ntime = cg_dgSolver(N, Q, nel, Np, ax, bx, integration_type,\
          method_type,icase,diss,u, CFL, time_final, kstages, initial_condition, ti_method, time_integration_type)

        toc = perf_counter()
        wall += toc - tic

        # Compute L2- norm
        num = 0
        denom = 0
        error = zeros(Np)
        for i in range(Np):
            num = num + (q[i]-qexact[i])**2
            error[i] = abs(q[i]-qexact[i])
            denom = denom + (qexact[i])**2

        e1 = sum(error)
        e3 = max(error)
        e2 = sqrt(num/denom)


        w = wall
        s = t_list[0]   # setup time
        I = t_list[1]   # integration time
        t = array((nel, w,s,I,ntime,dt,CFL,e1,e2,e3),dtype=dt_data)
        file.write(t)
        
        #Compute a gridpoint solution
        x_sol = zeros(Np);
        for ie in range(1,nel+1):
            for i in range(N+1):
                ip = int(intma[i,ie-1])
                x_sol[ip] = coord[i,ie-1]

    if(iplot == True):
        figure(iN)
        plot(x_sol, qexact, label = 'Exact')
        plot(x_sol, q, '--', label = 'Computed')
        legend()



# CG Simulation
output_file = 'cg_adv_data.dat_RK2'

order = [1,2]
integType = [1]

Nv = [8,16,32,64,128,256,512,1024]
        

# Simulation. Loop over integration type, and then loop over polynomial order. 
file = open(output_file,"wb")
file.write(array([len(Nv)]))

for integration_type in integType:
    
    for iN,N in enumerate(order):
        
        # Run over range of N = 32,64,128,256,...
        file.write(array((N,integration_type)))

        print("order = {:d}; Integration_type = {:d}".format(N,integration_type))
        dg_simulation(file,iN, N,integration_type,Nv)
        #print("")
        
file.close()

print("End")
#fout = open(output_file,"rb")
