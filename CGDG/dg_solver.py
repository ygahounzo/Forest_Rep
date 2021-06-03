#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 10 18:57:31 2021

@author: mathuser
"""

import matplotlib.pyplot as plt
from numpy import *
import pandas as pd
from time import perf_counter

#import sys
#sys.path.append("/Users/mathuser/Documents/Code/Forest_Rep/CG")

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


## Call of CG/DG solver

order = array([1,2])         # polynomial order
N_element = array([8,16,32,64,128,256,512,1024,2048])

kstages = 4               # RK2, RK3, RK4
dt = 1e-2                 # time-step, fraction of one revolution
CFL = 0.4                 # CFL number
time_final = 1.0          # final time in revolutions
integration_type = 2      # % = 1 is inexact and = 2 is exact
method_type = 'dg'        # CG or DG
iplot = False             # plot the solution
icase = 2                 # case number: 1 is a Gaussian and 2 is a sinusoidal
diss = 0 
u = 1
ax = 0
bx = 1
#ntime = time_final/dt;

Nv = N_element

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
                                method_type,icase,diss,u, CFL, time_final, kstages, initial_condition)

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
        I = t_list[1]   # integration time; 
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



# DG Simulation
output_file = 'dg_adv_data.dat'

order = [1,2,3,4,6]
integType = [1,2]

Nv = [8,16,32,64,128,256,512,1024]
        

# Simulation. Loop over integration type, and then loop over polynomial order. 
file = open(output_file,"wb")
file.write(array([len(Nv)]))

for integration_type in integType:
    
    for iN,N in enumerate(order):
        
        # Run over range of N = 8, 16, 32,64,128,256,...
        file.write(array((N,integration_type)))

        print("order = {:d}; Integration_type = {:d}".format(N,integration_type))
        dg_simulation(file,iN, N,integration_type,Nv)
        #print("")
        
file.close()
print("End")
#fout = open(output_file,"rb")
