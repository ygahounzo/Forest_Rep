
import numpy as np
import matplotlib.pyplot as plt
# set up nice tick marks for log data
def set_xticks(P):
    p0 = np.log2(P[0])
    p1 = np.log2(P[-1])
    plt.xlim([2**(p0-1), 2**(p1+1)])
    
    Pstr = (['{:d}'.format(int(p)) for p in P])
    plt.xticks(P,Pstr)
