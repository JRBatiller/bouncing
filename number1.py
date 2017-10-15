# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 10:46:16 2017

@author: Joby
"""

"""So i just realized t must be an array"""

"""
Use the Euler method to solve the ODE 



Make a function called solve_ode that uses solve_to and euler_step to generate a series of numerical solution estimates 
This should be similar to scipy’s odeint function.
"""



#solve ODE dxdt=x
import numpy as np
import scipy as sp

import matplotlib
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy import stats

def dxdt(x,t):
    dxdt=x
    return dxdt


#Ensure that you have a function called euler_step that does a single Euler step.
def euler_step(x1, t, dt=0.01):
    #get x and t and time step size
    x2=x1+dxdt(x1,t)*dt
    return x2 

#Also make a function called solve_to which solves from x1,t1 to x2,t2 in steps no bigger than deltat_max.
def solve_to(xstart,tstart,tend,deltat=0.01):
    # should return only last variable
    #steps=int((tend-tstart)/deltat)
    xnow=xstart
    tnow=tstart
    fsteps=(tend-tstart)/deltat
    steps=int((tend-tstart)/deltat)
    if(fsteps-steps>=0.000005):
        steps+=1
          
    deltat=(tend-tstart)/steps
    for step in range(steps):
        xnow =euler_step(xnow,tnow, deltat)
        tnow += deltat
    return xnow
    
    
    
def solve_ode(x0, t, deltat=0.01):
    # t is a bloddy array!!!
    x_list=np.array([solve_to(x0,t0,t1, deltat) for t1 in t ])
    #having t0 here is irritating me somehow
    #this will become really slow
    return x_list
    
x0=1
t0=0
"""
You should use this to estimate x(1) using different timesteps. 
Produce a (nicely formatted) plot with double logarithmic scale showing how the error depends on the size of the timestep 
"""

t=np.array([1])
#deltat_list=[(10**-(z+1)) for z in range(10)]
array1=np.array( [10**-(z+1) for z in range(7)])
array2=array1/2
deltat_array=np.append(array1,array2)
deltat_array.sort()

#deltat_array=np.array([0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001, 0.000005, 0.000001)
 
x1=np.array([])
for z in range(len(deltat_array)):
    newx1=solve_ode(x0,t, deltat_array[z])
    x1=np.append(x1,newx1)
    print("delta t is now {}".format(deltat_array[z]))

err=np.exp(1)-x1

plt.loglog(deltat_array,err)




#xode=odeint(dxdt,x0,t)




    
    