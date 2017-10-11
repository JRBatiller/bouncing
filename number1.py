# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 10:46:16 2017

@author: Joby
"""

"""So i just realized t must be an array"""

#solve ODE dxdt=x
import numpy as np
import scipy as sp
import matplotlib
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def dxdt(x,t):
    dxdt=x
    return dxdt

def euler_step(x, t, dt=0.01):
    #get x and t and time step size
    x2=x+dxdt(x,t)*dt
    t2=t+dt
    return [x2,t2] 

def solve_to(xstart,tstart,tend,deltat=0.01):
    # should return only last variable
    #must check that deltat will not overshoot tend
    #steps=int((tend-tstart)/deltat)
    """x_list=np.array([[xstart,tstart]])
    x_list.reshape(2,1)
    tnow=tstart
    while tnow<tend:
        #check size
        if tend-tnow<=deltat:
            deltat=tend-tnow
        [x, t]= x_list[-1]
        next_step=euler_step(x,t, deltat)
        x_list=np.vstack((x_list,next_step))
        tnow+=deltat
    return x_list
    """
    xnow=xstart
    tnow=tstart
    fsteps=(tend-tstart)/deltat
    steps=int((tend-tstart)/deltat)
    if(fsteps-steps>=0.000005):
        steps+=1
          
    deltat=(tend-tstart)/steps
    for step in range(steps):
        [xnow, tnow] =euler_step(xnow,tnow, deltat)
    return xnow
    
    """while tnow<tend:
        #check size
        if tend-tnow<=deltat:
            deltat=tend-tnow
        [xnow, tnow] =euler_step(xnow,tnow, deltat)
    return xnow"""
    
def solve_ode(x0, t, deltat=0.01):
    # t is a bloddy array!!!
    x_list=np.array([solve_to(x0,t0,t1) for t1 in t ])
    #having t0 here is irritating me somehow
    return x_list
    
x0=1
t0=0

t=np.arange(10)
x=solve_ode(x0,t)

plt.plot(t,x)

plt.plot(t, odeint(dxdt,x0,t))


plt.plot(t,x, deltat=0.00001)
    
    
    
    