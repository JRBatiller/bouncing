# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 10:46:16 2017

@author: Joby
"""

"""So i just realized t must be an array"""

#solve ODE dxdt=x
import numpy as np
import scipy as sp

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
    
    while tnow<tend:
        #check size
        if tend-tnow<=deltat:
            deltat=tend-tnow
        [xnow, tnow] =euler_step(xnow,tnow, deltat)
    return [xnow, tnow]
    
def solve_ode():
    # this returns plottable points
    x_list=np.array([[xstart,tstart]])
    x_list.reshape(2,1)
    



x0=1
t0=0


    
    
    
    