# -*- coding: utf-8 -*-
"""
Created on Sun Oct 15 23:13:54 2017

@author: Joby
"""

"""
Repeat part 1 using the 4th-order Runge-Kutta method.
Make it so that when calling solve_ode you can choose whether to use the Euler method or RK4.
How does the error depend on Δt now? 
How does this compare with the error for the Euler method (put this in the same plot)?
Find step-sizes for each method that give you the same error - how long does each method take? 
(you can use the time command when running your Python script)
"""

"""I can import stuff but"""

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

def RK4(x1,t, dt=0.01):
    k1=dxdt(x1,t)
    k2=dxdt(x1+k1*dt/2,t+dt/2)
    k3=dxdt(x1+k2*dt/2,t+dt/2)
    k4=dxdt(x+k3*dt,t+dt)
    
    x2=x1+(k1+2*k2+2*k3+k4)*dt/6
    
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