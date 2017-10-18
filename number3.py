# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 00:02:24 2017

@author: Joby
"""
"""
Consider the bouncing ball problem referred to above. 
Using the Euler method with the same step-size as in the provided script, 
show how the numerically generated trajectories appear in the phase-plane 
(plot y against v). 
How should the trajectories look? 
Considering the phase-plane argue that the Euler method is guaranteed to diverge 
in the way that it does for any initial conditions and any step-size.
"""

#FINE lets try importing
from number2 import euler_step
from number2 import RK4
from number2 import adjust_deltat
from number2 import solve_to
from number2 import solve_ode

def dxdt(x,t):
    dxdt=x
    return dxdt

if __name__ == "__main__":
    x0=1
    t0=0

    t=np.arange(10)
    t=t/2
    euler=solve_ode(euler_step, x0, t0 ,t, deltat=0.05)
    v=np.array([solve_to(euler_step,x0,t0,time,deltat=0.05) for time in t])
    fig1=plt.figure()
    ax1=fig1.add_subplot(1,1,1)
    ax1.plot(euler,v ,'bs')
    ax1.set_ylabel('y')
    ax1.set_xlabel('v')
    fig1.savefig('phase space.svg')
    