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
import time

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
    k4=dxdt(x1+k3*dt,t+dt)
    
    x2=x1+(k1+2*k2+2*k3+k4)*dt/6
    
    return x2

def adjust_deltat(tstart,tend,deltat):
    #this is not rk4 step adaptive step size
    fsteps=(tend-tstart)/deltat
    steps=int((tend-tstart)/deltat)
    if(fsteps-steps>=0.000005):
        steps+=1  
    deltat=(tend-tstart)/steps
    return deltat, steps

#Also make a function called solve_to which solves from x1,t1 to x2,t2 in steps no bigger than deltat_max.
def solve_to(f, xstart,tstart,tend,deltat=0.01):
    # should return only last variable
    #steps=int((tend-tstart)/deltat)
    xnow=xstart
    tnow=tstart
    deltat, steps=adjust_deltat(tstart, tend, deltat)
    
    for step in range(steps)[:-1]:
        xnow =f(xnow,tnow, deltat)
        tnow += deltat
    else: #last step to reach tend exactly
        deltat=tend-tnow
        xnow =f(xnow,tnow, deltat)      
            
    return xnow
    
    
def solve_ode(f, x0, t0, t, deltat=0.01):
    # t is a bloddy array!!!
    #x_list=np.array([solve_to(x0,t0,t1, deltat) for t1 in t ])
    #having t0 here is irritating me somehow
    #this will become really slow
    
    x_list=np.array([])
    tnow=t0 
    xnow=x0
    for time_step in t:
        new_x=solve_to(f, xnow, tnow, time_step, deltat)
        x_list=np.append(x_list,new_x)
        tnow=time_step
        xnow=new_x
        #print(xnow)
    
    return x_list



if __name__ == "__main__":
    x0=1
    t0=0

    t=np.array([1])
    #deltat_list=[(10**-(z+1)) for z in range(10)]
    array1=np.array( [10**-(z+1) for z in range(7)])
    array2=array1/2
    deltat_array=np.append(array1,array2)
    deltat_array.sort()
    
    #deltat_array=np.array([0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001, 0.000005, 0.000001)

 
    euler=np.array([])
    runge_katta=np.array([])
    euler_time=np.array([])
    rk_time=np.array([])
    for z in range(len(deltat_array)):
        #euler set
        start_time=time.time()
        neweuler=solve_ode(euler_step, x0,t0, t, deltat_array[z])
        elapsed_time = time.time() - start_time
        euler=np.append(euler,neweuler)
        euler_time=np.append(euler_time, elapsed_time)
    
        #rk4 set
        start_time=time.time()
        newrunge=solve_ode(RK4, x0,t0, t, deltat_array[z])
        elapsed_time = time.time() - start_time
        runge_katta=np.append(runge_katta,newrunge)
        rk_time=np.append(rk_time, elapsed_time)
    
        print("delta t is now {}".format(deltat_array[z]))

    err1=abs(np.exp(1)-euler)
    err2=abs(np.exp(1)-runge_katta)


    fig1=plt.figure()
    ax1=fig1.add_subplot(1,1,1)
    ax1.loglog(deltat_array,err2, 'ro', deltat_array, err1, 'bs')
    ax1.set_ylabel('Error')
    ax1.set_xlabel('Delta t')
    fig1.savefig('err vs deltat.svg')

    fig2=plt.figure()
    ax2=fig2.add_subplot(1,1,1)
    ax2.loglog(deltat_array,rk_time, 'ro' , deltat_array, euler_time, 'bs')
    ax2.set_ylabel('Time')
    ax2.set_xlabel('Delta t')
    fig2.savefig('time scale.svg')

    #spending too much time trying to make legends