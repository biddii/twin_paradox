#ODE SOLVER USING RUNGE-KUTTA METHOD
import numpy as np
import math
from sympy import *
import matplotlib.pyplot as plt
plt.show(block=True)

x, y, z, t, vx, vy, vz = symbols("x y z t vx vy vz")

#initial values
state0 = np.array([100., 100.]) #x0, vx0
t0 = 0.
dim = 2
h = 0.0001 #setting step size
n = 200000
labels = ["x(t)", "vx(t)"]

def dSdt(state, t): #where state is an array
    x, vx = state
    dxdt = (vx)
    dvxdt = (-9.8)
    val = np.array([dxdt, dvxdt])
    return val

def odesolver(t, n, h): #For number of iterations 'n' and stepsize 'h'
    #values = np.zeros([n, 3])
    state = state0
    t = t0
    tval = np.zeros([n, 1])
    values = np.zeros([n, dim])
    for j in range(0,n):
        #Runge-Kutta for x,y,z respect to t
        k1 = dSdt(state, t)
        k2 = dSdt((state+((h*0.5)*k1)), t+h*0.5)
        k3 = dSdt((state+((h*0.5)*k2)), t+h*0.5)
        k4 = dSdt((state+((h*0.5)*k3)), t+h)
        state = state + (h/6.)*(k1+2*k2+2*k3+k4)
        t += h

        #adding each position value to an array, so it can be plotted
        values[j] = state
        tval[j] = t 
    #plotting
    for i in range(dim):
         print(f"max value for {labels[i]} is {max(values[:,i])}.")
         abs_diff = np.abs(values-state0)
         closest_val = np.argmin(abs_diff, axis=0)
         print(f"time taken to reach initial {labels[i]} is {closest_val[i]}")
         plt.plot(tval, values[:,i], label = labels[i])
    plt.xlabel("Time (t)")
    plt.ylabel("Values")
    plt.title("Runge-Kutta Solution")
    plt.grid()
    plt.legend()
    plt.show()
#running the function
odesolver(t, n, h) 
