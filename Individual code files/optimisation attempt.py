#ODE SOLVER USING RUNGE-KUTTA METHOD
import numpy as np
import math
from sympy import * 
import matplotlib.pyplot as plt
plt.show(block=True)

x, y, z, t = symbols("x y z t")

#initial values
state0 = np.array([1., 1., 1.])
t0 = 0.
dim = 3
h = 0.001 #setting step size
n = 10000
labels = ["x(t)", "y(t)", "z(t)"]

def dSdt(state, t): #where state is an array
    x, y, z = state
    dxdt = (y) #defining dxdt as y, thus x should be the integration of y in the plot
    dydt = (0) #arbitrary function i picked for dy/dt 
    dzdt = (x)
    val = np.array([dxdt, dydt, dzdt])
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
         plt.plot(tval, values[:,i], label = labels[i])
    plt.xlabel("Time (t)")
    plt.ylabel("Values")
    plt.title("Runge-Kutta Solution")
    plt.grid()
    plt.legend()
    plt.show()
#running the function
odesolver(t, n, h) 
