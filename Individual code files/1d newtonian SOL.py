#ODE SOLVER USING RUNGE-KUTTA METHOD
import numpy as np
import math
from sympy import *
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
plt.show(block=True)

x, y, z, t, vx, vy, vz = symbols("x y z t vx vy vz")

#initial values
state0 = np.array([0., 0]) #x0, vx0 
t0 = 0.
dim = 2
h = 0.00001 #setting step size
n = 250000
labels = ["x(t)", "vx(t)"]

def dSdt(state, t): #where state is an array
    x, vx = state
    dxdt = (vx)
    dvxdt = (np.sin(t*10))
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
    xval = values[:,0]
    tval= tval
    t_ints = []
    for i in range(len(xval) - 1):
        if min(xval[i], xval[i + 1]) <= state0[0] <= max(xval[i], xval[i + 1]):
            local_interpolator = interp1d([xval[i], xval[i + 1]], [tval[i].flatten()[0], tval[i + 1].flatten()[0]], bounds_error=True)
            t_ints.append(local_interpolator(state0[0]))

    if not t_ints:
        print("No intersection found where x(t) = initial state")
    else:
        print(f"Intersection times where x(t) = initial state: {[float(time) for time in t_ints]}")
        print(f"Analytic result to compare: {[float(2.*state0[1]/9.8)]}")

    for i in range(dim):
        max_value = max(values[:, i])
        max_index = np.argmax(values[:, i])
        max_time = tval[max_index]
        print(f"Max value for {labels[i]} is {max_value} at time t = {max_time}.")

        plt.plot(tval, values[:,i], label = labels[i])
    plt.xlabel("Time (t)")
    plt.ylabel("Values")
    plt.title("Runge-Kutta Solution")
    plt.grid()
    plt.legend()
    plt.show()

intersection_times = odesolver(t0, n, h)
