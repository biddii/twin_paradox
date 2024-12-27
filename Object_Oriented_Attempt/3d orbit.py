#ODE SOLVER USING RUNGE-KUTTA METHOD
import numpy as np
import math
from sympy import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from scipy.interpolate import interp1d
from Objects import *

plt.show(block=True)

x, y, z, t, vx, vy, vz = symbols("x y z t vx vy vz")

#initial values
state0 = np.array([10000., 10000., 0., np.sqrt((((2*6.6743*10**-11)*2*10**15)/np.sqrt(2*10**8))/2), -np.sqrt((((2*6.6743*10**-11)*2*10**15)/np.sqrt(2*10**8))/2), 0.]) #x0, y0, z0, vx0, vy0, vz0
t0 = 0.
dim = len(state0) #dimensions of state0
h = 1 #setting step size 
n = 50000
mass1 = 2.*(10**15)
Gconst = 6.67430*(10**-11)
labels = ["x(t)", "y(t)", "z(t)", "vx(t)", "vy(t)", "vz(t)"]
#    velocity_mag = np.sqrt((state[3])**2+(state[4])**2+(state[5])**2) #velocity magnitude from target mass


def dSdt(state, t): #where state is an array
    x, y, z, vx, vy, vz = state
    pos_mag = np.sqrt((state[0])**2+(state[1])**2+(state[2])**2) #position magnitude from target mass
    normal_pos = state[0:3]/pos_mag #position direction but magnitude is 1.
    gravity_acceleration = -((Gconst*mass1)/(pos_mag**2))*normal_pos
    dxdt = (vx) 
    dydt = (vy)
    dzdt = (vz)
    dvxdt = (gravity_acceleration[0])
    dvydt = (gravity_acceleration[1])
    dvzdt = (gravity_acceleration[2]) 
    val = np.array([dxdt, dydt, dzdt, dvxdt, dvydt, dvzdt])
    return val

def odesolver(t, n, h): #For number of iterations 'n' and stepsize 'h'
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

    t_ints = []
    for j in range(int(dim*0.5-1)): #to include z get rid of -1
        val = values[:, j]
        for i in range(len(val) - 2): #to include z change -2 to -1
            if min(val[i], val[i + 1]) <= state0[j] <= max(val[i], val[i + 1]):
                local_interpolator = interp1d([val[i], val[i + 1]], [tval[i].flatten()[0], tval[i + 1].flatten()[0]], bounds_error=True)
                t_ints.append(local_interpolator(state0[j]))

        if not t_ints:
            print(f"No intersection found where {labels[j]} = initial state")
        else:
            print(f"Intersection times where {labels[j]} = initial state: {[float(time) for time in t_ints]}")
        t_ints = [] 
    for i in range(dim):
        max_value = max(values[:, i])
        max_index = np.argmax(values[:, i])
        max_time = tval[max_index]
        print(f"Max value for {labels[i]} is {max_value} at time t = {max_time}.")
      
    plotter = stationary_plot(values)
    plotter.plotf()

odesolver(t0, n, h)
