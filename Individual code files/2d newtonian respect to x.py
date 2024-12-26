#ODE SOLVER USING RUNGE-KUTTA METHOD
import numpy as np
import math
from sympy import *
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
plt.show(block=True)
x, y, z, t, vx, vy, vz = symbols("x y z t vx vy vz")
#initial values & conditions IF I WANT TO CHANGE STUFF IT WILL ONLY NEED TO BE CHANGED UP TO ODESOLVER


vmagnitude = 100 #i can change this also
launch_angle = np.pi*0.25 # i can change this value
v_x = np.cos(launch_angle)*vmagnitude
v_y = np.sin(launch_angle)*vmagnitude
print(f"vx0 = {v_x}, vy0 = {v_y}")
state0 = np.array([0., 10., v_x, v_y]) #x0, y0, z0, vx0, vy0, vz0
t0 = 0.
dim = 4
h = 0.0001 #setting step size
n = 145000
labels = ["x(t)", "y(t)", "vx(t)", "vy(t)"]
xlabels = ["x(x)", "y(x)", "vx(x)", "vy(x)"]


def dSdt(state, t): #where state is an array
    x, y, vx, vy = state
    dxdt = (vx)
    dydt = (vy)
    dvxdt = (0.)
    dvydt = (-9.8)
    val = np.array([dxdt, dydt, dvxdt, dvydt])
    return val
def odesolver(t, n, h): #For number of iterations 'n' and stepsize 'h' THE ACTUAL SOLVER 
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
    x_ints = []
    for j in range(int(dim*0.5)):
        val = values[:, j]
        xval = values[:,0]
        for i in range(len(val) - 1):
            if min(val[i], val[i + 1]) <= state0[j] <= max(val[i], val[i + 1]):
                local_interpolator = interp1d([val[i], val[i + 1]], [tval[i].flatten()[0], tval[i + 1].flatten()[0]], bounds_error=True)
                t_ints.append(local_interpolator(state0[j]))
                if j > 0:
                    local_interpolator = interp1d([val[i], val[i + 1]], [xval[i].flatten()[0], xval[i + 1].flatten()[0]], bounds_error=True)
                    x_ints.append(local_interpolator(state0[j]))       
        if not t_ints:
            print(f"No time found where {labels[j]} = initial state")
        else:
            print(f"Times where {labels[j]} = initial state: {[float(time) for time in t_ints]}")
        if not x_ints and j > 0:
            print(f"No x found where {labels[j]} = initial state")
        if x_ints and j > 0:
            print(f"x value where {labels[j]} = initial state: {[float(time) for time in x_ints]}")
        t_ints = []


    for i in range(dim):
        max_value = max(values[:, i])
        max_index = np.argmax(values[:, i])
        max_time = tval[max_index]
        print(f"Max value for {labels[i]} is {max_value} at time t = {max_time}.")
        if 0 < i <= dim*0.5-1:
            xmax = values[max_index, 0]
            print(f"Max value for {xlabels[i]} is {max_value} at x = {xmax}.")

    plt.plot(values[:,0], values[:,1], label = labels[1])
    plt.xlabel("x Position")
    plt.ylabel("y Position")
    plt.title("Displacement path of the particle")
    plt.grid()
    plt.legend()
    plt.show()
odesolver(t0, n, h) #running the code
