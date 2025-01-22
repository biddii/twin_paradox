#ODE SOLVER USING RUNGE-KUTTA METHOD
import numpy as np
import math
from sympy import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from scipy.interpolate import interp1d
from Objects import *

plt.show(block=True)

x, y, z, t, vx, vy, vz, r, th, ph, vr, vth, vph= symbols("x y z t vx vy vz r th ph vr vth vph")

#initial values
mass1 = 2.E15
Gconst = 6.67430E-11
state0 = np.array([1e7, 0.5*np.pi, 0, 0, 0, np.sqrt((Gconst*mass1)/10000000)]) #r0, th0, ph0, vr0, vth0, vph0
t0 = 0.
dim = len(state0) #dimensions of state0
h = 0.001 #setting step size 
n = 500000
labels = ["r(t)", "th(t)", "ph(t)", "vr(t)", "vth(t)", "vph(t)"]
#    velocity_mag = np.sqrt((state[3])**2+(state[4])**2+(state[5])**2) #velocity magnitude from target mass
#    pos_mag = np.sqrt((state[0])**2+(state[1])**2+(state[2])**2) #position magnitude from target mass
#    normal_pos = state[0:3]/pos_mag #position direction but magnitude is 1.

def dSdt(state, t):
    r, th, ph, vr, vth, vph = state
    drdt = vr
    dthdt = vth/r
    dphdt = vph/(r*np.sin(th))
    dvrdt = (-Gconst*mass1)/(r**2) + (r*(vth**2)) + (r*(np.sin(th))**2)*(vph**2)
    dvthdt = ((-2/r)*vr*vth) + ((np.sin(th)*np.cos(th))*(vph**2))
    dvphdt = ((-2/r)*vr*vph) - (((2*np.cos(th))/(np.sin(th)))*vph*vth)
    return np.array([drdt, dthdt, dphdt, dvrdt, dvthdt, dvphdt])

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

        values[j] = state
        tval[j] = t 
    r_vals = values[:, 0]
    th_vals = values[:, 1]
    ph_vals = values[:, 2]
    x_vals = r_vals * np.sin(th_vals) * np.cos(ph_vals)
    y_vals = r_vals * np.sin(th_vals) * np.sin(ph_vals)
    z_vals = r_vals * np.cos(th_vals)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(x_vals, y_vals, z_vals)
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.set_zlabel('z (m)')
    plt.show()

odesolver(t0, n, h)
