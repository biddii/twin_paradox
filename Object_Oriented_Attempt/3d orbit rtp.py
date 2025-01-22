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

#i set these values
r0 = np.sqrt(2)*10000
th0 = 0.5*np.pi
ph0 = 0.
v_r0 = 0
v_th0 = 0
v_ph0 = 1.3*np.sqrt(Gconst*mass1/r0)

#do not change these formulas
r_dot0 = v_r0
th_dot0 = v_th0/r0
ph_dot0 = v_ph0/(r0*np.sin(th0))

state0 = np.array([r0, th0, ph0, r_dot0, th_dot0, ph_dot0])
print(state0)


t0 = 0.
dim = len(state0) #dimensions of state0
h = 1 #setting step size 
n = 200000
labels = ["r(t)", "th(t)", "ph(t)", "vr(t)", "vth(t)", "vph(t)"]
#    velocity_mag = np.sqrt((state[3])**2+(state[4])**2+(state[5])**2) #velocity magnitude from target mass
#    pos_mag = np.sqrt((state[0])**2+(state[1])**2+(state[2])**2) #position magnitude from target mass
#    normal_pos = state[0:3]/pos_mag #position direction but magnitude is 1.



def dSdt(state, t):
    r, th, ph, vr, vth, vph = state
    drdt = vr
    dthdt = vth
    dphdt = vph
    dvrdt = (-Gconst*mass1)/(r**2) + (r*(vth**2)) + (r*(np.sin(th))**2)*(vph**2)
    dvthdt = ((-2/r)*vr*vth) + ((np.sin(th)*np.cos(th))*(vph**2))
    dvphdt = ((-2/r)*vr*vph) - (((2*np.cos(th))/(np.sin(th)))*vph*vth)
    return np.array([drdt, dthdt, dphdt, dvrdt, dvthdt, dvphdt])


def odesolver(t, n, h): #For number of iterations 'n' and stepsize 'h'
    state = state0
    t = t0
    tval = np.zeros([n, 1])
    values = np.zeros([n, dim])
    t_ints = []
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
    val = values[:, 2]
    n=1
    #testing if phi is increasing or decreasing to see when an orbit happens
    if val[0]-val[1] <= 0: 
        for i in range(np.shape(val)[0]-1):
            if min(val[i], val[i + 1]) <= 2*n*np.pi + state0[2] <= max(val[i], val[i + 1]):
                local_interpolator = interp1d([val[i], val[i + 1]], [tval[i].flatten()[0], tval[i + 1].flatten()[0]], bounds_error=True)
                t_ints.append(local_interpolator(2*n*np.pi+state0[2]))
                n+=1
        if not t_ints:
            print(f"Couldn't find a time where the orbit reaches initial state. Orbit is anticlockwise.")
        else:
            print(f"Times where the orbit = initial state: {[float(time) for time in t_ints]} seconds. Orbit is anticlockwise")
    else:
        for i in range(np.shape(val)[0]-1):
            if min(val[i], val[i + 1]) <= -2*n*np.pi + state0[2] <= max(val[i], val[i + 1]):
                local_interpolator = interp1d([val[i], val[i + 1]], [tval[i].flatten()[0], tval[i + 1].flatten()[0]], bounds_error=True)
                t_ints.append(local_interpolator(-2*n*np.pi+state0[2]))
                n+=1
        if not t_ints:
            print(f"Couldn't find a time where the orbit reaches initial state. Orbit is clockwise")
        else:
            print(f"Times where the orbit = initial state: {[float(time) for time in t_ints]} seconds. Orbit is clockwise.")
    #converting to cartesian to plot
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
    plt.axis('equal')
    plt.show()
odesolver(t0, n, h)