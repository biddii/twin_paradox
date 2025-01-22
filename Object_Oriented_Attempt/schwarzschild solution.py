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

#setting G constant and blackhole mass to unity. 
mass1 = 1.
Gconst = 1.
q_const = -1 #-1 for mass, 0 for massless

#initial conditions variables to change, ensure not greater than 1 for v_
r0 = 3.
th0 = 0.5*np.pi
ph0 = 0.
v_r0 = 0
v_th0 = 0.
v_ph0 = np.sqrt(Gconst*mass1/r0)

#do not change these formulas
r_dot0 = v_r0
th_dot0 = v_th0/r0
ph_dot0 = v_ph0/(r0*np.sin(th0))
f_dot0 = np.sqrt((q_const-(((1-(2*mass1/r0))**(-1))*(r_dot0**2))-((r0**2)*((th_dot0**2)+((np.sin(th0))**2)*(ph_dot0**2))))/(-1+(2*mass1/r0)))
#initial time velocity calc ^^^^^

state0 = np.array([r0, th0, ph0, f_dot0])


E_val = -1*(1-((2*mass1)/r0))*f_dot0
L_val = (r0**2)*((np.sin(th0))**2)*ph_dot0
K_val = (r0**4)*((th_dot0**2)+(((np.sin(th0))**2)*(ph_dot0**2)))


#important numbers
t0 = 0.
dim = len(state0) #dimensions of state0
h = 0.01 #setting step size 
n = 10000
labels = ["r(t)", "th(t)", "ph(t)", "r_dot(t)", "th_dot(t)", "ph_dot(t)"]

def dSdt(state, t):
    r, th, ph, f = state
    ph_dot = L_val/((r**2)*((np.sin(th))**2))
    r_dot2 = (E_val**2)+((1-(2*mass1/r))*(q_const-(K_val/(r**2))))
    th_dot2 = (1/(r**4))*(K_val-((L_val**2)/(np.sin(th))**2))
    f_dot = -E_val*((1-((2*mass1)/r))**-1)

    if th_dot2 < 0:
        th_dot2 = -1*th_dot2
        th_dot = -np.sqrt(th_dot2)
    else: 
        th_dot = np.sqrt(th_dot2)

    if r_dot2 < 0:
        r_dot2 = -1*r_dot2
        r_dot = -np.sqrt(r_dot2)
    else:
        r_dot = np.sqrt(r_dot2)
    

    return np.array([r_dot, th_dot, ph_dot, f_dot])


def odesolver(t, n, h): #For number of iterations 'n' and stepsize 'h'
    state = state0
    t = t0
    tval = np.zeros([n, 1])
    values = np.zeros([n, dim])
    t_ints = []
    for j in range(0,n):
        #Runge-Kutta for r, theta, phi, f respect to t
        k1 = dSdt(state, t)
        k2 = dSdt((state+((h*0.5)*k1)), t+h*0.5)
        k3 = dSdt((state+((h*0.5)*k2)), t+h*0.5)
        k4 = dSdt((state+((h*0.5)*k3)), t+h)
        state = state + (h/6.)*(k1+2*k2+2*k3+k4)
        t += h
        values[j] = state
        tval[j] = t 
    val = values[:, 2] #setting val to all of the phi values
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