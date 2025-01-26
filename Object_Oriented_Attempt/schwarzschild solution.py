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
r0 = 5
th0 = 0.5*np.pi
ph0 = 0.
f0 = 0.
v_r0 = 0.
v_th0 = 0.
v_ph0 = np.sqrt(Gconst*mass1/r0)

#do not change these formulas
r_dot0 = v_r0
th_dot0 = (v_th0)/(r0)
ph_dot0 = (v_ph0)/((r0)*(np.sin(th0)))
f_dot0 = np.sqrt((((r_dot0**2)/(1-2*mass1/r0)) + (r0**2)*((th_dot0**2)+(((np.sin(th0))**2)*(ph_dot0**2))) - q_const)/(1-(2*mass1/r0)))

#initial time velocity calc ^^^^^ above f0

state0 = np.array([r0, th0, ph0, f0])

#constants of motion
E_val = -(f_dot0)*(1-(2*mass1/r0))
L_val = (ph_dot0)*(r0**2)*(((np.sin(th0))**2))
K_val = (r0**4)*((th_dot0**2)+(((np.sin(th0))**2)*(ph_dot0**2)))

#important numbers
t0 = 0.
dim = len(state0) #dimensions of state0
h = 0.001 #setting step size 
n = 200000
labels = ["r(t)", "th(t)", "ph(t)", "r_dot(t)", "th_dot(t)", "ph_dot(t)"]


def dSdt(state, t, sign_r_dot, sign_th_dot):
    r, th, ph, f = state
    ph_dot = (L_val)/((r**2)*((np.sin(th))**2))
    r_dot2 = (E_val)**2 + ((1-(2*mass1/r))*(q_const - (K_val/(r**2))))
    th_dot2 = (K_val-((L_val**2)/((np.sin(th))**2)))/(r**4)
    f_dot = (-E_val)/(1-(2*mass1/r))

    #theta dot and r dot value fixing due to the square
    if r_dot2 < 0:
        r_dot2 = -r_dot2  #making sqrt positive
        sign_r_dot *= -1  #flipping the sign
    r_dot = sign_r_dot * np.sqrt(r_dot2)

    if th_dot2 < 0:
        th_dot2 = -th_dot2  #making sqrt positive
        sign_th_dot *= -1  #flipping the sign
    th_dot = sign_th_dot * np.sqrt(th_dot2)

    return np.array([r_dot, th_dot, ph_dot, f_dot]), sign_r_dot, sign_th_dot

def odesolver(t, n, h):
    state = np.copy(state0)
    t = t0
    tval = np.zeros([n, 1])
    values = np.zeros([n, dim])
    t_ints = []
    f_ints = []
    
    #setting sign values
    sign_r_dot = -1 if r_dot0 < 0 else 1
    sign_th_dot = -1 if th_dot0 < 0 else 1

    for j in range(n):
        #Runge-Kutta and sign tracking
        k1, sign_r_dot, sign_th_dot = dSdt(state, t, sign_r_dot, sign_th_dot)
        k2, sign_r_dot, sign_th_dot = dSdt(state + 0.5 * h * k1, t + 0.5 * h, sign_r_dot, sign_th_dot)
        k3, sign_r_dot, sign_th_dot = dSdt(state + 0.5 * h * k2, t + 0.5 * h, sign_r_dot, sign_th_dot)
        k4, sign_r_dot, sign_th_dot = dSdt(state + h * k3, t + h, sign_r_dot, sign_th_dot)

        #updating state and time
        state += (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
        t += h
        values[j] = state
        tval[j] = t
    val = values[:, 2]
    fval = values[:, 3]
    n=1
    #testing if phi is increasing or decreasing to see when an orbit happens
    if val[0]-val[1] <= 0: 
        for i in range(np.shape(val)[0]-1):
            if min(val[i], val[i + 1]) <= 2*n*np.pi + state0[2] <= max(val[i], val[i + 1]):
                tlocal_interpolator = interp1d([val[i], val[i + 1]], [tval[i].flatten()[0], tval[i + 1].flatten()[0]], bounds_error=True)
                t_ints.append(tlocal_interpolator(2*n*np.pi+state0[2]))

                flocal_interpolator = interp1d([val[i], val[i + 1]], [fval[i].flatten()[0], fval[i + 1].flatten()[0]], bounds_error=True)
                f_ints.append(flocal_interpolator(2*n*np.pi+state0[2]))
                n+=1
        if not t_ints:
            print(f"Couldn't find a time where the orbit reaches initial state. Orbit is anticlockwise.")
        else:
            print(f"Times where the orbit = initial state: {[float(time) for time in t_ints]} proper time units (orbiting). Orbit is anticlockwise")
            print(f"Times where the orbit = initial state: {[float(time) for time in f_ints]} time units (stationary). Orbit is anticlockwise")
    else:
        for i in range(np.shape(val)[0]-1):
            if min(val[i], val[i + 1]) <= -2*n*np.pi + state0[2] <= max(val[i], val[i + 1]):
                tlocal_interpolator = interp1d([val[i], val[i + 1]], [tval[i].flatten()[0], tval[i + 1].flatten()[0]], bounds_error=True)
                t_ints.append(tlocal_interpolator(2*n*np.pi+state0[2]))

                flocal_interpolator = interp1d([val[i], val[i + 1]], [fval[i].flatten()[0], fval[i + 1].flatten()[0]], bounds_error=True)
                f_ints.append(flocal_interpolator(2*n*np.pi+state0[2]))
                n+=1
        if not t_ints:
            print(f"Couldn't find a time where the orbit reaches initial state. Orbit is anticlockwise.")
        else:
            print(f"Times where the orbit = initial state: {[float(time) for time in t_ints]} proper time units (orbiting). Orbit is anticlockwise")
            print(f"Times where the orbit = initial state: {[float(time) for time in f_ints]} time units (stationary). Orbit is anticlockwise")
    #converting to xyz for plot
    r_vals = values[:, 0]
    th_vals = values[:, 1]
    ph_vals = values[:, 2]
    x_vals = r_vals * np.sin(th_vals) * np.cos(ph_vals)
    y_vals = r_vals * np.sin(th_vals) * np.sin(ph_vals)
    z_vals = r_vals * np.cos(th_vals)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(x_vals, y_vals, z_vals)
    ax.plot(0,0,0, 'ro')
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_zlabel("z (m)")
    plt.axis('equal')
    plt.show()


#running it
odesolver(t0, n, h)


#add in values the state0 as i = 0 of array
