#ODE SOLVER USING RUNGE-KUTTA METHOD
import numpy as np
import math
from mpmath import *
from sympy import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from scipy.interpolate import interp1d
from Objects import *

plt.show(block=True)

x, y, z, t, vx, vy, vz = symbols("x y z t vx vy vz")

#initial values
mp.dps = 15; mp.pretty = True
h = 1 #setting step size 
n = 60000
mass1 = 2.E15
Gconst = 6.67430E-11
state0_xyz = np.array([10000., 10000., 0., -np.sqrt((((Gconst*mass1))/np.sqrt(2*10**8))/2), np.sqrt(((Gconst*mass1)/np.sqrt(2*10**8))/2), 0.]) #x0, y0, z0, vx0, vy0, vz0
t0 = 0.
dim = len(state0_xyz) #dimensions of state0
labels = ["x(t)", "y(t)", "z(t)", "vx(t)", "vy(t)", "vz(t)"]
#    velocity_mag = np.sqrt((state[3])**2+(state[4])**2+(state[5])**2) #velocity magnitude from target mass

polar_radius = np.sqrt(state0_xyz[0]**2 + state0_xyz[1]**2 + state0_xyz[2]**2)
polar_theta = np.arccos(state0_xyz[2] / polar_radius)
polar_phi = np.arctan2(state0_xyz[1], state0_xyz[0])

x, y, z, vx, vy, vz = state0_xyz
v_r = (x * vx + y * vy + z * vz) / polar_radius
v_theta = (z * (x * vx + y * vy) - polar_radius**2 * vz) / (polar_radius**2 * np.sqrt(x**2 + y**2))
v_phi = (x * vy - y * vx) / (x**2 + y**2)

state0 = np.array([polar_radius, polar_theta, polar_phi, v_r, v_theta, v_phi])

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
    for i in range(np.shape(val)[0]-1):
        if min(val[i], val[i + 1]) <= 2*np.pi + state0[2] <= max(val[i], val[i + 1]):
            local_interpolator = interp1d([val[i], val[i + 1]], [tval[i].flatten()[0], tval[i + 1].flatten()[0]], bounds_error=True)
            t_ints.append(local_interpolator(2*np.pi+state0[2]))

    if not t_ints:
        print(f"Couldn't find a time where the orbit reaches initial state.")
    else:
        print(f"Times where the orbit = initial state: {[float(time) for time in t_ints]} seconds")


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
