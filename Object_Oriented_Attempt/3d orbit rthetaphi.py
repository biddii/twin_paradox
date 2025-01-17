import numpy as np
import matplotlib.pyplot as plt

# Constants
mass1 = 2.E15  # Mass of the central body
Gconst = 6.67430E-11  # Gravitational constant

# Initial conditions
r0 = 1e7  # Initial radial distance
th0 = np.pi / 2  # Initial angle in the equatorial plane
ph0 = 0  # Initial azimuthal angle
vr0 = 1000  # Initial radial velocity (experiment with this)
vth0 = 0  # Initial angular velocity in the theta direction
vph0 = np.sqrt((Gconst * mass1) / r0)  # Initial phi velocity

state0 = np.array([r0, th0, ph0, vr0, vth0, vph0])  # Initial state
t0 = 0  # Start time
dim = len(state0)  # Dimensions of the state vector
h = 0.1  # Step size
n = 10000  # Number of iterations

# Define the system of equations (derivatives)
def dSdt(state, t):
    r, th, ph, vr, vth, vph = state
    drdt = vr
    dthdt = vth
    dphdt = vph
    dvrdt = (-Gconst * mass1) / (r**2) + (vth**2 / r) + (vph**2 / r)
    dvthdt = - (2 * vr * vth) / (r**2) + (vph**2 * np.cos(th)) / (r**2 * np.sin(th))
    dvphdt = - (2 * vr * vph) / (r**2 * np.sin(th)) - (vth * vph * np.cos(th)) / (r**2 * np.sin(th)**2)
    
    return np.array([drdt, dthdt, dphdt, dvrdt, dvthdt, dvphdt])

# Apply periodic correction to keep the motion in the equatorial plane
def apply_correction(state):
    r, th, ph, vr, vth, vph = state
    z = r * np.cos(th)  # Calculate the z-coordinate
    if np.abs(z) > 1e-15:  # If the z-coordinate deviates too much
        th = np.pi / 2  # Correct the theta angle to stay in the equatorial plane
        state[1] = th
    return state

# Runge-Kutta ODE solver
def odesolver(t, n, h):
    state = state0
    t = t0
    tval = np.zeros([n, 1])
    values = np.zeros([n, dim])
    
    for j in range(0, n):
        # Runge-Kutta integration
        k1 = dSdt(state, t)
        k2 = dSdt((state + ((h * 0.5) * k1)), t + h * 0.5)
        k3 = dSdt((state + ((h * 0.5) * k2)), t + h * 0.5)
        k4 = dSdt((state + ((h * 0.5) * k3)), t + h)
        
        state = state + (h / 6.) * (k1 + 2 * k2 + 2 * k3 + k4)
        t += h

        # Apply periodic correction to keep motion in the equatorial plane
        state = apply_correction(state)
        
        # Store values
        values[j] = state
        tval[j] = t

    r_vals = values[:, 0]
    th_vals = values[:, 1]
    ph_vals = values[:, 2]
    
    # Convert to Cartesian coordinates
    x_vals = r_vals * np.sin(th_vals) * np.cos(ph_vals)
    y_vals = r_vals * np.sin(th_vals) * np.sin(ph_vals)
    z_vals = r_vals * np.cos(th_vals)

    # Plot the orbit in 3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(x_vals, y_vals, z_vals)
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.set_zlabel('z (m)')
    plt.show()

# Solve the ODE and plot the result
odesolver(t0, n, h)
