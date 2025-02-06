import numpy as np
import matplotlib.pyplot as plt



t0, x0, y0, z0 = -5, 0, 0, 0  
n = 10
h = 1

#Euler method solver
def odesolver(t, x, y, z, n, h):
    t_val, x_val, y_val, z_val = [t], [x], [y], [z]

    for _ in range(n):
        dx = t          
        dy = np.sin(t)     
        dz = t**2       

        x += dx * h
        y += dy * h
        z += dz * h
        t += h

        t_val.append(t)
        x_val.append(x)
        y_val.append(y)
        z_val.append(z)
    print(f"t = {t}, x = {x}, y = {y}, z = {z}")
    plt.plot(t_val, x_val, label="x")
    plt.plot(t_val, y_val, label="y")
    plt.plot(t_val, z_val, label="z")
    plt.legend()
    plt.show()

odesolver(t0, x0, y0, z0, n, h)
