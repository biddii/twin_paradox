import numpy as np
import matplotlib.pyplot as plt

def odesolver(f, x0, y0, n, h):
    x = x0
    y = y0
    x_val = [x]
    y_val = [y]
    
    for i in range(n):
        dydx = f(x)
        x += h
        y += dydx * h
        x_val.append(x)
        y_val.append(y)
    
    print(f"x = {x}, y = {y}")
    plt.plot(x_val, y_val)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Euler's Method ODE Solver")
    plt.show()

#defining dy/dx
def func(x):
    return np.sin(x)

#running the code
odesolver(func, x0=-5, y0=2, n=10, h=1)