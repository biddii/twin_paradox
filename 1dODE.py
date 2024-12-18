#ODE SOLVER USING EULER's METHOD
# Define an ODE solver function
# input dy/dx 
# input starting conditions on y0 and x0
# input n and h (step size)
# execute the function
import numpy as np
import math
from sympy import * 
import matplotlib.pyplot as plt

plt.show(block=True)

x, y, z, t = symbols("x y z t")

x = "-10"
y = "2"

x_val = []
y_val = []

def odesolver(x, y, n, h): #For number of iterations 'n' and stepsize 'h'   
    equation = input("What is dydx?")
    for j in range(0,n):
        dydx = eval(equation)
        x = round(x + h, 50)
        y = round(y + dydx*h, 50)
        x_val.append(x)
        y_val.append(y)
    print("x = ",x,"y = ",y)
    plt.plot(x_val, y_val)
    plt.show()

odesolver(eval(x), eval(y), 20000, 0.001)