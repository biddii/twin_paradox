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

#t = (input("What is t?"))
#x = (input("What is x?"))
#y = (input("What is y?"))
#z = (input("What is z?"))
t = "0"
x = "1"
y = "2"
z = "-1*pi/2"

x_val = []
y_val = []
z_val = []
t_val = []

k1 = input("What is ∂x/∂t?")
def fx(t, x):
    k1
    return k1

k2 = input("What is ∂y/∂t?")
def fy(t, y):
    k2
    return k2

k3 = input("What is ∂z/∂t?")
def fz(t, z):
    k3
    return k3

n = input("What is the desired iteration count?")
h = input("What is the desired delta t?")

#eq1 = "x + y/z - y**2"
#eq2 = "x + 2*y + z"
#eq3 = "sin(z)+y"


#i end up with 3 strings representing equations all definining a function.

def odesolver(t, x, y, z, n, h): #For number of iterations 'n' and stepsize 'h'    
    for j in range(0,n):
        dxdt1 = eval(fx(t, x))
        dxdt2 = (fx(t+h/2, x+(dxdt1*h/2)))
        print(dxdt2)

        dydt = eval(fy(t, y))
        dzdt = eval(fz(t, z))
        x = round(x + dxdt1*h, 50)
        y = round(y + dydt*h, 50)
        z = round(z + dzdt*h, 50)
        t = round(t + h, 50)
        x_val.append(x)
        y_val.append(y)
        z_val.append(z)
        t_val.append(t)
        print(dxdt, dydt, dzdt)
    print("t = ",t,"x = ",x,"y = ",y,"z = ",z)
    plt.plot(t_val, x_val)
    plt.plot(t_val, y_val)
    plt.plot(t_val, z_val)
    plt.show()

odesolver(eval(t), eval(x), eval(y), eval(z), eval(n), eval(h))