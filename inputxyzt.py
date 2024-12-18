#ODE SOLVER USING RUNGE-KUTTA METHOD
import numpy as np
import math
from sympy import * 
import matplotlib.pyplot as plt
import re

plt.show(block=True)

x, y, z, t = symbols("x y z t")
#initial values
t = "-5"
x = "0"
y = "0"
z = "0"
def string_extract(j):
    replace = {'x':'{x}', 'y':'{y}', 'z': '{z}', 't': '{t}'}
    exempt = {"sqrt", "tan"}
    pattern = r'(tan|sin|cos|sqrt|[txyz])'
    def replacement(match):
        word = match.group()
        if word in exempt:
            return word
        return replace.get(word, word)
    g = (re.sub(pattern, replacement, j))
    g = f"{g}"

    return(re.sub(pattern, replacement, j))

dx_input = string_extract(input("What is dx/dt"))
dy_input = string_extract(input("What is dy/dt"))
dz_input = string_extract(input("What is dz/dt"))

def fx(x, y, z, t): #xt partial
    d=eval(eval("f'{}'".format(dx_input)))
    return d

def fy(x, y, z, t): #yt partial
    d=eval(eval("f'{}'".format(dy_input)))
    return d

def fz(x, y, z, t): #yt partial
    d=eval(eval("f'{}'".format(dz_input)))
    return d

h = "0.001" #setting step size
def odesolver(t, x, y, z, n, h): #For number of iterations 'n' and stepsize 'h'
    x_val = [] #empty arrays to be filled with data to plot
    y_val = []
    z_val = []
    t_val = []
    for j in range(0,n):
        #Runge-Kutta for x,y,z respect to t
        dx1 = fx(x,y,z,t)
        dy1 = fy(x,y,z,t)
        dz1 = fz(x,y,z,t)

        dx2 = fx(x+h*dx1/2, y+h*dy1/2, z+h*dz1/2, t+h/2)
        dy2 = fy(x+h*dx1/2, y+h*dy1/2, z+h*dz1/2, t+h/2)
        dz2 = fz(x+h*dx1/2, y+h*dy1/2, z+h*dz1/2, t+h/2)

        dx3 = fx(x+h*dx2/2, y+h*dy2/2, z+h*dz2/2, t+h/2)
        dy3 = fy(x+h*dx2/2, y+h*dy2/2, z+h*dz2/2, t+h/2)
        dz3 = fz(x+h*dx2/2, y+h*dy2/2, z+h*dz2/2, t+h/2)
        
        dx4 = fx(x+h*dx3, y+h*dy3, z+h*dz3, t+h)
        dy4 = fy(x+h*dx3, y+h*dy3, z+h*dz3, t+h)
        dz4 = fz(x+h*dx3, y+h*dy3, z+h*dz3, t+h)
    
        #updating each term
        x = round(x + (h/6)*(dx1+2*dx2+2*dx3+dx4), 50)
        y = round(y + (h/6)*(dy1+2*dy2+2*dy3+dy4), 50)
        z = round(z + (h/6)*(dz1+2*dz2+2*dz3+dz4), 50)
        t = round(t + h, 50)

        #adding each position value to an array, so it can be plotted
        x_val.append(x)
        y_val.append(y)
        z_val.append(z)
        t_val.append(t)

    #plotting
    plt.plot(t_val, x_val, label="x(t)")
    plt.plot(t_val, y_val, label="y(t)")
    plt.plot(t_val, z_val, label="z(t)")
    plt.xlabel("Time (t)")
    plt.ylabel("Values")
    plt.title("Runge-Kutta Solution")
    plt.grid()
    plt.legend()
    plt.show()
#running the function
odesolver(eval(t), eval(x), eval(y), eval(z), 10000, eval(h)) 
