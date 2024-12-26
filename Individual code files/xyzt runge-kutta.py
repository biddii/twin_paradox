#ODE SOLVER USING RUNGE-KUTTA METHOD
import numpy as np
import math
from sympy import * 
import matplotlib.pyplot as plt

plt.show(block=True)

x, y, z, t = symbols("x y z t")
#initial values
t = "-5"
x = "0"
y = "0"
z = "0"
def fx(x, y, z, t): #xt
    k = (t)**2
    return k
def fy(x, y, z, t): #yt
    k = (t)
    return k
def fz(x, y, z, t): #zt
    k = (t+2)
    return k

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
        print(dx1, dy1, dz1)

        dx2 = fx(x+h*dx1*0.5, y+h*dy1*0.5, z+h*dz1*0.5, t+h*0.5)
        dy2 = fy(x+h*dx1*0.5, y+h*dy1*0.5, z+h*dz1*0.5, t+h*0.5)
        dz2 = fz(x+h*dx1*0.5, y+h*dy1*0.5, z+h*dz1*0.5, t+h*0.5)
        print(dx2, dy2, dz2)

        dx3 = fx(x+h*dx2*0.5, y+h*dy2*0.5, z+h*dz2*0.5, t+h*0.5)
        dy3 = fy(x+h*dx2*0.5, y+h*dy2*0.5, z+h*dz2*0.5, t+h*0.5)
        dz3 = fz(x+h*dx2*0.5, y+h*dy2*0.5, z+h*dz2*0.5, t+h*0.5)
        print(dx3, dy3, dz3)
        
        dx4 = fx(x+h*dx3, y+h*dy3, z+h*dz3, t+h)
        dy4 = fy(x+h*dx3, y+h*dy3, z+h*dz3, t+h)
        dz4 = fz(x+h*dx3, y+h*dy3, z+h*dz3, t+h)
        print(dx4, dy4, dz4)
    
        #updating each term
        x += (h/6.)*(dx1+2*dx2+2*dx3+dx4)
        y += (h/6.)*(dy1+2*dy2+2*dy3+dy4)
        z += (h/6.)*(dz1+2*dz2+2*dz3+dz4)
        t += h

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
