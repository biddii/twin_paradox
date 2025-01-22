import numpy as np
import math
from sympy import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from scipy.interpolate import interp1d

class Rotating_plot:
    def __init__(self, values):
         self.values = values
         
    def plotf(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot3D (self.values[:,0], self.values[:,1], self.values[:,2])
        ax.set_title('3d xyz plot')
        plt.xlabel("x Position")
        plt.ylabel("y Position")
        for angle in range(0, 360):
            ax.view_init(10, angle)
            plt.draw()
            plt.pause(.01)

class Stationary_plot:
    def __init__(self, values):
        self.values = values
    
    def plotf(self):
        fig = plt.figure()
        ax = plt.axes(projection = '3d')
        ax.plot3D(self.values[:,0], self.values[:,1], self.values[:,2], 'red')
        ax.set_title('3d spatial path of orbit')
        ax.set_xlabel('x position')
        ax.set_ylabel('y position')
        ax.set_zlabel('z position')
        plt.show()

class Polar_stationary_plot:
    def __init__(self, xyz_values):
        self.values = xyz_values

    def plotf(self):
        fig = plt.figure()
        ax = plt.axes(projection = '3d')
        ax.plot3D(self.values[:,0], self.values[:,1], self.values[:,2], 'red')
        ax.set_title('3d spatial path of orbit')
        ax.set_xlabel('x position')
        ax.set_ylabel('y position')
        ax.set_zlabel('z position')
        plt.axis('equal')
        plt.show()
    