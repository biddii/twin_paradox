import numpy as np
import math
from sympy import * 
import matplotlib.pyplot as plt
import re

plt.show(block=True)

x, y, z, t = symbols("x y z t")
#initial values
t = "0"
x = "1"
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

k = string_extract(input("What is dy/dx"))
print(k)

def fx(x, y, z, t): #xt partial
    d=eval(eval("f'{}'".format(k)))
    print(d)
    return d

(fx(2, 4, 2, 1))

