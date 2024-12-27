import numpy as np
import math
from sympy import * 
import matplotlib.pyplot as plt
plt.show(block=True)

x, y, z, t = symbols("x y z t")

y = 1.

print(np.sqrt(y))

dim = 4.
for i in range(int(dim*0.5)):
    print(i)


somearray = np.array([2. , 5. , 10. , 8.])

otherarray = somearray[0:3]/2
print(otherarray)