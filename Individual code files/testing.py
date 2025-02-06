import numpy as np
import math
from sympy import * 
import matplotlib.pyplot as plt
plt.show(block=True)




n = 2374725897235
s = (n - 2**(floor((np.log2(n)))))*2 + 1
print(s)
