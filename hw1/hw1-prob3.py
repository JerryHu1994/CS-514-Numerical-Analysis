# -*- coding: utf-8 -*-
#ＣＳ514 HW1 - Problem 2
# Import the necessary modeules
import numpy as np
from matplotlib import pyplot as plt

def prob3(x):
    f = x - x**3 - 4*x**2 + 10;
    return f

plt.figure(2)
x = 1.5
# perform function iteration with x0 = 1.5
for k in np.arange(6):
    print ('x%d = %.5f' %(k,x))
    x = prob3(x)
    plt.plot(k, x, 'b.')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.title('Function iteration with X0 = 1.5')
plt.show()