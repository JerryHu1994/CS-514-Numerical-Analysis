# -*- coding: utf-8 -*-
# CS/MATH 514 HW1 - Problem 2
# Jieru Hu
# Import the necessary modeules
import numpy as np
from matplotlib import pyplot as plt

# Function implements the RHS of the equation
def prob2(x, c=0.3):
    f = 1.0/2.0*(x**2 + c);
    return f

# Function implements the iterative calculation process and do the plotting in a single graph
def iterative(x0, iteNumber, pointType, legend):
    x = x0
    xlist = np.linspace(0,20,21)
    ylist = []
    for k in np.arange(21):
        ylist.append(x)
        x = prob2(x)
    print x
    plt.plot(xlist, ylist, pointType, label=legend)
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.legend(loc = "upper right")
    plt.title('Plot on function iteration with c = 0.5')


# part a (x0 = 0)
iterative(0.0, 101, 'k.', 'x0=0')
# part b (x0 = 1.83)
iterative(1.83, 101, 'r*', 'x0=1.83')
# part c (x0 = -1)
iterative(-1.0, 101, 'b^', 'x0=-1')
plt.show()
# part d (x0 = 1.84)
# The iterative process does not converge.
#iterative(1.84, 101, 3)