# Math/CS 514 HW3
# Jieru Hu  ID:9070194544
# -*- coding: utf-8 -*-
import numpy as np

# define the function
def f(x):
    return np.sin(np.exp(x))


# define the Composite Simpson method
def compositeSimpson(h, xstart, xend):
    n = (xend - xstart)/h
    sum = f(xstart) + f(xend)
    for i in range(1,int(n)):
        if i % 2:
            sum += 4*f(xstart+i*h)
        else:
            sum += 2*f(xstart+i*h)
    return sum*h/3


# Call the Composite Trapezoidal method
In = compositeSimpson(0.1, 0.0, 2.0)
I2n = compositeSimpson(0.05, 0.0, 2.0)
I4n = compositeSimpson(0.025, 0.0, 2.0)
print 'The integral using composite simpson rule with h = 1/10 is : ' + str(In)
print 'The integral using composite simpson rule with h = 1/20 is : ' + str(I2n)
print 'The integral using composite simpson rule with h = 1/40 is : ' + str(I4n)
print 'The ratio R using composite simpson rule is : ' + str((In - I2n)/(I2n - I4n))
