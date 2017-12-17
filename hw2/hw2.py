# Math/CS 514
# Jieru Hu  ID:9070194544
# -*- coding: utf-8 -*-
import numpy as np
from matplotlib import pyplot as plt

# global variables
tol = 1e-12
printouts = False
iteLimit = 200

# Define the function
def f(x):
    return np.sin(np.power(x,2)) + 1.02 - np.exp(-x)

# Define the derivative
def f_prime(x):
    return np.cos(np.power(x,2))*2*x + np.exp(-x)

# Define the Newton's method
def Newton(x0):
    print 'Start running Newton method...'
    xk = x0
    residual = abs(f(xk))
    xk_list=[xk]
    residual_list=[residual]


    while residual>tol:

        xk = xk - f(xk)/f_prime(xk)
        residual = abs(f(xk))

        xk_list.append(xk)
        residual_list.append(residual)

        if printouts:
            print xk, residual

        if len(xk_list)==iteLimit:
            print 'Did not converge after 200 iterations'
            break

    return xk_list, residual_list

# Define the secant method
def Secant(x0,x1):
    print 'Start running Secant method...'
    xkm1 = x0
    xk = x1

    residual = abs(f(xk))
    xk_list=[xk]
    residual_list=[residual]

    fkm1 = abs(f(xkm1))

    while residual>tol:

        fk = f(xk)
        xkp1 = xk - (xk-xkm1)/(fk-fkm1) * fk

        xkm1 = xk
        xk = xkp1
        fkm1 = np.copy(fk)

        residual = abs(fk)

        if printouts:
            print xk, residual

        xk_list.append(xk)
        residual_list.append(residual)

        if len(xk_list)==iteLimit:
            print 'Did not converge after 200 iterations'
            break

    return xk_list, residual_list

# Define the bisection method, starting from [xleft, xright]
def Bisection(x0, x1):
    print 'Start running Bisection method...'
    x_l = x0
    x_r = x1
    xk_list = []
    residual_list = []
    residual = float("inf")

    # check if the initial guess satisfy the f(x) = 0
    residual_l = abs(f(x_l))
    residual_r = abs(f(x_r))
    if residual_l <= tol:
        return [x_l], [residual_l]
    if residual_r <= tol:
        return [x_r], [residual_r]


    while residual>tol:
        xmid = (x_l + x_r)/2.0
        fk = f(xmid)
        residual = abs(fk)

        if printouts:
            print xmid, residual

        xk_list.append(xmid)
        residual_list.append(residual)

        if f(xmid)*f(x_l) > 0:
            x_l = xmid
        else:
            x_r = xmid

        if(len(xk_list) == iteLimit):
            print 'Did not converge after 200 iterations'
            break

    return xk_list, residual_list

# Define the Illinois algorithm, with starting point x0, x1
def Illinois(x0, x1):
    print 'Start running Illinois method...'
    xkm1 = x0
    xk = x1
    fkm1 = f(xkm1)

    residual = abs(f(xk))

    # check if the initial guess satisfy the f(x) = 0
    residual_m1 = abs(f(xkm1))
    residual = abs(f(xk))
    if residual_m1 <= tol:
        return [xkm1], [residual_m1]
    if residual <= tol:
        return [xk], [residual]

    xk_list = []
    residual_list = []

    while residual>tol:
        # calculate x_k+1
        xkp1 = xk - (xk-xkm1)/(f(xk)-fkm1)*f(xk)
        residual = abs(f(xkp1))

        if printouts:
            print xkp1, residual

        xk_list.append(xkp1)
        residual_list.append(residual)

        # update
        if f(xkp1) == 0:
            break
        if f(xkp1)*f(xk) < 0:
            xkm1, fkm1 = xk, f(xk)
        if f(xkp1)*f(xk) > 0:
            xkm1, fkm1 = xkm1, 1/2*fkm1

        if (len(xk_list) == iteLimit):
            print 'Did not converge after 200 iterations'
            break
        xk = xkp1

    return xk_list, residual_list
################################################

######## Running / plotting ##############
# plot iterative process for the function
def plotIteration(xlist, legend):
    # first plot the function
    x = np.linspace(-3, 3, 400)
    fig1 = plt.figure(1)
    plt.plot(x, f(x))
    plt.plot(x, np.zeros_like(x), 'k-', label='Original Function')
    plt.xlabel('x')
    plt.ylabel('f(x)')

    # plot the  iteration
    fig1 = plt.figure(1)
    plt.plot(np.array(xlist), f(np.array(xlist)), 'o', label=legend)
    plt.legend(loc="lower right")
    title = "Plot on the function iteration with " + legend + " method"
    plt.title(title)
    plt.show()

def plotResidual(residual_list, legend):
    fig = plt.figure(2)
    plt.plot(np.log10(np.array(residual_list)), 'o-')
    plt.xlabel('k')
    plt.ylabel('$\log_{10}(r_k)$')
    title = "Residual plot for " + legend + " method"
    plt.title(title)
    plt.show()

# run all four methods
# First Root (x1 = -0.0202026931)

xk_list_newton, residual_list_newton = Newton(x0=0.2)
print 'Newton method number of iterations = ', len(xk_list_newton)
xk_list_Secant, residual_list_Secant = Secant(x0=0.2, x1 =0.3)
print 'Secant method number of iterations = ', len(xk_list_Secant)
xk_list_bisect, residual_list_bisect = Bisection(x0=-0.2, x1 =0.2)
print 'Bisect method number of iterations = ', len(xk_list_bisect)
xk_list_illinois, residual_list_illinois = Illinois(x0=-0.2, x1 =0.2)
print 'Illinois method number of iterations = ', len(xk_list_illinois)
print xk_list_illinois
print residual_list_illinois

# Second Root (x2 = 2.0602482004)
# xk_list_newton, residual_list_newton = Newton(x0=2.0)
# print 'Newton method number of iterations = ', len(xk_list_newton)
# xk_list_Secant, residual_list_Secant = Secant(x0=2.0, x1 =1.9)
# print 'Secant method number of iterations = ', len(xk_list_Secant)
# xk_list_bisect, residual_list_bisect = Bisection(x0=2.0, x1 =2.1)
# print 'Bisect method number of iterations = ', len(xk_list_bisect)
# xk_list_illinois, residual_list_illinois = Illinois(x0=2.0, x1 =2.1)
# print 'Illinois method number of iterations = ', len(xk_list_illinois)

# Third Root (x3 = 2.26386025285)
# xk_list_newton, residual_list_newton = Newton(x0=2.4)
# print 'Newton method number of iterations = ', len(xk_list_newton)
# xk_list_Secant, residual_list_Secant = Secant(x0=2.4, x1 =2.5)
# print 'Secant method number of iterations = ', len(xk_list_Secant)
# xk_list_bisect, residual_list_bisect = Bisection(x0=2.2, x1 =2.4)
# print 'Bisect method number of iterations = ', len(xk_list_bisect)
# xk_list_illinois, residual_list_illinois = Illinois(x0=2.2, x1 =2.4)
# print 'Illinois method number of iterations = ', len(xk_list_illinois)

Plotting = True
if Plotting:

    # plot the iterations for different methods
    plotIteration(xk_list_newton, 'Newton')
    plotResidual(residual_list_newton, 'Newton')
    plotIteration(xk_list_Secant, 'Secant')
    plotResidual(residual_list_Secant, 'Secant')
    plotIteration(xk_list_bisect, 'Bisect')
    plotResidual(residual_list_bisect, 'Bisect')
    plotIteration(xk_list_illinois, 'Illinois')
    plotResidual(residual_list_illinois, 'Illinois')

# Analysis: Taking the root 1 as an example, the Newton's method only takes 5 steps, and has the fastest convergence to the root starting from a initial guess
# of 0.2, with quadratic convergence rate. The Secant method is a little bit slower than the Newton, which takes 8 iterations to converge. Because we are
# approximating the derivative by computing a difference quotient over two close points. However, even though the Newton's method takes less steps, it involves
# computing the derivative of the f'(x) which may sometimes take substantial amount of time. So overall, there is some tradeoff in between choosing these two methods.
# The bisect method takes 36 to find the root while the Illinois method only takes 17 steps to find the root, which shows us that
# the Illinois method is more efficient. This is because in the Illinois method, instead of setting the Xk+1 to a point between
# two boundary each time, the method set (Xk-1, fk-1) to (Xk-1, 1/2*fk-1) when fk+1fk > 0. In this case, the convergence for the next
# iteration would be faster, so that the overall convergence performance of Illinois method is better.
# The residual plots I made agrees with the discussion on the convergence rate above, for each three roots.
# The convergence on other two roots have the similar behavior on these four methods.