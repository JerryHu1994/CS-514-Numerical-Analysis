# -*- coding: utf-8 -*-
import numpy as np
from matplotlib import pyplot as plt

tol = 1e-12  # If the "residual" |f(x_k)| < tol, stop iteration
printouts = True

# Define the function, its derivative, and all the methods
def f(x):
# Basic idea, superlinear convergence, comparing Newton, Secant, and Newton2:
#    return x**2+np.sin(x)+0.1

# Convergence to a double root:
    return x**2*(1.-x)
#    return x**2
#    return x**3*np.cos(x)

# Stable two-cycle
#    return x**3-2*x+2

# Multiple roots and boundary behavior
#    return  x*(x**2-1)*(x-3)*np.exp(-0.5*(x-1)**2)

    
def f_prime(x):
#    return 2*x+np.cos(x)

     return 2*x-3*x**2
#    return 2*x
#    return 3*x**2*np.cos(x)-x**2*np.sin(x)

#    return 3*x**2-2

#    return  np.exp(-0.5*(x-1)**2)*(3+x-13*x**2+2*x**3+4*x**4-x**5)

######## Define iterative methods ##############

def Newton(x0):

    # Take an initial guess:    
    xk = x0
    
    xk_list=[xk]
    residual_list=[abs(f(xk))]
    
    residual = abs(f(xk))
    while residual>tol:
        
        xk = xk - f(xk)/f_prime(xk)
        residual = abs(f(xk))
        
        xk_list.append(xk)
        residual_list.append(residual)
    
        if printouts:
            print xk, residual
        
        if len(xk_list)==200:
            print 'Did not converge after 200 iterations'
            break

    return xk_list, residual_list
   
             
def Secant(x0,x1):
    
    # Take two initial guesses:    
    xkm1 = x0
    xk = x1
    
    xk_list=[xk]
    residual_list=[abs(f(xk))]
    
    residual = abs(f(xk))

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
        
        if len(xk_list)==200:
            print 'Did not converge after 200 iterations'
            break

    return xk_list, residual_list
    
def Relaxation(x0):

    # Relaxation parameter:
    l_ambda = .2
    
    # Take an initial guess:    
    xk = x0
    
    xk_list=[xk]
    residual_list=[abs(f(xk))]
    
    residual = abs(f(xk))
    while residual>tol:
        
        xk = xk - l_ambda*f(xk)
        residual = abs(f(xk))
        
        xk_list.append(xk)
        residual_list.append(residual)
    
        if printouts:
            print xk, residual
        
        if len(xk_list)==200:
            print 'Did not converge after 200 iterations'
            break
            
    return xk_list, residual_list
    

################################################ 


######## Running / plotting ##############

xk_list, residual_list = Newton(x0=1.01)

print 'Number of iterations = ', len(xk_list)

Plotting = True

if Plotting:
    # Plot f(x)
    x = np.linspace(-3,3,400)
    fig = plt.figure(1)
    plt.plot(x,f(x))
    plt.plot(x,np.zeros_like(x),'k-')
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.show()
            
    fig = plt.figure(1)
    plt.plot(np.array(xk_list),f(np.array(xk_list)),'o')

    fig = plt.figure(2)
    plt.plot(xk_list,'o')
    plt.xlabel('k')
    plt.ylabel('$x_k$')
    plt.show()

    fig = plt.figure(3)
    plt.plot(np.log10(np.array(residual_list)),'o-')
    plt.xlabel('k')
    plt.ylabel('$\log_{10}(r_k)$')
    plt.show()


# Other topics: stable two-cycle, complex variables, systems

#x0s = np.linspace(-4,4,1000)
#final_values = np.zeros(1000)
#
#for j, x0 in enumerate(x0s):
#    xk_list, error_list = Newton(x0)
#    final_values[j]=xk_list[-1]
#
#fig = plt.figure(10)
#plt.plot(x0s,final_values,'.')
#plt.xlabel('$x_0$')
#plt.ylabel('$ \\xi $')
#plt.show()
