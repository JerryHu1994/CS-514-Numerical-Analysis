# Math/CS 514 HW5
# Jieru Hu  ID:9070194544
import scipy.integrate as integrate
import numpy as np
from matplotlib import pyplot as plt

#define the function
def f(x):
    return np.exp(-x)*np.sin(10*x)

# the integration function of the function f(x)
def fintegration(x):
    return -10/float(101)*np.exp(-x)*np.cos(10*x) - 1/float(101)*np.exp(-x)*np.sin(10*x)

# define the Composite Trapezoidal method
def CompTrapRule(a,b,n):
    h = (b-a)/float(n)
    arr = np.linspace(a,b,n)
    sum = 2*np.sum(f(arr))
    sum = sum - f(a) - f(b)
    return sum*h/2

# define the Composite Simpson method
def CompSimpRule(a,b,n):
    h = (b - a)/float(n)
    arr = np.linspace(a,b,n)
    sum = 2*np.sum(f(arr)) - f(a) - f(b)
    oddarr = np.arange(a+h, b, 2*h)
    sum += 2*np.sum(f(oddarr))
    return sum*h/3

# define the Gaussian quadrature method
def Gauss(a,b,n):
    I = integrate.fixed_quad(f,a=a,b=b,n=n)[0]
    return I


# Question 0
# Plot the function
x0 = np.linspace(-1, 1, 100)
fig1 = plt.figure(0)
plt.plot(x0, f(x0))
plt.xlabel('x')
plt.ylabel('f(x)')
plt.title("Question 0: Plot of f(x) on [-1, 1]")
plt.show()

# Integration by hand
# The process is written with the paper after the homework problems.
Istar = fintegration(1) - fintegration(-1)
print "Istar is ",Istar
# Question 1
# CompTrapRule(a,b,n) as above

# Question 2
# CompSimpRule(a,b,n) as above

# Question 3
# Gauss(a,b,n) as above

# Question 4
# 4a) Composite trapezoidal rule with n=100 points.
print "The result and relative error calculated by Composite trapezoidal rule with 100 points is :", CompTrapRule(-1,1,100), (CompTrapRule(-1,1,100)-Istar)/Istar
# 4b) Composite trapezoidal rule with n=200 points.
print "The result and relative error calculated by Composite trapezoidal rule with 200 points is :", CompTrapRule(-1,1,200), (CompTrapRule(-1,1,200)-Istar)/Istar
# 4c) Composite Simpson's rule with n=100 points.
print "The result and relative error calculated by Composite Simpson rule with 100 points is :", CompSimpRule(-1,1,100), (CompSimpRule(-1,1,100)-Istar)/Istar
# 4d) Composite Simpson's rule with n=200 points.
print "The result and relative error calculated by Composite Simpson rule with 200 points is :", CompSimpRule(-1,1,200), (CompSimpRule(-1,1,200)-Istar)/Istar
# 4e) Gaussian quadrature with just n=10 points.
print "The result and relative error calculated by Gaussian quadrature rule with 10 points is :", Gauss(-1,1,10), (Gauss(-1,1,10)-Istar)/Istar
# 4f) Gaussian quadrature with just n=20 points.
print "The result and relative error calculated by Gaussian quadrature rule with 20 points is :", Gauss(-1,1,20), (Gauss(-1,1,20)-Istar)/Istar
# Question 5
print "Composite trapezoidal rule: The reduction ratio in error when doubling n is : ", np.abs(Istar - CompTrapRule(-1,1,100))/np.abs(Istar - CompTrapRule(-1,1,200))
print "Composite Simpson rule: The reduction ratio in error when doubling n is : ", np.abs(Istar - CompSimpRule(-1,1,100))/np.abs(Istar - CompSimpRule(-1,1,200))
print "Gaussian quadrature rule: The reduction ratio in error when doubling n is : ", np.abs(Istar - Gauss(-1,1,10))/np.abs(Istar - Gauss(-1,1,20))

# Question 6
nList = np.linspace(1,20,20)
logList = []
for i in range(1,21,1):
    logList.append(np.log10(np.abs((Gauss(-1,1,i) - Istar)/Istar)))
fig2 = plt.figure(1)
plt.plot(nList, logList)
plt.xlim([0,21])
plt.xlabel('x')
plt.ylabel('Log10 (Error)')
plt.title("The Log10 (Error) vs number of Nodes using Gaussian quadrature for Integration")
plt.show()