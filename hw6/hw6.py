# Math/CS 514 HW6
# Jieru Hu  ID:9070194544
import numpy as np
from matplotlib import pyplot as plt
import math

fiveDecimal = 0.00001
eulerNum = 2.718281828259
lowerLimit = 743300
upperLimit = 800000
MAXITE = 10000000

# derivative of y
def f(y):
    return y

# forward Euler method
def forward_Euler(N):
    yzero = 1
    yn = yzero
    h = 1/float(N)
    for i in range(1,N+1):
        ynp1 = yn + h*f(yn)
        yn = ynp1
    return yn

# improved euler method
def improved_Euler(N):
    yzero = 1
    yn = yzero
    h = 1/float(N)

    for i in range(1, N + 1):
        ynp1 = yn + 0.5 * h * (f(yn) + f(yn + h*f(yn)))
        yn = ynp1
    return yn

# fourth-order- Runge-Kutta method
def RK4(N):
    yzero = 1
    yn = yzero
    h = 1 / float(N)

    for i in range(1, N + 1):
        k1 = f(yn)
        k2 = f(yn + 0.5*h*k1)
        k3 = f(yn + 0.5*h*k2)
        k4 = f(yn + h*k3)
        ynp1 = yn + 1/float(6) * h * (k1 + 2*k2 + 2*k3 + k4)
        yn = ynp1
    return yn

def fiveDecimalAccurate(num1, num2):
    return math.floor(num1*100000) == math.floor(num2*100000)

# Question 1: Forward Euler Method
print "Problem 1)"
for i in range(743320,800000):
    if fiveDecimalAccurate(eulerNum, forward_Euler(i)):
        print "The Forward Euler Method requires", i, "steps to achieve five decimal accurate. Calculated e:", "{0:.12f}".format(
            forward_Euler(i)), "(rounded to 12 decimal digits)"
        break

print "Problem 2)"
# Question 2: Improved Euler Method
for i in range(400, 700000):
    if fiveDecimalAccurate(eulerNum, improved_Euler(i)):
        print "The Improved Euler Method requires", i, "steps to achieve five decimal accurate. Calculated e:", "{0:.12f}".format(
            improved_Euler(i)), "(rounded to 12 decimal digits)"
        break

print "Problem 3)"
# Question 3: RK4 method
for i in range(10, 500):
    if fiveDecimalAccurate(eulerNum, RK4(i)):
        print "The RK4 Method requires", i, "steps to achieve five decimal accurate. Calculated e:", "{0:.12f}".format(
            RK4(i)), "(rounded to 12 decimal digits)"
        break

# Question 4
print "Problem 4)"
V0 = 25
theta = math.pi/4
g = 9.8

def RK4_Soccer(N,a, full):
    x0 ,y0, u0, v0= 0, 0, V0 * math.sin(theta), V0 * math.cos(theta)
    h = 2 / float(N)
    x,y,u,v = x0 ,y0, u0, v0
    xlist, ylist = [x], [y]
    end = MAXITE if full else N
    for i in range(1, end + 1):
        k1x, k1y, k1u, k1v = x_prime(u), y_prime(v), u_prime(u,v,a), v_prime(u,v,a)
        k2x, k2y, k2u, k2v = x_prime(u+0.5*h*k1x), y_prime(v+0.5*h*k1y), u_prime(u+0.5*h*k1u,v,a), v_prime(u,v+0.5*h*k1v,a)
        k3x, k3y, k3u, k3v = x_prime(u+0.5*h*k2x), y_prime(v+0.5*h*k2y), u_prime(u+0.5*h*k2u,v,a), v_prime(u,v+0.5*h*k2v,a)
        k4x, k4y, k4u, k4v = x_prime(u+h*k3x), y_prime(v+h*k3y), u_prime(u+h*k3u,v,a), v_prime(u,v+h*k3v,a)
        xp1, yp1, up1, vp1 = computeRK4(k1x, k2x, k3x, k4x, h, x), computeRK4(k1y, k2y, k3y, k4y, h, y),computeRK4(k1u, k2u, k3u, k4u, h, u),computeRK4(k1v, k2v, k3v, k4v, h, v)
        x, y, u, v = xp1, yp1, up1, vp1
        if y < 0:
            break
        xlist.append(x)
        ylist.append(y)
    return xlist, ylist

# Given k1, k2, k3, k4 for RK4 method, and returns y_plus1
def computeRK4(k1, k2, k3, k4, h, yn):
    return yn + 1/float(6)*h*(k1+2*k2+2*k3+k4)
# x-velocity
def x_prime(u):
    return u
# y-velocity
def y_prime(v):
    return v
# x-acceleration
def u_prime(u,v,a):
    absv = math.sqrt(math.pow(u,2) + math.pow(v,2))
    return -a*absv*u
# y-acceleration
def v_prime(u,v,a):
    absv = math.sqrt(math.pow(u, 2) + math.pow(v, 2))
    return - g - a*absv*v

# Part a)
N = 100

x2N, y2N= RK4_Soccer(N, 0.024, False)
x22N, y22N = RK4_Soccer(2*N, 0.024, False)
x24N, y24N = RK4_Soccer(4*N, 0.024, False)
print "Question a)"
print "The R ratio is : ", (x2N[-1]-x22N[-1])/(x22N[-1]-x24N[-1])
print "The expected R ratio is 16."

# Part b)
# Added the y check in the RK4_Soccer() method

# Part c)
Nc = int(2/0.001)
xc, yc = RK4_Soccer(Nc, 0.024, True)
fig1 = plt.figure(1)
plt.plot(xc, yc, color="black",label='a = 0.024')
plt.title("Trajectory of the ball simulation")
plt.xlabel("x-position (unit:m)")
plt.ylabel("y-position (unit:m)")

# Part d)
xd, yd = RK4_Soccer(Nc, 0.0, True)
plt.plot(xd, yd, '--', color="blue",label='a = 0')
plt.legend(loc="best")
plt.show()

# Part e)
max_height_c = max(yc)
horizontal_distance_c = xc[-1] - xc[0]
max_height_d = max(yd)
horizontal_distance_d = xd[-1] - xd[0]
print "Question e)"
print "Part c: The maximum value of height attained", max_height_c, "m\nThe total horizontal distance travvelled before the ball hits the ground is", horizontal_distance_c,"m"
print "Part d: The maximum value of height attained", max_height_d, "m\nThe total horizontal distance travvelled before the ball hits the ground is", horizontal_distance_d,"m"
