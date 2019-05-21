#!/bin/python
#
# Author: Philipp Schuette
# Date: 09/04/2019
import matplotlib.pyplot as plt
import numpy as np
#import sympy as sym

from scipy.integrate import odeint


# Define the linearized pendulum ODE:
def pend_ODE1(X, t):
    x = X[0]
    y = X[1]

    dxdt = y
    dydt = x
    return([dxdt, dydt])

# Define the non-linear pendulum ODE:
def pend_ODE2(X, t):
    x = X[0]
    y = X[1]

    dxdt = y
    dydt = np.sin(x)
    return([dxdt, dydt])


# Set up initial values and solve the ODEs:
X0 = [0.25*np.pi, 0.45]
t = np.linspace(0, 30, 3000)
solution1 = odeint(pend_ODE1, X0, t)
x1 = solution1[:, 0] % (2*np.pi)
y1 = solution1[:, 1]
solution2 = odeint(pend_ODE2, X0, t)
x2 = solution2[:, 0]
y2 = solution2[:, 1]

# Set up solution plotting in the phase plane (x-y-plane):
plt.figure(1)
plt.subplot(111)
plt.plot(x1, y1, "g-", label="Linearized Pendulum")
plt.plot(x2, y2, "r-", label="Nonlinear Pendulum")
plt.plot(X0[0], X0[1], "r.")

plt.title("PID controlled Inverted Pendulum")
plt.xlabel("$\phi$", fontsize=15)
plt.ylabel("$d\phi/dt$", fontsize=15)
plt.tick_params(labelsize=15)

plt.xlim(0.0, 6.0)
plt.ylim(-2.0, 6.0)

plt.legend(loc="upper left", shadow=True)
plt.grid(True)

# Plot vector field of free inverted pendulum:
Y, Z = np.mgrid[0.0:6.0:30j, -2.0:6.0:40j]
u = 0.0*Y + 1.0*Z
v = 1.0*Y + 0.0*Z
plt.quiver(Y, Z, u, v, color='b')

# Show plots:
plt.show()

# Plot time-dependent pendulum angle for both cases:
plt.figure(1)
plt.subplot(211)
plt.plot(t, x1, "g-")

plt.title("PID controlled Inverted Pendulum")
plt.xlabel("$t$", fontsize=15)
plt.ylabel("$\phi$", fontsize=15)
plt.tick_params(labelsize=15)

plt.xlim(0.0, 30.0)
plt.ylim(0.0, 2*np.pi)

plt.grid(True)

plt.subplot(212)
plt.plot(t, x2, "r-")

plt.xlabel("$t$", fontsize=15)
plt.ylabel("$\phi$", fontsize=15)
plt.tick_params(labelsize=15)

plt.xlim(0.0, 30.0)
plt.ylim(0.0, 2*np.pi)

plt.grid(True)
plt.subplots_adjust(hspace=0.3)

# Show plots:
plt.show()


if __name__ == '__main__':
	print("\n" + "0")