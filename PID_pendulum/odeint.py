#!/bin/python
#
# Author: Philipp Schuette
# Date: 09/04/2019
import matplotlib.pyplot as plt
import numpy as np
#import sympy as sym

from scipy.integrate import odeint


# Set PID control parameters and maximal control value:
alpha = 4.4
beta = 2.0
mu = 0.9
MAX_CONTROL = 2.0

# Define the control/input function.  The maximal control output is
# bounded by MAX_CONTROL and the relative application frequency of the
# control is bounded by the variable frequency with standard value 2:
def control(X, t, frequency=1):
    control_value = -alpha*X[1] - beta*X[2] - mu*X[0]

    if (int(100*t) % frequency) == 0:
        if (control_value <= MAX_CONTROL):
            if (control_value >= -MAX_CONTROL):
                return(control_value)
            else:
                return(-MAX_CONTROL)
        else:
            return(MAX_CONTROL)
    else:
        return(0.0)


# Define the pendulum ODE:
def pend_ODE(X, t):
    x = X[0]
    y = X[1]
    z = X[2]

    dxdt = y
    dydt = z
    dzdt = y + control(X, t)
    return([dxdt, dydt, dzdt])


# Set up initial values and solve the ODE:
X0 = [0.0, 0.25*np.pi, 0.45]
t = np.linspace(0, 30, 3000)
solution = odeint(pend_ODE, X0, t)
x = solution[:, 0]
y = solution[:, 1]
z = solution[:, 2]

control_output = []
for i in range(0, 3000):
    control_output.append(-control([x[i], y[i], z[i]], t[i]))


# Set up solution plotting in the phase plane (y-z-plane):
plt.figure(1)
plt.subplot(211)
plt.plot(y, z, "g-")
plt.plot(X0[1], X0[2], "r.")

plt.title("PID controlled Inverted Pendulum")
plt.xlabel("$\phi$", fontsize=15)
plt.ylabel("$d\phi/dt$", fontsize=15)
plt.tick_params(labelsize=15)

plt.xlim(-1.0, 1.0)
plt.ylim(-1.2, 1.2)

plt.grid(True)

# Plot vector field of free inverted pendulum:
Y, Z = np.mgrid[-1.0:1.0:25j, -1.2:1.2:20j]
u = 0.0*Y + 1.0*Z
v = 1.0*Y + 0.0*Z
plt.quiver(Y, Z, u, v, color='b')

# Set up plotting of control values over time:
plt.subplot(212)
plt.plot(t, mu*x, "b--", label="Integral Control")
plt.plot(t, alpha*y, "y--", label="Proportional Control")
plt.plot(t, beta*z, "r--", label="Derivative Control")
plt.plot(t, control_output, "g--", label="Total Control Output")
plt.title("Control Values")
plt.xlabel("t", fontsize=15)
plt.ylabel("Controls", fontsize=15)
plt.tick_params(labelsize=15)

plt.xlim(0, 30)
plt.ylim(-3.5, 3.5)

plt.grid(True)
plt.legend(loc='upper right', shadow=True)
plt.subplots_adjust(hspace=0.3)

# Show plots:
plt.show()

# Plot time-dependent pendulum angle:
plt.figure()
plt.plot(t, y, "g-")

plt.title("PID controlled Inverted Pendulum")
plt.xlabel("$t$", fontsize=15)
plt.ylabel("$\phi$", fontsize=15)
plt.tick_params(labelsize=15)

plt.xlim(0.0, 30.0)
plt.ylim(-1.0, 1.0)
plt.grid(True)

# Show plots:
plt.show()


if __name__ == '__main__':
	print("\n" + "0")