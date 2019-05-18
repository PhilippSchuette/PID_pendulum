#!/bin/python
#
# This file is intended for debugging of the classes and methods in PID_control
# only!
#
# Author: Philipp Schuette
# License: GPL-3.0
# Date: 14/05/2019
import numpy as np

from PID_control import AnimatedPendulum, Pendulum

if __name__ == '__main__':
    # Set PID control parameters:
    ALPHA = 4.4
    BETA = 2.0
    MU = 1.2
    MAX_CONTROL = 2.6
    FREQUENCY = 30
    DEADBAND = 0.01
    SET_POINT = -0.0 * np.pi
    PRECISION = 5

    t_start = 0.0
    t_end = 45.0
    N = 9000
    LENGTH = 0.1

    # Perturbation could perhaps be randomized; something like (0.01*np.pi)
    # seems to be a good value for this particular parameter set.
    PERTURBATION = 0.0 * np.pi

    f1 = np.sin

    def f2(x):
        return(x + PERTURBATION)

    # The following is an unused prototype for a perturbed nonlinear pendulum:
    def f3(x):
        return(np.sin(x) + PERTURBATION)

    phi0 = 0.5 * np.pi
    phi0_dot = 0.3 * np.pi

    # After specifying all necessary data, the Pendulum class solves the
    # ODE within three statements:  creation of an appropriate Pendulum
    # instance, a call to the solve() method and a call to the plot()
    # method:
    ode1 = Pendulum(t_start, t_end, N, f1, L=LENGTH)
    ode1.solve(phi0, phi0_dot, ALPHA, BETA, MU, MAX_CONTROL, FREQUENCY,
               DEADBAND, SET_POINT, PRECISION)
    ode1.plot("nonlinearPID", parameter=True)

    ode2 = Pendulum(t_start, t_end, N, f2, L=LENGTH)
    ode2.solve(phi0, phi0_dot, ALPHA, BETA, MU, MAX_CONTROL, FREQUENCY,
               DEADBAND, SET_POINT, PRECISION)
    ode2.plot("linearPID", parameter=True)

    # Test pendulum animation features:
    animatedpendulum = AnimatedPendulum(
        phi0, phi0_dot, ALPHA, BETA, MU, MAX_CONTROL, FREQUENCY, DEADBAND,
        SET_POINT, PRECISION, t_start=0.0, t_end=6.5, N=650, L=LENGTH, f=f1
    )
    animatedpendulum.animate("pendulum1")

    pendulum1 = Pendulum(t_start=0, t_end=5, N=500, f=np.sin, L=LENGTH)
    pendulum1.solve(phi0, phi0_dot, ALPHA, BETA, MU, MAX_CONTROL, FREQUENCY,
                    DEADBAND, SET_POINT, PRECISION)
    pendulum1.animate("pendulum2")
