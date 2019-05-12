# PID control for inverted pendulum.  Control is implemented as a class,
# such that other systems/ODE can be controlled without much additional
# effort.
# Several phenomena present in real world control applications are also
# implemented, such as bounded maximal control, limited controller
# speed and limited access to system variables.  A random, but constant
# noise perturbes the system.
#
# Further features include the possibility to switch between linearized
# and non-linear pendulum and an optimization algorithm for optimal
# control parameters with respect to multiple optimality criteria.
#
# Author: Philipp Schuette
# Date: 20.04.2019


import matplotlib.pyplot as plt
import numpy as np

##################
# PID Controller #
##################


class PIDControl():
    """
    Class implementation of a PID controller.
    """

    def __init__(self, alpha, beta, mu, frequency,
                 max_control, set_point, deadband):
        """
        Initialize PID controller from given data like parameters for P, I and
        D, controller speed with respect to the modelled system and a bound for
        controller ouput.  Especially the latter two are quite important in
        practical applications.  The set point, that the controller tries to
        reach, can also be adjusted manually.  A parameter deadband reduces
        strain on possible machinery behind the controller; the control output
        is only changed, if the newly calculated output differs from the
        previous one by at least the amount specified by deadband.

        INPUT:

        - alpha: Proportional control parameter (float > 0).

        - beta: Derivative control parameter (float > 0).

        - mu: Integral control parameter (float > 0).

        - frequency: Controller speed parameter (int >= 1).

        - max_control: Controller output bound (float > 0).

        - set_point: Desired value of controlled system (float).

        - deadband: Minimum difference between calculated control outputs
                    (float).

        OUTPUT:

        - Representation of a PID controller with parameters alpha, beta
          and mu.
        """
        self.alpha = alpha
        self.beta = beta
        self.mu = mu
        self.frequency = frequency
        self.max_control = max_control
        self.set_point = set_point
        self.deadband = deadband

        # Attribute to memorize the previous or latest change in controller
        # output:
        self.output = 0

        # Attribute to memorize the integral control value up to the current
        # time:
        self.integral = 0

        # Attribute to keep track of controller speed:
        self.action = 0

    def total_output(self, x1, x2, t1, t2, precision=4):
        """
        Method returning the total controller output in response to system
        value x at time t.  The output is bounded by the max_control attribute
        and controller adjustment speed is bounded by the frequency attribute.
        Controller does not adjust, if successive outputs differ by less than
        the deadband attribute. Limited measurement precision is included by
        rounding controller input.

        INPUT:

        - x1: system value at time t1 (float).

        - x2: system value at time t2 (float).

        - t1: last time point (float).

        - t2: current time point (float).

        - precision: system measurement precision (int > 0).
        """
        if (self.action == 0):
            control_value = (
                self.proportional_output(np.round(x2, precision), t2)
                + self.derivative_output(np.round(x1, precision),
                                         np.round(x2, precision), t1, t2)
                + self.integral_output(np.round(x1, precision),
                                       np.round(x2, precision), t1, t2)
            )
            if (np.abs(control_value - self.output) >= self.deadband):
                if (control_value > self.max_control):
                    self.output = self.max_control
                elif (control_value < -self.max_control):
                    self.output = -self.max_control
                else:
                    self.output = control_value

        self.action = (self.action + 1) % self.frequency

    def proportional_output(self, x, t):
        """
        Method returning the proportional or I controller output, depending on
        the controller attribute alpha.

        INPUT:

        - x: system value at time t (float).

        - t: current time point (float).
        """
        control_value = -self.alpha * (x - self.set_point)

        return(control_value)

    def derivative_output(self, x1, x2, t1, t2):
        """
        Method returning the derivative or D controller output, depending on
        the attribute beta.  A numerical approximation of the derivative value
        is necessary to compute the ouput.  The trapezoid rule was chosen for
        that.

        INPUT:

        - x1: system value at time t1 (float).

        - x2: system value at time t2 (float).

        - t1: last time point (float).

        - t2: current time point (float).
        """
        control_value = -self.beta*(x2 - x1)/(t2 - t1)

        return(control_value)

    def integral_output(self, x1, x2, t1, t2):
        """
        Method returning the integral or I controller output, depending on the
        controller attribute mu.  A numerical approximation of the integral
        value is necessary to compute the output.

        INPUT:

        - x1: system value at time t1 (float).

        - x2: system value at time t2 (float).

        - t1: last time point (float).

        - t2: current time point (float).
        """
        control_value = (self.integral - self.mu*(t2 - t1)
                         * (x1 - self.set_point + x2 - self.set_point)/2.0)

        self.integral = control_value
        return(control_value)


##################################
# Controlled Pendulum ODE Solver #
##################################

class ODESolver():
    """
    Class implementation of a simple ODE solver for the inverted pendulum.
    Intended to be used in conjunction with the PIDControl class.  The solved
    ODE is of the form x''(t) = f(x(t)) + u(t).
    """

    def __init__(self, t_start, t_end, N, f):
        """
        Initialize second order ODE solver for given time parameters t_start,
        t_end, number of support points N.  Initial values are given to the
        solve method for flexibility in calculating solutions for different
        initial conditions.

        INPUT:

        - t_start: Starting time for ODE solving (float).

        - t_end: Final time for ODE solving (float > t_start).

        - N: Number of numerical support points for ODE solving (int > 0).

        - f: Right hand side for the second order pendulum ODE (funtion).

        OUTPUT:

        - Object representing the second order ODE with right hand side f.
        """
        self.t_start = t_start
        self.t_end = t_end
        self.N = N
        self.f = f

        # Calculate time step width:
        self.h = (t_end - t_start)/N
        # Define array with time support points:
        self.t = [(self.t_start + i*self.h) for i in range(self.N + 1)]

    def solve(self, phi0, phi0_dot, alpha, beta, mu, max_control, frequency,
              deadband, set_point, precision):
        """
        Method solving the ODE for given initial conditions and with PID
        controller, that has the given parameters.

        INPUT:

        - phi0: Initial value (float).

        - phi0_dot: Initial velocity (float).

        - alpha: Proportional control parameter (float > 0).

        - beta: Derivative control parameter (float > 0).

        - mu: Integral control parameter (float > 0).

        - max_control: Controller output bound (float > 0).

        - frequency: Controller speed parameter (int >= 1).

        - deadband: Minimum difference between calculated control outputs
                    (float).

        - set_point: Desired value of controlled system (float).

        - precision: Measurement precision of controller input (int > 0).

        OUTPUT:

        - Numerically calculates the attributes (float arrays) phi,
          output_array, P_array, I_array and D_array.
        """
        self.phi0 = phi0
        self.phi0_dot = phi0_dot

        # Initialize solution array with given initial data:
        self.phi = [self.phi0, self.phi0 + self.phi0_dot * self.h]
        self.output_array = [0, 0]
        self.P_array = [0, 0]
        self.I_array = [0, 0]
        self.D_array = [0, 0]

        # Initialize PID controller:
        self.alpha = alpha
        self.beta = beta
        self.mu = mu
        self.max_control = max_control
        self.frequency = frequency
        self.deadband = deadband
        self.set_point = set_point
        self.precision = precision

        # Create an instance of a PIDControl:
        self.controller = PIDControl(self.alpha, self.beta, self.mu,
                                     self.frequency, self.max_control,
                                     self.set_point, self.deadband)

        # Integrate the controlled ODE:
        for n in range(1, N):
            # Calculate the current control value:
            self.controller.total_output(self.phi[n-1], self.phi[n],
                                         self.t[n-1], self.t[n],
                                         precision=precision)
            u_n = self.controller.output

            # Save control values for control value plot:
            self.output_array.append(u_n)
            self.P_array.append(
                self.controller.proportional_output(self.phi[n], self.t[n]))
            self.D_array.append(
                self.controller.derivative_output(
                    self.phi[n-1], self.phi[n], self.t[n-1], self.t[n]
                )
            )
            self.I_array.append(self.controller.integral)

            # After the following calculation, tmp contains the entry
            # phi[n + 1]:
            tmp = (2.0*self.phi[n] + self.f(self.phi[n])*self.h**2
                   - self.phi[n-1] + u_n*self.h**2)
            self.phi.append(tmp)

    def plot(self, file_name, parameter=False):
        """
        Method to plot solutions generated with ODESolver class. PID Parameters
        can be written in the filename.  This might be useful for numerical
        experiments with several distinct sets of parameters.

        INPUT:

        - file_name: Name for the .png file, in which the solution is stored
                     (string).

        - parameter: If True, controller parameters are written in the filename
                     (Bool).

        TODO:

        - Fix naming scheme (Linear vs Nonlinear; testing for
          (lambda x: x) == self.f is not implemented very nicely).
        """
        # Plot time dependencies of system, i.e. angle phi, and control output:
        plt.figure()
        plt.plot(self.t, self.phi, "g-", label="Pendulum Angle")

        # Try to fit the correct name to the kind of pendulum integrated.  This
        # fails to make sense, if self.f is anything else than np.sin or
        # (lambda x: x)
        if self.f == np.sin:
            name = " Nonlinear "
        elif (self.f(np.e) == np.e) and (self.f(np.sqrt(2)) == np.sqrt(2)):
            name = " Linear "
        else:
            name = " Disturbed "
        plt.plot(self.t, self.output_array, "r-", label="Controller Output")
        plt.title(
            "PID controlled" + name + "Inverted Pendulum with Parameters\n"
            "alpha = {}, beta = {}, mu = {}, max_output = {}, frequency = {}, "
            "deadband = {}, precision = {}".format(
                self.alpha, self.beta, self.mu,
                self.max_control, self.frequency,
                self.deadband, self.precision
            )
        )

        plt.xlabel("$t$ (s)", fontsize=15)
        plt.ylabel("$\phi$ ($2\pi$)", fontsize=15)
        plt.tick_params(labelsize=15)

        plt.xlim(self.t_start, self.t_end)
        plt.ylim(-self.max_control - 0.05, self.max_control + 0.05)
        plt.legend(loc="upper right", shadow=True)

        plt.grid(True)
        try:
            if parameter:
                plt.savefig("pics/" + file_name + "_{}_{}_{}_1.png".format(
                    self.alpha, self.beta, self.mu
                ))
            else:
                plt.savefig("pics/" + file_name + "_1.png")
        except FileNotFoundError:
            print("Warning: Please create directory `pics' to save plots to.")
            exit(1)

        # Plot P, I, D and total control outputs together in one diagram:
        plt.figure()
        plt.plot(self.t, self.output_array, "b--",
                 label="Total Control Output")
        plt.plot(self.t, self.P_array, "y--",
                 label="Proportional Control Output")
        plt.plot(self.t, self.I_array, "r--", label="Integral Control Output")
        plt.plot(self.t, self.D_array, "g--",
                 label="Derivative Control Output")

        plt.title("Control Values")
        plt.xlabel("$t$", fontsize=15)
        plt.ylabel("Controls", fontsize=15)

        plt.xlim(self.t_start, self.t_end)
        plt.ylim(-self.max_control - 0.05, self.max_control + 0.05)

        plt.legend(loc="upper right", shadow=True)

        plt.grid(True)
        if parameter:
            plt.savefig("pics/" + file_name
                        + "_{}_{}_{}_2.png".format(
                            self.alpha, self.beta, self.mu)
                        )
        else:
            plt.savefig("pics/" + file_name + "_2.png")
        plt.show()


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

    # Perturbation could perhaps be randomized;  something like (0.01*np.pi)
    # seems to be a good value for this particular parameter set.
    PERTURBATION = 0.0 * np.pi
    f1 = np.sin

    def f2(x): return x + PERTURBATION
    # The following is am unused prototype for a perturbed nonlinear pendulum:

    def f3(x): return np.sin(x) + PERTURBATION

    phi0 = 0.5 * np.pi
    phi0_dot = 0.3 * np.pi

    # After specifying all necessary data, the ODESolver class solves the ODE
    # within three statements:  creation of an appropriate ODESolver instance,
    # a call to the solve() method and a call to the plot() method:
    ode1 = ODESolver(t_start, t_end, N, f1)
    ode1.solve(phi0, phi0_dot, ALPHA, BETA, MU, MAX_CONTROL, FREQUENCY,
               DEADBAND, SET_POINT, PRECISION)
    ode1.plot("nonlinearPID", parameter=True)

    ode2 = ODESolver(t_start, t_end, N, f2)
    ode2.solve(phi0, phi0_dot, ALPHA, BETA, MU, MAX_CONTROL, FREQUENCY,
               DEADBAND, SET_POINT, PRECISION)
    ode2.plot("linearPID", parameter=True)


###############################################################################
# TODO: Implement and investigate (random) noise; use PERTURBATION parameter or
#       maybe even a stochastic function of time
# TODO: Implement optimality criteria (at least speed of first achieving the
#       set point and overswing width around the set point).
# TODO: Implement functionality for variable set point.
# TODO: Implement integral anti-windup (deactivate integrator, if difference
#       between system value and set point is too small)?
# TODO: Implement (in ODE) friction.
# TODO: Include physical parameters pendulum length (L), gravitational constant
#       (G)!
###############################################################################
