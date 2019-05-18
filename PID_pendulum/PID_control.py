#!/bin/python
#
# PID control for inverted pendulum. Control is implemented as a class,
# such that other systems/ODE can be controlled without much additional
# effort.
# Several phenomena present in real world control applications are also
# implemented, such as bounded maximal control, limited controller
# speed and limited access to system variables. A random, but constant
# noise perturbes the system.
#
# Further features include the possibility to switch between linearized
# and non-linear pendulum and an optimization algorithm for optimal
# control parameters with respect to multiple optimality criteria.
#
# For debugging and testing purposes, use the file debug.py!
#
# Author: Philipp Schuette
# License: GPL-3.0
# Date: 20/04/2019
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation


##################
# PID Controller #
##################
class PIDControl():
    """
    Class implementation of a PID controller.  This class is initialize
    from given data like parameters for P, I and D, controller speed
    with respect to the modelled system and a bound for controller
    ouput.  Especially the latter two are quite important in practical
    applications.  The set point, that the controller tries to reach,
    can also be adjusted manually.  A parameter deadband reduces strain
    on possible machinery behind the controller; the control output is
    only changed, if the newly calculated output differs from the
    previous one by at least the amount specified by deadband.
    """

    def __init__(self, alpha, beta, mu, frequency,
                 max_control, set_point, deadband):
        """
        Initialize PIDControl class.

        :type alpha: float > 0
        :param alpha: proportional control parameter

        :type beta: float > 0
        :param beta: derivative control parameter

        :type mu: float > 0
        :param mu: integral control parameter

        :type frequency: int >= 1
        :param frequency: controller speed parameter

        :type max_control: float > 0
        :param max_control: controller output bound

        :type set_point: float between 0 and 2*pi
        :param set_point: desired valued of controlled system

        :type deadband: float > 0
        :param deadband: minimum difference for controller adjustment

        >>> import numpy as np; from PID_control import *
        >>> ALPHA = 4.4; BETA = 2.0; MU = 1.2; MAX_CONTROL = 2.6
        >>> FREQUENCY = 30; DEADBAND = 0.01; SET_POINT = -0.0*np.pi
        >>> controller = PIDControl(
        ... ALPHA, BETA, MU, FREQUENCY, MAX_CONTROL, SET_POINT, DEADBAND
        ... )
        >>> print(controller)
        PID Controller with alpha = 4.4, beta = 2.0, mu = 1.2
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

    def __repr__(self):
        """
        Return string representation of PID controller.
        """
        return("PID Controller with alpha = {}, beta = {}, mu = {}".format(
            self.alpha, self.beta, self.mu
        ))

    def total_output(self, x1, x2, t1, t2, precision=4):
        """
        Method returning the total controller output in response to system
        value x at time t.  The output is bounded by the max_control attribute
        and controller adjustment speed is bounded by the frequency attribute.
        Controller does not adjust, if successive outputs differ by less than
        the deadband attribute. Limited measurement precision is included by
        rounding controller input.

        :type x1: float
        :param x1: system value at time t1

        :type x2: float
        :param x2: system value at time t2

        :type t1: float
        :param t1: last time point

        :type t2: float
        :param t2: current time point

        :type precision: int > 0
        :param precision: system measurement precision

        :output: total control value
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
        return(self.output)

    def proportional_output(self, x, t):
        """
        Method returning the proportional or I controller output, depending on
        the controller attribute alpha.

        :type x: float
        :param x: system value at time t

        :type t: float
        :param t: current time point

        :output: Proportional control value
        """
        control_value = -self.alpha * (x - self.set_point)

        return(control_value)

    def derivative_output(self, x1, x2, t1, t2):
        """
        Method returning the derivative or D controller output, depending
        on the attribute beta.  A numerical approximation of the
        derivative value is necessary to compute the ouput.  The
        trapezoid rule was chosen for that.

        :type x1: float
        :param x1: system value at time t1

        :type x2: float
        :param x2: system value at time t2

        :type t1: float
        :param t1: last time point

        :type t2: float
        :param t2: current time point

        :output: derivative control value
        """
        control_value = -self.beta*(x2 - x1)/(t2 - t1)

        return(control_value)

    def integral_output(self, x1, x2, t1, t2):
        """
        Method returning the integral or I controller output, depending on the
        controller attribute mu.  A numerical approximation of the integral
        value is necessary to compute the output.

        :type x1: float
        :param x1: system value at time t1

        :type x2: float
        :param x2: system value at time t2

        :type t1: float
        :param t1: last time point

        :type t2: float
        :param t2: current time point

        :output: integral control value
        """
        control_value = (self.integral - self.mu*(t2 - t1)
                         * (x1 - self.set_point + x2 - self.set_point)/2.0)

        self.integral = control_value
        return(control_value)


##################################
# Controlled Pendulum ODE Solver #
##################################
class Pendulum():
    """
    Class implementation of a simple ODE solver for the inverted pendulum.
    Intended to be used in conjunction with the PIDControl class.  The solved
    ODE is of the form x''(t) = f(x(t)) + u(t).

    Initialize with time parameters t_start, t_end, number of support
    points N.  Initial values are given to the solve method for
    flexibility in calculating solutions for different initial
    conditions.  Right hand side can be adjusted, even though sine and
    identity are the most reasonable choices.
    """

    def __init__(self, t_start, t_end, N, f, L=0.1, G=9.81):
        """
        Initialize linear or nonlinear pendulum.

        :type t_start: float
        :param t_start: starting time for ODE solving

        :type t_end: float > t_start
        :param t_end: final time for ODE solving

        :type N: int > 0
        :param N: number of numerical support points for ODE solving

        :type f: function
        :param f: right hand side for the second order pendulum ODE

        :type L: float
        :param L: pendulum length parameter in [m]

        :type G: float > 0
        :param G: gravitational constant in [m/s^2]

        :output: object representing the second order pendulum ODE with
                 right hand side f
        """
        self.t_start = t_start
        self.t_end = t_end
        self.N = N
        self.f = f
        self.L = L
        self.G = G

        # Calculate time step width:
        self.h = (self.t_end - self.t_start)/self.N
        # Define array with time support points:
        self.t = [(self.t_start + i*self.h) for i in range(self.N + 1)]

    def __repr__(self):
        """
        Return string representation of inverted pendulum.
        """
        return("Inverted Pendulum of Length {}".format(self.L))

    def solve(self, phi0, phi0_dot, alpha, beta, mu, max_control, frequency,
              deadband, set_point, precision):
        """
        Method solving the ODE for given physical initial conditions,
        i.e. initial angle and velocity, and with PID controller, that
        has the given parameters.

        :type phi0: float
        :param phi0: initial value

        :type phi0_dot: float
        :param phi0_dot: initial angular velocity

        :type alpha: float > 0
        :param alpha: proportional control parameter

        :type beta: float > 0
        :param beta: derivative control parameter

        :type mu: float > 0
        :param mu: integral control parameter

        :type max_control: float > 0
        :param max_control: controller output bound

        :type frequency: int >= 1
        :param frequency: controller speed parameter

        :type deadband: float
        :param deadband: minimum difference between calculated control
                         outputs

        :type set_point: float
        :param set_point: desired value of controlled system

        :type precision: int > 0
        :param precision: measurement precision of controller input

        :output: numerically calculates the attributes (float arrays)
                 phi, output_array, P_array, I_array and D_array
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
        for n in range(1, self.N):
            # Calculate the current control value:
            u_n = self.controller.total_output(self.phi[n-1], self.phi[n],
                                               self.t[n-1], self.t[n],
                                               precision=self.precision)

            # Save control values for control value plot:
            self.output_array.append(u_n)
            self.P_array.append(
                self.controller.proportional_output(self.phi[n], self.t[n])
            )
            self.D_array.append(
                self.controller.derivative_output(
                    self.phi[n-1], self.phi[n], self.t[n-1], self.t[n]
                )
            )
            self.I_array.append(self.controller.integral)

            # After this calculation, `tmp' contains the angle value phi[n+1]:
            tmp = (
                2.0*self.phi[n] + self.L*self.G*self.f(self.phi[n])*self.h**2
                - self.phi[n-1] + u_n*self.h**2
            )
            self.phi.append(tmp)

    def solve_from_angles(self, phi0, phi1, alpha, beta, mu, max_control,
                          frequency, deadband, set_point, precision):
        """
        Method solving the ODE for given numerical initial conditions,
        i.e. two initial angles.  Works exactly like the solve() method.

        :type phi0: float
        :param phi0: first initial value

        :type phi1: float
        :param phi1: second initial value

        :type alpha: float > 0
        :param alpha: proportional control parameter

        :type beta: float > 0
        :param beta: derivative control parameter

        :type mu: float > 0
        :param mu: integral control parameter

        :type max_control: float > 0
        :param max_control: controller output bound

        :type frequency: int >= 1
        :param frequency: controller speed parameter

        :type deadband: float
        :param deadband: minimum difference between calculated control
                         outputs

        :type set_point: float
        :param set_point: desired value of controlled system

        :type precision: int > 0
        :param precision: measurement precision of controller input

        :output: Numerically calculates the attributes (float arrays)
                 phi, output_array, P_array, I_array and D_array

        >>> import numpy as np; from PID_control import *
        >>> ALPHA = 4.4; BETA = 2.0; MU = 1.2; MAX_CONTROL = 2.6
        >>> FREQUENCY = 30; DEADBAND = 0.01; SET_POINT = -0.0*np.pi
        >>> PRECISION = 5; t_start = 0.0; t_end = 45.0; N = 9000; LENGTH = 0.1
        >>> phi0 = 0.5 * np.pi; phi0_dot = 0.3 * np.pi
        >>> pendulum = Pendulum(t_start, t_end, N, np.sin, L=LENGTH)
        >>> pendulum
        Inverted Pendulum of Length 0.1
        >>> phi1 = phi0 + phi0_dot*pendulum.h
        >>> pendulum.solve(
        ... phi0, phi0_dot, ALPHA, BETA, MU, MAX_CONTROL, FREQUENCY, DEADBAND,
        ... SET_POINT, PRECISION
        ... )
        >>> x = pendulum.phi[10]
        >>> pendulum.solve_from_angles(
        ... phi0, phi1, ALPHA, BETA, MU, MAX_CONTROL, FREQUENCY, DEADBAND,
        ... SET_POINT, PRECISION
        ... )
        >>> y = pendulum.phi[10]
        >>> print(np.abs(x - y) < 5e-5)
        True
        """
        self.phi0 = phi0
        self.phi1 = phi1

        # Initialize solution array with given initial angles:
        self.phi = [self.phi0, self.phi1]
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
        for n in range(1, self.N):
            # Calculate the current control value:
            u_n = self.controller.total_output(self.phi[n-1], self.phi[n],
                                               self.t[n-1], self.t[n],
                                               precision=self.precision)

            # Save control values for control value plot:
            self.output_array.append(u_n)
            self.P_array.append(
                self.controller.proportional_output(self.phi[n], self.t[n])
            )
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

    def get_func_values(self):
        """
        A convenience method that returns the calculated function values
        or an empty list if `solve` was never called on this pendulum.

        :output: array containing angle solution values (floats)
        """
        if hasattr(self, "phi"):
            return self.phi
        else:
            return []

    def get_support_values(self):
        """
        A convenience method that returns the support values (== time
        points)or an empty list if `solve` was never called on this
        pendulum.

        :output: array containing support values (floats)
        """
        if hasattr(self, "t"):
            return self.t
        else:
            return []

    def plot(self, file_name, parameter=False):
        """
        Method to plot solutions generated with Pendulum class.  One needs to
        call the solve() method, before plot() can be called.  PID Parameters
        can be written in the filename.  This might be useful for numerical
        experiments with several distinct sets of parameters.

        :type file_name: string
        :param file_name: name for the .png file, in which the solution
                          gets stored

        :type parameter: bool
        :param parameter: if true, controller parameters are written in
                          the filename

        :todo: fix naming scheme (Linear vs Nonlinear; testing for
               (lambda x: x) == self.f is not implemented very nicely).
        """
        # Plot time dependencies of system, i.e. angle phi, and control
        # output:
        plt.figure()
        plt.plot(self.t, self.phi, "g-", label="Pendulum Angle")

        # Try to fit the correct name to the kind of pendulum
        # integrated.  This fails to make sense, if self.f is anything
        # else than np.sin or (lambda x: x)!
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
        plt.ylabel("$\phi$ ($2\pi$)", fontsize=15)  # noqa: LaTeX escape seq
        plt.tick_params(labelsize=15)

        plt.xlim(self.t_start, self.t_end)
        plt.ylim(-self.max_control - 0.05, self.max_control + 0.05)
        plt.legend(loc="upper right", shadow=True)

        plt.grid(True)

        # Catch exception, if directory pics/ does not exist:
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

    def animate(self, anim_name):
        """
        Method to animate solutions generated with Pendulum class.  One
        needs to call the solve() method, before plot() makes sense.

        :type anim_name: string
        :param anim_name: name for the mp4 file, in which the animation
                          gets stored
        """
        # Set up formatting for the movie files:
        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=35, bitrate=1800)

        fig = plt.figure()
        ax = plt.axes(xlim=(self.t_start, self.t_end), ylim=(
            -self.max_control - 0.01, self.max_control + 0.01
        ))
        plt.grid(True)
        line, = ax.plot([], [], lw=2)

        plt.xlabel('$t$')
        plt.ylabel('$\phi$')  # noqa: LaTeX escape seq

        # Initialize the animation:
        def init():
            line.set_data([], [])
            return line,

        # Define the function, that is animated, from pendulum solution data:
        def animate(i):
            x = self.t[:i]
            y = self.phi[:i]
            line.set_data(x, y)
            return line,

        anim = animation.FuncAnimation(
            fig, animate, init_func=init, frames=self.N, interval=2, blit=True
        )

        # Saves with the formatting set up earlier:
        try:
            anim.save("pics/" + anim_name + ".mp4", writer=writer)
        except FileNotFoundError:
            print("Warning: Please create directory `pics' to save animation.")
            exit(1)


################################
# Animated Controlled Pendulum #
################################
class AnimatedPendulum():
    """
    Class especially tailored to creating, solving and finally animating
    a controlled inverted pendulum. While the Pendulum class method
    animate() can be called after solving the pendulum ODE with solve(),
    the AnimatedPendulum class wraps these steps. Simply create an
    instance and call the animate() method.
    """

    def __init__(self, phi0, phi0_dot, alpha, beta, mu, max_control, frequency,
                 deadband, set_point, precision, t_start, t_end, N, L, f):
        """
        :type phi0: float
        :param phi0: initial angle value

        :type phi0_dot: float
        :param phi0_dot: initial angle velocity value

        :type alpha: float > 0
        :param alpha: proportional control parameter

        :type beta: float > 0
        :param beta: derivative control parameter

        :type mu: float > 0
        :param mu: integral control parameter

        :type max_control: float > 0
        :param max_control: controller output bound

        :type frequency: int >= 1
        :param frequency: controller speed parameter

        :type deadband: float
        :param deadband: minimum difference between calculated control
                         outputs

        :type set_point: float
        :param set_point: desired value of controlled system

        :type precision: int > 0
        :param precision: measurement precision of controller input

        :type t_start: float
        :param t_start: starting time for pendulum dynamics

        :type t_end: float > t_start
        :param t_end: ending time for pendulum dynamics

        :type N: int > 0
        :param N: number of support points

        :type L: float > 0
        :param L: pendulum length

        :type f: function
        :param f: right hand side for pendulum ODE

        :output: Creates an instance of a pendulum, waiting to be
                 animated
        """
        self.phi0 = phi0
        self.phi0_dot = phi0_dot
        self.alpha = alpha
        self.beta = beta
        self.mu = mu
        self.max_control = max_control
        self.frequency = frequency
        self.deadband = deadband
        self.set_point = set_point
        self.precision = precision
        self.t_start = t_start
        self.t_end = t_end
        self.N = N
        self.L = L
        self.f = f

    def animate(self, anim_name):
        """
        Animates an instance of AnimatedPendulum.

        :type anim_name: string
        :param anim_name: name for the mp4 file, in which the animation
                          gets stored
        """
        # Set up formatting for the movie files:
        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=35, bitrate=1800)

        fig = plt.figure()
        ax = plt.axes(xlim=(self.t_start, self.t_end), ylim=(
            -self.max_control - 0.01, self.max_control + 0.01
        ))
        plt.grid(True)
        line, = ax.plot([], [], lw=2)

        plt.xlabel('$t$')
        plt.ylabel('$\phi$')  # noqa: LaTeX escape seq

        # Initialize the animation:
        def init():
            line.set_data([], [])
            return line,

        # Initialize pendulum, this differentiations AnimatedPendulum
        # method animate() from Pendulum method animate():
        pendulum = Pendulum(self.t_start, self.t_end, self.N, self.f, self.L)
        pendulum.solve(
            self.phi0, self.phi0_dot, self.alpha, self.beta, self.mu,
            self.max_control, self.frequency, self.deadband, self.set_point,
            self.precision
        )

        # Define the function, that is animated, from pendulum data:
        def animate(i):
            x = pendulum.get_support_values()[:i]
            y = pendulum.get_func_values()[:i]
            line.set_data(x, y)
            return line,

        anim = animation.FuncAnimation(
            fig, animate, init_func=init, frames=self.N, interval=2, blit=True
        )

        # Saves with the formatting set up earlier:
        try:
            anim.save("pics/" + anim_name + ".mp4", writer=writer)
        except FileNotFoundError:
            print("Warning: Please create directory `pics' to save animation.")
            exit(1)


###############################################################################
# :todo: Implement and investigate (random) noise; use PERTURBATION parameter
#        or maybe even a stochastic function of time
# :todo: Implement optimality criteria (at least speed of first achieving the
#        set point and overswing width around the set point).
# :todo: Implement functionality for variable set point.
# :todo: Implement integral anti-windup (deactivate integrator, if difference
#        between system value and set point is too small)?
# :todo: Implement (in ODE) friction.
# :todo: Raise exceptions, when incorrect parameters are supplied!
# :todo: Improve doctest coverage!
###############################################################################
