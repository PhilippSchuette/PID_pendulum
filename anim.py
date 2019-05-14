# This code uses the PIDControl and Pendulum classes, to animate a
# controlled pendulum. The animation is generated with matplotlib.
#
# Author: Philipp Schuette
# Date: 13.05.2019


import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from PID_control import Pendulum


# Set up formatting for the movie files:
Writer = animation.writers['ffmpeg']
writer = Writer(fps=35, bitrate=1800)

fig = plt.figure()
ax = plt.axes(xlim=(0, 20), ylim=(-2.5, 2.5))
plt.grid(True)
line, = ax.plot([], [], lw=2)

plt.xlabel('$t$')
plt.ylabel('$\phi$')


# Initialize the animation:
def init():
    line.set_data([], [])
    return line,


# Set PID control and pendulum parameters:
ALPHA = 4.4
BETA = 2.0
MU = 1.2
MAX_CONTROL = 2.6
FREQUENCY = 30
DEADBAND = 0.01
SET_POINT = -0.0 * np.pi
PRECISION = 5

t_start = 0.0
t_end = 20.0
N = 4000
LENGTH = 0.1

f1 = np.sin

phi0 = 0.6 * np.pi
phi0_dot = 0.4 * np.pi

# Initialize pendulum:
pendulum = Pendulum(t_start, t_end, N, f1, L=LENGTH)
pendulum.solve(phi0, phi0_dot, ALPHA, BETA, MU, MAX_CONTROL, FREQUENCY,
               DEADBAND, SET_POINT, PRECISION)


# Define the function, that is animated, from pendulum data:
def animate(i):
    x = pendulum.get_support_values()[:i]
    y = pendulum.get_func_values()[:i]
    line.set_data(x, y)
    return line,


anim = animation.FuncAnimation(fig, animate, init_func=init, frames=4000,
                               interval=2, blit=True)

# Saves with the formatting set up earlier:
anim.save("pics/animation.mp4", writer=writer)

plt.show()
