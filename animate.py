#! /usr/bin/python3

"""Skript to animate all time steps written by quantum_dynamics.py"""

import sys
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

filename = "plot.dat"

# Delay between each frame (ms)
if len(sys.argv) != 2:
  print("")
  print(" No suitable number of arguments given")
  print(" Usage: animate.py [delay between frames (ms)]")
  sys.exit(1)
delay = float(sys.argv[1])

fig = plt.figure("Quantum Dynamics")
ax = fig.add_subplot(1,1,1)

# plot objects that should be updated in each frame
psi, = ax.plot([], [], 'r')
time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)

print("Reading in", filename, "...")
# Read in generated values from file plot.dat
data = np.genfromtxt(filename, comments='#')

# Determine number of x values = number of gridpoints
plot_file = open(filename, 'r')
plot_file.readline()
ngridpoints = 0
for line in plot_file:
  if re.search('#', line):
    break
  ngridpoints += 1
ngridpoints -= 2

# Determine time steps (frames) that should be displayed
plot_file.seek(0)
plot_file.readline()
frames = []
for line in plot_file:
  if re.search('#', line):
    val = float(line.rstrip().strip('#'))
    frames.append(val)

dt = frames[1] - frames[0]

plot_file.close()
print("done")

x = data[0:ngridpoints,0]
potential = data[0:ngridpoints,1]

###############################################################################
############################### Animation #####################################
###############################################################################
def update(i):
  # Set up new wave function for ith frame
  y = data[(i+1)*ngridpoints:(i+2)*ngridpoints,1]

  time_text.set_text(f"time = {i*dt:.2f}")
  psi.set_data(x,y)

  return psi, time_text

anim = FuncAnimation(fig, update, frames=len(frames), blit=True, \
                     interval=delay, repeat_delay=1000)

ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$|\Psi(x)|^2$')
ax.set_ylim(0.0, 1.5)
ax.grid()

ax2 = ax.twinx()
ax2.plot(x, potential, 'b')
ax2.set_ylabel(r'$V(x)$')
ax2.set_ylim(0.0, potential[-1] + 1.0)

plt.show()
