#! /mingw64/bin/python

# Copyright: Maximilian Bechtel

"""Skript to animate all time steps written by quantum_dynamics.py"""

import sys
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

filename = "plot.dat"
energy_file = "energy.dat"
potential_file = "potential.dat"

# Delay between each frame (ms)
if len(sys.argv) != 2:
  print("")
  print(" No suitable number of arguments given")
  print(" Usage: animate.py [delay between frames (ms)]")
  sys.exit(1)

try:
  delay = float(sys.argv[1])
except:
  print("")
  print(" No suitable delay given")
  print(" Usage: animate.py [delay between frames (ms)]")
  sys.exit(2)

fig = plt.figure("Quantum Dynamics")
ax = fig.add_subplot(1,1,1)

# plot objects that should be updated in each frame
psi, = ax.plot([], [], 'r-')
loc, = ax.plot([], [], 'go')
time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)
ekin_text = ax.text(0.22, 0.95, '', transform=ax.transAxes)
epot_text = ax.text(0.48, 0.95, '', transform=ax.transAxes)
etot_text = ax.text(0.74, 0.95, '', transform=ax.transAxes)

# Read in values from energy.dat, potential.dat, plot.dat
print("Reading in files to be plotted ...")
start = time.time()
try:
  energy_data = np.genfromtxt(energy_file, comments='#')
  potential_data = np.genfromtxt(potential_file, comments='#')
  wavefunction_data = np.genfromtxt(filename, comments='#')
except:
  print("")
  print("Either plot.dat, energy.dat or potential.dat was not found ...")
  print("Aborting ...")
  sys.exit(3)
end = time.time()
print(f"Finished in {(end - start):.3f} s")

# Extract needed values
time_steps = energy_data[:,0]
norm = energy_data[:,1]
epot = energy_data[:,2]
ekin = energy_data[:,3]
etot = energy_data[:,4]
x_loc = energy_data[:,6]

# Determine length of time step
dt = time_steps[1] - time_steps[0]
# Determine x and potential values to be plotted
ngridpoints = len(potential_data)
x = potential_data[:,0]
v = potential_data[:,1]

###############################################################################
############################### Animation #####################################
###############################################################################
def update(i):
  # Set up new wave function for ith frame
  y = wavefunction_data[i*ngridpoints:(i+1)*ngridpoints,1]

  time_text.set_text(f"time = {i*dt:.2f}")
  ekin_text.set_text(f"<Ekin> = {ekin[i]:.4f}")
  epot_text.set_text(f"<Epot> = {epot[i]:.4f}")
  etot_text.set_text(f"<Etot> = {etot[i]:.4f}")
  psi.set_data(x,y)

  y_init, y_end = ax.get_ylim()
  loc.set_data(x_loc[i], y_end/2)

  return psi, loc, time_text, ekin_text, epot_text, etot_text

anim = FuncAnimation(fig, update, frames=len(time_steps), blit=True, \
                     interval=delay, repeat_delay=1000)

ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$|\Psi(x)|^2$')
ax.set_ylim(0.0, 1.5)
ax.grid()

ax2 = ax.twinx()
ax2.plot(x, v, 'b')
ax2.set_ylabel(r'$V(x)$')
ax2.set_ylim(0.0, v[-1] + 2.0)

plt.show()
