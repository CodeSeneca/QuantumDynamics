#! /mingw64/bin/python

"""Skript to animate all time steps written by quantum_dynamics.py"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

fig = plt.figure("Quantum Dynamics")
ax = fig.add_subplot(1,1,1)

# plot object that should be updated in each frame
psi, = ax.plot([], [], 'ro')

# Read in generated values from file plot.dat
values = np.genfromtxt("plot.dat", comments='#')
#print(values[990:1010,:])

def update():
  return psi,

frames = np.linspace(0, nsteps*dt, dt)
anim = FuncAnimation(fig, update, frames=frames, blit=True, intervals=20)

plt.show()
