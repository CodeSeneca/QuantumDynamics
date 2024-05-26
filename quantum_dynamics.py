#! /mingw64/bin/python

"""Main program to drive everything"""

import time
import sys
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from pkg.functions import harmonic_potential, gaussian
from pkg.simulation import calc_norm, calc_b, calc_d, solve_les_thomas

###############################################################################
############################## Input parameters ###############################
###############################################################################

"""Definition of simulation parameters
Atomic units (a.u.) will be used in the following

dt                      time step
nsteps                  number of time steps (int)
                        -> [dt, nsteps*dt]
dx                      grid spacing (float)
ngridpoints             total number of grid points (int)
x0                      start location of initial Gaussain
p0                      start momentum of initial Gaussian
sigma                   standard deviation of initial Gaussain
k                       harmonic force constant
mass                    particle mass
output_mode             extent of output written
output_step             write out energy for only each nth time step
"""

# Default values
dt = 0.005
nsteps = 100
dx = 0.01
ngridpoints = 1000
x0 = 0.0
p0 = 1.0
sigma = 1.0
k = 5.0
mass = 1.0
output_mode = 1
output_step = 10

# Get input file name via command line argument
argc = len(sys.argv)
if argc != 2:
  print("\n No suitable number of arguments given")
  print(" Usage: quantum_dynamics.py <input name>")
  sys.exit(1)
filename = sys.argv[1]

try:
  input_file = open(filename, 'r')
except:
  print("File", filename, "could not be found!")
  sys.exit(2)

for line in input_file:
  if re.search("=", line):
    line = line.rstrip().split("=")
    param, val = line
    param = param.rstrip()
    if param == "dt":
      dt = float(val)
    if param == "nsteps":
      nsteps = int(val)
    if param == "dx":
      dx = float(val)
    if param == "ngridpoints":
      ngridpoints = int(val)
    if param == "x0":
      x0 = float(val)
    if param == "p0":
      p0 = float(val)
    if param == "sigma":
      sigma = float(val)
    if param == "k":
      k = float(val)
    if param == "mass":
      mass = float(val)
    if param == "output_mode":
      output_mode = int(val)
    if param == "output_step":
      output_step = int(val)

input_file.close()

print("""
Project 5: Simulation of Quantum Dynamics
=========================================
Author: Maximilian Bechtel <maxi.bechtel@fau.de>
Project of the course Scientific Programming at FAU Erlangen-Nuernberg
""")

print("Input will be read from", filename)

print(f"""
-------------SIMULATION PARAMETERS---------------
Atomic units (a.u.) will be used

Time step dt = {dt}
Number of time steps nsteps = {nsteps}
Grid spacing dx = {dx}
Number of grid points ngridpoints = {ngridpoints}
Start location of particle x0 = {x0}
Start momentum of particle p0 = {p0}
Standard deviation of Gaussian sigma = {sigma}
Particle mass = {mass}

Force constant for harmonic potential = {k}

Energies will be written for each {output_step}th time step
-------------------------------------------------
""")

# File handler for writing output
output_name = "energy.dat"
plot_name = "plot.dat"
output_file = open(output_name, 'w')
plot_file = open(plot_name, 'w')
output_file.write("#Time step Norm Total Energy Kinetic Energy " \
                + "Potential Energy\n")

# Create an array with the x-values of the grid
# The grid points (ngridpoints) are equally distributed around x0
print(" Creating simulation grid ...", end='')
x_values = np.linspace(x0 - 0.5*ngridpoints*dx, x0 + 0.5*ngridpoints*dx, \
ngridpoints)
print("done")

# Create an array with the values of the potential on the grid
print("\n Creating potential ...", end='')
v_values = harmonic_potential(x_values, k, x0)
print("done")
# Generate start configuration of psi (t = 0)
print("\n Creating initial Gaussian wave packet ...", end='')
psi = gaussian(x_values, x0, sigma, p0)
norm = calc_norm(psi, dx)
psi *= norm
print("done")

# Test if psi is really normalized
psi_2 = psi*np.conjugate(psi)
norm = sum(psi_2*dx).real
output_file.write(f"0.0000 {norm:.5f}\n")
plot_file.write(f"#0.0000\n {psi_2}\n")
#norm = 0.0
#for i in range(len(psi)):
  #norm += abs(psi[i])**2 * dx
#print(f"\nNormalization: {norm:.5f}")

# Allocate arrays a,b,c,d for Thomas algorithm
# ngridpoints = number of equations in LES
a = np.ones(ngridpoints - 1)                    # subdiagonal
c = np.ones(ngridpoints - 1)                    # supradiagonal
b = calc_b(v_values, ngridpoints, dx, dt,mass)  # main diagonal
d = np.zeros(ngridpoints, dtype=complex)        # right hand side vector

###############################################################################
############################### Main Loop #####################################
###############################################################################

print("\nEntering main loop ...")

fig, ax = plt.subplots()
ax.grid()
ax.set_xlabel("x")

# Loop over all time steps
ims = []
start = time.time()
for i in range(1, nsteps+1):
  ##### STEP 1: With current psi calculate new vector d
  calc_d(d, v_values, psi, dt, dx)
  ##### STEP 2: Solve LES with Thomas algorithm -> new psi
  solve_les_thomas(a, b, c, d, psi)
  ##### STEP 3: Write output
  psi_2 = abs(psi)**2

  im = ax.plot(x_values, psi_2, 'b')
  ims.append(im)
  norm = sum(psi_2*dx).real

  plot_file.write(f"#{i*dt:.4f}\n")
  if i%output_step == 0:
      output_file.write(f"{i*dt:.4f} {norm:.5f}\n")

end = time.time()
diff = end - start
print(f"Finished in {diff:0.3f} s")

output_file.close()
plot_file.close()
print("\n Written Norm, Total Energy, Kinetic Energy and Potential Energy to ",\
      output_name, "...")
print("\n Written wave function to", plot_name, "for plotting ...")

# Plot all created images
ani = animation.ArtistAnimation(fig, ims, interval=20, blit=True)
plt.show()
