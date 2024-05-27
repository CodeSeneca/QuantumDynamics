#! /mingw64/bin/python

"""Main program to drive everything"""

import time
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from pkg.functions import harmonic_potential, gaussian
from pkg.simulation import calc_norm, calc_b, calc_d, solve_les_thomas
from pkg.input import read_input

###############################################################################
############################## Input parameters ###############################
###############################################################################

# Get input file name via command line argument
argc = len(sys.argv)
if argc != 2:
  print("\n No suitable number of arguments given")
  print(" Usage: quantum_dynamics.py <input name>")
  sys.exit(1)
input_filename = sys.argv[1]

# Set simulation parameters
dt, nsteps, dx, ngridpoints, x0, p0, sigma, k, mass, output_mode, \
output_step = read_input(input_filename)

print("""
Project 5: Simulation of Quantum Dynamics
=========================================
Author: Maximilian Bechtel <maxi.bechtel@fau.de>
Project of the course Scientific Programming at FAU Erlangen-Nuernberg
""")

print("Input will be read from", input_filename)

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

# File handlers for writing output
output_name = "energy.dat"
plot_name = "plot.dat"
output_file = open(output_name, 'w')
plot_file = open(plot_name, 'w')
output_file.write("#Time step Norm Total Energy Kinetic Energy " \
                + "Potential Energy\n")

###############################################################################
###################### Initializations for t = 0 ##############################
###############################################################################

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

# Generate start configuration of psi (gaussian wave packet)
print("\n Creating initial Gaussian wave packet ...", end='')
psi = gaussian(x_values, x0, sigma, p0)
norm = calc_norm(psi, dx)
psi *= norm
print("done")

# Write out start configuration
plot_file.write("#Potential\n")
for i in range(len(x_values)):
  plot_file.write(f"{x_values[i]:.4f}    {v_values[i]:.5f}\n")
plot_file.write("\n\n")

psi_2 = np.abs(psi)**2
norm = np.sum(psi_2*dx)
output_file.write(f"0.0000 {norm:.5f}\n")
plot_file.write(f"#0.0000\n")
for i in range(len(x_values)):
  plot_file.write(f" {x_values[i]:.4f}    {psi_2[i]:.5f}\n")
plot_file.write("\n\n")

###############################################################################
############################### Main Loop #####################################
###############################################################################

# Allocate arrays a,b,c,d for Thomas algorithm
# ngridpoints = number of equations in LES
a = np.ones(ngridpoints - 1)                    # subdiagonal
c = np.ones(ngridpoints - 1)                    # supradiagonal
b = calc_b(v_values, ngridpoints, dx, dt,mass)  # main diagonal
d = np.zeros(ngridpoints, dtype=complex)        # right hand side vector

print("\nEntering main loop ...")

fig, ax = plt.subplots()
ax.grid()
ax.set_xlabel("x")
ax.set_ylim([0, 1.5])

# Loop over all time steps
ims = []
start = time.time()
for i in range(1, nsteps+1):
  ##### STEP 1: With current psi calculate new vector d
  calc_d(d, v_values, psi, dt, dx)
  ##### STEP 2: Solve LES with Thomas algorithm -> new psi
  solve_les_thomas(a, b, c, d, psi)
  ##### STEP 3: Write output
  psi_2 = np.abs(psi)**2
  norm = np.sum(psi_2*dx)

  im = ax.plot(x_values, psi_2, 'b')
  ims.append(im)

  plot_file.write(f"#{i*dt:.4f}\n")
  for j in range(len(x_values)):
    plot_file.write(f" {x_values[j]:.4f}    {psi_2[j]:.5f}\n")
  plot_file.write("\n\n")

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
plt.plot(x_values, v_values)
plt.show()
