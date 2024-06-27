#! /usr/bin/python3
"""Main program to drive everything"""

import time
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from pkg.input import read_input, write_output
from pkg.functions import harmonic_potential, morse_potential, gaussian
from pkg.simulation import calc_norm, calc_b, calc_Epot, calc_Ekin
from pkg.les import solve_les, calc_d, calc_ekin

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
dt, nsteps, dx, ngridpoints, x0, p0, sigma, k, mass, potential, output_mode, \
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

Results will be written for each {output_step}th time step
""")

if output_mode == 0:
    print(f"""
No output will be written
""")

elif output_mode == 1:
    print(f"""
Wave function and energies will be written out
""")

if potential == 1:
  print(f"""
A free particle will be simulated
-------------------------------------------------
""")

elif potential == 2:
  print(f"""
A harmonic potential will be simulated
Force constant for harmonic potential = {k}
-------------------------------------------------
""")

elif potential == 3:
  print(f"""
A Morse potential will be simulated
-------------------------------------------------
""")

# File handlers for writing output
if output_mode != 0:
  output_name = "energy.dat"
  plot_name = "plot.dat"
  output_file = open(output_name, 'w')
  plot_file = open(plot_name, 'w')
  output_file.write("#Time step Norm Potential Energy Kinetic Energy " \
                + "Total Energy\n")

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
# Free particle
if potential == 1:
  v_values = np.zeros(ngridpoints)
# Harmonic potential
elif potential == 2:
  v_values = harmonic_potential(x_values, k, x0)
# Morse potential
elif potential == 3:
  v_values = morse_potential(x_values, 0.5, x0)
elif potential == 4:
  v_values = np.zeros(ngridpoints)
  v_values[ngridpoints//2 + 200:ngridpoints//2 + 300] = 0.5
print("done")

# Generate start configuration of psi (gaussian wave packet)
print("\n Creating initial Gaussian wave packet ...", end='')
#psi = gaussian(x_values, x0, sigma, p0)
b = np.pi/(2*(x0 - 0.5*ngridpoints*dx))
psi = np.cos(b * x_values)
norm = calc_norm(psi, dx)
psi *= norm
print("done")

# Write out start configuration
# First block: Potential
if output_mode != 0:
  plot_file.write("#Potential\n")
  for i in range(len(x_values)):
    plot_file.write(f"{x_values[i]:.4f}    {v_values[i]:.5f}\n")
  plot_file.write("\n\n")
  # Next blocks: |Psi|^2
  epot = calc_Epot(psi, v_values, dx)
  ekin = calc_Ekin(psi, dx, mass)
  etot = ekin + epot
  write_output(0, plot_file, output_file, psi, x_values, dx, dt, epot, \
               ekin, etot)

###############################################################################
############################### Main Loop #####################################
###############################################################################

# Allocate arrays a,b,c,d for Thomas algorithm
# ngridpoints = number of equations in LES
#a = np.ones(ngridpoints - 1)                    # subdiagonal
#c = np.ones(ngridpoints - 1)                    # supradiagonal
d = np.zeros(ngridpoints, dtype=complex)         # right hand side vector

print("\nEntering main loop ...")

# Loop over all time steps
start = time.time()
for i in range(1, nsteps+1):
  # print(f"Loop #{i}")
  ##### STEP 1: With current psi calculate new vector d and overwritten vector b
  b = calc_b(v_values, ngridpoints, dx, dt, mass)
  d = calc_d(v_values, psi, dt, dx, mass)
  ##### STEP 2: Solve LES with Thomas algorithm -> new psi
  psi = solve_les(b, d)
  # print(type(psi[0]))
  ##### STEP 3: Calculate energies and write output for each nth step
  if output_mode != 0 and i%output_step == 0:
    epot = calc_Epot(psi, v_values, dx)
    #ekin = calc_Ekin(psi, dx, mass)
    ekin = calc_ekin(psi, dx, mass)
    etot = ekin + epot
    write_output(i, plot_file, output_file, psi, x_values, dx, dt, epot, \
                 ekin, etot)

  perc = 100 * i/nsteps
  if perc%10 == 0:
    print(f" ... {int(perc)} % done")

end = time.time()
diff = end - start
print(f"Finished in {diff:0.3f} s")

if output_mode != 0:
  output_file.close()
  plot_file.close()
  print("\n Written Norm, Total Energy, Kinetic Energy and Potential Energy",\
        "to", output_name, "...")
  print("\n Written wave function to", plot_name, "for plotting ...")
