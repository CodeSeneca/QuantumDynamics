#! /usr/bin/python3

# Copyright: Maximilian Bechtel

"""Main program to drive everything"""

import time
import sys
import numpy as np
from pkg.IO import read_input, write_output
from pkg.functions import harmonic_potential, morse_potential, gaussian
from pkg.simulation import calc_b
from pkg.les import solve_les, calc_norm, calc_d, calc_epot, calc_ekin, calc_p, calc_x

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
dt, nsteps, dx, ngridpoints, x0, p0, sigma, k, alpha, box, mass, potential, wavefunction, \
output_mode, output_step = read_input(input_filename)

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
""", end='')

if output_mode == 0:
    print(f"""
No output will be written
""", end='')

elif output_mode == 1:
    print(f"""
Wave function and energies will be written out
""", end='')

if potential == 1:
  print(f"""
A particle in a box will be simulated
box heigt = {box}
""", end='')

elif potential == 2:
  print(f"""
A harmonic potential will be simulated
Force constant for harmonic potential = {k}
""", end='')

elif potential == 3:
  print(f"""
A Morse potential will be simulated
Stiffness of potential alpha = {alpha}
""", end='')

elif potential == 4:
  print(f"""
A potential wall will be simulated
potential wall height = {box}
""", end='')

if wavefunction == 1:
  print(f"""
A gaussian wave packet was chosen
-------------------------------------------------
""")

elif wavefunction == 2:
  print(f"""
A sinus function was chosen (eigenfunction of
particle in the box
-------------------------------------------------
""")

# File handlers for writing output
if output_mode != 0:
  output_name = "energy.dat"
  plot_name = "plot.dat"
  potential_name = "potential.dat"
  output_file = open(output_name, 'w')
  plot_file = open(plot_name, 'w')
  potential_file = open(potential_name, 'w')
  output_file.write("#Time step Norm Potential Energy Kinetic Energy " \
                + "Total Energy momentum p location x\n")
  potential_file.write("#Potential\n")

###############################################################################
###################### Initializations for t = 0 ##############################
###############################################################################

# Create an array with the x-values of the grid
# The grid points are equally distributed around x0
print(" Creating simulation grid ...")
x_values = np.linspace(x0 - 0.5*ngridpoints*dx, x0 + 0.5*ngridpoints*dx, \
ngridpoints)

# Create an array with the values of the potential on the grid
print(" Creating potential ...")
# Particle in the box
if potential == 1:
  v_values = np.zeros(ngridpoints)
  v_values[0] = box
  v_values[-1] = box
# Harmonic potential
elif potential == 2:
  v_values = harmonic_potential(x_values, k, x0)
# Morse potential
elif potential == 3:
  v_values = morse_potential(x_values, alpha, x0)
# Potential wall
elif potential == 4:
  v_values = np.zeros(ngridpoints)
  start = int(0.70*ngridpoints)
  end = int(0.80*ngridpoints)
  v_values [start:end] = box

# Generate start configuration of psi
print(" Creating initial wave function ...")

if wavefunction == 1:
  psi = gaussian(x_values, x0, sigma, p0)
elif wavefunction == 2:
  b = np.pi/(2*(x0 - 0.5*ngridpoints*dx))
  psi = np.cos(b * x_values)

norm = calc_norm(psi, dx)
psi *= norm

# Write out start configuration
if output_mode != 0:
  # Write out potential
  for i in range(ngridpoints):
    potential_file.write(f"{x_values[i]:.4f}    {v_values[i]:.5f}\n")

  # Write out expectation values
  epot = calc_epot(psi, v_values, dx)
  ekin = calc_ekin(psi, dx, mass)
  etot = ekin + epot
  p = calc_p(psi)
  x = calc_x(psi, x_values, dx)
  write_output(0, plot_file, output_file, psi, x_values, dx, dt, epot, \
               ekin, etot, p, x)

###############################################################################
############################### Main Loop #####################################
###############################################################################

d = np.zeros(ngridpoints, dtype=complex)         # right hand side vector

print("\nEntering main loop ...")

# Loop over all time steps
start = time.time()
for i in range(1, nsteps+1):
  ##### STEP 1: With current psi calculate new vector d and overwritten vector b
  b = calc_b(v_values, ngridpoints, dx, dt, mass)
  d = calc_d(v_values, psi, dt, dx, mass)
  ##### STEP 2: Solve LES with Thomas algorithm -> new psi
  psi = solve_les(b, d)
  norm = calc_norm(psi, dx)
  psi *= norm
  ##### STEP 3: Calculate energies and write output for each nth step
  if output_mode != 0 and i%output_step == 0:
    epot = calc_epot(psi, v_values, dx)
    ekin = calc_ekin(psi, dx, mass)
    etot = ekin + epot
    p = calc_p(psi)
    x = calc_x(psi, x_values, dx)
    write_output(i, plot_file, output_file, psi, x_values, dx, dt, epot, \
                 ekin, etot, p, x)

  perc = 100 * i/nsteps
  if perc%10 == 0:
    print(f" ... {int(perc)} % done")

end = time.time()
diff = end - start
print(f"Finished in {diff:0.3f} s")

if output_mode != 0:
  output_file.close()
  plot_file.close()
  potential_file.close()
  print("\n Written Potential to", potential_name, "for plotting ...")
  print(" Written wave function to", plot_name, "for plotting ...")
  print(" Written energies to", output_name, "...")
