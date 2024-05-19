#! /mingw64/bin/python

import numpy as np
from pkg.input import *
from pkg.functions import harmonic_potential, gaussian
from pkg.simulation import norm, calc_b, calc_lescoeff, solve_les_thomas
from pkg.plot import plot

print("""
Project 5: Simulation of Quantum Dynamics
=========================================
Author: Maximilian Bechtel
Project of the course Scientific Programming at FAU Erlangen-Nuernberg
""")

print(f"""Simulation parameters:
Atomic units (a.u.) will be used
dt = {dt}
nsteps = {nsteps}
dx = {dx}
ngridpoints = {ngridpoints}
x0 = {x0}
p0 = {p0}
sigma = {sigma}
mass = {mass}""")

# x-values of the grid
x_values = np.linspace(x0 - 0.5*ngridpoints*dx, x0 + 0.5*ngridpoints*dx, \
ngridpoints)
print("\nCreated simulation grid")

# values of the potential
v_values = harmonic_potential(x_values, k, x0)
print("Created harmonic potential")
# Generate start configuration of psi (t = 0)
psi = gaussian(x_values, x0, sigma, p0)
norm = norm(psi, dx)
psi *= norm
print("Created initial Gaussian wave packet (psi at t = 0)")

# Test normalization
norm = 0.0
for i in range(len(psi)):
  norm += abs(psi[i])**2 * dx
print(f"\nNormalization: {norm}")

# Allocate arrays a,b,c,d for Thomas algorithm
# ngridpoints = number of equations in LES
a = np.ones(ngridpoints - 1)               # subdiagonal
c = np.ones(ngridpoints - 1)               # supradiagonal
b = calc_b(v_values, ngridpoints, dx, dt)  # main diagonal
d = np.zeros(ngridpoints, dtype=complex)   # right hand side vector

###############################################################################
############################### Main Loop #####################################
###############################################################################

print("\nEntering main loop ...")

# Loop over all time steps
t_start = 0.0
t_end = t_start + (nsteps+1)*dt
time_steps = np.arange(t_start + dt, t_end, dt)
for t in time_steps:
  ##### STEP 1: With current psi calculate new vector d
  calc_lescoeff(d, v_values, psi, dt, dx)
  ##### STEP 2: Solve LES with Thomas algorithm -> new psi
  solve_les_thomas(a, b, c, d, psi)

# Plot the harmonic potential and initial Gaussian
#plot(x_values, v_values, psi)
