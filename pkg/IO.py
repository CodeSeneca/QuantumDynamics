import re
import sys
import numpy as np

"""Definition of simulation parameters and I/O
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
alpha                   stiffness of the Morse potential
box                     height of the box for particle in a box
mass                    particle mass
potential               1 = particle in the box, 2 = harmonic, 3 = Morse potential
                        4 = potential wall
wavefunction            1 = gaussian wave packet, 2 = sinus (eigenfucntion
                        of particle in box)
output_mode             extent of output written
output_step             write out energy for only each nth time step
"""

def read_input(filename:str) -> list:
  """Read simulation parameters from input file"""

  # Default values
  dt = 0.005
  nsteps = 100
  dx = 0.01
  ngridpoints = 1000
  x0 = 0.0
  p0 = 1.0
  sigma = 1.0
  k = 5.0
  alpha = 0.5
  box = 1000
  mass = 1.0
  potential = 2
  wavefunction = 1
  output_mode = 1
  output_step = 10

  try:
    input_file = open(filename, 'r')
  except:
    print("File", filename, "could not be found!")
    sys.exit(2)

# Change default parameters with values from input file
  for line in input_file:
    if re.search("=", line):
      line = line.rstrip().split("=")
      param, val = line
      param = param.rstrip()
      if param == "dt":
        dt = complex(val)
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
      if param == "alpha":
        alpha = float(val)
      if param == "box":
        box = float(val)
      if param == "mass":
        mass = float(val)
      if param == "potential":
        potential = int(val)
      if param == "wavefunction":
        wavefunction = int(val)
        print("Wavefunction: ", wavefunction)
      if param == "output_mode":
        output_mode = int(val)
      if param == "output_step":
        output_step = int(val)

  input_file.close()

  return dt, nsteps, dx, ngridpoints, x0, p0, sigma, k, alpha, box, mass, potential, \
         wavefunction, output_mode, output_step

def write_output(step, plot_file, output_file, psi, x_values, dx, dt, epot, ekin, etot, p, x):
  """Write output to files plot_file and output_file"""

  if dt.real == 0.0:
    dt = dt.imag
  elif dt.imag == 0.0:
    dt = dt.real

  # Calculate |Psi|^2 and norm for writing
  psi_2 = np.abs(psi)**2
  norm = np.sum(psi_2*dx)

  plot_file.write(f"#{step*dt:.4f}\n")
  for i in range(len(x_values)):
    plot_file.write(f" {x_values[i]:.4f}    {psi_2[i]:.5f}\n")
  plot_file.write("\n\n")

  output_file.write(f"{step*dt:.5f}    {norm:.5f}    {epot:.5f}    {ekin:.5f}    {etot:.5f}    {p:.5f}    {x:.5f}\n")
