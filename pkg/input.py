import re
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
mass                    particle mass
output_mode             extent of output written
output_step             write out energy for only each nth time step
"""


def read_input(filename:str) -> list:
  """Read simulation parameters from input file

  filename           name of the input file
  return             list of simulation parameters
                     -> if parameter is not given the default value will be
                        used
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

  return dt, nsteps, dx, ngridpoints, x0, p0, sigma, k, mass, output_mode, \
         output_step
