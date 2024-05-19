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
mass                    mass of the particle (1 a.u.)
"""

dt = 0.005
nsteps = 100
dx = 0.01
ngridpoints = 1000
x0 = 0
p0 = 1
sigma = 1
k = 5
mass = 1
