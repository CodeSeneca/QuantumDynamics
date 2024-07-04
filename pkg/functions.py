# Copyright: Maximilian Bechtel

"""Definition of all needed mathematical functions"""

import numpy as np

def harmonic_potential(x: float, k: float, x0:float) -> float:
  """Harmonic potential

  x                  location of Gaussian
  k                  harmonic force constant
  x0                 start location of initial Gaussian
  """

  potential = 0.5*k*(x - x0)**2
  return potential

def morse_potential(x: float, alpha:float, x0:float) -> float:
  """ Morse potential

  x                   location of Gaussian
  alpha               stifness of potential
  x0                  start location of initial Gaussian
  """

  morse = (1 - np.exp(-alpha*(x-x0)))**2

  return morse

def gaussian(x:float, x0:float, sigma:float, p0:float) -> complex:
  """Gaussian wave packet

  h = 1 (for a.u.)
  x0                 start location of initial gaussian
  p0                 start momentum of initial gaussian
  sigma              standard deviation of initial gaussian
  """

  gauss = np.exp(-0.5*((x-x0)/sigma)**2) * np.exp(1j*p0*x)
  return gauss
