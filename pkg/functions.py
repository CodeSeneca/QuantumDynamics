"""Definition of all needed mathematical functions"""

from numpy import exp

def harmonic_potential(x: float, k: float, x0:float) -> float:
  """Harmonic potential

  0.5k*(x - x0)
  x                  location of Gaussian
  k                  harmonic force constant
  x0                 start location of initial Gaussian
  """

  potential = 0.5*k*(x - x0)**2
  return potential

def gaussian(x:float, x0:float, sigma:float, p0:float) -> complex:
  """Gaussian wave packet

  exp(-0.5*((x-x0)/sigma)**2) * exp(i*p0*x/h)
  h = 1 (for a.u.)
  x0                 start location of initial gaussian
  p0                 start momentum of initial gaussian
  sigma              standard deviation of initial gaussian
  CAUTION: return value is complex
  """

  gauss = exp(-0.5*((x-x0)/sigma)**2) * exp(1j*p0*x)
  return gauss
