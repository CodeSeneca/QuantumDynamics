# Copyright: Maximilian Bechtel

"""Definition of simulation related functions

slow implementation in Python
faster subroutines are used written in FORTRAN -> file les.f90
by using f2py for compilation
"""

import numpy as np

def calc_norm(psi:'complex ndarray', dx:float) -> float:
  """Calculate the norm of given gaussian wave packet"""

  psi_2 = np.abs(psi)**2
  norm = np.sum(psi_2*dx)

  norm = 1/np.sqrt(norm)

  return norm

def calc_b(V:'ndarray', n:int, dx:float, dt:float, mass:float) -> 'ndarray':
  """Calculate the elements bi of the main diagonal for Thomas algorithm"""

  b = np.zeros(n, dtype=complex)
  b = mass*dx**2 * (2j/dt - V) - 1

  return 2*b

def calc_d(V:'ndarray', psi:'complex ndarray', dt, dx, mass) -> 'ndarray':
  """Calculate the right hand side vector d for Thomas algorithm

  BOUNDARY CONDITIONS:
  d_0:   no d_i-1
  d_n-1: no d_i+1
  """

  n = len(psi)
  d = np.empty(n, dtype=complex)

  # Loop over elements [1,n-2]
  for i in range(1,n-1):
    d[i] = -psi[i+1] + 2*( mass*dx**2 * (2j/dt + V[i]) + 1 )*psi[i] \
           -psi[i-1]

  # Most left element
  d[0] = -psi[1] + 2*( mass*dx**2 * (2j/dt + V[0]) + 1 )*psi[0]

  # Most right element
  d[-1] = 2*( mass*dx**2 * (2j/dt + V[-1]) + 1 )*psi[-1] - psi[-2]

  return d

##############################
##### FAST IMPLEMENTATION
##############################

def solve_les(b, d):
  """Solve the Linear Equation System (LES) by using the Thomas algorithm

  return psi       solution vector of LES = wavefunction on grid
  CAUTION: b and d are overwritten

  REFERENCE: https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
  """

  n = len(d)
  psi = np.empty(n, dtype=complex)

  for i in range(1,n):
    omega = 1/b[i-1]
    b[i] = b[i] - omega
    d[i] = d[i] - omega*d[i-1]

  # Back substitution
  psi[-1] = d[-1]/b[-1]
  for i in range(n-2,-1,-1):
    psi[i] = ( d[i] - psi[i+1] ) / b[i]

  return psi

################################
##### VERY SLOW IMPLEMENTATION
################################

def solve_les_thomas(a:'real', b:'complex', c:'real', d:'complex', x) -> None:
  """Solve a Linear Equation System (LES) by using the Thomas algorithm

  Ax = d with A being a tridiagonal matrix

  a[]              subdiagonal, indexed from [0,...,n-2]
  b[]              main diagonal, indexed from [0,...,n-1]
  c[]              supradiagonal, indexed from [0,...,n-2]
  d[]              right hand side vector, indexed from [0,...,n-1]
  SOLUTION VECTOR x will be overwritten -> complex array = psi
  type: ndarray for all vectors (dytpe = complex)

  REFERENCE: https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
  """

  n = len(d)  # Number of equations in LES

  ##### STEP 1: Transform vector c to c_trans and vector d to d_trans
  c_trans = np.zeros(n-1, dtype=complex)
  d_trans = np.zeros(n, dtype=complex)
  c_trans[0] = c[0]/b[0]
  d_trans[0] = d[0]/b[0]

  # Loop from 1 to n-2 for c_trans and from 1 to n-1 for d_trans
  for i in range(1, n-1):
    c_trans[i] = c[i] / (b[i] - a[i-1]*c_trans[i-1])
  #c_trans[1:] = c[1:] / ( b[1:-1] - a[:-1]*c_trans[:-1] )
 
  for i in range(1, n):
    d_trans[i] = ( d[i] - a[i-1]*d_trans[i-1] ) \
               / ( b[i] - a[i-1]*c_trans[i-1] )
  #d_trans[1:] = ( d[1:] - a[:]*d_trans[:-1] ) / (b[1:] - a[:]*c_trans[:] )

  ##### STEP 2: Find the coefficients xi by back substitution
  x[-1] = d_trans[-1]
  for i in range(n-2,-1,-1):
    # print(i, end=' ')
    x[i] = d_trans[i] - c_trans[i]*x[i+1]

def calc_Epot(psi, V, dx) -> float:
  """Calculate <V(t)>: the expectation value of the potential energy"""

  e_pot = np.sum(np.abs(psi)**2 * V) * dx

  return e_pot

def calc_Ekin(psi, dx, m) -> float:
  """Calculate <Ekin(t)>"""

  n = len(psi)
  # Calculate second derivative of psi = psi_diff
  psi_diff = np.empty(n, dtype=complex)

  # Middle elements [2, n-2]
  for i in range(1,n-1):
    psi_diff[i] = psi[i+1] - 2*psi[i] + psi[i-1]
  # Elements at the left and right border
  psi_diff[0] = psi[1] - 2*psi[0]
  psi_diff[-1] = -2*psi[-1] + psi[-2]

  psi_diff /= dx**2

  e_kin = np.sum(np.conjugate(psi) * psi_diff) * dx
  e_kin = -e_kin/(2*m)

  return e_kin.real
