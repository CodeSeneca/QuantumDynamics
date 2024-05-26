"""Definition of simulation related functions"""

from numpy import sqrt, zeros, conjugate

def calc_norm(psi:'complex ndarray', dx:float) -> float:
  """Calculate the norm of given gaussian wave packet

  psi               numpy array with complex values of psi on grid
  dx                grid spacing
  return norm       norm of gaussian wave packet -> float
  """

  # Create an array with the complex conjugate elements of psi
  psi_conj = conjugate(psi)
  norm = sum(psi*psi_conj*dx)
  
  # This is the vectorized version of:

  #norm = 0.0
  #for i in range(len(psi)):
    #norm += (psi[i].real**2 + psi[i].imag**2) * dx

  norm = 1/sqrt(norm)

  return norm

def calc_b(V:'ndarray', n:int, dx:float, dt:float, mass:float) -> 'ndarray':
  """Calculate the elements bi of the main diagonal for Thomas algorithm

  V                 array with elements Vi of the potential on the grid
  n                 number of equations in LES = ngridpoints
  dx                grid spacing
  dt                time step
  mass              particle mass
  return b          array with elements of main diagonal
                    -> complex elements
  """

  b = zeros(n, dtype=complex)
  b = mass*dx**2 * (2j/dt - V) - 1

  return 2*b

def calc_d(d:'complex ndarray', V:'ndarray', psi:'complex ndarray', dt, dx) -> None:
  """Calculate the right hand side vector d for Thomas algorithm

  d                 array with coefficients di -> will be overwritten
                    -> indexed from [0,...,n-1]
                    -> complex elements
  potential         array with potential on the grid
  psi               array with current wave function on the grid
  dt                time step
  dx                grid spacing

  BOUNDARY CONDITIONS:
  d_0:   no d_i-1
  d_n-1: no d_i+1
  """

  n = len(d)
  # Loop over elements [1,n-2]
  for i in range(1,n-1):
    d[i] = -psi[i+1] + 2*( dx**2 * (2j/dt + V[i]) + 1 )*psi[i] \
           -psi[i-1]

  # Most left element
  d[0] = -psi[i+1] + 2*( dx**2 * (2j/dt + V[i]) + 1 )*psi[i] \

  # Most right element
  d[n-1] = 2*( dx**2 * (2j/dt + V[i]) + 1 )*psi[i] - psi[i-1]

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
  c_trans = zeros(n-1, dtype=complex)
  d_trans = zeros(n, dtype=complex)
  c_trans[0] = c[0]/b[0]
  d_trans[0] = d[0]/b[0]

  # Loop from 1 to n-2
  #for i in range(1, n-1):

  # Loop from 1 to n-2 for c_trans and from 1 to n-1 for d_trans
  for i in range(1, n):
    if i < n - 1:
      c_trans[i] = c[i] / (b[i] - a[i-1]*c_trans[i-1])

    d_trans[i] = ( d[i] - a[i-1]*d_trans[i-1] ) \
               / ( b[i] - a[i-1]*c_trans[i-1] )

  ##### STEP 2: Find the coefficients xi by back substitution
  x[n-1] = d_trans[n-1]
  for i in range(n-2,-1,-1):
    x[i] = d_trans[i] - c_trans[i]*x[i+1]
