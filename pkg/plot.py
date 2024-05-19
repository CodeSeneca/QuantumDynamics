"""Plot the given potential and wavefunction on the grid"""

import matplotlib.pyplot as plt

def plot(x_values:'ndarray', v_values:'ndarray', psi:'ndarray') -> None:
  """Plot everything

  x_values         values of the grid
  v_values         values of the potential on the grid
  psi              normalized values of Gaussian on the grid
                   plotted values: abs(psi(x))**2
  """

  fig = plt.figure("Quantum Dynamics")
  ax = fig.add_subplot(1,1,1)
  ax.set_xlabel("x")
  ax.set_ylim(bottom=0, top=1.5*max(abs(psi)**2))

  ax.plot(x_values, v_values, label="V(x)")
  ax.plot(x_values, abs(psi)**2, label=r"$|\psi|^2$(x)")
  #ax.plot(x_values, psi.real, label="Re($\psi$(x))")
  ax.legend()
  #ax.grid()

  plt.show()
