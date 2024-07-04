! All simulation related subroutines
! Much faster implementation as in Python by using f2py

!!!!! THOMAS ALGORITHM
subroutine solve_les(b, d, psi, n)
  implicit none

  ! Dummy variables
  integer :: n
  double complex, dimension(n) :: b, d, psi
  !f2py intent(inout) :: b,d
  !f2py intent(out) :: psi
  !f2py intent(hide) :: n

  ! Local variables
  double complex :: omega
  integer :: i

  do i = 2, n
    omega = (1.0d0, 0.0d0)/b(i-1)
    b(i) = b(i) - omega
    d(i) = d(i) - omega*d(i-1)
  end do

  ! Back substitution
  psi(n) = d(n)/b(n)
  do i = n-1,1,-1
    psi(i) = ( d(i) - psi(i+1) ) / b(i)
  end do
end subroutine solve_les

subroutine calc_d(V, psi, dt, dx, mass, n, d)
  implicit none

  ! Dummy variables
  integer :: n
  double precision :: dx, mass
  double complex :: dt
  double precision, dimension(n) :: V
  double complex, dimension(n) :: psi, d

  !f2py intent(hide) :: n
  !f2py intent(in) :: dt, dx, mass, V, psi
  !f2py intent(out) :: d

  ! Local variables
  integer :: i

  ! First element of d
  d(1) = -psi(2) + 2*( mass*dx**2 * ((0.0d0, 2.0d0)/dt + V(1)) + 1 )*psi(1)
  ! Last element of d
  d(n) = 2*( mass*dx**2 * ((0.0d0, 2.0d0)/dt + V(n)) + 1 )*psi(n) - psi(n-1)
  ! Loop over the rest elements [2,n-1]
  do i = 2,n-1
    d(i) = -psi(i+1) + 2*( mass*dx**2 * ((0.0d0, 2.0d0)/dt + V(i)) + 1 )*psi(i) &
         & -psi(i-1)
  end do

end subroutine

!!!!! Expectation value of Ekin
subroutine calc_ekin(psi, n, dx, m, ekin)
  implicit none

  ! Dummy variables
  double precision :: dx, m, ekin
  integer :: n
  double complex, dimension(n) :: psi

  !f2py intent(hide) :: n
  !f2py intent(in) :: dx, m, psi
  !f2py intent(out) :: ekin

  ! Local variables
  double complex, dimension(n) :: psi_2  ! Second derivative of psi
  integer :: i

  ! Calculate second derivative of psi
  ! Middle elements [2, n-1]
  do i = 2, n-1
    psi_2(i) = psi(i+1) - 2*psi(i) + psi(i-1)
  end do
  ! First and last element
  psi_2(1) = psi(2) - 2*psi(1)
  psi_2(n) = -2*psi(n) + psi(n-1)

  psi_2 = psi_2 / (dx**2)

  ekin = 0.0d0
  do i = 1, n
    ekin = ekin + real(( conjg(psi(i)) * psi_2(i) ) * dx)
  end do
  ekin = -ekin/(2.0d0*m)
end subroutine calc_ekin

!!!!! Expectation value of p
subroutine calc_p(psi, n, p_real)
  implicit none

  ! Dummy variables
  integer :: n
  double complex, dimension(n) :: psi
  double precision :: p_real

  !f2py intent(hide) :: n
  !f2py intent(in) :: psi
  !f2py intent(out) :: p_real

  ! Local variables
  double complex, dimension(n) :: dpsi
  double complex :: p
  integer :: i

  ! Calculate first derivative of psi
  do i = 2, n-1
    dpsi(i) = 0.5d0*(psi(i+1) - psi(i-1))
  end do
  ! Last element
  dpsi(n) = (psi(n) - psi(n-1))
  ! First element
  dpsi(1) = psi(2) -  psi(1)

  p = 0.0d0
  do i = 1, n
    p = p + conjg(psi(i)) * dpsi(i)
  end do
  p = p * (0.0d0, -1.0d0)
  p_real = real(p)
end subroutine calc_p

!!!!! Expectation value of x
subroutine calc_x(psi, x_values, n, dx, x_real)
  ! Dummy variables
  integer :: n
  double complex, dimension(n) :: psi
  double precision, dimension(n) :: x_values
  double precision :: dx, x_real

  !f2py intent(hide) :: n
  !f2py intent(in) :: psi, x_values, dx
  !f2py intent(out) :: x_real

  ! Local varaibles
  double complex :: x
  integer :: i

  x = 0.0d0
  do i = 1, n
    x = x + conjg(psi(i)) * x_values(i) * psi(i)
  end do
  x = x * dx
  x_real = real(x)

end subroutine calc_x
