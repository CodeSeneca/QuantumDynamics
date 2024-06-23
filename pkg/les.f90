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
  double precision :: dt, dx, mass
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

