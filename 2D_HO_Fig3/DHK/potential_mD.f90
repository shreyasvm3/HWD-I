!########################################
!########################################

! This subroutine computes V, grad(V),
! and hess(V) for an anharmonic
! oscillator coupled to a harmonic
! oscillator.

!########################################
!########################################

!##################
! Anharmonic
!##################
subroutine vdv_AnHO(x,v,dv,d2v)
use parameters, only : Mass, Ndof, coupling
implicit none

real*8, intent(in) ::  x(Ndof)
integer :: i, j
real*8 :: ky, wbath
real*8, intent(out) :: v, dv(Ndof), d2v(Ndof,Ndof)

v   = 0.d0
dv  = 0.d0
d2v = 0.d0

wbath    = 0.333d0 

ky = Mass(2,2)*wbath**2

v = x(1)**2 + 0.5d0*ky*x(2)**2 + coupling*x(1)*x(2)

dv(1) = 2.d0*x(1) + coupling*x(2)
dv(2) = ky*x(2) + coupling*x(1)

d2v(1,1) = 2.d0
d2v(1,2) = coupling
d2v(2,1) = coupling
d2v(2,2) = ky

end subroutine vdv_AnHO
!########################
