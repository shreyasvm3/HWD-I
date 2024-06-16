!########################################
!########################################

! This subroutine computes V, grad(V),
! and hess(V) for an anharmonic
! oscillator coupled to a harmonic
! oscillator.

!########################################
!########################################
!Last Updated 10/11/2018 M.Church msc336@cornell.edu

!##################
! Anharmonic
!##################
subroutine vdv_AnHO(x,v,dv,d2v)
use parameters, only : Mass, Ndof
implicit none

real*8, intent(in) ::  x(Ndof)
integer :: i, j
real*8 :: ky, wbath, coupling
real*8, intent(out) :: v, dv(Ndof), d2v(Ndof,Ndof)

v   = 0.d0
dv  = 0.d0
d2v = 0.d0

coupling = input-couplingd0
wbath    = 1.d0/3.d0 

ky = Mass(2,2)*wbath**2

v = x(1)**2-0.1d0*x(1)**3+0.1d0*x(1)**4+0.5d0*ky*x(2)**2+coupling*x(1)*x(2)

dv(1) = 2.d0*x(1)-0.3d0*x(1)**2+0.4d0*x(1)**3+coupling*x(2)
dv(2) = ky*x(2) + coupling*x(1)

d2v(1,1) = 2.d0-0.6d0*x(1)+1.2d0*x(1)**2
d2v(1,2) = coupling
d2v(2,1) = coupling
d2v(2,2) = ky

end subroutine vdv_AnHO
!########################
