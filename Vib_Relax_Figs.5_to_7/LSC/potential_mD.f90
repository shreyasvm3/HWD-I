!########################################
!########################################

! This subroutine computes V, grad(V),
! and hess(V) for Morse oscilaltor linearly 
! coupled to a bath of harmonic oscillators

!########################################
!########################################

!##################
! Morse + HO bath
!##################
subroutine vdv_Morse(x,v,dv,d2v)
use parameters, only : Ndof, Mass, De, Re, alpha, coupling, wbath
implicit none

real*8, intent(in) ::  x(Ndof)
integer :: i
real*8  :: f
real*8, intent(out) :: v, dv(Ndof), d2v(Ndof,Ndof)

v   = 0.d0
dv  = 0.d0
d2v = 0.d0

f = dexp(-alpha*(x(1)-Re))

v        = De*(1.d0-f)**2               ! Morse dof 
dv(1)    = 2.d0*alpha*De*(f-f**2)            ! Derv. of Morse wrt R
d2v(1,1) = 2.d0*alpha**2*De*(2.d0*f**2-f)    ! 2nd derv. of Morse wrt R

do i = 2, Ndof
  v        = v + 0.5d0*wbath(i)**2*(x(i) + coupling(i)*(x(1)-Re)/wbath(i)**2)**2   ! Add in V_b and V_sb
  
  dv(1)    = dv(1) + coupling(i)*(x(i) + coupling(i)*(x(1)-Re)/wbath(i)**2) ! Add in bath part for dV/ds
  dv(i)    = wbath(i)**2*x(i) + coupling(i)*(x(1)-Re)          ! dV/dQ_j
   
  d2v(1,1) = d2v(1,1) + (coupling(i)/wbath(i))**2 !Add in bath part of d2V/ds2
  d2v(1,i) = coupling(i)  !d2V/dsdQ_j
  d2v(i,1) = coupling(i) 
  d2v(i,i) = wbath(i)**2  !d2V/dQ_j^2

end do
end subroutine vdv_Morse
!########################
