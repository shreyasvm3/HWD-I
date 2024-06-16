! ############################################## !
! ### wigner transform of  nuclear               ### !
! ### delta(R-Rf)                                ### !
! ############################################## !
subroutine operBW(Rf,qt,dR,wf)
use parameters, only : pi
implicit none

real*8, intent(in)      :: Rf, qt, dR
complex*16, intent(out) :: wf

wf = dexp(-(qt-Rf)**2/2.d0/dR**2)/dsqrt(2.d0*pi*dR**2)
if (dble(wf).lt.1.d-10) then
   wf = 0.d0
end if

end subroutine
!######################################################!
!######################################################!
!               Symplectic Integrator                  !       
!######################################################!
subroutine integrator(q,p,dt,s,M11,M12,M21,M22)
use parameters, only : Mass, Ndof
implicit none

real*8, intent(in)    :: dt
real*8, intent(inout) :: p(Ndof), q(Ndof), s
real*8, intent(inout) :: M11(Ndof,Ndof), M12(Ndof,Ndof), M21(Ndof,Ndof), M22(Ndof,Ndof)

integer               :: j, i, k

real*8                :: ke, v, const
real*8                :: temp(Ndof), dv(Ndof), d2v(Ndof,Ndof), a(4), b(4)

const = 1.d0/dsqrt(3.d0)
a(1)  = 0.5d0*(1.d0-const)*dt
a(2)  =  const*dt
a(3)  = -const*dt
a(4)  = 0.5d0*(1.d0+const)*dt
b(1)  = 0.d0
b(2)  = 0.5d0*(0.5d0+const)*dt
b(3)  = 0.5d0*dt
b(4)  = 0.5d0*(0.5d0-const)*dt

do j = 1, 4

   if (j.gt.1) then
      call vdv_Morse(q,v,dv,d2v)
      p   = p - b(j)*dv
      M22 = M22 - b(j)*matmul(d2v,M12)
      M21 = M21 - b(j)*matmul(d2v,M11)
      s   = s - b(j)*v
   end if

   ke = 0.d0
   do i = 1, Ndof
      temp(i) = p(i)/Mass(i,i)
      ke = ke + 0.5d0*p(i)**2/Mass(i,i)
   enddo

   q   = q + a(j)*temp
   M12 = M12 + a(j)*divm(M22)
   M11 = M11 + a(j)*divm(M21)

   s = s + a(j)*ke

enddo

contains

function divm(M)
use parameters, only : Ndof, Mass
implicit none

integer :: i
real*8  :: divm(Ndof,Ndof), M(Ndof,Ndof)

do i = 1, Ndof
   divm(i,:) = M(i,:)/Mass(i,i)
enddo

end function divm

end subroutine integrator
!######################################################!

!######################################################!
!             Compute MPI iteration range              !       
!######################################################!
subroutine para_range(n1,n2,nprocs,irank,ista,iend)
implicit none
integer, intent(in)  :: n1, n2, nprocs, irank
integer, intent(out) :: ista, iend
integer              :: iwork1, iwork2

iwork1 = (n2 - n1 + 1)/nprocs
iwork2 = mod(n2 - n1 + 1, nprocs)
ista   = irank*iwork1 + n1 + min(irank, iwork2)
iend   = ista + iwork1 -1

if (iwork2 .gt. irank) then
        iend = iend + 1
endif

return

end subroutine para_range
!######################################################!
