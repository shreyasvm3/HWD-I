!Last Updated 02/08/2024 Shreyas M. sm2699@cornell.edu
!###################################
!           DHK/LSC
!###################################
subroutine sample_HWD(PScoord)
use parameters, only: qIn, pIn, Width0, Ndof, NSC, Beta, pi, wbath
implicit none

!Sampling distribution:    w = |<p0,q0|pIn,qIn><pIn,qIn|p0',q0'>| for SC dofs
!                          w = [exp(-beta*H_B)]_W for Cl dofs
!Sample in zb and dZ variables. Set dz_cl = 0.

real*8, intent(inout) :: PScoord(4,Ndof)                     ! (q0,p0,q0',p0')
real*8                :: gauss,sigmaP, sigmaQ                ! GRN gen.
integer               :: i 

!Sample Morse dof
   PScoord(1,1) = gauss(dsqrt(1.d0/Width0(1,1))) + qIn     !q0b
   PScoord(2,1) = gauss(dsqrt(1.d0*Width0(1,1))) + pIn     !p0b
   PScoord(3,1) = gauss(dsqrt(4.d0/Width0(1,1)))           !dq0
   PScoord(4,1) = gauss(dsqrt(4.d0*Width0(1,1)))           !dp0

!Sample bath dofs
do i = NSC+1, Ndof 
      ! For Cl dofs, diff. variable is zero
      ! Mean is sampled from W.T.
   sigmaP = dsqrt(Width0(i,i)/2.d0/dtanh(Beta*wbath(i)/2.d0))
   sigmaQ = dsqrt(0.5d0/width0(i,i)/dtanh(Beta*wbath(i)/2.d0))
   PScoord(1,i) = gauss(sigmaQ)                                !q0b
   PScoord(2,i) = gauss(sigmaP)                                !p0b
   PScoord(3,i) = 0.d0                                         !dq0
   PScoord(4,i) = 0.d0                                         !dp0
enddo

end subroutine
!###################################
!###################################
!       Initialize Random Seed
!###################################
subroutine init_random_seed(myrank)

integer, intent(in)  :: myrank
integer              :: i, n, clock
integer, allocatable :: seed(:)

call random_seed(size = n)
allocate(seed(n))

call system_clock(count=clock)
seed = clock + 37 * (/ (i - 1, i = 1, n) /) + clock*myrank

call random_seed(put = seed)

deallocate(seed)

end subroutine init_random_seed
!##############################################

!##############################################
!     Gaussian random number generator
!##############################################
double precision function gauss(sigma)

implicit none
real*8, intent(in)      :: sigma
real*8                  :: w, v1, v2, l
real*8                  :: s1, s2

w = 2.d0

do
        call random_number(s1)
        call random_number(s2)

        v1 = 2.d0*s1 - 1.d0
        v2 = 2.d0*s2 - 1.d0
        w = v1*v1 + v2*v2

        if (w.lt.1.d0) exit
end do

l     = v1*sqrt(-2.d0*log(w)/(w))
l     = sigma*l
gauss = l

end function gauss
!##############################################
