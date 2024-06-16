!Last Updated 05/03/2023 Shreyas M. sm2699@cornell.edu
!###################################
!           DHK/LSC
!###################################
subroutine sample_HWD(PScoord)
use parameters, only: qIn, pIn, Width0, Ndof, NSC
implicit none

!Sampling distribution:    w = |<p0,q0|pIn,qIn><pIn,qIn|p0',q0'>| for SC dofs
!                          w = [|pIn,qIn><pIn,qIn|_W for Cl dofs
!Sample in zb and dZ variables. Set dz_cl = 0.

real*8, intent(inout) :: PScoord(4,Ndof)                     ! (q0,p0,q0',p0')
real*8                :: gauss                               ! GRN gen.
integer               :: j 

do j = 1, NSC                                          ! Loop over degrees of freedom
   PScoord(1,j) = gauss(dsqrt(1.d0/Width0(j,j))) + qIn(j)  !q0b
   PScoord(2,j) = gauss(dsqrt(1.d0*Width0(j,j))) + pIn(j)  !p0b
   PScoord(3,j) = gauss(dsqrt(4.d0/Width0(j,j)))           !dq0
   PScoord(4,j) = gauss(dsqrt(4.d0*Width0(j,j)))           !dp0
enddo
do j = NSC+1, Ndof 
      ! For Cl dofs, diff. variable is zero
      ! Mean is sampled from W.T.
   PScoord(1,j) = gauss(1.d0/dsqrt(2.d0*Width0(j,j))) + qIn(j) !q0b 
   PScoord(2,j) = gauss(dsqrt(Width0(j,j)/2.d0))      + pIn(j) !p0b
   PScoord(3,j) = 0.d0                                         !dq0
   PScoord(4,j) = 0.d0                                         !dp0
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
