!###################################
!            LSC-IVR 
!###################################
subroutine sample_LSC(PScoord)
use parameters, only: Ndof, Beta, pi, Width, qIn, pIn, wbath
implicit none

!Sampling distribution: Wigner transform of |zi><zi| x exp(-beta*H_b)

real*8, intent(out)       :: PScoord(2,Ndof)                  ! (q0,p0)
real*8                    :: gauss, sigmaP, sigmaQ            ! GRN gen.
integer                   :: i

!Sample Morse from [|zi><zi|]_W
sigmaQ = dsqrt(0.5d0/Width(1,1))
sigmaP = dsqrt(0.5d0*Width(1,1))
PScoord(1,1) =gauss(sigmaQ) + qIn  
PScoord(2,1) =gauss(sigmaP) + pIn
 
! bath dofs
do i = 2, Ndof
   sigmaP = dsqrt(Width(i,i)/2.d0/dtanh(Beta*wbath(i)/2.d0))
   sigmaQ = dsqrt(0.5d0/width(i,i)/dtanh(Beta*wbath(i)/2.d0))
   PScoord(1,i) = gauss(sigmaQ) 
   PScoord(2,i) = gauss(sigmaP) 
end do
!}

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
