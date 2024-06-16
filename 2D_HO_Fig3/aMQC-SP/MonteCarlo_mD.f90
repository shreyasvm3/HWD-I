!#########################################################
!#########################################################

! This file contains subroutines that generate arrays of 
! initial phase space points for different SC-IVR 
! implementations using a Gaussian RNG.

! In each case the initial state for each dof is a
! coherent state and operator A is taken to be

!       A = |pIn,qIn><pIn,qIn|

! File contents:

! Subroutine                           Label               Line    
! ----------                           ---                 ----

! Husimi-IVR                           sample_HIVR         42

! LSC-IVR                              sample_LSC          72    

! DHK-IVR                              sample_DHK          102

! MQC-IVR
!  -Forward-Backward (B=B(q))          sample_FBMQC_Bq     134
!  -Forward-Backward (B=B(p))          sample_FBMQC_Bp     166
!  -Forward-Backward (General B)       sample_FBMQC        198
!  -Forward-Forward                    sample_DFMQC        231
!
! aMQC-IVR                             sample_aMQC         272 
!                               
! Initialize Random Seed               init_random_seed    310
!
! Gaussian Random Number Generator     guass               331

!#############################################################################
!Last Updated 07/07/2023 S.Nandi and S.Malpathak

!##################################
!          Husimi-IVR
!##################################
subroutine sample_HIVR(PScoord)
use parameters, only: qIn, pIn, Width0, NumberMCsteps, Ndof
implicit none

!Sampling distribution:         w = |<p0,q0|pIn,qIn>|^2

real*8, intent(inout)     :: PScoord(2,Ndof)   ! (q0,p0)
real*8                    :: gauss                           ! GRN gen.
integer                   :: i, j, k

  do j = 1, Ndof                                             ! Loop over degrees of freedom
    PScoord(1,j) = gauss(1.d0/dsqrt(Width0(j,j))) + qIn(j)
    PScoord(2,j) = gauss(dsqrt(Width0(j,j))) + pIn(j)
  enddo

end subroutine
!###################################

!###################################
!            LSC-IVR 
!###################################
subroutine sample_LSC(PScoord)
use parameters, only: qIn, pIn, Width0, NumberMCsteps, Ndof
implicit none

!Sampling distribution: Wigner transform of A

real*8, intent(inout)     :: PScoord(2,Ndof)    ! (q0,p0)
real*8                    :: gauss                            ! GRN gen.
integer                   :: i, j, k

  do j = 1, Ndof                                               ! Loop over degrees of freedom
    PScoord(1,j) = gauss(1.d0/dsqrt(Width0(j,j))) + qIn(j)
    PScoord(2,j) = gauss(dsqrt(Width0(j,j))) + pIn(j)
  enddo

end subroutine
!###################################

!###################################
!           DHK-IVR
!###################################
subroutine sample_DHK(PScoord)
use parameters, only: qIn, pIn, Width0, NumberMCsteps, Ndof
implicit none

!Sampling distribution:    w = |<p0,q0|pIn,qIn><pIn,qIn|p0',q0'>|

real*8, intent(inout) :: PScoord(4,Ndof)       ! (q0,p0,q0',p0')
real*8                :: gauss                               ! GRN gen.
integer               :: i, j, k

  do j = 1, Ndof                                             ! Loop over degrees of freedom
    PScoord(1,j) = gauss(dsqrt(2.d0/Width0(j,j))) + qIn(j) !q0
    PScoord(2,j) = gauss(dsqrt(2.d0*Width0(j,j))) + pIn(j) !p0
    PScoord(3,j) = gauss(dsqrt(2.d0/Width0(j,j))) + qIn(j) !q0'
    PScoord(4,j) = gauss(dsqrt(2.d0*Width0(j,j))) + pIn(j) !p0'
  enddo

end subroutine
!###################################

!###################################
!      FB-MQC-IVR: B = B(q)
!###################################
subroutine sample_FBMQC_Bq(PScoord)
use parameters, only: qIn, pIn, Width0, TuningP, NumberMCsteps, Ndof
implicit none

!       Sampling distribution:
!       w = |<p0,q0|pIn,qIn(j)>| exp(-cp Dp**2 /2)

real*8, intent(inout) :: PScoord(3,Ndof)       ! (q0,p0,Dp)
integer               :: i, j, k
real*8                :: gauss                               ! GRN gen.

  do j = 1, Ndof                                             ! Loop over degrees of freedom
    PScoord(1,j) = gauss(dsqrt(2.d0/Width0(j,j))) + qIn(j) ! q0
    PScoord(2,j) = gauss(dsqrt(2.d0*Width0(j,j))) + pIn(j) ! p0
    PScoord(3,j) = gauss(1.d0/dsqrt(TuningP(j,j)))         ! Dp = pt'-pt
  enddo

end subroutine
!###################################

!###################################
!       FB-MQC-IVR: B = B(p)
!###################################
subroutine sample_FBMQC_Bp(PScoord)
use parameters, only: qIn, pIn, Width0, TuningQ, NumberMCsteps, Ndof
implicit none

!       Sampling distribution:
!       w = |<p0,q0|pIn,qIn(j)>| exp(-cq Dq**2 /2)

real*8, intent(inout) :: PScoord(3,Ndof)       ! (q0,p0,Dq)
integer               :: i, j, k
real*8                :: gauss                               ! GRN gen.

  do j = 1, Ndof                                             ! Loop over degrees of freedom
    PScoord(1,j) = gauss(1.d0/dsqrt(Width0(j,j))) + qIn(j) ! q0
    PScoord(2,j) = gauss(dsqrt(Width0(j,j))) + pIn(j)      ! p0
    PScoord(3,j) = gauss(1.d0/dsqrt(TuningQ(j,j)))         ! Dq = qt'-qt
  enddo

end subroutine
!###################################

!###################################
!      FB-MQC-IVR: General B
!###################################
subroutine sample_FBMQC(PScoord)
use parameters, only: qIn, pIn, Width0, TuningQ, TuningP, NumberMCsteps, Ndof
implicit none

!       Sampling distribution:
!       w = |<p0,q0|pIn,qIn>| exp(-cq Dq**2 /2) exp(-cp Dp**2 /2)

real*8, intent(inout) :: PScoord(4,Ndof)
integer               :: i, j, k
real*8                :: gauss                                ! GRN gen.

  do j = 1, Ndof                                              ! Loop over degrees of freedom
    PScoord(1,j) = gauss(dsqrt(2.d0/Width0(j,j))) + qIn(j)  !q0
    PScoord(2,j) = gauss(dsqrt(2.d0*Width0(j,j))) + pIn(j)  ! p0
    PScoord(3,j) = gauss(1.d0/dsqrt(TuningQ(j,j)))          ! Dq = qt'-qt
    PScoord(4,j) = gauss(1.d0/dsqrt(TuningP(j,j)))          ! Dp = pt'-pt
  enddo

end subroutine
!###################################

!###################################
!           DF-MQC-IVR
!###################################
subroutine sample_DFMQC(PScoord)
use parameters, only: qIn, pIn, Width0, TuningQ, TuningP, NumberMCsteps, Ndof
implicit none

!       Sampling distribution:
!       w = |<pb,qb|pIn,qIn(j)>| exp(-cq Dq0**2 /2) exp(-cp Dp0**2 /2)
!       qb  = 0.5 ( q0' + q0 )      pb  = 0.5 ( p0' + p0 )
!       Dq0 = q0' - q0              Dp0 = p0' - p0

real*8, intent(inout) :: PScoord(4,Ndof)
integer               :: i, j, k
real*8                :: sample(4)
real*8                :: gauss                            ! GRN gen.

  do j = 1, Ndof                                          ! Loop over degrees of freedom
    sample(1) = gauss(dsqrt(2.d0/Width0(j,j))) + qIn(j)   !qav
    sample(2) = gauss(dsqrt(2.d0*Width0(j,j))) + pIn(j)   !pav
    sample(3) = gauss(1.d0/dsqrt(TuningQ(j,j)))           !Dq0
    sample(4) = gauss(1.d0/dsqrt(TuningP(j,j)))           !Dp0

    PScoord(1,j) = sample(1)-0.5d0*sample(3)            !q0
    PScoord(2,j) = sample(2)-0.5d0*sample(4)            !p0
    PScoord(3,j) = sample(1)+0.5d0*sample(3)            !q0'
    PScoord(4,j) = sample(2)+0.5d0*sample(4)            !p0'
  enddo

end subroutine
!###################################

!###################################
!           DF-MQC-IVR
!###################################
subroutine sample_aMQC(PScoord)
use parameters, only: qIn, pIn, Width0, NumberMCsteps, Ndof
implicit none

!       Sampling distribution:
!       w(anharmonic) = |<p0,q0|pIn,qIn><pIn,qIn|p0',q0'>|
!       w(harmonic)   = |<p0,q0|pIn,qIn>|^2

real*8, intent(inout) :: PScoord(4,Ndof)
integer               :: i, j, k
real*8                :: sample(4)
real*8                :: gauss                             ! GRN gen.

  PScoord(1,1) = gauss(dsqrt(2.d0/Width0(1,1))) + qIn(1)
  PScoord(2,1) = gauss(dsqrt(2.d0*Width0(1,1))) + pIn(1)
  PScoord(3,1) = gauss(dsqrt(2.d0/Width0(1,1))) + qIn(1)
  PScoord(4,1) = gauss(dsqrt(2.d0*Width0(1,1))) + pIn(1)

  PScoord(1,2) = gauss(dsqrt(1.d0/Width0(2,2))) + qIn(2)
  PScoord(2,2) = gauss(dsqrt(1.d0*Width0(2,2))) + pIn(2)
  PScoord(3,2) = PScoord(1,2)
  PScoord(4,2) = PScoord(2,2)

end subroutine
!###################################

!###################################
!       Initialize Random Seed
!###################################
subroutine init_random_seed(myrank)

integer              :: i, n, clock
integer, allocatable :: seed(:)

call random_seed(size = n)
allocate(seed(n))

call system_clock(count=clock)
seed = clock + 37 * (/ (i - 1, i = 1, n) /) + myrank*clock

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
