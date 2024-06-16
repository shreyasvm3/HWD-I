!#################################################
!#################################################          

! The module defines the global variables for a
! correlation function of a multidimensional system.                  

! The additional subroutine, 'input',  reads the 
! input file and allocates arrays accordingly.

!#################################################
!#################################################
!Last Updated 10/13/2018 M.Church msc336@cornell.edu
module parameters
implicit none

! Define pi and the imaginary unit
real*8, parameter     :: pi=dacos(-1.d0)
complex*16, parameter :: Iu=(0.d0,1.d0)

! Degrees of freedom
integer               :: Ndof

! Number of trajectories
integer               ::  NumberMCsteps

! Number of timesteps
integer               ::  Ntime

! Timestep, threshold of energy conservation
real*8                ::  Timestep, EnergyTolerance, coupling

! Masses
real*8, allocatable   :: Mass(:,:)

! Monodromy matrix elements
real*8, allocatable   :: Mqq(:,:), Mqp(:,:), Mpq(:,:), Mpp(:,:)

! Coherent state parameters
real*8, allocatable   :: Width0(:,:), InverseWidth0(:,:) 
real*8, allocatable   :: WidthT(:,:), InverseWidthT(:,:)  
real*8, allocatable   :: qIn(:), pIn(:)    

! Tuning matrices
real*8, allocatable   :: TuningQ(:,:), TuningP(:,:)
real*8, allocatable   :: InverseTuningQ(:,:), InverseTuningP(:,:)

end module
!######################################################

!######################################################
!               Input Subroutine
!######################################################
subroutine input
use parameters
implicit none

integer       :: i, j, k, dummy
character*75  :: infostr

open(555,file='input_mD',status='old')

read(555,'(a75)') infostr
read(555,*) Ndof

allocate(Mass(Ndof,Ndof),Width0(Ndof,Ndof),WidthT(Ndof,Ndof))
allocate(InverseWidth0(Ndof,Ndof),InverseWidthT(Ndof,Ndof))
allocate(qIn(Ndof),pIn(Ndof),TuningP(Ndof,Ndof),TuningQ(Ndof,Ndof))
allocate(InverseTuningP(Ndof,Ndof),InverseTuningQ(Ndof,Ndof))
allocate(Mqq(Ndof,Ndof),Mqp(Ndof,Ndof),Mpq(Ndof,Ndof),Mpp(Ndof,Ndof))

Mass            = 0.d0
Width0          = 0.d0
WidthT          = 0.d0
InverseWidth0   = 0.d0
InverseWidthT   = 0.d0
TuningQ         = 0.d0
TuningP         = 0.d0
InverseTuningQ  = 0.d0
InverseTuningP  = 0.d0

do i = 1, 4
  read(555,'(a75)') infostr
  read(555,*)       dummy
enddo

read(555,'(a75)') infostr
read(555,*) (Mass(j,j),j=1,Ndof)

read(555,'(a75)') infostr
read(555,*) coupling

read(555,'(a75)') infostr
do k = 1, Ndof
  read(555,*) qIn(k), pIn(k), Width0(k,k), WidthT(k,k), TuningQ(k,k), TuningP(k,k)
enddo

read(555,'(a75)') infostr
read(555,*) TimeStep, Ntime, EnergyTolerance

read(555,'(a75)') infostr
read(555,*) NumberMCsteps

close(555)

! Inverses
do k = 1, Ndof
  InverseWidth0(k,k)  = 1.d0/Width0(k,k)
  InverseWidthT(k,k)  = 1.d0/WidthT(k,k)
  InverseTuningQ(k,k) = 1.d0/TuningQ(k,k)
  InverseTuningP(k,k) = 1.d0/TuningP(k,k)
enddo

end subroutine input
!######################################################
