!#################################################
!#################################################          

! The module defines the global variables for a
! correlation function of a multidimensional system.                  

! The additional subroutine, 'input',  reads the 
! input file and allocates arrays accordingly.

!#################################################
!#################################################
module parameters
implicit none

! Define pi and the imaginary unit
real*8, parameter     :: pi=dacos(-1.d0)
complex*16, parameter :: Iu=(0.d0,1.d0)

! Degrees of freedom: total, SC and NCl
integer               :: Ndof, NSC, NCl
!1:NSC are SC dofs, NSC+1:NSC+NCl=Ndof are Cl dofs.

! Number of trajectories
integer               ::  NumberMCsteps

! Number of timesteps
integer               ::  Ntime

! Timestep, threshold of energy conservation
real*8                ::  Timestep, EnergyTolerance, SympTolerance

!Coupling b/w oscillators
real*8                :: coupling
! Masses
real*8, allocatable   :: Mass(:,:), MInv(:,:)

! Coherent state parameters
real*8, allocatable   :: Width0(:,:), InverseWidth0(:,:) 
real*8, allocatable   :: qIn(:), pIn(:)    

!J Matrix, Gamma Matrix, J.G matrix
real*8, allocatable   :: JMat(:,:), GMat(:,:), JGMat(:,:)

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
read(555,*) Ndof, NSC, NCl

allocate(Mass(Ndof,Ndof),Width0(Ndof,Ndof))
allocate(InverseWidth0(Ndof,Ndof))
allocate(qIn(Ndof),pIn(Ndof),MInv(Ndof,Ndof))
allocate(JMat(2*NSC,2*NSC),Gmat(2*NSC,2*NSC))
allocate(JGMat(2*NSC,2*NSC))

Mass            = 0.d0
MInv            = 0.d0
Width0          = 0.d0
InverseWidth0   = 0.d0
Jmat            = 0.d0
Gmat            = 0.d0
JGmat           = 0.d0

read(555,'(a75)') infostr
read(555,*) (Mass(j,j),j=1,Ndof)

read(555,'(a75)') infostr
do k = 1, Ndof
  read(555,*) qIn(k), pIn(k), Width0(k,k)
enddo

read(555,'(a75)') infostr
read(555,*) coupling

read(555,'(a75)') infostr
read(555,*) TimeStep, Ntime, EnergyTolerance, SympTolerance

read(555,'(a75)') infostr
read(555,*) NumberMCsteps

close(555)

! Inverses and other matrices
do k = 1, Ndof
  InverseWidth0(k,k)  = 1.d0/Width0(k,k)
  MInv(k,k)           = 1.d0/Mass(k,k)
end do
do k = 1, NSC
  JMat(k,k+NSC)      = 1.d0
  JMat(k+NSC,k)      =-1.d0
  GMat(k,k)          = Width0(k,k)
  Gmat(k+NSC,k+NSC)  = 1.d0/Width0(k,k)
enddo
JGMat = matmul(Jmat,Gmat)

end subroutine input
!######################################################
