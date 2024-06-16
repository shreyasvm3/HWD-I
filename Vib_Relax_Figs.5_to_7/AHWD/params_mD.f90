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
real*8, parameter     :: pi = dacos(-1.d0), hbar = 1.d0, fstoau=4.134d1, Angtoau=1.8897d0 
real*8, parameter     :: au2cm_1 = 219474.63d0, kb = 8.6173324d-5/27.2114d0 !<--- au/K
complex*16, parameter :: Iu = (0.d0,1.d0) 

! Degrees of freedom: total, SC and NCl. NR: # of grid points for density
integer               :: Ndof, NSC, NCl, NR
!1:NSC are SC dofs, NSC+1:NSC+NCl=Ndof are Cl dofs.

! Number of trajectories
integer               ::  NumberMCsteps

! Number of timesteps
integer               ::  Ntime

! Timestep, threshold of energy conservation
real*8                ::  Timestep, EnergyTolerance, SympTolerance

! Masses
real*8, allocatable   :: Mass(:,:), MInv(:,:)

! Coherent state parameters
real*8, allocatable   :: Width0(:,:), InverseWidth0(:,:) 
real*8                :: qIn, pIn    

!J Matrix, Gamma Matrix, J.G matrix
real*8, allocatable   :: JMat(:,:), GMat(:,:), JGMat(:,:)

!grid paramteres
real*8                :: Rmin, Rmax

!Morse dof parameters
real*8                :: mdw, De, Re, alpha

!Bath parameters
real*8                :: friction, wcutoff, Tk, kbT, Beta
real*8, allocatable   :: wbath(:), coupling(:)


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
allocate(MInv(Ndof,Ndof))
allocate(JMat(2*NSC,2*NSC),Gmat(2*NSC,2*NSC))
allocate(JGMat(2*NSC,2*NSC))
allocate(wbath(Ndof),coupling(NSc+1:Ndof))

Mass            = 0.d0
MInv            = 0.d0
Width0          = 0.d0
InverseWidth0   = 0.d0
Jmat            = 0.d0
Gmat            = 0.d0
JGmat           = 0.d0
wbath           = 0.d0
coupling        = 0.d0

read(555,'(a75)') infostr
read(555,*) mdw, De, Re, Alpha 

read(555,'(a75)') infostr
read(555,*) Tk, wcutoff, friction

read(555,'(a75)') infostr
read(555,*) qIn, pIn

read(555,'(a75)') infostr
read(555,*) TimeStep, Ntime, EnergyTolerance, SympTolerance

read(555,'(a75)') infostr
read(555,*) Rmin, Rmax, NR

Rmin = Rmin*Angtoau
Rmax = Rmax*Angtoau

read(555,'(a75)') infostr
read(555,*) NumberMCsteps

close(555)

kbT       = kb*Tk
Beta      = 1.d0/kbT
De        = De/au2cm_1
Re        = Re*Angtoau
alpha     = alpha/Angtoau
wcutoff   = wcutoff/au2cm_1
Timestep  = Timestep*fstoau
qIn       = qIn*Angtoau

Mass(1,1) = mdw !Mass of Morse dof
wbath(1)  = alpha*dsqrt(2.d0*De/mdw) !Harmonic freq. of morse dof at bottom
friction  = friction*mdw*wbath(1)

do i = NSC+1, Ndof
   Mass(i,i)   =  1.d0
   wbath(i)    = -wcutoff*dlog((dfloat(i)-1.5d0)/dble(Ndof-1))
   coupling(i)  = wbath(i)*dsqrt(2.d0*friction*Mass(i,i)*wcutoff/(pi*dble(Ndof-1)))
enddo

do i = 1, Ndof
   Width0(i,i) = Mass(i,i)*wbath(i)
end do

open(111,file="Frequency.out",status='unknown')

write(111,*)
write(111,*) "Discretized bath frequencies: au, cm-1"
do i = 1, Ndof
  write(111,*) wbath(i), wbath(i)*au2cm_1
enddo
write(111,*)

write(111,*)
write(111,*) 'Coupling'
write(111,*)
do i = NSc+1, Ndof
  write(111,*) wbath(i), Coupling(i)
enddo 
write(111,*)


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
