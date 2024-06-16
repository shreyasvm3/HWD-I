module parameters
implicit none

real*8, parameter     :: pi = dacos(-1.d0), hbar = 1.d0, fstoau=4.134d1, Angtoau=1.8897d0 
real*8, parameter     :: au2cm_1 = 219474.63d0, kb = 8.6173324d-5/27.2114d0 !<--- au/K
complex*16, parameter :: Iu = (0.d0,1.d0) 

integer               :: NumberMCsteps, Ntime, Ndof, NR
real*8                :: TimeStep, EnergyTolerance, Beta

real*8                :: pIn, qIn, Rmin, Rmax
real*8                :: mdw, De, Re, alpha
real*8                :: friction, wcutoff, Tk, kbT
real*8, allocatable   :: Mass(:,:), wbath(:), coupling(:), Width(:,:)

real*8, allocatable   :: Mpp(:,:), Mqq(:,:), Mpq(:,:), Mqp(:,:)

end module

subroutine input
use parameters
implicit none

integer         :: i, j, k
character*75    :: infostr
real*8          :: u

open(555,file='input_mD',status='old')

read(555,'(a75)') infostr
read(555,*) mdw, De, Re, Alpha 

read(555,'(a75)') infostr
read(555,*) Ndof, Tk, wcutoff, friction

read(555,'(a75)') infostr
read(555,*) qIn, pIn

Ndof = Ndof + 1

allocate(Mass(Ndof,Ndof),wbath(Ndof),coupling(2:Ndof))
allocate(Width(Ndof,Ndof))!,Tuning(Ndof,Ndof),InverseWidth(Ndof,Ndof),InverseTuning(Ndof,Ndof))
allocate(Mpp(Ndof,Ndof),Mqq(Ndof,Ndof),Mpq(Ndof,Ndof),Mqp(Ndof,Ndof))

Mass    = 0.d0
wbath   = 0.d0
coupling= 0.d0
Width   = 0.d0
Mpp     = 0.d0
Mqq     = 0.d0
Mpq     = 0.d0
Mqp     = 0.d0

read(555,'(a75)') infostr
read(555,*) TimeStep, Ntime, EnergyTolerance

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

do i = 2, Ndof
   Mass(i,i)   =  1.d0
   wbath(i)    = -wcutoff*dlog((dfloat(i)-1.5d0)/dble(Ndof-1))
   coupling(i)  = wbath(i)*dsqrt(2.d0*friction*Mass(i,i)*wcutoff/(pi*dble(Ndof-1)))
enddo

do i = 1, Ndof
   Width(i,i) = Mass(i,i)*wbath(i)
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
do i = 2, Ndof
  write(111,*) wbath(i), Coupling(i)
enddo 
write(111,*)

end subroutine input


