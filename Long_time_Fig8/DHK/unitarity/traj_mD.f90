!###########################################
!###########################################

! This file contains subroutines that 
! propagate trajectories and the monodromy
! matrix for any 1D SC-IVR

! Trajectory flags are different for 
! Forward-Backward and Forward-Forward
! implementations so they use separate
! subroutines

!###########################################
!###########################################  
!Last Updated 09/2/2016 M.Church msc336@cornell.edu

!######################################################!
!       FB: Propagate forward trajectory               !
!######################################################!
subroutine PropagateFwd_FB(Initialq,Initialp,coord,Sf,Monodromy,flag)
use parameters, only: TimeStep, Ntime, Mass, EnergyTolerance,&
Mqq, Mqp, Mpq, Mpp, Ndof
implicit none

real*8, intent(in)     :: Initialq(Ndof), Initialp(Ndof)
real*8, intent (out)   :: coord(2,Ndof,0:Ntime), Sf(0:Ntime), Monodromy(Ndof,Ndof,4,0:Ntime)
integer, intent(inout) :: flag                          

integer :: i, j, k, icount
real*8  :: p(Ndof), q(Ndof), Mono(Ndof,Ndof), v, dv(Ndof), d2v(Ndof,Ndof)
real*8  :: InitialEnergy, Energy, s, det

s = 0.d0
q = Initialq
p = Initialp

Mpp = 0.d0
Mqq = 0.d0
Mpq = 0.d0
Mqp = 0.d0

do i = 1, Ndof
        Mpp(i,i) = 1.d0
        Mqq(i,i) = 1.d0
enddo

call vdv_AnHO(q,v,dv,d2v)

InitialEnergy = v
do j = 1, Ndof
   InitialEnergy = InitialEnergy + p(j)**2/2.d0/Mass(j,j)
enddo

do i = 0, Ntime        ! Loop over timesteps

        if (i.eq.0) goto 111

        ! Propagate for one timestep
        call Integrator(q,p,TimeStep,s,Mqq,Mqp,Mpq,Mpp)

        111 continue

        ! Check energy conservation
        call vdv_AnHO(q,v,dv,d2v)

        Energy = v
        do j = 1, Ndof
                Energy = Energy + p(j)**2/2.d0/Mass(j,j)
        enddo

        ! Save phase space coordinates, action and monodromy matrix
        ! at each timestep
        Sf(i)          = s
        coord(1,:,i)     = q
        coord(2,:,i)     = p
        Monodromy(:,:,1,i) = Mqq
        Monodromy(:,:,2,i) = Mqp
        Monodromy(:,:,3,i) = Mpq
        Monodromy(:,:,4,i) = Mpp
   
        ! Check symplecticity 
        Mono = matmul(transpose(Mqq),Mpp)-matmul(transpose(Mpq),Mqp)
        det  = Mono(1,1)*Mono(2,2)-Mono(1,2)*Mono(2,1)

        ! Flag trajectory if energy cons or sympecticity breaks
        if (dabs(1.d0-Energy/InitialEnergy).ge.EnergyTolerance.or.dabs(det-1.d0).ge.1.d-5) then
                flag = 1
                goto 112 
        endif
enddo

112 continue

return

end subroutine
!######################################################!

!######################################################!
!       FB: Propagate backward trajectory              !
!######################################################!
subroutine PropagateBck_FB(Nt,Initialq,Initialp,coord,Sb,Monodromy,flag)
use parameters, only: TimeStep, Ntime, EnergyTolerance, Mass,&
Mqq, Mqp, Mpq, Mpp, Ndof
implicit none

integer, intent(in)    :: Nt
real*8,  intent(in)    :: Initialq(Ndof), Initialp(Ndof)
real*8,  intent(out)   :: coord(2,Ndof,0:Ntime), Sb(0:Ntime), Monodromy(Ndof,Ndof,4,0:Ntime)
integer, intent(inout) :: flag

integer :: i, j, k, icount, z, zp
real*8  :: p(Ndof), q(Ndof), v, dv(Ndof), d2v(Ndof,Ndof), s
real*8  :: InitialEnergy, Energy, dt
real*8  :: Mono(Ndof,Ndof), det

dt = -TimeStep

s = 0.d0
q = initialq
p = initialp

Mqq = 0.d0
Mpq = 0.d0
Mqp = 0.d0
Mpp = 0.d0

do i = 1, Ndof
        Mpp(i,i) = 1.d0
        Mqq(i,i) = 1.d0
enddo

Sb        = 0.d0
coord     = 0.d0
Monodromy = 0.d0

call vdv_AnHO(q,v,dv,d2v)

InitialEnergy = v
do j = 1, Ndof
   InitialEnergy = InitialEnergy + p(j)**2/2.d0/Mass(j,j)
enddo

! Loop over timesteps
do i = 0, Nt            ! Propagate with a negative timestep

        if (i.eq.0) goto 111

        ! Propagate back for a single timestep
        call Integrator(q,p,dt,s,Mqq,Mqp,Mpq,Mpp)

        111 continue

        ! Check energy conservation
        call vdv_AnHO(q,v,dv,d2v)

        Energy = v
        do j = 1, Ndof
                Energy = Energy+p(j)**2/2.d0/Mass(j,j)
        enddo

        ! Save phase space coordinates, action and monodromy
        ! matrix at each timestep
        Sb(i)=s
        coord(1,:,i)=q
        coord(2,:,i)=p
        Monodromy(:,:,1,i)=Mqq
        Monodromy(:,:,2,i)=Mqp
        Monodromy(:,:,3,i)=Mpq
        Monodromy(:,:,4,i)=Mpp

        ! Check symplecticity
        Mono = matmul(transpose(Mqq),Mpp)-matmul(transpose(Mpq),Mqp)
        det  = Mono(1,1)*Mono(2,2)-Mono(1,2)*Mono(2,1)

        ! Flag trajectory if energy cons or symplecticity breaks

        if (dabs(1.d0-Energy/InitialEnergy).ge.EnergyTolerance .or. dabs(det-1.d0).ge.1.d-5) then
                flag = 1
                goto 112
        endif

enddo

112 continue

return

end subroutine
!######################################################!

!######################################################!
!       DF: Propagate trajectory forward               !
!######################################################!
subroutine PropagateFwd_DF(InitialQ,InitialP,coord,Sf,Monodromy,flag)
use parameters, only: TimeStep, Ntime, Mqq, Mqp, Mpq, Mpp, Mass,&
EnergyTolerance, Ndof
implicit none

real*8, intent(in)     :: InitialQ(Ndof), InitialP(Ndof)
real*8, intent (out)   :: coord(2,Ndof,0:Ntime), Sf(0:Ntime), Monodromy(Ndof,Ndof,4,0:Ntime)
integer, intent(inout) :: flag

integer :: i, j, k, icount
real*8  :: p(Ndof), q(Ndof), Mono(Ndof,Ndof), v, dv(Ndof), d2v(Ndof,Ndof)
real*8  :: InitialEnergy, Energy, s, det

s = 0.d0
q = InitialQ
p = InitialP

Mpp = 0.d0
Mqq = 0.d0
Mpq = 0.d0
Mqp = 0.d0

do i = 1, Ndof
        Mpp(i,i) = 1.d0
        Mqq(i,i) = 1.d0
enddo

call vdv_AnHO(q,v,dv,d2v)

InitialEnergy = v
do j = 1, Ndof
   InitialEnergy = InitialEnergy + p(j)**2/2.d0/Mass(j,j)
enddo

do  i  =  0,  Ntime     ! Loop over timesteps

        if (i.eq.0) goto 111

        ! Propagate for one timestep
        call Integrator(q,p,TimeStep,s,Mqq,Mqp,Mpq,Mpp)

  111 continue

        ! Check energy conservation
        call vdv_AnHO(q,v,dv,d2v)

        Energy = v
        do j = 1, Ndof
                Energy = Energy + p(j)**2/2.d0/Mass(j,j)
        enddo

        ! Save phase space coordinates, action and monodromy matrix
        ! at each timestep
        Sf(i)              = s
        coord(1,:,i)       = q
        coord(2,:,i)       = p
        Monodromy(:,:,1,i) = Mqq
        Monodromy(:,:,2,i) = Mqp
        Monodromy(:,:,3,i) = Mpq
        Monodromy(:,:,4,i) = Mpp

        Mono = matmul(transpose(Mqq),Mpp)-matmul(transpose(Mpq),Mqp)
        det  = Mono(1,1)*Mono(2,2)-Mono(1,2)*Mono(2,1)

        if (dabs(1.d0-Energy/InitialEnergy).ge.EnergyTolerance .or. dabs(det-1.d0).ge.1.d-5) then
                flag = 1
                goto 112
        endif

enddo

112 continue

return
end subroutine
!######################################################!
