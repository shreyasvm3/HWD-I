!###########################################
!###########################################

! This file contains subroutines that 
! propagate trajectories and the monodromy
! matrix 
!###########################################
!###########################################  

!######################################################!
!       Propagate trajector for HWD                    !
!######################################################!
subroutine Propagate_HT(z,qp,Sfwd,Mb,dM,flagE,flagS)
use parameters, only: NSC, NCl, TimeStep, Ntime, EnergyTolerance, Ndof, SympTolerance, MInv, Jmat
implicit none

real*8, intent(in)     :: z(4,Ndof)
real*8, intent(out)    :: qp(4,Ndof,0:Ntime)
real*8, intent(out)    :: Sfwd(0:Ntime)
real*8, intent(out)    :: Mb(2*NSC,2*NSC,0:Ntime), dM(2*NSC,2*NSC,0:Ntime)
integer, intent(inout) :: flagE, flagS
integer                :: i, j, k, icount
real*8                 :: pb(Ndof), qb(Ndof), dq(NSC), dp(NSC) 
real*8                 :: Mono(4*NSC,4*NSC), Maxdiff, temp(4*NSC,4*NSC)
real*8                 :: Vp, dVp(Ndof), d2Vp(Ndof,Ndof), dMt(2*NSC,2*NSC)
real*8                 :: Vm, dVm(Ndof), d2Vm(Ndof,Ndof), Mbt(2*NSC,2*NSC)
real*8                 :: InitialEnergy, Energy, Sf, JJ(4*NSC,4*NSC) 

Sf = 0.d0

!Assign initial phase space coordinates.
qb              = z(1,:)
dq              = z(3,1:NSC)
pb              = z(2,:)
dp              = z(4,1:NSC)

Mbt   = 0.d0
dMt   = 0.d0

Mono  = 0.d0
JJ    = 0.d0
temp  = 0.d0

!Mb(0) = identity, dM(0) = 0.
do i = 1, 2*NSC
        Mbt(i,i)                = 1.d0
enddo

JJ(1:2*NSC,1:2*NSC)                   = 2.d0*Jmat
JJ(2*NSC+1:4*NSC,2*NSC+1:4*NSC)       = 0.5d0*Jmat

! Get V(q+) and V(q-)
call vdv_Morse(qb(1:NSC)+dq/2.d0,qb(NSC+1:Ndof),Vp,dVp,d2Vp)
call vdv_Morse(qb(1:NSC)-dq/2.d0,qb(NSC+1:Ndof),Vm,dVm,d2Vm)

!InitialEnergy =  dot_product(pb(1:NSC),matmul(MInv(1:NSC,1:NSC),dp(1:NSC))) + Vp - Vm
!Get energy convservation using H+ + H-
InitialEnergy = dot_product(pb,matmul(MInv,pb)) &
              + dot_product(dp(1:NSC),matmul(MInv(1:NSC,1:NSC),dp(1:NSC)))/4.d0 &
              + Vp + Vm

!print*, "E_i =", InitialEnergy

do  i  =  0,  Ntime     ! Loop over timesteps

        if (i.eq.0) goto 111

        ! Propagate for one timestep
        call Integrator(qb,pb,dq,dp,Sf,Mbt,dMt)

        111 continue

        ! Get V(q+) and V(q-)
        call vdv_Morse(qb(1:NSC)+dq/2.d0,qb(NSC+1:Ndof),Vp,dVp,d2Vp)
        !print*, "Vp", Vp
        !print*, "dVp", dVp
        !print*, "d2Vp", d2Vp
        call vdv_Morse(qb(1:NSC)-dq/2.d0,qb(NSC+1:Ndof),Vm,dVm,d2Vm)
        !print*, "Vm", Vm
        !print*, "dVm", dVm
        !print*, "d2Vm", d2Vm
       Energy = dot_product(pb,matmul(MInv,pb)) &
              + dot_product(dp(1:NSC),matmul(MInv(1:NSC,1:NSC),dp(1:NSC)))/4.d0 &
              + Vp + Vm
       !print*, "E_t for t=", i, ":", Energy
     
       ! Define full M matrix
        Mono(1:2*NSC,1:2*NSC)               = Mbt
        Mono(1:2*NSC,2*NSC+1:4*NSC)         = dMt
        Mono(2*NSC+1:4*NSC,1:2*NSC)         = dMt*4.d0
        Mono(2*NSC+1:4*NSC,2*NSC+1:4*NSC)   = Mbt
     
       ! symplecticity conservation check
        temp    = matmul(transpose(Mono),matmul(JJ,Mono))
        MaxDiff = maxval(dabs(temp-JJ))
        
        if (dabs(1.d0-Energy/InitialEnergy).ge.EnergyTolerance) then
           flagE = 1
           goto 112
        elseif (MaxDiff.ge.SympTolerance) then
           flagS = 1
       !    print*, "breaks at time:", i
           goto 112
        endif
      
        ! Save phase space coordinates, action and monodromy matrix
        ! at each timestep
        Sfwd(i)                 = Sf
        qp(1,:,i)               = qb
        qp(2,:,i)               = pb
        qp(3,1:NSC,i)           = dq
        qp(4,1:NSC,i)           = dp
        qp(3,NSC+1:Ndof,i)      = 0.d0  ! dq_cl = 0
        qp(4,NSC+1:Ndof,i)      = 0.d0  ! dp_cl = 0 
        Mb(:,:,i)               = Mbt
        dM(:,:,i)               = dMt

enddo

112 continue

return
end subroutine
!######################################################!
