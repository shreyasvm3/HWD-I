!######################################################!
!            Matrix element of B = delt(x-x_0)         !
!           x of all dofs is evaluated                 !
!######################################################!
subroutine DF_B(z,grid,Bq)
use parameters, only : Iu, Width0, NR, pi
implicit none

real*8, intent(in)      :: z(4), grid(NR)
complex*16, intent(out) :: Bq(NR)
integer                 :: i
real*8                  :: qb, pb, dq, dp !Only morse dof
 
qb = z(1)
pb = z(2)
dq = z(3) 
dp = z(4) 
   
!<zt'|delta(x-x_0)|zt> for SC dofs, 1_W for Cl dofs
!Already includes -pb_t.dq_t from S^tilda.

do i = 1, NR
   Bq(i) = dsqrt(Width0(1,1)/pi)*dexp(-Width0(1,1)*(grid(i)-qb)**2 - 2.5d-1*Width0(1,1)*dq**2) &
         * cdexp(Iu*(dp*(grid(i)-qb)-dq*pb))
end do 

end subroutine
!#######################################################!
!#######################################################!
!       Prefactor for HWD with SP                       !    
!#######################################################!
subroutine HWDpref(Mb,dM,pref)
use parameters, only : Iu, JMat, GMat, JGMat, NSC
implicit none

real*8, intent(in)      :: Mb(2*NSC,2*NSC), dM(2*NSC,2*NSC)
complex*16, intent(out) :: pref

integer                 :: i, info, sgn
integer, allocatable    :: ipiv(:)

complex*16              :: pref1(2*NSC,2*NSC)

! HWD pref written using Mb and dM

pref1 = 0.5d0*(matmul(GMat,Mb) + 2.d0*Iu*matmul(Gmat,matmul(dM,JGMat)) &
      + 2.d0*Iu*matmul(Jmat,dM) - matmul(JMat,matmul(Mb,JGMat)))

!Calculate determinant using LU decomposition
allocate(ipiv(2*NSC))

call zgetrf(2*NSC,2*NSC,pref1,2*NSC,ipiv,info)

sgn   = 1
pref  = 1.d0

do i = 1, 2*NSC
        pref = pref*pref1(i,i)
        if (ipiv(i).ne.i) sgn = -sgn
enddo

pref = pref*sgn

end subroutine
!######################################################!

!######################################################!
!               Symplectic Integrator                  !       
!######################################################!
subroutine integrator(qb,pb,dq,dp,Sf,Mbt,dMt)
use parameters, only : TimeStep, Ndof, NSC, MInv
implicit none

real*8, intent(inout) :: pb(Ndof), qb(Ndof), dq(NSC), dp(NSC), Sf
real*8, intent(inout) :: Mbt(2*NSC,2*NSC), dMt(2*NSC,2*NSC)
integer               :: j
real*8                :: V, const, dt
real*8                :: a(4), b(4)
real*8                :: Vp, dVp(Ndof), d2Vp(Ndof,Ndof)
real*8                :: Vm, dVm(Ndof), d2Vm(Ndof,Ndof)
real*8,dimension(Ndof,Ndof) :: d2Vs, d2Vd
real*8,dimension(NSC,NSC)   :: Mb11, Mb12, Mb21, Mb22
real*8,dimension(NSC,NSC)   :: dM11, dM12, dM21, dM22

dt    = TimeStep
const = 1.d0/dsqrt(3.d0)
a(1)  = 0.5d0*(1.d0-const)*dt
a(2)  =  const*dt
a(3)  = -const*dt
a(4)  = 0.5d0*(1.d0+const)*dt
b(1)  = 0.d0
b(2)  = 0.5d0*(0.5d0+const)*dt
b(3)  = 0.5d0*dt
b(4)  = 0.5d0*(0.5d0-const)*dt

d2Vs  = 0.d0
d2Vd  = 0.d0

Mb11  = Mbt(1:NSC,1:NSC)
Mb12  = Mbt(1:NSC,NSC+1:2*NSC)
Mb21  = Mbt(NSC+1:2*NSC,1:NSC)
Mb22  = Mbt(NSC+1:2*NSC,NSC+1:2*NSC)

dM11  = dMt(1:NSC,1:NSC)
dM12  = dMt(1:NSC,NSC+1:2*NSC)
dM21  = dMt(NSC+1:2*NSC,1:NSC)
dM22  = dMt(NSC+1:2*NSC,NSC+1:2*NSC)

! pb,qp,dq,dp, Sf and all 4 M matrix elements are propagated under the full extended
! Hamiltonian

do j = 1, 4

   if (j.gt.1) then

      ! Calculate V+, V- and Vtilda
      call vdv_Morse(qb(1:NSC)+dq/2.d0,qb(NSC+1:Ndof),Vp,dVp,d2Vp)
      call vdv_Morse(qb(1:NSC)-dq/2.d0,qb(NSC+1:Ndof),Vm,dVm,d2Vm)
      d2Vs = d2Vp + d2Vm
      d2Vd = d2Vp - d2Vm
      !Propagate pb and dp
      pb   = pb  - b(j)*(dVp+dVm)/2.d0
      dp   = dp  - b(j)*(dVp(1:NSC)-dVm(1:NSC))
      ! Propagate action under H^tilda
      Sf   = Sf - b(j)*(Vp-Vm)
      !Propagate M_pp and M_pq for Mb and dM
      Mb22 = Mb22 - b(j)*(matmul(d2Vs(1:NSC,1:NSC)/2.d0,Mb12) + matmul(d2Vd(1:NSC,1:NSC),dM12))
      Mb21 = Mb21 - b(j)*(matmul(d2Vs(1:NSC,1:NSC)/2.d0,Mb11) + matmul(d2Vd(1:NSC,1:NSC),dM11))
      dM22 = dM22 - b(j)*(matmul(d2Vs(1:NSC,1:NSC)/2.d0,dM12) + matmul(d2Vd(1:NSC,1:NSC)/4.d0,Mb12))
      dM21 = dM21 - b(j)*(matmul(d2Vs(1:NSC,1:NSC)/2.d0,dM11) + matmul(d2Vd(1:NSC,1:NSC)/4.d0,Mb11))

   end if

   !Propagate qb and dq
   qb  = qb + a(j)*matmul(MInv,pb)
   dq  = dq + a(j)*matmul(MInv(1:NSc,1:NSc),dp)
   !Propagate St
   Sf  = Sf + a(j)*dot_product(pb(1:NSc),matmul(MInv(1:Nsc,1:NSc),dp(1:NSc)))
   !Propagate M11 nnd M12 for dM and Mb
   Mb12 = Mb12 + a(j)*matmul(MInv(1:NSC,1:NSC),Mb22)
   Mb11 = Mb11 + a(j)*matmul(MInv(1:NSC,1:NSC),Mb21)
   dM12 = dM12 + a(j)*matmul(MInv(1:NSC,1:NSC),dM22)
   dM11 = dM11 + a(j)*matmul(MInv(1:NSC,1:NSC),dM21)

enddo

Mbt(1:NSC,1:NSC)             = Mb11
Mbt(1:NSC,NSC+1:2*NSC)       = Mb12
Mbt(NSC+1:2*NSC,1:NSC)       = Mb21
Mbt(NSC+1:2*NSC,NSC+1:2*NSC) = Mb22

dMt(1:NSC,1:NSC)             = dM11
dMt(1:NSC,NSC+1:2*NSC)       = dM12
dMt(NSC+1:2*NSC,1:NSC)       = dM21
dMt(NSC+1:2*NSC,NSC+1:2*NSC) = dM22

end subroutine integrator
!######################################################!

!######################################################!
!             Compute MPI iteration range              !       
!######################################################!
subroutine para_range(n1,n2,nprocs,irank,ista,iend)
implicit none
integer, intent(in)  :: n1, n2, nprocs, irank
integer, intent(out) :: ista, iend
integer              :: iwork1, iwork2

iwork1 = (n2 - n1 + 1)/nprocs
iwork2 = mod(n2 - n1 + 1, nprocs)
ista   = irank*iwork1 + n1 + min(irank, iwork2)
iend   = ista + iwork1 -1

if (iwork2 .gt. irank) then
        iend = iend + 1
endif

return

end subroutine para_range
!######################################################!
