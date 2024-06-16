!Last Updated 10/18/23 SHreyas Malpathak
!######################################################!
!            Matrix element of B = x                   !
!           x of all dofs is evaluated                 !
!######################################################!
subroutine DF_B(z,Bq)
use parameters, only : Iu, Width0, InverseWidth0, Ndof, NSC
implicit none

real*8, intent(in)      :: z(4,Ndof)
complex*16, intent(out) :: Bq(Ndof)
integer                 :: i
real*8,dimension(Ndof)  :: qb, pb, dq, dp  !dp_cl and dq_cl are 0.
complex*16              :: csoverlap
 
qb = z(1,:)
pb = z(2,:)
dq = z(3,:) !dq_cl is set to zero in z
dp = z(4,:) !dp_cl is set to zero in z
   
!<zt'|zt> for SC dofs. Since dq_cl and dp_cl are 0, full Width0 and InvWidth0
!matrices can be used. Same logic for last dotproduct.
!Already includes -pb_t.dq_t from S^tilda.
csoverlap = dexp(-0.25d0*dot_product(dq,matmul(Width0,dq)) &
                 -0.25d0*dot_product(dp,matmul(InverseWidth0,dp))) &
          * cdexp(-Iu*dot_product(pb,dq))

! For SC dofs, we have <zt'|B|zt> * 1_W
Bq(1:NSC) = (qb(1:NSC)+0.5d0*Iu*matmul(InverseWidth0(1:NSC,1:NSC),dp(1:NSC)))*csoverlap
! For classical dofs, we have (B_Cl)_W * <zt'|zt>
Bq(NSC+1:Ndof) = qb(NSC+1:Ndof)*csoverlap

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
      call vdv_AnHO(qb(1:NSC)+dq/2.d0,qb(NSC+1:Ndof),Vp,dVp,d2Vp)
      call vdv_AnHO(qb(1:NSC)-dq/2.d0,qb(NSC+1:Ndof),Vm,dVm,d2Vm)
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
