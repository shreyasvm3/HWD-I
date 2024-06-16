!########################################################
!########################################################                    

! This file contains subroutines that 
! compute coherent state matrix elements
! and SC prefactors for multi-body sysmtes                                                                                               
! Subroutine                       Label         Line #                          
! ----------                       -----         ------                          
!                                                                                  
! Init. Coherent state overlap     CSoverlap     41                                   
!                                                                                  
! Operator B  matrix elements:                                             
!   - Forward-Backward: B = B(q)   FB_Bq         65                                               
!   - Forward-Backward: B = B(p)   FB_Bp         82                                    
!   - Double-Forward  : B = B(q)   DF_Bq         99                                   
!   - Double-Forward  : B = B(p)   DF_Bp         126                                     
!   - Husimi-IVR      : B = B(q)   Hus_Bq        153
!   - Husimi-IVR      : B = B(p)   Hus_Bp        170
!   - LSC-IVR         : B = B(q)   LSC_Bq        187                                
!
! MQC-IVR Prefactors:                                                     
!   - Forward-Backward: B = B(q)   MQCprefFB_Bq  203                                    
!   - Forward-Backward: B = B(p)   MQCprefFB_Bp  253                               
!   - Forward-Backward: General B  MQCprefFB     303                               
!   - Double-Forward  : General B  MQCprefDF     370 
!
! aMQC-IVR Prefactor               aMQCpref      431 
!                               
! DHK-IVR Prefactor                DHKpref       527 
!
! Symplectic Integrator            integrator    607
!
! MPI Iteration Range              para_range    675                                                     
!                                                                                  
!######################################################
!Last Updated 10/13/2018 M.Church msc336@cornell.edu

!######################################################!
!   Init. Coherent State Overlap = <p0,q0|pIn,qIn>     !
!######################################################!
subroutine CSoverlap(q,p,WF)
use parameters, only : Width0, InverseWidth0, qIn, pIn, Iu, Ndof
implicit none

integer                 :: i
real*8, intent(in)      :: q(Ndof), p(Ndof)
real*8                  :: qp(Ndof), pp(Ndof)
complex*16, intent(out) :: WF

do i = 1, Ndof
   qp(i) = sum(Width0(i,:)*(q-qIn))
   pp(i) = sum(InverseWidth0(i,:)*(p-pIn))
enddo

WF =  dexp(-dot_product(q-qIn,qp)/4.d0-dot_product(p-pIn,pp)/4.d0)*&
     cdexp(Iu*dot_product(p+pIn,q-qIn)/2.d0)

end subroutine
!######################################################!

!######################################################!
!          Matrix element of B = B(q) = q              !
!                  FB MQC-IVR                          !
!######################################################!
subroutine FB_Bq(lab,q,Bq)
use parameters, only : Ndof
implicit none

real*8, intent(in)      ::      q(Ndof)
integer, intent(in)     ::      lab
real*8, intent(out)     ::      Bq

Bq = q(lab)

end subroutine
!######################################################!

!######################################################!
!       Matrix element of B = B(p) = p                 !
!               FB MQC-IVR                             !
!######################################################!
subroutine FB_Bp(lab,p,Bp)
use parameters, only : Ndof
implicit none

real*8, intent(in)      ::      p(Ndof)
integer, intent(in)     ::      lab
real*8, intent(out)     ::      Bp

Bp = p(lab)

end subroutine
!######################################################!

!######################################################!
!            Matrix element of B = HO Energy           !
!                   DF MQC-IVR                         !
!######################################################!
subroutine DF_B(q,p,qp,pp,Et)
use parameters, only : Iu, WidthT, InverseWidthT, Ndof, Mass, coupling
implicit none

real*8, intent(in)      :: q(Ndof), p(Ndof), qp(Ndof), pp(Ndof)
complex*16, intent(out) :: Et
integer                 :: i
real*8                  :: q1(Ndof), p1(Ndof), wbath, ky
complex*16              :: xi(Ndof), eta(Ndof)
    
Et      = 0.d0
xi      = ((q+qp) + Iu*matmul(InverseWidthT,(p-pp)))/2.d0
eta     = ((p+pp) - Iu*matmul(WidthT,(q-qp)))/2.d0
wbath   = 0.333d0 
ky      = Mass(2,2)*wbath**2

Et = 0.5d0*(eta(1)**2 + WidthT(1,1)/2.d0)/Mass(1,1) & ! px^2/2mx
   + 0.5d0*(eta(2)**2 + WidthT(2,2)/2.d0)/Mass(2,2) & ! py^2/2my
   + (xi(1)**2 + 0.5d0/WidthT(1,1)) - 0.1d0*(xi(1)**3 + 1.5d0*xi(1)/WidthT(1,1)) & !x^2 - 0.1*x^3
   + 0.1d0*(xi(1)**4 + 3.d0*xi(1)**2/WidthT(1,1) + 7.5d-1/WidthT(1,1)**2)     &    ! + 0.1 x^4
   + coupling*xi(1)*xi(2) + 0.5d0*ky*(xi(2)**2 + 0.5d0/WidthT(2,2))                ! kxy + 0.5*my*wy^2*y^2

do i = 1, Ndof
   q1(i) = sum(WidthT(i,:)*(q-qp))
   p1(i) = sum(InverseWidthT(i,:)*(p-pp))
enddo

Et = Et*dexp(-0.25d0*dot_product(q-qp,q1)-0.25d0*dot_product(p-pp,p1))&
       *cdexp(-0.5d0*Iu*dot_product(p+pp,q-qp))

end subroutine
!######################################################!
!             Matrix element of B = p                  !
!                   DF MQC-IVR                         !
!######################################################!
subroutine DF_Bp(lab,q,p,qp,pp,Bp)
use parameters, only : Iu, WidthT, InverseWidthT, Ndof
implicit none

real*8, intent(in)      ::      q(Ndof), p(Ndof), qp(Ndof), pp(Ndof)
integer, intent(in)     ::      lab
complex*16, intent(out) ::      Bp

integer                 ::      i
real*8                  ::      q1(Ndof), p1(Ndof)

do i = 1, Ndof
   q1(i) = sum(WidthT(i,:)*(q-qp))
   p1(i) = sum(InverseWidthT(i,:)*(p-pp))
enddo

Bp = 0.5d0*((pp(lab)+p(lab))-Iu*WidthT(lab,lab)*(q(lab)-qp(lab)))&
    *dexp(-0.25d0*dot_product(q-qp,q1)-0.25d0*dot_product(p-pp,p1))&
   *cdexp(-0.5d0*Iu*dot_product(p+pp,q-qp))

end subroutine
!######################################################!

!######################################################!
!              Matrix element of B = q                 !
!                      Husimi-IVR                      !
!######################################################!
subroutine Hus_Bq(lab,q,Bq)
use parameters, only : Ndof
implicit none

real*8, intent(in)      ::      q(Ndof)
integer, intent(in)     ::      lab
real*8, intent(out)     ::      Bq

Bq = q(lab)

end subroutine
!######################################################!

!######################################################!
!              Matrix element of B = q                 !
!                      Husimi-IVR                      !
!######################################################!
subroutine Hus_Bp(lab,p,Bp)
use parameters, only : Ndof
implicit none

real*8, intent(in)      ::      p(Ndof)
integer, intent(in)     ::      lab
real*8, intent(out)     ::      Bp

Bp = p(lab)

end subroutine
!######################################################!

!######################################################!
!              Matrix element of B = q                 !
!                      LSC-IVR                         !
!######################################################!
subroutine LSC_Bq(lab,q,Bq)
use parameters, only : Ndof
implicit none

real*8, intent(in)      ::      q(Ndof)
integer, intent(in)     ::      lab
real*8, intent(out)     ::      Bq

Bq = q(lab)

end subroutine
!######################################################!

!######################################################!
!            FB MQC-IVR Prefactor with B=B(q)          !
!######################################################!
subroutine MQCprefFB_Bq(Mfwd,Mbck,pref)
use parameters, only : Width0, InverseWidth0, TuningP, InverseTuningP, Iu, Ndof
implicit none

real*8, intent(in)      :: Mfwd(Ndof,Ndof,4), Mbck(Ndof,Ndof,4) 
complex*16, intent(out) :: pref

complex*16              :: A(Ndof,Ndof), B(Ndof,Ndof), C(Ndof,Ndof)
complex*16              :: D(Ndof,Ndof), E(Ndof,Ndof), F(Ndof,Ndof)
complex*16              :: G(Ndof,Ndof), pref1(Ndof,Ndof)

integer                 :: i, info, sgn
integer, allocatable    :: ipiv(:)

! Monodromy matrix arrays defined as
! Mbck = (MqqB,MqpB,MpqB,MppB)
! Mfwd = (MqqF,MqpF,MpqF,MppF)

A = Mbck(:,:,4)-Iu*matmul(Width0,Mbck(:,:,2))
B = matmul(Mfwd(:,:,4),Width0)+Iu*Mfwd(:,:,3)
C = matmul(Width0,Mbck(:,:,1))+Iu*Mbck(:,:,3)
D = Mfwd(:,:,1)-Iu*matmul(Mfwd(:,:,2),Width0)
E = matmul(InverseTuningP,D)
F = matmul(A,B+E)
G = matmul(C,D)

pref1 = matmul(TuningP,F+G)
pref1 = matmul(InverseWidth0,pref1)

!Calculate determinant using LU decomposition
allocate(ipiv(Ndof))

call zgetrf(Ndof,Ndof,pref1,Ndof,ipiv,info)

sgn   = 1
pref  = 1.d0

do i = 1, Ndof
   pref = pref*pref1(i,i)
   if (ipiv(i).ne.i) sgn = -sgn
enddo

pref = pref*sgn

end subroutine
!######################################################!

!######################################################!
!            FB MQC-IVR Prefactor with B=B(p)          !
!######################################################!
subroutine MQCprefFB_Bp(Mfwd,Mbck,pref)
use parameters, only : Width0, InverseWidth0, TuningQ, InverseTuningQ, Iu, Ndof
implicit none

real*8, intent(in)      :: Mfwd(Ndof,Ndof,4), Mbck(Ndof,Ndof,4)
complex*16, intent(out) :: pref

complex*16              :: A(Ndof,Ndof), B(Ndof,Ndof), C(Ndof,Ndof)
complex*16              :: D(Ndof,Ndof), E(Ndof,Ndof), F(Ndof,Ndof)
complex*16              :: G(Ndof,Ndof), pref1(Ndof,Ndof)

integer                 :: i, info, sgn
integer, allocatable    :: ipiv(:)

! Monodromy matrix arrays defined as
! Mbck = (MqqB,MqpB,MpqB,MppB)
! Mfwd = (MqqF,MqpF,MpqF,MppF)

A = Mbck(:,:,4)-Iu*matmul(Width0,Mbck(:,:,2))
B = matmul(Mfwd(:,:,4),Width0)+Iu*Mfwd(:,:,3)
C = matmul(Width0,Mbck(:,:,1))+Iu*Mbck(:,:,3)
D = Mfwd(:,:,1)-Iu*matmul(Mfwd(:,:,2),Width0)
E = matmul(InverseTuningQ,B)
F = matmul(C,D+E)
G = matmul(A,B)

pref1 = matmul(TuningQ,F+G)
pref1 = matmul(InverseWidth0,pref1)

!Calculate determinant using LU decomposition
allocate(ipiv(Ndof))

call zgetrf(Ndof,Ndof,pref1,Ndof,ipiv,info)

sgn   = 1
pref  = 1.d0

do i = 1, Ndof
   pref = pref*pref1(i,i)
   if (ipiv(i).ne.i) sgn = -sgn
enddo

pref = pref*sgn

end subroutine
!######################################################!

!######################################################!
!            FB MQC-IVR Prefactor for general B        !
!######################################################!
subroutine MQCprefFB(Mfwd,Mbck,pref)
use parameters, only : Iu, Width0, InverseWidth0, WidthT, InverseWidthT,&
TuningQ, InverseTuningQ, TuningP, InverseTuningP, Ndof
implicit none

real*8, intent(in)      :: Mfwd(Ndof,Ndof,4), Mbck(Ndof,Ndof,4)
complex*16, intent(out) :: pref

real*8                  :: Id(Ndof,Ndof), G(Ndof,Ndof), GInv(Ndof,Ndof)
real*8                  :: G1(Ndof,Ndof), Gp(Ndof,Ndof), Gq(Ndof,Ndof)

complex*16              :: A(Ndof,Ndof), B(Ndof,Ndof), C(Ndof,Ndof)
complex*16              :: D(Ndof,Ndof), pref1(Ndof,Ndof)

integer                 :: i, info, sgn
integer, allocatable    :: ipiv(:)

! Monodromy matrix arrays defined as
! Mbck = (MqqB,MqpB,MpqB,MppB)
! Mfwd = (MqqF,MqpF,MpqF,MppF)

G    = matmul(TuningQ+WidthT,TuningP) + matmul(TuningQ,InverseWidthT+TuningP)
GInv = 0.d0
Id   = 0.d0

do i = 1, Ndof
        Ginv(i,i) = 1.d0/G(i,i)
        Id(i,i)   = 1.d0
enddo

G1   = Ginv+Id
Gp   = matmul(0.5d0*InverseWidthT + TuningP,GInv)
Gq   = matmul(0.5d0*WidthT + TuningQ,GInv)

A  = Mbck(:,:,4)-Iu*matmul(Width0,Mbck(:,:,2))
B  = matmul(Mfwd(:,:,4),Width0)+Iu*Mfwd(:,:,3)
C  = matmul(Width0,Mbck(:,:,1))+Iu*Mbck(:,:,3)
D  = Mfwd(:,:,1)-Iu*matmul(Mfwd(:,:,2),Width0)

pref1 = 0.5d0*matmul(matmul(A,G1),B) + &
              matmul(matmul(C,Gp),B) + &
        0.5d0*matmul(matmul(C,G1),D) + &
              matmul(matmul(A,Gq),D)

pref1 = matmul(matmul(0.50d0*InverseWidth0,G),pref1)

!Calculate determinant using LU decomposition
allocate(ipiv(Ndof))

call zgetrf(Ndof,Ndof,pref1,Ndof,ipiv,info)

sgn   = 1
pref  = 1.d0

do i = 1, Ndof
   pref = pref*pref1(i,i)
   if (ipiv(i).ne.i) sgn = -sgn
enddo

pref = pref*sgn

end subroutine
!######################################################!

!######################################################!
!            DF MQC-IVR Prefactor                      !
!######################################################!
subroutine MQCprefDF(Mfwd,Mbck,coeff)
use parameters, only : Width0, TuningQ, TuningP, InverseWidth0, InverseTuningQ, Iu, Ndof,&
WidthT, InverseWidthT, InverseTuningP
implicit none

real*8, intent(in)               :: Mfwd(Ndof,Ndof,4), Mbck(Ndof,Ndof,4)
complex*16, intent(out)          :: coeff
integer                          :: i, info, sgn
complex*16, dimension(Ndof,Ndof) :: t1, t2, t3, t4, sq, sp, temp, temp2
complex*16, dimension(Ndof,Ndof) :: A, B, C, D, G, InverseG, coeff2, Um
complex*16, allocatable          :: W(:), VL(:,:), VR(:,:), WORK(:)
real*8, allocatable              :: RWORK(:)
integer, allocatable             :: ipiv(:)

A = Mfwd(:,:,4)-Iu*matmul(WidthT,Mfwd(:,:,2))
B = matmul(Mbck(:,:,4),WidthT)+Iu*Mbck(:,:,3)
C = matmul(WidthT,Mfwd(:,:,1))+Iu*Mfwd(:,:,3)
D = Mbck(:,:,1)-Iu*matmul(Mbck(:,:,2),WidthT)
G = matmul(TuningP,TuningQ+Width0)+matmul(TuningQ,TuningP+InverseWidth0)

Um=0.d0
InverseG=0.d0

do i=1,Ndof
   Um(i,i)=1.d0
   InverseG(i,i)=1.d0/G(i,i)
enddo

sq=matmul(0.5d0*Width0+TuningQ,InverseG)
sp=matmul(0.5d0*InverseWidth0+TuningP,InverseG)

t1=matmul(InverseG+Um,B)
t2=matmul(sp,B)
t3=Matmul(InverseG+Um,D)
t4=matmul(sq,D)

temp=0.5d0*matmul(A,t1)+matmul(C,t2)+0.5d0*matmul(C,t3)+matmul(A,t4)
temp2=matmul(G,temp)
coeff2=matmul(0.5d0*InverseWidthT,temp2)

!Calculate determinant using LU decomposition
allocate(ipiv(Ndof))

call zgetrf(Ndof,Ndof,coeff2,Ndof,ipiv,info)

sgn=1
coeff=1.d0

do i=1, Ndof
   coeff=coeff*coeff2(i,i)
   if (ipiv(i).ne.i) sgn=-sgn
enddo

coeff=coeff*sgn

end subroutine
!######################################################!

!######################################################!
!               aMQC-IVR Prefactor   
!######################################################!
subroutine aMQCpref(M,Mp,coeff)
use parameters, only : Width0, WidthT, InverseWidth0, InverseWidthT, Iu, Ndof
implicit none

real*8, intent(in)      :: M(Ndof,Ndof,4), Mp(Ndof,Ndof,4)
complex*16, intent(out) :: coeff

real*8                  :: Mf(2*Ndof,2*Ndof), Mfp(2*Ndof,2*Ndof)
complex*16              :: Jmat(4*Ndof,4*Ndof)
real*8                  :: Imat1(Ndof,Ndof)
real*8                  :: Imat2(2*Ndof,2*Ndof), Imat4(4*Ndof,4*ndof)
complex*16              :: Xmat(2*Ndof,2*Ndof), Ymat(2*Ndof,2*Ndof)
complex*16              :: Kmat(4*Ndof,4*Ndof)
complex*16, allocatable :: work(:)

integer                 :: i, info, sgn, lwork
integer, allocatable    :: ipiv(:)

lwork=4*Ndof
allocate(work(lwork))

! Identity matrices
Imat1 = 0.d0
Imat2 = 0.d0
Imat4 = 0.d0
do i = 1, Ndof
 Imat1(i,i)           = 1.d0
 Imat2(i,i)           = 1.d0
 Imat2(i+Ndof,i+Ndof) = 1.d0
enddo
Imat4 = 0.d0
do i = 1, 4*Ndof
 Imat4(i,i) = 1.d0
enddo

Mf(1:Ndof,1:Ndof)                = M(:,:,1)
Mf(1:Ndof,Ndof+1:2*Ndof)         = M(:,:,2)
Mf(Ndof+1:2*Ndof,1:Ndof)         = M(:,:,3)
Mf(Ndof+1:2*Ndof,Ndof+1:2*Ndof)  = M(:,:,4)
Mfp(1:Ndof,1:Ndof)               = Mp(:,:,1)
Mfp(1:Ndof,Ndof+1:2*Ndof)        = Mp(:,:,2)
Mfp(Ndof+1:2*Ndof,1:Ndof)        = Mp(:,:,3)
Mfp(Ndof+1:2*Ndof,Ndof+1:2*Ndof) = Mp(:,:,4)

Xmat = 0.d0
Xmat(1:Ndof,1:Ndof)               = +0.5d0*Iu*Width0
Xmat(1:Ndof,Ndof+1:2*Ndof)        = -0.5d0*Imat1
Xmat(Ndof+1:2*Ndof,1:Ndof)        = +0.5d0*Imat1
Xmat(Ndof+1:2*Ndof,Ndof+1:2*Ndof) = +0.5d0*Iu*InverseWidth0

Ymat = 0.d0
Ymat(1:Ndof,1:Ndof)               = +0.5d0*Iu*WidthT
Ymat(1:Ndof,Ndof+1:2*Ndof)        = +0.5d0*Imat1
Ymat(Ndof+1:2*Ndof,1:Ndof)        = -0.5d0*Imat1
Ymat(Ndof+1:2*Ndof,Ndof+1:2*Ndof) = +0.5d0*Iu*InverseWidthT

Kmat = 0.d0
Kmat(1:2*Ndof,1:2*Ndof)               = Xmat
Kmat(1:2*Ndof,2*Ndof+1:4*Ndof)        = dconjg(Xmat)
Kmat(2*Ndof+1:4*Ndof,1:2*Ndof)        = matmul(Ymat,Mfp)
Kmat(2*Ndof+1:4*Ndof,2*Ndof+1:4*Ndof) = matmul(dconjg(Ymat),Mf)

Jmat = 0.d0
Jmat(1:2*Ndof,1:2*Ndof)               = +Imat2*Iu
Jmat(1:2*Ndof,2*Ndof+1:4*Ndof)        = -Imat2*Iu
Jmat(2*Ndof+1:4*Ndof,1:2*Ndof)        =  Mfp
Jmat(2*Ndof+1:4*Ndof,2*Ndof+1:4*Ndof) = -Mf

! KT IN MIXED LIMIT = Ktil
do i = 2, Ndof
  Kmat(i,:) = Jmat(i,:)
enddo
do i = Ndof+2, 2*Ndof
  Kmat(i,:) = Jmat(i,:)
enddo

allocate(ipiv(4*Ndof))
call zgetrf(4*Ndof,4*Ndof,Kmat,4*Ndof,ipiv,info)

sgn   = 1
coeff = 1.d0

do i = 1, 4*Ndof
   coeff = coeff*Kmat(i,i)
   if (ipiv(i).ne.i) sgn = -sgn
enddo
deallocate(ipiv)

coeff = coeff*sgn*(-1.d0)**Ndof

end subroutine
!######################################################!

!######################################################!
!       Herman-Kluk Prefactors for DHK-IVR             !    
!######################################################!
subroutine DHKpref(Mfwd,Mbck,pref,prefp)
use parameters, only : Width0, WidthT, InverseWidth0, InverseWidthT, Iu, Ndof
implicit none

real*8, intent(in)      :: Mfwd(Ndof,Ndof,4), Mbck(Ndof,Ndof,4)
complex*16, intent(out) :: pref, prefp

integer                 :: i, info, sgn, info2, sgn2
integer, allocatable    :: ipiv(:), ipiv2(:)

real*8                  :: A(Ndof,Ndof), B(Ndof,Ndof), C(Ndof,Ndof)
real*8                  :: D(Ndof,Ndof), G(Ndof,Ndof), H(Ndof,Ndof), J(Ndof,Ndof)
real*8                  :: K(Ndof,Ndof)

complex*16              :: pref1(Ndof,Ndof), pref1p(Ndof,Ndof)

A  =  matmul(dsqrt(Width0),Mfwd(:,:,1))
B  =  matmul(A,dsqrt(InverseWidthT))

C  =  matmul(dsqrt(InverseWidth0),Mfwd(:,:,4))
D  =  matmul(C,dsqrt(WidthT))

G  =  matmul(dsqrt(Width0),Mfwd(:,:,2))
H  =  matmul(G,dsqrt(WidthT))

J  =  matmul(dsqrt(InverseWidth0),Mfwd(:,:,3))
K  =  matmul(J,dsqrt(InverseWidthT))

pref1  =  0.5d0*(B + D -Iu*H + Iu*K)

!Calculate determinant using LU decomposition
allocate(ipiv(Ndof))

call zgetrf(Ndof,Ndof,pref1,Ndof,ipiv,info)

sgn   = 1
pref  = 1.d0

do i = 1, Ndof
        pref = pref*pref1(i,i)
        if (ipiv(i).ne.i) sgn = -sgn
enddo

pref = pref*sgn

A  =  matmul(dsqrt(WidthT),Mbck(:,:,1))
B  =  matmul(A,dsqrt(InverseWidth0))

C  =  matmul(dsqrt(InverseWidthT),Mbck(:,:,4))
D  =  matmul(C,dsqrt(Width0))

G  =  matmul(dsqrt(WidthT),Mbck(:,:,2))
H  =  matmul(G,dsqrt(Width0))

J  =  matmul(dsqrt(InverseWidthT),Mbck(:,:,3))
K  =  matmul(J,dsqrt(InverseWidth0))

pref1p  =  0.5d0*(B + D - Iu*H + Iu*K)

!Calculate determinant using LU decomposition
allocate(ipiv2(Ndof))

call zgetrf(Ndof,Ndof,pref1p,Ndof,ipiv2,info)

sgn    = 1
prefp  = 1.d0

do i = 1, Ndof
        prefp = prefp*pref1p(i,i)
        if (ipiv2(i).ne.i) sgn = -sgn
enddo

prefp = prefp*sgn

end subroutine
!######################################################!

!######################################################!
!               Symplectic Integrator                  !       
!######################################################!
subroutine integrator(q,p,dt,s,M11,M12,M21,M22)
use parameters, only : Mass, Ndof
implicit none

real*8, intent(in)    :: dt
real*8, intent(inout) :: p(Ndof), q(Ndof), s
real*8, intent(inout) :: M11(Ndof,Ndof), M12(Ndof,Ndof), M21(Ndof,Ndof), M22(Ndof,Ndof)

integer               :: j, i, k

real*8                :: ke, v, const
real*8                :: temp(Ndof), dv(Ndof), d2v(Ndof,Ndof), a(4), b(4)

const = 1.d0/dsqrt(3.d0)
a(1)  = 0.5d0*(1.d0-const)*dt
a(2)  =  const*dt
a(3)  = -const*dt
a(4)  = 0.5d0*(1.d0+const)*dt
b(1)  = 0.d0
b(2)  = 0.5d0*(0.5d0+const)*dt
b(3)  = 0.5d0*dt
b(4)  = 0.5d0*(0.5d0-const)*dt

do j = 1, 4

   if (j.gt.1) then
      call vdv_AnHO(q,v,dv,d2v)
      p   = p - b(j)*dv
      M22 = M22 - b(j)*matmul(d2v,M12)
      M21 = M21 - b(j)*matmul(d2v,M11)
      s   = s - b(j)*v
   end if

   ke = 0.d0
   do i = 1, Ndof
      temp(i) = p(i)/Mass(i,i)
      ke = ke + 0.5d0*p(i)**2/Mass(i,i)
   enddo

   q   = q + a(j)*temp
   M12 = M12 + a(j)*divm(M22)
   M11 = M11 + a(j)*divm(M21)

   s = s + a(j)*ke

enddo

contains

function divm(M)
use parameters, only : Ndof, Mass
implicit none

integer :: i
real*8  :: divm(Ndof,Ndof), M(Ndof,Ndof)

do i = 1, Ndof
   divm(i,:) = M(i,:)/Mass(i,i)
enddo

end function divm

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
