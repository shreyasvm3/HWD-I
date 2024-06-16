!#######################################################################
!#######################################################################

! This program uses aMQC-IVR to compute the TCF of a multi-body system

! The default setting is to quantize the anharmonic mode and treat
! the bath mode in the classical limit

! The output file  'TCF.out'  contains the real and imaginary
! parts of the TCF as a function of time.

!#######################################################################
!#######################################################################   
!Last Updated 07/07/2023 S.Nandi and S.Malpathak
program computeTCF
use parameters
implicit none

include 'mpif.h'

integer                 :: i, j, imc, flag, Maslov, bad, badtot, Maslovp
integer                 :: ierr, myrank, nprocs, istart, iend, doflabel

real*8                  :: normC
real*8, allocatable     :: PScoord(:,:), qpfwd(:,:,:)
real*8, allocatable     :: qppfwd(:,:,:), Spfwd(:), iTCF(:)
real*8, allocatable     :: InitialP(:), InitialQ(:), Sfwd(:)
real*8, allocatable     :: InitialPp(:), InitialQp(:), qav(:), pav(:)
real*8, allocatable     :: Mfwd(:,:,:,:), Mpfwd(:,:,:,:), rTCF(:)
real*8, allocatable     :: MonodromyFwd(:,:,:), MonodromyBck(:,:,:)

complex*16              :: TCF, OverlapIn, OverlapS, pref, posn, prev
complex*16              :: Overlap, OverlapP, prefp, prevp
complex*16, allocatable :: TCFsum(:)

! Read input file
call input

allocate(PScoord(4,Ndof),qpfwd(2,Ndof,0:Ntime),TCFsum(0:Ntime))
allocate(qppfwd(2,Ndof,0:Ntime),Spfwd(0:Ntime),InitialPp(Ndof),InitialQp(Ndof))
allocate(InitialP(Ndof),InitialQ(Ndof),Sfwd(0:Ntime),Mfwd(Ndof,Ndof,4,0:Ntime))
allocate(Mpfwd(Ndof,Ndof,4,0:Ntime),rTCF(0:Ntime),MonodromyFwd(Ndof,Ndof,4))
allocate(MonodromyBck(Ndof,Ndof,4),iTCF(0:Ntime))
allocate(qav(Ndof),pav(Ndof))

bad     = 0
badtot  = 0
TCFsum  = 0.d0
rTCF    = 0.d0
iTCF    = 0.d0

! Normalization constant from sampling
normC = 1.d0/16.d0/pi**2/2.d0/pi


! Set up parallelization
call mpi_init(ierr)                                             ! Create child processes
call mpi_comm_size(mpi_comm_world,nprocs,ierr)                  ! Find the total number of processes
call mpi_comm_rank(mpi_comm_world,myrank,ierr)                  ! Find the rank of each process
call para_range(1,NumberMCsteps,nprocs,myrank,istart,iend)      ! Divide configurations between processes
call mpi_barrier(mpi_comm_world,ierr)

call init_random_seed(myrank)

do imc = istart, iend ! Loop over initial conditions

  ! Sampling initial conditions
  call sample_aMQC(PScoord)

  Maslov  = 0
  Maslovp = 0
  prev    = 1.d0
  prevp   = 1.d0

  InitialQ  = PScoord(1,:)
  InitialP  = PScoord(2,:)
  InitialQp = PScoord(3,:)
  InitialPp = PScoord(4,:)

  ! Compute trajectories
  flag = 0
  call PropagateFwd_DF(InitialQ,InitialP,qpfwd,Sfwd,Mfwd,flag)
  if (flag.eq.1) then
    bad = bad + 1
    goto 222
  endif

  flag = 0
  call PropagateFwd_DF(InitialQp,InitialPp,qppfwd,Spfwd,Mpfwd,flag)
  if (flag.eq.1) then
    bad = bad + 1
    goto 222
  endif

  ! CS matrix element of A, corrected for sampling
  OverlapIn = cdexp(+0.50d0*Iu*(Initialp(1)+pIn(1))*(Initialq(1)-qIn(1)))*&
              cdexp(-0.50d0*Iu*(Initialpp(1)+pIn(1))*(Initialqp(1)-qIn(1)))

  do j = 0, Ntime

    ! CS matrix element of B
    call DF_Bq(1,qpfwd(1,:,j),qpfwd(2,:,j),qppfwd(1,:,j),qppfwd(2,:,j),posn)

    ! Monodromy matrices for prefactor
    MonodromyFwd(:,:,:) = Mfwd(:,:,:,j)

    MonodromyBck(:,:,1) =  transpose(Mpfwd(:,:,4,j))
    MonodromyBck(:,:,2) = -transpose(Mpfwd(:,:,2,j))
    MonodromyBck(:,:,3) = -transpose(Mpfwd(:,:,3,j))
    MonodromyBck(:,:,4) =  transpose(Mpfwd(:,:,1,j))

    ! Compute prefactor (only for system M block)
    call DHKpref(MonodromyFwd,MonodromyBck,pref,prefp)

    !Track branch cut
    if (dble(pref).lt.0.d0 .and. dimag(pref)*dimag(prev).lt.0d0) Maslov = Maslov + 1
    prev  = pref

    if (dble(prefp).lt.0.d0 .and. dimag(prefp)*dimag(prevp).lt.0.d0) Maslovp = Maslovp + 1
    prevp = prefp        

 
    ! TCF estimator
    TCF = OverlapIn*(-1)**Maslov*cdsqrt(pref)*(-1)**Maslovp*cdsqrt(prefp)*posn*&
          cdexp(Iu*(Sfwd(j)-Spfwd(j)))/normC/(2.d0*pi)**(Ndof+1)

    TCFsum(j) = TCFsum(j) + TCF

  enddo

  222 continue

enddo

call mpi_barrier(mpi_comm_world,ierr)

! Collect data
call mpi_reduce(dble(TCFsum),rTCF,Ntime+1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
call mpi_reduce(dimag(TCFsum),iTCF,Ntime+1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
call mpi_reduce(bad,badtot,1,mpi_integer,mpi_sum,0,mpi_comm_world,ierr)

if (myrank.eq.0) then

open(555,file='TCF.out',status='unknown')

do j = 0, Ntime
  rTCF(j) = rTCF(j)/dble(NumberMCsteps-badtot)
  iTCF(j) = iTCF(j)/dble(NumberMCsteps-badtot)
  write(555,'(f14.6 :2e25.8)') j*TimeStep, rTCF(j), iTCF(j)
enddo

close(555)

endif

call mpi_finalize(ierr)

end program
