!#######################################################################
!#######################################################################

! This program uses DHK-IVR to compute the TCF of a multi-body system

! The output file  'TCF.out'  contains the real and imaginary
! parts of the TCF as a function of time.

!#######################################################################
!#######################################################################   
program computeTCF
use parameters
implicit none

include 'mpif.h'

integer                 :: i, j, imc, flag, Maslov, Maslovp
integer                 :: ierr, myrank, nprocs, istart, iend, doflabel
integer                 :: bad, badtot

real*8                  :: normC
real*8, allocatable     :: PScoord(:,:), qpfwd(:,:,:)
real*8, allocatable     :: qppfwd(:,:,:), Spfwd(:)
real*8, allocatable     :: InitialP(:), InitialQ(:), Sfwd(:)
real*8, allocatable     :: InitialPp(:), InitialQp(:)
real*8, allocatable     :: Mfwd(:,:,:,:), Mpfwd(:,:,:,:)
real*8, allocatable     :: MonodromyFwd(:,:,:), MonodromyBck(:,:,:)

complex*16              :: OverlapIn, pref, prefp, prev, prevP
complex*16, allocatable :: TCFsum(:,:), posn(:), TCF(:), TCFT(:,:)

! Read input file
call input

allocate(PScoord(4,Ndof),qpfwd(2,Ndof,0:Ntime),TCFsum(Ndof,0:Ntime))
allocate(qppfwd(2,Ndof,0:Ntime),Spfwd(0:Ntime),InitialPp(Ndof),InitialQp(Ndof))
allocate(InitialP(Ndof),InitialQ(Ndof),Sfwd(0:Ntime),Mfwd(Ndof,Ndof,4,0:Ntime))
allocate(Mpfwd(Ndof,Ndof,4,0:Ntime),TCFT(Ndof,0:Ntime),MonodromyFwd(Ndof,Ndof,4))
allocate(MonodromyBck(Ndof,Ndof,4),posn(Ndof),TCF(Ndof))

bad     = 0
badtot  = 0
TCFsum  = 0.d0
TCFT    = 0.d0
PScoord = 0.d0
posn    = 0.d0
TCF     = 0.d0

! Normalization constant from sampling
normC   = (1.d0/16.d0/pi**2)**Ndof

! Set up parallelization
call mpi_init(ierr)                                             ! Create child processes
call mpi_comm_size(mpi_comm_world,nprocs,ierr)                  ! Find the total number of processes
call mpi_comm_rank(mpi_comm_world,myrank,ierr)                  ! Find the rank of each process
call para_range(1,NumberMCsteps,nprocs,myrank,istart,iend)      ! Divide configurations between processes
call mpi_barrier(mpi_comm_world,ierr)

call init_random_seed()                                      ! Initialize random seed

do imc = istart, iend  ! Loop over initial conditions

  ! Sample initial conditions
  call sample_DHK(PScoord)

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
  OverlapIn = 1.d0
  do i = 1, Ndof
    OverlapIn = OverlapIn*cdexp(+0.50d0*Iu*(Initialp(i) +pIn(i))*(Initialq(i) -qIn(i)))*&
                          cdexp(-0.50d0*Iu*(Initialpp(i)+pIn(i))*(Initialqp(i)-qIn(i)))
  enddo

  do j = 0, Ntime

    ! CS matrix element of B
    call DF_B(qpfwd(1,:,j),qpfwd(2,:,j),qppfwd(1,:,j),qppfwd(2,:,j),posn(1))

    ! Monodromy matrices for prefactor
    MonodromyFwd(:,:,:) = Mfwd(:,:,:,j)

    MonodromyBck(:,:,1) =  transpose(Mpfwd(:,:,4,j))
    MonodromyBck(:,:,2) = -transpose(Mpfwd(:,:,2,j))
    MonodromyBck(:,:,3) = -transpose(Mpfwd(:,:,3,j))
    MonodromyBck(:,:,4) =  transpose(Mpfwd(:,:,1,j))

    ! Compute prefactor
    call DHKpref(MonodromyFwd,MonodromyBck,pref,prefp)

    !Track branch cut
    if (dble(pref).lt.0.d0 .and. dimag(pref)*dimag(prev).lt.0d0) Maslov = Maslov + 1
    prev  = pref

    if (dble(prefp).lt.0.d0 .and. dimag(prefp)*dimag(prevp).lt.0.d0) Maslovp = Maslovp + 1
    prevp = prefp        

    ! TCF estimator
    TCF(1) = OverlapIn*(-1)**Maslov*(-1)**Maslovp*cdsqrt(pref)*cdsqrt(prefp)*&
             posn(1)*cdexp(Iu*(Sfwd(j)-Spfwd(j)))/normC/(2.d0*pi)**(2*Ndof)

    TCFsum(1,j) = TCFsum(1,j) + TCF(1)

  enddo

  222 continue

enddo

call mpi_barrier(mpi_comm_world,ierr)

!COLLECTING DATA FROM DIFFERENT PROCESSES
call mpi_reduce(TCFsum(1,:),TCFT(1,:),(Ntime+1),mpi_double_complex,mpi_sum,0,mpi_comm_world,ierr)
call mpi_reduce(bad,badtot,1,mpi_integer,mpi_sum,0,mpi_comm_world,ierr)

if (myrank.eq.0) then

 open(555,file='TCF.out',status='unknown')
 TCFT(1,:) = TCFT(1,:)/dble(NumberMCsteps-badtot)

 do j = 0, Ntime

   write(555,'(E15.6)',advance='no'), dble(j*TimeStep)
   write(555,'(E15.6)',advance='no'), real(TCFT(1,j))
   !write(555,'(E15.6)',advance='no'), dimag(TCFT(1,j))
   !write(555,'(E15.6)',advance='no'), real(TCFT(2,j))
   write(555,'(E15.6)'), dimag(TCFT(1,j))

  enddo

 close(555)

 print*, "# of Bad trajectories", badtot, dble(badtot)/dble(NumberMCsteps)*1.d2
endif

call mpi_finalize(ierr)

end program
