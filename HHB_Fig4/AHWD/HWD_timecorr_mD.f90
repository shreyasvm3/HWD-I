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

integer                 :: i, j, imc, flagE, flagS, Maslov
integer                 :: ierr, myrank, nprocs, istart, iend
integer                 :: badE, badS, badtotE, badtotS

real*8                  :: normC
real*8, allocatable     :: PScoord(:,:), qp(:,:,:)
real*8, allocatable     :: Sfwd(:)
real*8, allocatable     :: Mbt(:,:,:), dMt(:,:,:)
real*8, allocatable     :: Mb(:,:), dM(:,:)

complex*16              :: OverlapIn, pref, prev
complex*16, allocatable :: TCFsum(:,:), posn(:), TCF(:), TCFT(:,:)

! Read input file
call input

allocate(PScoord(4,Ndof),qp(4,Ndof,0:Ntime),TCFsum(Ndof,0:Ntime))
allocate(Sfwd(0:Ntime))
allocate(Mbt(2*NSC,2*NSC,0:Ntime),dMt(2*NSC,2*NSC,0:Ntime))
allocate(Mb(2*NSC,2*NSC),dM(2*NSC,2*NSC))
allocate(TCFT(Ndof,0:Ntime))
allocate(posn(Ndof),TCF(Ndof))

badE    = 0
badtotE = 0
badS    = 0
badtotS = 0
TCFsum  = 0.d0
TCFT    = 0.d0
PScoord = 0.d0
posn    = 0.d0
TCF     = 0.d0
Mb      = 0.d0
dM      = 0.d0
Mbt     = 0.d0
dMt     = 0.d0
qp      = 0.d0
Sfwd    = 0.d0

! Normalization constant from sampling
normC   = (1.d0/16.d0/pi**2)**NSC*(1.d0/pi)**NCl
! Set up parallelization
call mpi_init(ierr)                                             ! Create child processes
call mpi_comm_size(mpi_comm_world,nprocs,ierr)                  ! Find the total number of processes
call mpi_comm_rank(mpi_comm_world,myrank,ierr)                  ! Find the rank of each process
call para_range(1,NumberMCsteps,nprocs,myrank,istart,iend)      ! Divide configurations between processes
call mpi_barrier(mpi_comm_world,ierr)

call init_random_seed(myrank)                                   ! Initialize random seed
do imc = istart, iend  ! Loop over initial conditions

  ! Sample initial conditions
  call sample_HWD(PScoord)
  Maslov  = 0
  prev    = 1.d0

  ! Compute trajectories
  flagE = 0
  flagS = 0

  call Propagate_HT(PScoord,qp,Sfwd,Mbt,dMt,flagE,flagS)
  if (flagE.eq.1) then
    badE = badE + 1
    goto 222
  elseif (flagS.eq.1) then
    badS = badS + 1
    goto 222
  end if
  ! rho_A^tilda, corrected for sampling along with the t=0 portion of S^tilda
  OverlapIn = cdexp(+0.50d0*Iu*dot_product(PSCoord(2,1:NSC)+pIn(1:NSC),PSCoord(3,1:NSC))) &
            * cdexp(+0.50d0*Iu*dot_product(PSCoord(1,1:NSC)-qIn(1:NSC),PSCoord(4,1:NSC)))
  do j = 0, Ntime
    
    ! CS matrix element of B
    call DF_B(qp(:,:,j),posn)
    ! Monodromy matrices for prefactor
    Mb(:,:) = Mbt(:,:,j)
    dM(:,:) = dMt(:,:,j)
    
    ! Compute prefactor
    call HWDpref(Mb,dM,pref)

    !Track branch cut
    if (dble(pref).lt.0.d0 .and. dimag(pref)*dimag(prev).lt.0d0) Maslov = Maslov + 1
    prev  = pref

    ! TCF estimator
    TCF = OverlapIn*(-1)**Maslov*cdsqrt(pref)&
        * posn*cdexp(Iu*Sfwd(j))*(2.d0)**NCl&  !2 *NCl from rho wig
        / (2.d0*pi)**(2*NSC)/(2.d0*pi)**(NCl)/normC
       
    TCFsum(:,j) = TCFsum(:,j) + TCF
  
  enddo

  222 continue

enddo

call mpi_barrier(mpi_comm_world,ierr)
!COLLECTING DATA FROM DIFFERENT PROCESSES
call mpi_reduce(TCFsum,TCFT,Ndof*(Ntime+1),mpi_double_complex,mpi_sum,0,mpi_comm_world,ierr)
call mpi_reduce(badE,badtotE,1,mpi_integer,mpi_sum,0,mpi_comm_world,ierr)
call mpi_reduce(badS,badtotS,1,mpi_integer,mpi_sum,0,mpi_comm_world,ierr)

if (myrank.eq.0) then

 open(555,file='TCF.out',status='unknown')
 TCFT = TCFT/dble(NumberMCsteps-badtotE-badtotS)

 do j = 0, Ntime

   write(555,'(E15.6)',advance='no'), dble(j*TimeStep)
   do i = 1, Ndof-1
      write(555,'(E15.6)',advance='no'), real(TCFT(i,j))
      write(555,'(E15.6)',advance='no'), dimag(TCFT(i,j))
   end do
   write(555,'(E15.6)',advance='no'), real(TCFT(Ndof,j))
   write(555,'(E15.6)'), dimag(TCFT(Ndof,j))

  enddo

 close(555)

 print*, "% of broken energy trajectories", badtotE, dble(badtotE)/dble(NumberMCsteps)*1.d2
 print*, "% of broken symplecticity trajectories", badtotS, dble(badtotS)/dble(NumberMCsteps)*1.d2

endif
call mpi_finalize(ierr)

end program
