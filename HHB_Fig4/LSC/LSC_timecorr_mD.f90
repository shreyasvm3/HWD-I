!#######################################################################
!#######################################################################

! This program uses LSC-IVR to compute the TCF of a multi-body system

! The output file  'TCF.out'  contains the real and imaginary
! parts of the TCF as a function of time.

!#######################################################################
!#######################################################################   
!Last Updated 10/13/2018 M.Church msc336@cornell.edu
program computeTCF
use parameters
implicit none

include 'mpif.h'

integer                 :: i, j, imc, flag, doflabel, bad, badtot
integer                 :: ierr, myrank, nprocs, istart, iend

real*8                  :: TCF, pos, normC
real*8, allocatable     :: PScoord(:,:,:), qpfwd(:,:,:), TCFsum(:)
real*8, allocatable     :: InitialP(:), InitialQ(:), Sfwd(:)
real*8, allocatable     :: Mfwd(:,:,:,:), Mpfwd(:,:,:,:), rTCF(:)
real*8, allocatable     :: MonodromyFwd(:,:,:), MonodromyBck(:,:,:)

! Read input file
call input

allocate(PScoord(2,Ndof,NumberMCsteps),qpfwd(2,Ndof,0:Ntime),TCFsum(0:Ntime))
allocate(InitialP(Ndof),InitialQ(Ndof),Sfwd(0:Ntime),Mfwd(Ndof,Ndof,4,0:Ntime))
allocate(Mpfwd(Ndof,Ndof,4,0:Ntime),rTCF(0:Ntime),MonodromyFwd(Ndof,Ndof,4))
allocate(MonodromyBck(Ndof,Ndof,4))

! Indicate which degree of freedom to observe
doflabel = 1

bad     = 0
badtot  = 0
TCFsum  = 0.d0
rTCF    = 0.d0

! Normalization constant from sampling
normC   = 1.d0
do i = 1, Ndof
  normC = normC/pi
enddo

! Sample initial conditions
call sample_LSC(PScoord) 

! Set up parallelization
call mpi_init(ierr)                                             ! Create child processes
call mpi_comm_size(mpi_comm_world,nprocs,ierr)                  ! Find the total number of processes
call mpi_comm_rank(mpi_comm_world,myrank,ierr)                  ! Find the rank of each process
call para_range(1,NumberMCsteps,nprocs,myrank,istart,iend)      ! Divide configurations between processes
call mpi_barrier(mpi_comm_world,ierr)

do imc = istart, iend  ! Loop over initial conditions

   InitialQ = PScoord(1,:,imc)
   InitialP = PScoord(2,:,imc)

   ! Compute trajectory
   flag = 0
   call PropagateFwd_DF(InitialQ,InitialP,qpfwd,Sfwd,Mfwd,flag)
   if (flag.eq.1) then
     bad = bad + 1
     goto 222
   endif

   do j = 0, Ntime

        ! Wigner transform of B
        call LSC_Bq(doflabel,qpfwd(1,:,j),pos)

        ! TCF estimator
        TCFsum(j) = TCFsum(j) + 2.d0**Ndof*pos/normC/(2.d0*pi)**Ndof
   
   enddo

   222 continue

enddo

call mpi_barrier(mpi_comm_world,ierr)

! Collect data
call mpi_reduce(TCFsum,rTCF,Ntime+1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
call mpi_reduce(bad,badtot,1,mpi_integer,mpi_sum,0,mpi_comm_world,ierr)

if (myrank.eq.0) then

 open(555,file='TCF.out',status='unknown')

 do j = 0, Ntime
   rTCF(j) = rTCF(j)/(NumberMCsteps-badtot)
   write(555,'(f14.6 :2e25.8)') j*TimeStep, rTCF(j)
 enddo

 close(555)

endif

call mpi_finalize(ierr)

end program computeTCF

