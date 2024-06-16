!#######################################################################
!#######################################################################

! This program uses LSC-IVR to compute the TCF of a multi-body system

! The output file  'TCF.out'  contains the real and imaginary
! parts of the TCF as a function of time.

!#######################################################################
!#######################################################################   
program computeTCF
use parameters
implicit none

include 'mpif.h'

integer                 :: i, j, imc, flag, doflabel, bad, badtot
integer                 :: ierr, myrank, nprocs, istart, iend

real*8                  :: TCF, normC, init, sigmaP, sigmaQ, dR
real*8, allocatable     :: PScoord(:,:), qpfwd(:,:,:), TCFsum(:,:), summ(:)
real*8, allocatable     :: InitialP(:), InitialQ(:), Sfwd(:), pos(:), grid(:)
real*8, allocatable     :: Mfwd(:,:,:,:), Mpfwd(:,:,:,:), rTCF(:,:)
real*8, allocatable     :: MonodromyFwd(:,:,:), MonodromyBck(:,:,:)

! Read input file
call input

allocate(PScoord(2,Ndof),qpfwd(2,Ndof,0:Ntime),TCFsum(NR,0:Ntime))
allocate(InitialP(Ndof),InitialQ(Ndof),Sfwd(0:Ntime),Mfwd(Ndof,Ndof,4,0:Ntime))
allocate(Mpfwd(Ndof,Ndof,4,0:Ntime),rTCF(NR,0:Ntime),MonodromyFwd(Ndof,Ndof,4))
allocate(MonodromyBck(Ndof,Ndof,4),grid(NR),pos(NR),summ(0:Ntime))

bad     = 0
badtot  = 0
grid    = 0.d0
TCFsum  = 0.d0
rTCF    = 0.d0
summ    = 0.d0

! grid of final nuclear R dist
dR = (Rmax-Rmin)/NR
do i = 1, NR
   grid(i) = Rmin + dble(i-1)*dR
enddo


! Set up parallelization
call mpi_init(ierr)                                             ! Create child processes
call mpi_comm_size(mpi_comm_world,nprocs,ierr)                  ! Find the total number of processes
call mpi_comm_rank(mpi_comm_world,myrank,ierr)                  ! Find the rank of each process
call para_range(1,NumberMCsteps,nprocs,myrank,istart,iend)      ! Divide configurations between processes
call mpi_barrier(mpi_comm_world,ierr)
call init_random_seed(myrank)                                   ! Initialize random seed

do imc = istart, iend  ! Loop over initial conditions

   ! Sample initial conditions
   call sample_LSC(PScoord) 
   
   InitialQ = PScoord(1,:)
   InitialP = PScoord(2,:)

   ! Compute trajectory
   flag = 0
   call PropagateFwd_DF(InitialQ,InitialP,qpfwd,Sfwd,Mfwd,flag)
   if (flag.eq.1) then
     bad = bad + 1
     goto 222
   endif

   do j = 0, Ntime

       do i = 1, NR
        ! B matrix element: delta(R_t-Rf)
        call operBW(grid(i),qpfwd(1,1,j),dR,pos(i))

        ! TCF estimator
        TCFsum(i,j) = TCFsum(i,j) + pos(i)

      enddo
   enddo

   222 continue

enddo

call mpi_barrier(mpi_comm_world,ierr)

! Collect data
call mpi_reduce(TCFsum,rTCF,NR*(Ntime+1),mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
call mpi_reduce(bad,badtot,1,mpi_integer,mpi_sum,0,mpi_comm_world,ierr)

if (myrank.eq.0) then

 open(555,file='TCF.out',status='unknown')

 do j = 0, Ntime
   do i = 1, NR
      rTCF(i,j) = rTCF(i,j)/(NumberMCsteps-badtot)
   end do
 end do
   ! normalize nuclear distribution for each time
   do j = 0, Ntime
     summ(j) = 0.d0
     do i = 1, NR
       if (i.eq.1) then
         summ(j) = rTCF(i,j)
       elseif(i.eq.NR) then
         summ(j) = summ(j) + rTCF(i,j)
       else
         summ(j) = summ(j) + 2.d0*rTCF(i,j)
       endif
     enddo
     summ(j) = summ(j)*dR/2.d0
     ! write nuclear P distribution at each m
     do i = 1, NR
        write(555,'(e15.6)',advance='no'), dble(j)*Timestep/fstoau
        write(555,'(e15.6)',advance='no'), grid(i)/Angtoau
        write(555,'(e15.6)',advance='no'), rTCF(i,j)/summ(j)
        write(555,'(e15.6)'), summ(j)
     enddo
   !Also write nuclear distribution at final time in a separate file:
     if (j.eq.96) then  ! t = 192 fs
        do i = 1, NR
           write(192,'(e15.6)',advance='no'), grid(i)/Angtoau
           write(192,'(e15.6)',advance='no'), rTCF(i,j)/summ(j)
           write(192,'(e15.6)'), summ(j)
        enddo
      elseif (j.eq.320) then !t = 640 fs
        do i = 1, NR
           write(640,'(e15.6)',advance='no'), grid(i)/Angtoau
           write(640,'(e15.6)',advance='no'), rTCF(i,j)/summ(j)
           write(640,'(e15.6)'), summ(j)
        enddo
      elseif (j.eq.800) then !t = 1600 fs 
        do i = 1, NR
           write(1600,'(e15.6)',advance='no'), grid(i)/Angtoau
           write(1600,'(e15.6)',advance='no'), rTCF(i,j)/summ(j)
           write(1600,'(e15.6)'), summ(j)
        enddo
     end if
   end do

 close(555)

write(*,*), "% of broken trajs.", dble(badtot)/dble(NumberMCsteps)*1.d2

endif


call mpi_finalize(ierr)

end program computeTCF

