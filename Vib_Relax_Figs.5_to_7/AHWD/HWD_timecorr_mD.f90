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

real*8                  :: normC, dR
real*8, allocatable     :: PScoord(:,:), qp(:,:,:)
real*8, allocatable     :: Sfwd(:), grid(:), summ(:)
real*8, allocatable     :: Mbt(:,:,:), dMt(:,:,:)
real*8, allocatable     :: Mb(:,:), dM(:,:)

complex*16              :: OverlapIn, pref, prev
complex*16, allocatable :: TCFsum(:,:), posn(:), TCF(:), TCFT(:,:)

! Read input file
call input

allocate(PScoord(4,Ndof),qp(4,Ndof,0:Ntime),TCFsum(NR,0:Ntime))
allocate(Sfwd(0:Ntime),grid(NR))
allocate(Mbt(2*NSC,2*NSC,0:Ntime),dMt(2*NSC,2*NSC,0:Ntime))
allocate(Mb(2*NSC,2*NSC),dM(2*NSC,2*NSC))
allocate(TCFT(NR,0:Ntime),summ(0:Ntime))
allocate(posn(NR),TCF(NR))

badE    = 0
badtotE = 0
badS    = 0
badtotS = 0
summ    = 0.d0
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
grid    = 0.d0

! grid of final nuclear R dist
dR = (Rmax-Rmin)/NR
do i = 1, NR
   grid(i) = Rmin + dble(i-1)*dR
enddo

! Normalization constant from sampling
normC   = (1.d0/16.d0/pi**2)**NSC
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
  OverlapIn = cdexp(+0.50d0*Iu*((PSCoord(2,1)+pIn)*PSCoord(3,1))) &
            * cdexp(+0.50d0*Iu*((PSCoord(1,1)-qIn)*PSCoord(4,1)))
  do j = 0, Ntime
    
    ! Monodromy matrices for prefactor
    Mb(:,:) = Mbt(:,:,j)
    dM(:,:) = dMt(:,:,j)
    ! Compute prefactor
    call HWDpref(Mb,dM,pref)
    !Track branch cut
    if (dble(pref).lt.0.d0 .and. dimag(pref)*dimag(prev).lt.0d0) Maslov = Maslov + 1
    prev  = pref
    ! CS matrix element of B
    call DF_B(qp(:,1,j),grid,posn)
    ! TCF estimator
    TCF = OverlapIn*(-1)**Maslov*cdsqrt(pref)&
        * posn*cdexp(Iu*Sfwd(j))/(2.d0*pi)**(2*NSC)/normC
       
    TCFsum(:,j) = TCFsum(:,j) + TCF
  enddo

  222 continue

enddo

call mpi_barrier(mpi_comm_world,ierr)
!COLLECTING DATA FROM DIFFERENT PROCESSES
call mpi_reduce(TCFsum,TCFT,NR*(Ntime+1),mpi_double_complex,mpi_sum,0,mpi_comm_world,ierr)
call mpi_reduce(badE,badtotE,1,mpi_integer,mpi_sum,0,mpi_comm_world,ierr)
call mpi_reduce(badS,badtotS,1,mpi_integer,mpi_sum,0,mpi_comm_world,ierr)

if (myrank.eq.0) then

 open(555,file='TCF.out',status='unknown')
   TCFT = TCFT/(NumberMCsteps-badtotE-badtotS)
   ! normalize nuclear distribution for each time
   do j = 0, Ntime
     summ(j) = 0.d0
     do i = 1, NR
       if (i.eq.1) then
         summ(j) = dble(TCFT(i,j))
       elseif(i.eq.NR) then
         summ(j) = summ(j) + dble(TCFT(i,j))
       else
         summ(j) = summ(j) + 2.d0*dble(TCFT(i,j))
       endif
     enddo
     summ(j) = summ(j)*dR/2.d0 
     ! write nuclear P distribution at each m
     do i = 1, NR
        write(555,'(e15.6)',advance='no'), dble(j)*Timestep/fstoau
        write(555,'(e15.6)',advance='no'), grid(i)/Angtoau
        write(555,'(e15.6)',advance='no'), dble(TCFT(i,j))/summ(j)
        write(555,'(e15.6)'), dimag(TCFT(i,j))/summ(j)
     enddo
   !Also write nuclear distribution at final time in a separate file:
     if (j.eq.96) then  ! t = 192 fs
        do i = 1, NR
           write(192,'(e15.6)',advance='no'), grid(i)/Angtoau
           write(192,'(e15.6)',advance='no'), dble(TCFT(i,j))/summ(j)
           write(192,'(e15.6)'), dimag(TCFT(i,j))/summ(j)
        enddo
      elseif (j.eq.320) then !t = 640 fs
        do i = 1, NR
           write(640,'(e15.6)',advance='no'), grid(i)/Angtoau
           write(640,'(e15.6)',advance='no'), dble(TCFT(i,j))/summ(j)
           write(640,'(e15.6)'), dimag(TCFT(i,j))/summ(j)
        enddo
      elseif (j.eq.800) then !t = 1600 fs 
        do i = 1, NR
           write(1600,'(e15.6)',advance='no'), grid(i)/Angtoau
           write(1600,'(e15.6)',advance='no'), dble(TCFT(i,j))/summ(j)
           write(1600,'(e15.6)'), dimag(TCFT(i,j))/summ(j)
        enddo
     end if
   end do

 print*, "% of broken energy trajectories", badtotE, dble(badtotE)/dble(NumberMCsteps)*1.d2
 print*, "% of broken symplecticity trajectories", badtotS, dble(badtotS)/dble(NumberMCsteps)*1.d2

endif
call mpi_finalize(ierr)

end program
