program mpi_mp

use mpi
use omp_lib

implicit none

integer myrank, from_id
integer nslave, cluster_size, nrecv
integer mpierr
integer mpistatus(MPI_STATUS_SIZE)
integer, parameter :: npowers = 4
integer thread_check(npowers)
integer i

logical master


real(8) myx, myx_powers(npowers)
real(8), allocatable :: x(:), x_powers(:,:)

call mpi_init(mpierr)                             ! Introduce yourself to the MPI daemon
call mpi_comm_rank(MPI_COMM_WORLD,myrank,mpierr)  ! Get your rank number, store in myrank.  Master is rank 0.
if(myrank .eq. 0) then
  master=.true.
else
  master=.false.
endif

if(master) then
  write(*,*)
  write(*,'(a)') "This MPI program is a demonstration of using OpenMP and OpenMPI concurrently."
  write(*,'(a)') "It is run using: mpirun -n # <PATH_TO>/mpi_mp"
  write(*,'(a)') "mpirun will spawn # processes.  The first process is called the master, and"
  write(*,'(a)') "the remaining #-1 processes are called workers."
  write(*,'(a)') "Via MPI, the master will send each worker a real number."
  write(*,'(a)') "Via OpenMP, each worker will spawn 4 threads to calculate the first 4 powers of"
  write(*,'(a)') "the number it was assigned.  Each worker then returns those powers"
  write(*,'(a)') "to the master via MPI, and the master prints the results."
  write(*,*)
endif

if(master) then
  !Check that cluster has at least two nodes
  call mpi_comm_size(MPI_COMM_WORLD,cluster_size,mpierr)
  nslave=cluster_size-1
  if(nslave .eq. 0) then
    write(*,*) "ERROR: no slaves found in cluster.  At least two nodes"
    write(*,*) "must be available to run this program."
    call mpi_finalize(mpierr) !tell the MPI daemon that you are leaving.
    error stop
  else
    write(*,*) "number of MPI slaves: ", nslave
  endif
endif

if(master) then
  allocate(x(nslave))
  allocate(x_powers(nslave,npowers))
  do i=1, nslave
    x(i) = i*1.11111d0
  enddo
endif

if(master) then
  do i=1, nslave
    call mpi_send(x(i), 1, MPI_DOUBLE_PRECISION, i, 1, MPI_COMM_WORLD, mpierr)
  enddo

  nrecv = 0
  do while(nrecv .lt. nslave)
    call mpi_probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, mpistatus, mpierr)  !blocking
    from_id = mpistatus(MPI_SOURCE)
    call mpi_recv(x_powers(from_id,:), npowers, MPI_DOUBLE_PRECISION, from_id, 2, MPI_COMM_WORLD, mpistatus, mpierr)
    nrecv = nrecv + 1
  enddo

  write(*,'(10a14)') "slave id", "x", "x**1", "x**2", "x**3", "x**4"
  do i=1,nslave
    write(*,'(i14,20f14.5)') i, x(i), x_powers(i,:)
  enddo
else
  call mpi_recv(myx, 1, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, mpistatus, mpierr)

  thread_check = -8088
  !$OMP PARALLEL DO &
  !$OMP DEFAULT(PRIVATE), &
  !$OMP SHARED(myx,myx_powers,thread_check)
  do i=1,npowers
    thread_check(omp_get_thread_num()+1) = omp_get_thread_num()
    myx_powers(i) = myx**i
  enddo
  !$OMP END PARALLEL DO
  write(*,'(a,i3,a,20i3)') "mpi slave: ", myrank, " has omp thread numbers: ", thread_check

  call mpi_send(myx_powers, npowers, MPI_DOUBLE_PRECISION, 0, 2, MPI_COMM_WORLD, mpierr)
endif

call mpi_finalize(mpierr)
write(*,*) "mpi process ", myrank, " made it!"

end program
