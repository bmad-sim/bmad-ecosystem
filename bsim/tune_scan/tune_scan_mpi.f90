program tune_scan_mpi

use ts_mod
use mpi

implicit none

type (ts_params_struct) ts
type (ts_com_struct) ts_com
type (ts_data_struct), allocatable, target :: ts_dat(:,:,:), dat_arr(:)
type (ts_data_struct), pointer :: t

integer ja, jb, jz
integer num_slaves, slave_rank, stat(MPI_STATUS_SIZE)
integer i, j0, rc, ierr, mpi_n_proc, leng, ix, num_a, num_b, num_z
integer n_run, n_per_job, n_job, n_track, dat_size, ix_job, ix_job_last, is
integer, allocatable :: slave_ix_job(:)

character(MPI_MAX_PROCESSOR_NAME) name

! Initialize MPI

call mpi_init(ierr)
if (ierr /= MPI_SUCCESS) then
  print *,'Error starting MPI program. Terminating.'
  call mpi_abort(MPI_COMM_WORLD, rc, ierr)
end if

! Get the number of processors this job is using:
call mpi_comm_size(MPI_COMM_WORLD, mpi_n_proc, ierr)

! Get the rank of the processor this thread is running on.
! Each processor has a unique rank.
call mpi_comm_rank(MPI_COMM_WORLD, ts_com%mpi_rank, ierr)

! Get the name of this processor (usually the hostname)
call mpi_get_processor_name(name, leng, ierr)
if (ierr /= MPI_SUCCESS) then
  print *,'Error getting processor name. Terminating.'
  call mpi_abort(MPI_COMM_WORLD, rc, ierr)
end if

num_slaves = mpi_n_proc - 1
if (num_slaves /= 0) ts_com%using_mpi = .true.

!---------------------------------------------
! Init. Only the master should create a digested lat file

call run_timer ('ABS', ts_com%time_start)

if (ts_com%mpi_rank == master_rank$) then
  call ts_init_params (ts, ts_com)
  call mpi_Bcast (0, 1, MPI_INTEGER, master_rank$, MPI_COMM_WORLD, ierr)
else
  call mpi_Bcast (ix, 1, MPI_INTEGER, master_rank$, MPI_COMM_WORLD, ierr)
  call ts_init_params (ts, ts_com)
endif

!---------------------------------------------
! If only one thread!

if (.not. ts_com%using_mpi) then
  print '(a, i0)', 'Note! Number of threads is one!'
  allocate (ts_dat(0:ts_com%n_a, 0:ts_com%n_b, 0:ts_com%n_z))
  do ja = 0, ts_com%n_a
  do jb = 0, ts_com%n_b
  do jz = 0, ts_com%n_z
    call ts_track_particle (ts, ts_com, [ja, jb, jz], ts_dat(ja,jb,jz))
  enddo
  enddo
  enddo
  call ts_write_results (ts, ts_com, ts_dat)
  stop
endif

!---------------------------------------------
! Common init

num_a = ts_com%n_a + 1
num_b = ts_com%n_b + 1
num_z = ts_com%n_z + 1

n_run = num_a * num_b * num_z
n_per_job = 1 + n_run / (2*num_slaves)
n_job = ceiling(0.9999999_rp * n_run / n_per_job)
allocate (dat_arr(n_per_job))

print '(a, 3(2x, i0))', 'n_run, n_job, n_per_job:', n_run, n_job, n_per_job
!---------------------------------------------
! Master

if (ts_com%mpi_rank == master_rank$) then
  print '(a, i0)', 'Number of processes (including Master): ', mpi_n_proc
  call ts_print_mpi_info (ts, ts_com, 'Master: Starting...', .true.)

  allocate (ts_dat(0:ts_com%n_a, 0:ts_com%n_b, 0:ts_com%n_z))
  allocate (slave_ix_job(num_slaves))
  slave_ix_job = 0    ! Job index

  ! Init slaves
  ix_job_last = 0
  do is = 1, num_slaves
    ix_job_last = ix_job_last + 1
    if (ix_job_last > n_job) exit
    call mpi_send(ix_job_last, 1, MPI_INTEGER, is, job_tag$, MPI_COMM_WORLD, ierr)
    slave_ix_job(is) = ix_job_last
    call ts_print_mpi_info (ts, ts_com, 'Master: Initial positions to slave: ' // int_str(is), .true.)
  enddo

  do
    ! Wait for data
    call ts_print_mpi_info (ts, ts_com, 'Master: Waiting for data from a Slave... ')
    slave_rank = MPI_ANY_SOURCE
    call mpi_recv (n_track, 1, MPI_INTEGER, slave_rank, have_data_tag$, MPI_COMM_WORLD, stat, ierr)
    slave_rank = stat(MPI_SOURCE)
    call ts_print_mpi_info (ts, ts_com, 'Master: Slave ' // int_str(slave_rank) // ' signals have data: ' // int_str(n_track))
    dat_size = n_track * storage_size(dat_arr(1)) / 8
    call mpi_recv(dat_arr(1:n_track), dat_size, MPI_BYTE, slave_rank, results_tag$, MPI_COMM_WORLD, stat, ierr)
    j0 = slave_ix_job(slave_rank)
    slave_ix_job(slave_rank) = 0
    do i = 1, n_track
      ix = (j0-1)*n_per_job + i - 1
      ja = ix / (num_b*num_z)
      ix = ix - ja * num_b * num_z
      jb = ix / (ts_com%n_z+1)
      jz = ix - jb * num_z
      ts_dat(ja,jb,jz) = dat_arr(i)
    enddo
    call ts_print_mpi_info (ts, ts_com, 'Master: Gathered data from Slave: ' // int_str(slave_rank) // '  For job: ' // int_str(j0), .true.)

    ! Send job info to slave
    ix_job_last = ix_job_last + 1
    if (ix_job_last > n_job) then
      call ts_print_mpi_info (ts, ts_com, 'Master: Signal slave is done: ' // int_str(slave_rank))
      call mpi_send(0, 1, MPI_INTEGER, slave_rank, job_tag$, MPI_COMM_WORLD, ierr)
    else
      call ts_print_mpi_info (ts, ts_com, 'Master: Commanding Slave: ' // int_str(slave_rank) // &
                                              '  For job: ' // int_str(ix_job_last))
      call mpi_send(ix_job_last, 1, MPI_INTEGER, slave_rank, job_tag$, MPI_COMM_WORLD, ierr)
      slave_ix_job(slave_rank) = ix_job_last
    endif

    if (all(slave_ix_job == 0)) exit
  enddo

  call ts_write_results (ts, ts_com, ts_dat)
  call ts_print_mpi_info (ts, ts_com, 'Master: All done!', .true.)
  call mpi_finalize(ierr)

!---------------------------------------------
! Slave

else
  call ts_print_mpi_info (ts, ts_com, 'Slave Starting...')

  do
    call ts_print_mpi_info (ts, ts_com, 'Slave: Waiting for master.')
    call mpi_recv (ix_job, 1, MPI_INTEGER, master_rank$, job_tag$, MPI_COMM_WORLD, stat, ierr)
    call ts_print_mpi_info (ts, ts_com, 'Slave: Starting job: ' // int_str(ix_job))
    if (ix_job == 0) exit

    do i = 1, n_per_job
      ix = (ix_job-1) * n_per_job + i - 1
      ja = ix / (num_b * num_z)
      if (ja > ts_com%n_a) exit
      ix = ix - ja * num_b * num_z
      jb = ix / num_z
      jz = ix - jb * num_z
      call ts_track_particle (ts, ts_com, ja, jb, jz, dat_arr(i))
      call ts_print_mpi_info (ts, ts_com, 'Slave: Tracked particle at index: ' // &
                                    int_str(ja) // ', ' // int_str(jb) // ', ' // int_str(jz))
    enddo
    n_track = i - 1

    call ts_print_mpi_info (ts, ts_com, 'Slave: Have data for job: ' // int_str(ix_job))
    call mpi_send (n_track, 1, MPI_INTEGER, master_rank$, have_data_tag$, MPI_COMM_WORLD, ierr)
    call ts_print_mpi_info (ts, ts_com, 'Slave: Sending data for job: ' // int_str(ix_job))
    dat_size = n_track * storage_size(dat_arr(1)) / 8
    call mpi_send (dat_arr(1:n_track), dat_size, MPI_BYTE, master_rank$, results_tag$, MPI_COMM_WORLD, ierr)
  enddo

  call ts_print_mpi_info (ts, ts_com, 'Slave: All done!', .true.)
  call mpi_finalize(ierr)
endif

end program
