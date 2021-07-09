program lux_mpi

use lux_module
use mpi

implicit none

type (lux_param_struct) lux_param
type (lux_common_struct), target :: lux_com
type (lux_output_data_struct) lux_data
type (surface_grid_struct) :: slave_grid
type (surface_grid_struct), pointer :: detec_grid

integer master_rank, ierr, rc, leng, i, stat(MPI_STATUS_SIZE)
integer data_size, num_slaves
integer results_tag, is_done_tag, slave_rank
integer(8) num_photons_left

logical am_i_done
logical, allocatable :: slave_is_done(:)

character(MPI_MAX_PROCESSOR_NAME) name

! Initialize MPI 

call mpi_init(ierr)
if (ierr /= MPI_SUCCESS) then
  print *,'Error starting MPI program. Terminating.'
  call mpi_abort(MPI_COMM_WORLD, rc, ierr)
end if

! Get the number of processors this job is using:
call mpi_comm_size(MPI_COMM_WORLD, lux_com%mpi_n_proc, ierr)

! Get the rank of the processor this thread is running on.  (Each
! processor has a unique rank.)
call mpi_comm_rank(MPI_COMM_WORLD, lux_com%mpi_rank, ierr)

! Get the name of this processor (usually the hostname)
call mpi_get_processor_name(name, leng, ierr)
if (ierr /= MPI_SUCCESS) then
  print *,'Error getting processor name. Terminating.'
  call mpi_abort(MPI_COMM_WORLD, rc, ierr)
end if

master_rank = 0
num_slaves = lux_com%mpi_n_proc - 1

if (num_slaves < 1) then
  print *, 'ONLY ONE PROCESS EXISTS! WILL RUN SERIAL...'
  call lux_run_serial()
  call mpi_finalize(ierr)
  stop
endif

!-------------------------------
! Init Lux

lux_com%using_mpi = .true.
call lux_init (lux_param, lux_com)

detec_grid => lux_com%detec_ele%photon%surface%grid
allocate (slave_grid%pt(size(detec_grid%pt, 1), size(detec_grid%pt, 1)))

! storage_size returns size in bytes per point
data_size = size(detec_grid%pt, 1) * size(detec_grid%pt, 2) * storage_size(detec_grid%pt) / 8

results_tag = 1000
is_done_tag = 1001

!-------------------------------
! Master collects the work of the slaves

if (lux_com%mpi_rank == master_rank) then
  print '(a, i0)', 'Number of processes (including Master): ', lux_com%mpi_n_proc
  call print_mpi_info ('Master: Starting...')

  allocate (slave_is_done(num_slaves))
  slave_is_done = .false.

  ! Init the output arrays.
  call lux_init_data (lux_param, lux_com, lux_data)

  ! Slaves automatically start one round of tracking
  num_photons_left = lux_param%stop_num_photons - num_slaves * lux_com%n_photon_stop1

  do
    ! Get data from a slave
    call print_mpi_info ('Master: Waiting for a Slave...')
    call mpi_recv (slave_grid%pt, data_size, MPI_REAL8, MPI_ANY_SOURCE, results_tag, MPI_COMM_WORLD, stat, ierr)
    call lux_add_in_mpi_slave_data (slave_grid%pt, lux_param, lux_com, lux_data)
    slave_rank = stat(MPI_SOURCE)
    call print_mpi_info ('Master: Gathered data from Slave: ', slave_rank)

    ! Tell slave if more tracking needed
    call print_mpi_info ('Master: Commanding this Slave. Photons left:', num8 = num_photons_left)
    if (num_photons_left < 1) slave_is_done(slave_rank) = .true.
    call mpi_send (slave_is_done(slave_rank), 1, MPI_LOGICAL, slave_rank, is_done_tag, MPI_COMM_WORLD, ierr)
    if (.not. slave_is_done(slave_rank)) num_photons_left = num_photons_left - lux_com%n_photon_stop1

    ! All done?
    if (all(slave_is_done)) exit
  enddo

  call print_mpi_info ('Master: All done!')

  ! write results and quit

  call lux_write_data (lux_param, lux_com, lux_data)

  call mpi_finalize(ierr)
  stop

endif

!-------------------------------
! A slave process tracks photons

call print_mpi_info ('Slave: Starting...')

do
  ! Init the output arrays
  call print_mpi_info ('Slave: Tracking Photons...')
  call lux_init_data (lux_param, lux_com, lux_data)
  call lux_track_photons (lux_param, lux_com, lux_data)

  ! Send results to the Master
  call print_mpi_info ('Slave: Sending Data...')
  call mpi_send (detec_grid%pt, data_size, MPI_BYTE, master_rank, results_tag, MPI_COMM_WORLD, ierr)

  ! Query Master if more tracking needed
  call print_mpi_info ('Slave: Query to master...')
  call mpi_recv (am_i_done, 1, MPI_LOGICAL, master_rank, is_done_tag, MPI_COMM_WORLD, stat, ierr)
  if (am_i_done) exit

enddo

if (lux_param%photon1_out_file /= '') close (lux_com%iu_photon1_out)

call print_mpi_info ('Slave: All done!')

! And end

call mpi_finalize(ierr)


!---------------------------------------------------------------------
contains

subroutine print_mpi_info (line, num, num8)

character(*) line
character(20) dtime
integer, optional :: num
integer(8), optional :: num8

!

call date_and_time_stamp (dtime)
if (present(num)) then
  print '(a, 2x, i0, 2a, 1x, i0)', dtime, lux_com%mpi_rank, ': ', trim(line), num
elseif (present(num8)) then
  print '(a, 2x, i0, 2a, 1x, i0)', dtime, lux_com%mpi_rank, ': ', trim(line), num8
else
  print '(a, 2x, i0, 2a)', dtime, lux_com%mpi_rank, ': ', trim(line)
endif

end subroutine print_mpi_info

end program
