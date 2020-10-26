program long_term_tracking

use lt_tracking_mod
use mpi

implicit none

type (ltt_params_struct) lttp
type (ltt_com_struct) ltt_com
type (beam_init_struct) beam_init
type (ltt_sum_data_struct), allocatable, target :: sum_data_arr(:), sd_arr(:)
type (ltt_sum_data_struct) sum_data
type (ltt_sum_data_struct), pointer :: sd

real(rp) del_time

integer num_slaves, slave_rank, stat(MPI_STATUS_SIZE)
integer n, ix, ierr, rc, leng, data_size, num_particles_left, storage_size

logical am_i_done
logical, allocatable :: slave_is_done(:)

character(80) line
character(MPI_MAX_PROCESSOR_NAME) name

! Initialize MPI

call mpi_init(ierr)
if (ierr /= MPI_SUCCESS) then
  print *,'Error starting MPI program. Terminating.'
  call mpi_abort(MPI_COMM_WORLD, rc, ierr)
end if

! Get the number of processors this job is using:
call mpi_comm_size(MPI_COMM_WORLD, lttp%mpi_n_proc, ierr)

! Get the rank of the processor this thread is running on.
! Each processor has a unique rank.
call mpi_comm_rank(MPI_COMM_WORLD, lttp%mpi_rank, ierr)

! Get the name of this processor (usually the hostname)
call mpi_get_processor_name(name, leng, ierr)
if (ierr /= MPI_SUCCESS) then
  print *,'Error getting processor name. Terminating.'
  call mpi_abort(MPI_COMM_WORLD, rc, ierr)
end if

num_slaves = lttp%mpi_n_proc - 1
if (num_slaves /= 0) lttp%using_mpi = .true.

! If not doing BUNCH tracking then slaves have nothing to do.

call ltt_init_params(lttp, ltt_com, beam_init)

if (lttp%simulation_mode /= 'BUNCH' .and. lttp%mpi_rank /= master_rank$) then
  call mpi_finalize(ierr)
  stop
endif

! Only the master should create a map file if a file is to be created.

if (lttp%mpi_rank == master_rank$) then
  call ltt_init_tracking (lttp, ltt_com)
  call mpi_Bcast (0, 1, MPI_INTEGER, master_rank$, MPI_COMM_WORLD, ierr)
  call ltt_print_inital_info (lttp, ltt_com)

else
  call mpi_Bcast (ix, 1, MPI_INTEGER, master_rank$, MPI_COMM_WORLD, ierr)
  call ltt_init_tracking (lttp, ltt_com)
endif

! Calculation start.

call run_timer ('START')

select case (lttp%simulation_mode)
case ('CHECK');  call ltt_run_check_mode(lttp, ltt_com, beam_init)  ! A single turn tracking check
case ('SINGLE'); call ltt_run_single_mode(lttp, ltt_com, beam_init) ! Single particle tracking
case ('STAT');   call ltt_run_stat_mode(lttp, ltt_com)              ! Lattice statistics (radiation integrals, etc.).
case default;    print *, 'BAD SIMULATION_MODE: ' // lttp%simulation_mode

!-------------------------
! Only the BUNCH simulation mode uses mpi

case ('BUNCH')

  if (lttp%using_mpi .and. lttp%output_every_n_turns < 1) then
    if (lttp%mpi_rank == 0) print *, 'OUTPUT_EVERY_N_TURNS MUST BE POSITIVE WHEN USING MPI!'
    stop
  endif

  if (lttp%using_mpi .and. lttp%averages_output_file == '') then
    if (lttp%mpi_rank == 0) print *, 'AVERAGES_OUPUT_FILE MUST BE SET WHEN USING MPI!'
    stop
  endif

  n = lttp%n_turns / lttp%output_every_n_turns
  allocate (sum_data_arr(0:n), sd_arr(0:n))
  lttp%mpi_n_particles_per_run = nint(real(beam_init%n_particle, rp) / (lttp%mpi_num_runs * (lttp%mpi_n_proc - 1)))
  data_size = size(sum_data_arr) * storage_size(sum_data_arr(0)) / 8

  if (.not. lttp%using_mpi) then
    call ltt_run_bunch_mode(lttp, ltt_com, beam_init)  ! Beam tracking

  !-----------------------------------------
  elseif (lttp%mpi_rank == master_rank$) then

    print '(a, i0)', 'Number of processes (including Master): ', lttp%mpi_n_proc
    print '(a, i0, 2x, i0)', 'Number of particles per run: ', lttp%mpi_n_particles_per_run
    call print_mpi_info (lttp, 'Master: Starting...')

    allocate (slave_is_done(num_slaves))
    slave_is_done = .false.

    ! Slaves automatically start one round of tracking
    num_particles_left = beam_init%n_particle - num_slaves * lttp%mpi_n_particles_per_run

    do
      ! Get data from a slave
      call print_mpi_info (lttp, 'Master: Waiting for a Slave...')
      call mpi_recv (sd_arr, data_size, MPI_BYTE, MPI_ANY_SOURCE, results_tag$, MPI_COMM_WORLD, stat, ierr)

      ! Add to data
      do ix = lbound(sum_data_arr, 1), ubound(sum_data_arr, 1)
        sd => sum_data_arr(ix)
        sd%i_turn   = ix * lttp%output_every_n_turns
        sd%n_live   = sd%n_live + sd_arr(ix)%n_live
        sd%orb_sum  = sd%orb_sum + sd_arr(ix)%orb_sum
        sd%orb2_sum = sd%orb2_sum + sd_arr(ix)%orb2_sum
        sd%spin_sum = sd%spin_sum + sd_arr(ix)%spin_sum
      enddo

      slave_rank = stat(MPI_SOURCE)
      write (line, '(a, i0)') 'Master: Gathered data from Slave: ', slave_rank
      call print_mpi_info (lttp, line)

      ! Tell slave if more tracking needed
      write (line, '(a, i0, a, i0)') 'Master: Commanding slave: ', slave_rank, '. Particles left:', num_particles_left
      call print_mpi_info (lttp, line, .true.)
      if (num_particles_left < 1) slave_is_done(slave_rank) = .true.
      call mpi_send (slave_is_done(slave_rank), 1, MPI_LOGICAL, slave_rank, is_done_tag$, MPI_COMM_WORLD, ierr)
      if (.not. slave_is_done(slave_rank)) num_particles_left = num_particles_left - lttp%mpi_n_particles_per_run

      ! All done?
      if (all(slave_is_done)) exit
    enddo

    ! write results and quit

    call ltt_write_bunch_averages (lttp, -1, sum_data_arr)
    call ltt_write_sigma_matrix (lttp, -1, sum_data_arr)
    call print_mpi_info (lttp, 'Master: All done!', .true.)
    call mpi_finalize(ierr)

    call run_timer ('READ', del_time)
    print '(a, f8.2)', 'Tracking time (min):', del_time/60

  !-----------------------------------------
  else  ! Is a slave

    do
      ! Init the output arrays
      call print_mpi_info (lttp, 'Slave: Tracking Particles...')

      ! Run
      call ltt_run_bunch_mode(lttp, ltt_com, beam_init, sd_arr)  ! Beam tracking
      data_size = size(sd_arr) * storage_size(sd_arr(1)) / 8
      call print_mpi_info (lttp, 'Slave: Sending Data...')
      call mpi_send (sd_arr, data_size, MPI_BYTE, master_rank$, results_tag$, MPI_COMM_WORLD, ierr)

      ! Query Master if more tracking needed
      call print_mpi_info (lttp, 'Slave: Query to master...')
      call mpi_recv (am_i_done, 1, MPI_LOGICAL, master_rank$, is_done_tag$, MPI_COMM_WORLD, stat, ierr)
      if (am_i_done) exit
    enddo

    call print_mpi_info (lttp, 'Slave: All done!')
    call mpi_finalize(ierr)

  endif

end select

!

end program

