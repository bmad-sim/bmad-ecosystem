program long_term_tracking

use lt_tracking_mod
use mpi
use directory_mod

implicit none

type (ltt_params_struct) lttp
type (ltt_com_struct), target :: ltt_com
type (beam_init_struct) beam_init
type (ltt_sum_data_struct), allocatable, target :: sum_data_arr(:), sd_arr(:)
type (ltt_sum_data_struct) sum_data
type (ltt_sum_data_struct), pointer :: sd
type (ele_struct), pointer :: ele_start
type (lat_struct), pointer :: lat
type (beam_struct), target :: beam, beam2
type (bunch_struct), pointer :: bunch
type (bunch_struct) :: bunch0
type (coord_struct) orb

real(rp) time_now, mpi_particles_per_run

integer num_slaves, slave_rank, stat(MPI_STATUS_SIZE)
integer i, n, nn, ix, ierr, rc, leng, sd_arr_dat_size, storage_size, dat_size
integer ix0_p, ix1_p, mpi_n_proc, n_pass, ix_path, n_out, iu, ios, i_turn
integer, allocatable :: turn(:)

logical am_i_done, err_flag, ok
logical, allocatable :: slave_is_done(:)

character(100) line, path, basename
character(40) file
character(40), allocatable :: file_list(:)
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
call mpi_comm_rank(MPI_COMM_WORLD, ltt_com%mpi_rank, ierr)

! Get the name of this processor (usually the hostname)
call mpi_get_processor_name(name, leng, ierr)
if (ierr /= MPI_SUCCESS) then
  print *,'Error getting processor name. Terminating.'
  call mpi_abort(MPI_COMM_WORLD, rc, ierr)
end if

num_slaves = mpi_n_proc - 1
if (num_slaves /= 0) ltt_com%using_mpi = .true.

! If not doing BUNCH tracking then slaves have nothing to do.

call ltt_init_params(lttp, ltt_com, beam_init)

if (lttp%simulation_mode /= 'BUNCH' .and. ltt_com%mpi_rank /= master_rank$) then
  call mpi_finalize(ierr)
  stop
endif

! Only the master should create a map file if a file is to be created.

if (ltt_com%mpi_rank == master_rank$) then
  call ltt_init_tracking (lttp, ltt_com)
  call mpi_Bcast (0, 1, MPI_INTEGER, master_rank$, MPI_COMM_WORLD, ierr)
  call ltt_print_inital_info (lttp, ltt_com)

else
  call mpi_Bcast (ix, 1, MPI_INTEGER, master_rank$, MPI_COMM_WORLD, ierr)
  call ltt_init_tracking (lttp, ltt_com)
endif

! Calculation start.

call run_timer ('ABS', ltt_com%time_start)

select case (lttp%simulation_mode)
case ('CHECK');  call ltt_run_check_mode(lttp, ltt_com, beam_init)  ! A single turn tracking check
case ('SINGLE'); call ltt_run_single_mode(lttp, ltt_com, beam_init) ! Single particle tracking
case ('STAT');   call ltt_run_stat_mode(lttp, ltt_com)              ! Lattice statistics (radiation integrals, etc.).
case ('BUNCH');  ! Handled below. Only the BUNCH simulation mode uses mpi.
case default;    print *, 'BAD SIMULATION_MODE: ' // lttp%simulation_mode
end select

! Not using mpi if there is only one thread.

if (.not. ltt_com%using_mpi) then
  print '(a, i0)', 'Number of threads is one!'
  call ltt_run_bunch_mode(lttp, ltt_com, beam_init)  ! Beam tracking
  stop
endif

!-----------------------------------------
! MPI BUNCH simulation

ix_path = splitfilename (lttp%particle_output_file, path, basename)
ltt_com%mpi_data_dir = trim(path) // 'mpi_temp_dir/'

!

if (ltt_com%mpi_rank == master_rank$) then
  call system_command ('mkdir ' // trim(ltt_com%mpi_data_dir))

  allocate (slave_is_done(num_slaves))
  slave_is_done = .false.

  lat => ltt_com%tracking_lat
  call ltt_pointer_to_map_ends(lttp, lat, ele_start)
  call init_beam_distribution (ele_start, lat%param, beam_init, beam, err_flag, modes = ltt_com%modes)
  bunch => beam%bunch(1)
  mpi_particles_per_run = real(size(bunch%particle), rp) / real((lttp%mpi_runs_per_subprocess * (mpi_n_proc - 1)), rp)

  print '(a, i0)', 'Number of processes (including Master): ', mpi_n_proc
  print '(a, i0, 2x, i0)', 'Nominal number of particles per pass: ', nint(mpi_particles_per_run)
  call ltt_print_mpi_info (lttp, ltt_com, 'Master: Starting...', .true.)

  call ltt_allocate_sum_array(lttp, sd_arr)
  call ltt_allocate_sum_array(lttp, sum_data_arr)
  sd_arr_dat_size = size(sd_arr) * storage_size(sd_arr(1)) / 8

  n_pass = 0
  ix0_p = 0
  do i = 1, mpi_n_proc-1
    n_pass = n_pass + 1
    ix1_p = min(nint(n_pass * mpi_particles_per_run), size(bunch%particle))
    n = ix1_p-ix0_p
    dat_size = n * storage_size(bunch%particle(1)) / 8
    call ltt_print_mpi_info (lttp, ltt_com, 'Master: Init position data size to slave: ' // int_str(i))
    call mpi_send (n, 1, MPI_INTEGER, i, num_tag$, MPI_COMM_WORLD, ierr)
    call mpi_send (ix0_p, 1, MPI_INTEGER, i, num_tag$, MPI_COMM_WORLD, ierr)
    call ltt_print_mpi_info (lttp, ltt_com, 'Master: Initial positions to slave: ' // int_str(i) // &
                                                '  For particles: [' // int_str(ix0_p+1) // ':' // int_str(ix1_p) // ']', .true.)
    call mpi_send (bunch%particle(ix0_p+1:ix1_p), dat_size, MPI_BYTE, i, particle_tag$, MPI_COMM_WORLD, ierr)
    ix0_p = ix1_p
  enddo

  !

  do
    ! Get data from a slave
    call ltt_print_mpi_info (lttp, ltt_com, 'Master: Waiting for data from a Slave... ' // int_str(sd_arr_dat_size))
    call mpi_recv (sd_arr, sd_arr_dat_size, MPI_BYTE, MPI_ANY_SOURCE, results_tag$, MPI_COMM_WORLD, stat, ierr)

    slave_rank = stat(MPI_SOURCE)
    call ltt_print_mpi_info (lttp, ltt_com, 'Master: Gathered data from Slave: ' // int_str(slave_rank))

    ! Add to data
    do ix = lbound(sum_data_arr, 1), ubound(sum_data_arr, 1)
      sd => sum_data_arr(ix)
      sd%i_turn   = sd_arr(ix)%i_turn
      sd%n_live   = sd%n_live + sd_arr(ix)%n_live
      sd%n_count  = sd%n_count + sd_arr(ix)%n_count
      sd%orb_sum  = sd%orb_sum + sd_arr(ix)%orb_sum
      sd%orb2_sum = sd%orb2_sum + sd_arr(ix)%orb2_sum
      sd%spin_sum = sd%spin_sum + sd_arr(ix)%spin_sum
      sd%p0c_sum  = sd%p0c_sum + sd_arr(ix)%p0c_sum
      sd%time_sum = sd%time_sum + sd_arr(ix)%time_sum
      if (sd_arr(ix)%status == valid$) sd%status = valid$
    enddo

    ! Tell slave if more tracking needed

    if (ix0_p == size(bunch%particle)) slave_is_done(slave_rank) = .true.
    call mpi_send (slave_is_done(slave_rank), 1, MPI_LOGICAL, slave_rank, is_done_tag$, MPI_COMM_WORLD, ierr)
    if (all(slave_is_done)) exit       ! All done?
    if (ix0_p == size(bunch%particle)) cycle
    
    ! Give slave particle positions

    n_pass = n_pass + 1
    ix1_p = min(nint(n_pass*mpi_particles_per_run), size(bunch%particle))

    n = ix1_p-ix0_p
    dat_size = n * storage_size(bunch%particle(1)) / 8
    call ltt_print_mpi_info (lttp, ltt_com, 'Master: Position data size to slave: ' // int_str(slave_rank))
    call mpi_send (n, 1, MPI_INTEGER, slave_rank, num_tag$, MPI_COMM_WORLD, ierr)
    call mpi_send (ix0_p, 1, MPI_INTEGER, slave_rank, num_tag$, MPI_COMM_WORLD, ierr)
    call ltt_print_mpi_info (lttp, ltt_com, 'Master: Initial positions to slave: ' // int_str(slave_rank) // &
                                                 '  For particles: [' // int_str(ix0_p+1) // ':' // int_str(ix1_p) // ']', .true.)
    call mpi_send (bunch%particle(ix0_p+1:ix1_p), dat_size, MPI_BYTE, slave_rank, particle_tag$, MPI_COMM_WORLD, ierr)

    ix0_p = ix1_p
  enddo

  ! Write results

  call ltt_write_bunch_averages (lttp, sum_data_arr)
  call ltt_write_sigma_matrix (lttp, sum_data_arr)
  call ltt_print_mpi_info (lttp, ltt_com, 'Master: All done!', .true.)
  call mpi_finalize(ierr)

  ! Gather particle data files

  ok = dir_list (ltt_com%mpi_data_dir, file_list)
  iu = lunget()

  if (lttp%particle_output_file /= '') then
    n_out = 0   ! Number of output files to generage
    allocate (turn(size(file_list)))
    do i = 1, size(file_list)
      file = file_list(i)
      if (file(1:1) == 'b') cycle
      ix = index(file, '_')
      read (file(1:ix-1), *) n
      if (any(n == turn(1:n_out))) cycle
      n_out = n_out + 1
      turn(n_out) = n
    enddo

    do nn = 1, n_out
      do i = 1, size(file_list)
        file = file_list(i)
        if (file(1:1) == 'b') cycle
        ix = index(file, '_')
        read (file(1:ix-1), *) n
        if (n /= turn(nn)) cycle
        open (iu, file = trim(ltt_com%mpi_data_dir) // file, status = 'OLD')
        do
          read (iu, *, iostat = ios) ix, i_turn, orb%vec, orb%spin, orb%state
          if (ios /= 0) exit
          bunch%particle(ix) = orb
        enddo
        close (iu)
      enddo
      print *, 'Here:', i_turn, turn(nn)
      call ltt_write_particle_data (lttp, ltt_com, i_turn, bunch, .false.)
    enddo
  endif

  if (lttp%bunch_binary_output_file /= '') then
    do i = 1, size(file_list)
      file = file_list(i)
      if (file(1:1) /= 'b') cycle
      call hdf5_read_beam(trim(ltt_com%mpi_data_dir) // file, beam2, err_flag)
      read (file(2:), *) ix
      n = size(beam2%bunch(1)%particle)
      bunch%particle(ix+1:ix+n) = beam2%bunch(1)%particle
    enddo
    call write_beam_file (lttp%bunch_binary_output_file, beam)
  endif

  call system_command ('rm -rf ' // trim(ltt_com%mpi_data_dir))

  ! And end

  call run_timer ('ABS', time_now)
  print '(a, f8.2)', 'Tracking time (min):', (time_now - ltt_com%time_start) / 60

!-----------------------------------------
else  ! Is a slave

  do
    ! Init positions
    call ltt_print_mpi_info (lttp, ltt_com, 'Slave: Waiting for position size info.')
    call mpi_recv (n, 1, MPI_INTEGER, master_rank$, num_tag$, MPI_COMM_WORLD, stat, ierr)
    call mpi_recv (ltt_com%mpi_ix0_particle, 1, MPI_INTEGER, master_rank$, num_tag$, MPI_COMM_WORLD, stat, ierr)
    if (allocated(bunch0%particle)) then
      if (size(bunch0%particle) /= n) deallocate(bunch0%particle)
    endif
    if (.not. allocated(bunch0%particle)) allocate(bunch0%particle(n))
    call ltt_print_mpi_info (lttp, ltt_com, 'Slave: Waiting for position info for ' // int_str(n) // ' particles.')
    dat_size = n * storage_size(bunch0%particle(1)) / 8
    call mpi_recv (bunch0%particle, dat_size, MPI_BYTE, MPI_ANY_SOURCE, particle_tag$, MPI_COMM_WORLD, stat, ierr)

    ! Run
    call ltt_run_bunch_mode(lttp, ltt_com, beam_init, sd_arr, bunch0)  ! Beam tracking
    sd_arr_dat_size = size(sd_arr) * storage_size(sd_arr(1)) / 8
    call ltt_print_mpi_info (lttp, ltt_com, 'Slave: Sending Data... ' // int_str(sd_arr_dat_size))
    call mpi_send (sd_arr, sd_arr_dat_size, MPI_BYTE, master_rank$, results_tag$, MPI_COMM_WORLD, ierr)

    ! Query Master if more tracking needed
    call ltt_print_mpi_info (lttp, ltt_com, 'Slave: Query am-i-done to master...')
    call mpi_recv (am_i_done, 1, MPI_LOGICAL, master_rank$, is_done_tag$, MPI_COMM_WORLD, stat, ierr)
    if (am_i_done) exit
  enddo

  call ltt_print_mpi_info (lttp, ltt_com, 'Slave: All done!')
  call mpi_finalize(ierr)

endif

end program

