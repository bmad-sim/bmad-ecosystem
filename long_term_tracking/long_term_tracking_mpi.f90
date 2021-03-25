program long_term_tracking

use lt_tracking_mod
use mpi
use directory_mod

implicit none

type (ltt_params_struct) lttp
type (ltt_com_struct), target :: ltt_com
type (beam_init_struct) beam_init
type (ltt_beam_data_struct), target :: beam_data, beam_data_sum
type (ltt_bunch_data_struct), pointer :: bd
type (ele_struct), pointer :: ele_start
type (lat_struct), pointer :: lat
type (beam_struct), target :: beam, beam2
type (bunch_struct), pointer :: bunch
type (coord_struct) orb

real(rp) time_now, mpi_particles_per_run

integer num_slaves, slave_rank, stat(MPI_STATUS_SIZE)
integer i, n, nn, nb, ib, ix, ierr, rc, leng, bd_size, storage_size, dat_size
integer ix0_p, ix1_p, mpi_n_proc, n_pass, ix_path, n_out, iu, ios, i_turn, iarr(2)
integer, allocatable :: turn(:)

logical am_i_done, err_flag, ok
logical, allocatable :: slave_is_done(:)

character(100) line, path, basename
character(40) file, str
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

! If not doing BEAM tracking then slaves have nothing to do.

call ltt_init_params(lttp, ltt_com, beam_init)

if (lttp%simulation_mode /= 'BEAM' .and. ltt_com%mpi_rank /= master_rank$) then
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
case ('BEAM');   ! Handled below. Only the BEAM simulation mode uses mpi.
case default;    print *, 'BAD SIMULATION_MODE: ' // lttp%simulation_mode
end select

! Not using mpi if there is only one thread.

if (.not. ltt_com%using_mpi) then
  print '(a, i0)', 'Number of threads is one!'
  call ltt_run_beam_mode(lttp, ltt_com, beam_init)  ! Beam tracking
  stop
endif

!-----------------------------------------
! MPI BEAM simulation

ix_path = splitfilename (lttp%particle_output_file, path, basename)
ltt_com%mpi_data_dir = trim(path) // 'mpi_temp_dir/'

! Master:

if (ltt_com%mpi_rank == master_rank$) then
  call system_command ('mkdir ' // trim(ltt_com%mpi_data_dir))

  allocate (slave_is_done(num_slaves))
  slave_is_done = .false.

  lat => ltt_com%tracking_lat
  call ltt_pointer_to_map_ends(lttp, lat, ele_start)
  call init_beam_distribution (ele_start, lat%param, beam_init, beam, err_flag, modes = ltt_com%modes)
  mpi_particles_per_run = real(size(beam%bunch(1)%particle), rp) / real((lttp%mpi_runs_per_subprocess * (mpi_n_proc - 1)), rp)

  print '(a, i0)', 'Number of processes (including Master): ', mpi_n_proc
  print '(a, i0, 2x, i0)', 'Nominal number of particles per pass: ', nint(mpi_particles_per_run)
  call ltt_print_mpi_info (lttp, ltt_com, 'Master: Starting...', .true.)

  call ltt_allocate_beam_data(lttp, beam_data_sum, size(beam%bunch))
  call ltt_allocate_beam_data(lttp, beam_data, size(beam%bunch))
  bd_size = size(beam_data%turn(1)%bunch) * storage_size(beam_data%turn(1)%bunch(1)) / 8

  ! Init positions to slaves

  n_pass = 0
  ix0_p = 0
  do i = 1, mpi_n_proc-1
    n_pass = n_pass + 1
    ix1_p = min(nint(n_pass * mpi_particles_per_run), size(beam%bunch(1)%particle))
    n = ix1_p-ix0_p
    dat_size = n * storage_size(beam%bunch(1)%particle(1)) / 8
    call ltt_print_mpi_info (lttp, ltt_com, 'Master: Init data size to slave: ' // int_str(i))
    call mpi_send (size(beam%bunch(1)%particle), 1, MPI_INTEGER, i, num_tag$, MPI_COMM_WORLD, ierr)
    call mpi_send (n, 1, MPI_INTEGER, i, num_tag$, MPI_COMM_WORLD, ierr)       ! Number of particles to track.
    call mpi_send (ix0_p, 1, MPI_INTEGER, i, num_tag$, MPI_COMM_WORLD, ierr)   ! Index of 1st particle - 1.
    call ltt_print_mpi_info (lttp, ltt_com, 'Master: Initial positions to slave: ' // int_str(i) // &
                                                '  For particles: [' // int_str(ix0_p+1) // ':' // int_str(ix1_p) // ']', .true.)
    do nb = 1, size(beam%bunch)
      call ltt_print_mpi_info (lttp, ltt_com, 'Master: Initial send bunch: ' // int_str(nb))
      call mpi_send (beam%bunch(nb)%particle(ix0_p+1:ix1_p), dat_size, MPI_BYTE, i, particle_tag$, MPI_COMM_WORLD, ierr)
    enddo

    ix0_p = ix1_p
  enddo

  ! Get data and initiate more tracking loop

  do
    ! Get data from slave. Since arrays of beam_data are allocatable, need to communicate in pieces.
    call ltt_print_mpi_info (lttp, ltt_com, 'Master: Waiting for data from a Slave... ' // int_str(bd_size))
    slave_rank = MPI_ANY_SOURCE
    do ix = lbound(beam_data%turn, 1), ubound(beam_data%turn, 1)
      call mpi_recv (iarr, 2, MPI_INTEGER, slave_rank, results_tag$, MPI_COMM_WORLD, stat, ierr)
      beam_data%turn(ix)%status = iarr(1);  beam_data%turn(ix)%i_turn = iarr(2)
      call mpi_recv (beam_data%turn(ix)%bunch, bd_size, MPI_BYTE, slave_rank, results_tag$, MPI_COMM_WORLD, stat, ierr)
      slave_rank = stat(MPI_SOURCE)
    enddo

    ! Merge with existing data
    do ix = lbound(beam_data%turn, 1), ubound(beam_data%turn, 1)
      if (beam_data%turn(ix)%status /= valid$) cycle
      do ib = 1, size(beam%bunch)
        bd => beam_data_sum%turn(ix)%bunch(ib)
        beam_data_sum%turn(ix)%i_turn   = beam_data%turn(ix)%i_turn
        bd%n_live   = bd%n_live + beam_data%turn(ix)%bunch(ib)%n_live
        bd%n_count  = bd%n_count + beam_data%turn(ix)%bunch(ib)%n_count
        bd%orb_sum  = bd%orb_sum + beam_data%turn(ix)%bunch(ib)%orb_sum
        bd%orb2_sum = bd%orb2_sum + beam_data%turn(ix)%bunch(ib)%orb2_sum
        bd%spin_sum = bd%spin_sum + beam_data%turn(ix)%bunch(ib)%spin_sum
        bd%p0c_sum  = bd%p0c_sum + beam_data%turn(ix)%bunch(ib)%p0c_sum
        bd%time_sum = bd%time_sum + beam_data%turn(ix)%bunch(ib)%time_sum
        bd%species  = beam%bunch(1)%particle(1)%species
      enddo
      beam_data_sum%turn(ix)%i_turn = beam_data%turn(ix)%i_turn
      beam_data_sum%turn(ix)%status = valid$
    enddo

    call ltt_print_mpi_info (lttp, ltt_com, 'Master: Gathered data from Slave: ' // int_str(slave_rank))

    ! Tell slave if more tracking needed
    if (ix0_p == size(beam%bunch(1)%particle)) slave_is_done(slave_rank) = .true.
    call mpi_send (slave_is_done(slave_rank), 1, MPI_LOGICAL, slave_rank, is_done_tag$, MPI_COMM_WORLD, ierr)
    if (all(slave_is_done)) exit       ! All done?
    if (ix0_p == size(beam%bunch(1)%particle)) cycle
    
    ! Give slave particle positions
    n_pass = n_pass + 1
    ix1_p = min(nint(n_pass*mpi_particles_per_run), size(beam%bunch(1)%particle))

    n = ix1_p-ix0_p
    dat_size = n * storage_size(beam%bunch(1)%particle(1)) / 8
    call ltt_print_mpi_info (lttp, ltt_com, 'Master: Position data size to slave: ' // int_str(slave_rank))
    call mpi_send (n, 1, MPI_INTEGER, slave_rank, num_tag$, MPI_COMM_WORLD, ierr)
    call mpi_send (ix0_p, 1, MPI_INTEGER, slave_rank, num_tag$, MPI_COMM_WORLD, ierr)
    call ltt_print_mpi_info (lttp, ltt_com, 'Master: More Initial positions to slave: ' // int_str(slave_rank) // &
                                                 '  For particles: [' // int_str(ix0_p+1) // ':' // int_str(ix1_p) // ']', .true.)
    do nb = 1, size(beam%bunch)
      call mpi_send (beam%bunch(nb)%particle(ix0_p+1:ix1_p), dat_size, MPI_BYTE, slave_rank, particle_tag$, MPI_COMM_WORLD, ierr)
    enddo

    ix0_p = ix1_p
  enddo

  ! Write results

  call ltt_write_beam_averages (lttp, beam_data_sum)
  call ltt_write_sigma_matrix (lttp, beam_data_sum)
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
      if (file(1:1) == 'b') cycle  ! Skip binary files
      ix = index(file, '_')
      read (file(1:ix-1), *) n
      if (any(n == turn(1:n_out))) cycle
      n_out = n_out + 1
      turn(n_out) = n
    enddo

    do nn = 1, n_out
      do i = 1, size(file_list)
        file = file_list(i)
        if (file(1:1) == 'b') cycle   ! Skip binary files
        ix = index(file, '_')
        read (file(1:ix-1), *) n
        if (n /= turn(nn)) cycle

        str = file(ix+1:)
        ix = index(str, '_')
        read (str(1:ix-1), *) ib
        
        open (iu, file = trim(ltt_com%mpi_data_dir) // file, status = 'OLD')
        do
          read (iu, *, iostat = ios) ix, i_turn, orb%vec, orb%spin, orb%state
          if (ios /= 0) exit
          beam%bunch(ib)%particle(ix) = orb
        enddo
        close (iu)
      enddo
      call ltt_write_particle_data (lttp, ltt_com, i_turn, beam, .false.)
    enddo
  endif

  if (lttp%beam_binary_output_file /= '') then
    do i = 1, size(file_list)
      file = file_list(i)
      if (file(1:1) /= 'b') cycle
      call hdf5_read_beam(trim(ltt_com%mpi_data_dir) // file, beam2, err_flag)
      read (file(2:), *) ix
      do ib = 1, size(beam2%bunch)
        n = size(beam2%bunch(ib)%particle)
        beam%bunch(ib)%particle(ix+1:ix+n) = beam2%bunch(ib)%particle
      enddo
    enddo
    call write_beam_file (lttp%beam_binary_output_file, beam)
  endif

  call system_command ('rm -rf ' // trim(ltt_com%mpi_data_dir))

  ! And end

  call run_timer ('ABS', time_now)
  print '(a, f8.2)', 'Tracking time (min):', (time_now - ltt_com%time_start) / 60

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
else  ! Is a slave

  call ltt_print_mpi_info (lttp, ltt_com, 'Slave Starting...')
  call mpi_recv (ltt_com%n_particle, 1, MPI_INTEGER, master_rank$, num_tag$, MPI_COMM_WORLD, stat, ierr)

  do
    ! Init positions from master
    call ltt_print_mpi_info (lttp, ltt_com, 'Slave: Waiting for data size info.')
    call mpi_recv (n, 1, MPI_INTEGER, master_rank$, num_tag$, MPI_COMM_WORLD, stat, ierr)
    call mpi_recv (ltt_com%mpi_ix0_particle, 1, MPI_INTEGER, master_rank$, num_tag$, MPI_COMM_WORLD, stat, ierr)
    call ltt_print_mpi_info (lttp, ltt_com, 'Slave: Waiting for position info for ' // int_str(n) // ' particles.')
    call reallocate_beam(beam2, max(1, beam_init%n_bunch), n)
    dat_size = n * storage_size(beam2%bunch(1)%particle(1)) / 8
    do nb = 1, size(beam2%bunch)
      call ltt_print_mpi_info (lttp, ltt_com, 'Slave: Waiting for bunch: ' // int_str(nb))
      call mpi_recv (beam2%bunch(nb)%particle, dat_size, MPI_BYTE, master_rank$, particle_tag$, MPI_COMM_WORLD, stat, ierr)
    enddo

    ! Run
    call ltt_print_mpi_info (lttp, ltt_com, 'Slave: Track beam.')
    call ltt_run_beam_mode(lttp, ltt_com, beam_init, beam_data, beam2)  ! Beam tracking
    bd_size = size(beam_data%turn(1)%bunch) * storage_size(beam_data%turn(1)%bunch(1)) / 8

    ! Send data to master
    call ltt_print_mpi_info (lttp, ltt_com, 'Slave: Sending Beam Size Data... ' // int_str(bd_size))
    do ix = lbound(beam_data%turn, 1), ubound(beam_data%turn, 1)
      call mpi_send ([beam_data%turn(ix)%status, beam_data%turn(ix)%i_turn], 2, MPI_INTEGER, &
                                                            master_rank$, results_tag$, MPI_COMM_WORLD, ierr)
      call mpi_send (beam_data%turn(ix)%bunch, bd_size, MPI_BYTE, master_rank$, results_tag$, MPI_COMM_WORLD, ierr)
    enddo

    ! Query Master if more tracking needed
    call ltt_print_mpi_info (lttp, ltt_com, 'Slave: Query am-i-done to master...')
    call mpi_recv (am_i_done, 1, MPI_LOGICAL, master_rank$, is_done_tag$, MPI_COMM_WORLD, stat, ierr)
    if (am_i_done) exit
  enddo

  call ltt_print_mpi_info (lttp, ltt_com, 'Slave: All done!')
  call mpi_finalize(ierr)

endif

end program

