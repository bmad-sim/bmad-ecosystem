program long_term_tracking

use lt_tracking_mod
use mpi
use directory_mod

implicit none

type (ltt_params_struct) lttp
type (ltt_com_struct), target :: ltt_com
type (beam_struct), target :: beam, beam2
type (bunch_struct), pointer :: bunch
type (coord_struct) orb
type (coord_struct), pointer :: particle

real(rp) particles_per_thread, now_time

integer num_slaves, slave_rank, stat(MPI_STATUS_SIZE)
integer i, n, nn, nb, ib, ix, ierr, rc, leng, bd_size, storage_size, dat_size, ix_stage
integer ip, ip0, ip1, nt0, nt1, mpi_n_proc, iu, ios, n_particle, i_turn, n_part_tot, seed
integer, allocatable :: ix_stop_turn(:), ixp_slave(:)

integer, parameter :: base_tag$  = 1000

logical am_i_done, err_flag, ok, too_many_dead
logical, allocatable :: stop_here(:), slave_working(:)

character(200) pwd
character(100) line, path, basename
character(40) file, str
character(40), allocatable :: file_list(:)
character(MPI_MAX_PROCESSOR_NAME) name

!

track1_preprocess_ptr => ltt_track1_preprocess
track1_bunch_hook_ptr => ltt_track1_bunch_hook

! Initialize MPI

call run_timer ('ABS', ltt_com%time_start)

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

call ltt_read_params(lttp, ltt_com)

if (lttp%simulation_mode /= 'BEAM' .and. ltt_com%mpi_rank /= master_rank$) then
  call mpi_finalize(ierr)
  stop
endif

! Synchronize ramper ran state.
! The problem being solved here is due to the fact that each thread tracks a portion of the beam. If rampers 
! are used to, for example, simulate something like RF noise, all particles of a given bunch going through a 
! given element on a given turn need to see the RF waveform shifted by the same noise vector (it is always
! assumed by the long_term_tracking program that a given random process simulated by a ramper has variation
! that is small on the time scale of a bunch passage). That is, all threads must use the same random number
! sequence when dealing with rampers. However, at the same time, the random number sequence for radiation
! excitation must be different for all particles. This problem is solved by using a common random state
! for the rampers while each thread has its own unique random state for everything else 


if (ltt_com%mpi_rank == master_rank$) then
  seed = lttp%random_seed 
  if (seed /= 0) seed = seed + 1234
endif

call mpi_Bcast (seed, 1, MPI_INTEGER, master_rank$, MPI_COMM_WORLD, ierr)
if (ierr /= MPI_SUCCESS) call ltt_print_mpi_info (lttp, ltt_com, 'MPI_BCAST RAMPER_RAN_STATE ERROR!', .true.)

call ran_seed_put(seed)
ltt_com%ramper_ran_state = pointer_to_ran_state()

call mpi_barrier (MPI_COMM_WORLD, ierr)
if (ierr /= MPI_SUCCESS) call ltt_print_mpi_info (lttp, ltt_com, 'MPI ERROR!', .true.)

! Only the master should create a map file if a file is to be created.

call ltt_init_params(lttp, ltt_com)

call mpi_Bcast (ltt_com%time_start, 1, MPI_DOUBLE_PRECISION, master_rank$, MPI_COMM_WORLD, ierr)
if (ierr /= MPI_SUCCESS) call ltt_print_mpi_info (lttp, ltt_com, 'MPI_BCAST TIME_START ERROR!', .true.)

if (ltt_com%mpi_rank == master_rank$) then
  call ltt_print_mpi_info (lttp, ltt_com, 'Master: Init tracking', .true.)
  call ltt_init_tracking (lttp, ltt_com, beam)
  call mpi_barrier (MPI_COMM_WORLD, ierr)
  if (ierr /= MPI_SUCCESS) call ltt_print_mpi_info (lttp, ltt_com, 'MPI ERROR!', .true.)
  call ltt_print_inital_info (lttp, ltt_com)

else
  call mpi_barrier (MPI_COMM_WORLD, ierr)
  if (ierr /= MPI_SUCCESS) call ltt_print_mpi_info (lttp, ltt_com, 'MPI ERROR!', .true.)
  call ltt_init_tracking (lttp, ltt_com, beam)
endif

! Calculation start.

select case (lttp%simulation_mode)
case ('CHECK');  call ltt_run_check_mode(lttp, ltt_com)  ! A single turn tracking check
case ('SINGLE'); call ltt_run_single_mode(lttp, ltt_com) ! Single particle tracking
case ('STAT');   call ltt_run_stat_mode(lttp, ltt_com)              ! Lattice statistics (radiation integrals, etc.).
case ('BEAM')
  if (.not. ltt_com%using_mpi) then
    print '(a, i0)', 'Number of threads is one! (Need to use mpirun or mpiexec if on a single machine.)'
    call ltt_run_beam_mode(lttp, ltt_com, lttp%ix_turn_start, lttp%ix_turn_stop, beam)
    stop
  endif

case ('INDIVIDUAL')
  if (.not. ltt_com%using_mpi) then
    print '(a, i0)', 'Number of threads is one! (Need to use mpirun or mpiexec if on a single machine.)'
    call ltt_run_individual_mode(lttp, ltt_com)
    stop
  endif

case default
  print *, 'BAD SIMULATION_MODE: ' // lttp%simulation_mode
end select

! Not using mpi if there is only one thread.


!-----------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------
! INDIVIDUAL simulation

if (lttp%simulation_mode == 'INDIVIDUAL') then
  dat_size = storage_size(beam%bunch(1)%particle(1)) / 8

  if (ltt_com%mpi_rank == master_rank$) then
    print '(a, i0)', 'Number of processes (including Master): ', mpi_n_proc
    call ltt_print_mpi_info (lttp, ltt_com, 'Master: Initial Ramper Ran State: ' // int_str(ltt_com%ramper_ran_state%ix))

    n_particle = size(beam%bunch(1)%particle)
    allocate (slave_working(num_slaves), ixp_slave(num_slaves))
    slave_working = .false.

    do ib = 1, size(beam%bunch)
      do ip = 1, size(beam%bunch(ib)%particle)
        particle => beam%bunch(ib)%particle(ip)

        if (all(slave_working)) then
          slave_rank = MPI_ANY_SOURCE
          call ltt_print_mpi_info (lttp, ltt_com, 'Master: Waiting for data.')
          call mpi_recv (orb, dat_size, MPI_BYTE, slave_rank, base_tag$+3, MPI_COMM_WORLD, stat, ierr)
          if (ierr /= MPI_SUCCESS) call ltt_print_mpi_info (lttp, ltt_com, 'MPI ERROR!', .true.)
          slave_rank = stat(MPI_SOURCE)  ! Slave rank
          nn = ixp_slave(slave_rank)
          beam%bunch(ib)%particle(nn) = orb
          slave_working(slave_rank) = .false.
          call ltt_print_mpi_info (lttp, ltt_com, 'Master: Got data from slave: ' // int_str(slave_rank))
        endif

        do ix = 1, num_slaves
          if (.not. slave_working(ix)) exit
        enddo

        call ltt_print_mpi_info (lttp, ltt_com, 'Master: Tell slave ' // int_str(ix) // ' to be ready to track using Slave')
        call mpi_send (1, 1, MPI_INTEGER, ix, base_tag$+1, MPI_COMM_WORLD, ierr)
        if (ierr /= MPI_SUCCESS) call ltt_print_mpi_info (lttp, ltt_com, 'MPI ERROR #2!', .true.)             ! Tell slave more tracking needed.

        call ltt_print_mpi_info (lttp, ltt_com, 'Master: Starting particle ' // int_str(ip) // ' using Slave: ' // int_str(ix))
        call mpi_send (particle, dat_size, MPI_BYTE, ix, base_tag$+2, MPI_COMM_WORLD, ierr)
        if (ierr /= MPI_SUCCESS) call ltt_print_mpi_info (lttp, ltt_com, 'MPI ERROR #3!', .true.)

        ixp_slave(ix) = ip
        slave_working(ix) = .true.
      enddo
    enddo

    ! Finish getting data

    call ltt_print_mpi_info (lttp, ltt_com, 'Master: Finished broadcasting particle runs. Now collecting final data.')
    do
      if (all(.not. slave_working)) exit
      call mpi_recv (orb, dat_size, MPI_BYTE, slave_rank, base_tag$+3, MPI_COMM_WORLD, stat, ierr)
      if (ierr /= MPI_SUCCESS) call ltt_print_mpi_info (lttp, ltt_com, 'MPI ERROR!', .true.)
      slave_rank = stat(MPI_SOURCE)  ! Slave rank
      nn = ixp_slave(slave_rank)
      beam%bunch(ib)%particle(nn) = orb
      slave_working(slave_rank) = .false.
    enddo

    ! Tell slaves we are done.

    do i = 1, num_slaves
      if (.not. slave_working(i)) cycle
      call ltt_print_mpi_info (lttp, ltt_com, 'Slave: Waiting for init position info.')
      call mpi_send (0, 1, MPI_INTEGER, ix, base_tag$+2, MPI_COMM_WORLD, ierr)       ! No more tracking.
    enddo

    ! And write data.

    if (lttp%per_particle_output_file /= '') call write_beam_file (lttp%per_particle_output_file, beam, .true., ascii4$)
    if (lttp%beam_binary_output_file /= '')  call write_beam_file (lttp%beam_binary_output_file, beam, .true., hdf5$)
    call mpi_finalize(ierr)

  !---------------------------------------------------------
  ! INDIVIDUAL Slave

  else
    call ltt_print_mpi_info (lttp, ltt_com, 'Slave Starting...', .true.)
    call ltt_print_mpi_info (lttp, ltt_com, 'Slave: Initial Ramper Ran State: ' // int_str(ltt_com%ramper_ran_state%ix))

    do
      call ltt_print_mpi_info (lttp, ltt_com, 'Slave: Waiting for Master command.')
      call mpi_recv (nn, 1, MPI_INTEGER, master_rank$, base_tag$+1, MPI_COMM_WORLD, stat, ierr)
      if (nn == 0) exit
      call ltt_print_mpi_info (lttp, ltt_com, 'Slave: Waiting for particle info...')
      call mpi_recv (particle, dat_size, MPI_BYTE, master_rank$, base_tag$+2, MPI_COMM_WORLD, stat, ierr)
      call ltt_print_mpi_info (lttp, ltt_com, 'Slave: Starting tracking.')
      call ltt_run_single_mode(lttp, ltt_com, particle)
      call mpi_send (particle, dat_size, MPI_BYTE, master_rank$, base_tag$+3, MPI_COMM_WORLD, ierr)
    enddo

    call ltt_print_mpi_info (lttp, ltt_com, 'Slave: Master says all done!')
    call mpi_finalize(ierr)
    stop
  endif
endif

!-----------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------
! BEAM simulation

! Init beam distribution in master

if (ltt_com%mpi_rank == master_rank$) then
  n_particle = size(beam%bunch(1)%particle)

  do i = 1, num_slaves
    call mpi_send (n_particle, 1, MPI_INTEGER, i, base_tag$-1, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) call ltt_print_mpi_info (lttp, ltt_com, 'MPI ERROR!', .true.)
  enddo

else
  call mpi_recv(n_particle, 1, MPI_INTEGER, master_rank$, base_tag$-1, MPI_COMM_WORLD, stat, ierr)
  if (ierr /= MPI_SUCCESS) call ltt_print_mpi_info (lttp, ltt_com, 'MPI ERROR!', .true.)
endif

! Calculate which particles go to which slaves and 
! at what turns the slaves have to report back the beam distribution to the master.

! A complication is that a custom data file may need to print lattice parameters and if there is ramping
! then the master lattice needs to get tracked through.
! So allocate a 1/2 share of particles to the Master for tracking.

particles_per_thread = n_particle / (num_slaves + 0.5)
allocate (ixp_slave(0:num_slaves), stop_here(lttp%n_turns))

ixp_slave(0) = nint(0.5_rp * particles_per_thread) ! Share for Master
if (lttp%debug .and. ltt_com%mpi_rank == master_rank$) print *, 'ixp:', 0, ixp_slave(0)

do i = 1, num_slaves
  ixp_slave(i) = nint((i+0.5_rp)*particles_per_thread)
  if (lttp%debug .and. ltt_com%mpi_rank == master_rank$) print *, 'ixp:', i, ixp_slave(i)
enddo

! Calculate the  stopping points.

stop_here = .false.

select case (lttp%averages_output_every_n_turns)
case (-1); stop_here(lttp%n_turns) = .true.
case (0);  stop_here(lttp%n_turns) = .true.
case default
  do i = 1, lttp%n_turns/lttp%averages_output_every_n_turns + 1
    n = i * lttp%averages_output_every_n_turns 
    if (n > lttp%n_turns) exit
    stop_here(n) = .true.
  enddo
end select

select case (lttp%particle_output_every_n_turns)
case (-1); stop_here(lttp%n_turns) = .true.
case (0);  stop_here(lttp%n_turns) = .true.
case default
  do i = 1, lttp%n_turns/lttp%particle_output_every_n_turns + 1
    n = i * lttp%particle_output_every_n_turns 
    if (n > lttp%n_turns) exit
    stop_here(n) = .true.
  enddo
end select

n = count(stop_here)
allocate (ix_stop_turn(0:n))
ix_stop_turn(0) = lttp%ix_turn_start
n = 0
do i = 1, lttp%n_turns
  if (.not. stop_here(i)) cycle
  n = n + 1
  ix_stop_turn(n) = i + lttp%ix_turn_start
  if (lttp%debug .and. ltt_com%mpi_rank == master_rank$) print *, 'ix_stop:', n, ix_stop_turn(n)
enddo

!------------------------------------------------------------------------------------------
! BEAM Master:

if (ltt_com%mpi_rank == master_rank$) then
  print '(a, i0)', 'Number of processes (including Master): ', mpi_n_proc
  call ltt_print_mpi_info (lttp, ltt_com, 'Master: Initial Ramper Ran State: ' // int_str(ltt_com%ramper_ran_state%ix))
  call reallocate_beam(beam2, max(1, ltt_com%beam_init%n_bunch), ixp_slave(0))

  ! Init positions to slaves

  call ltt_print_mpi_info (lttp, ltt_com, 'Master: Master will track particles: [1:' // int_str(ixp_slave(0)) // ']')
  do i = 1, num_slaves
    ip0 = ixp_slave(i-1)+1; ip1 = ixp_slave(i)
    n = ip1 + 1 - ip0
    dat_size = n * storage_size(beam%bunch(1)%particle(1)) / 8
    call ltt_print_mpi_info (lttp, ltt_com, 'Master: Initial positions to slave: ' // int_str(i) // &
                                                '  For particles: [' // int_str(ip0) // ':' // int_str(ip1) // ']')
    do nb = 1, size(beam%bunch)
      call mpi_send (beam%bunch(nb)%particle(ip0:ip1), dat_size, MPI_BYTE, i, base_tag$, MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) call ltt_print_mpi_info (lttp, ltt_com, 'MPI ERROR!', .true.)
    enddo
  enddo

  call ltt_print_mpi_info (lttp, ltt_com, 'Master: Writing initial data.', .true.)
  call ltt_write_particle_data (lttp, ltt_com, 0, beam)
  call ltt_write_averages_data (lttp, 0, beam)
  call ltt_write_custom (lttp, ltt_com, 0, beam = beam)

  ! Loop over tracking states

  do ix_stage = 1, ubound(ix_stop_turn, 1)
    call mpi_barrier(MPI_COMM_WORLD, ierr)

    n_part_tot = 0
    do nb = 1, size(beam%bunch)
      beam%bunch(nb)%n_live = count(beam%bunch(nb)%particle%state == alive$)
      n_part_tot = n_part_tot + size(beam%bunch(nb)%particle)
    enddo
    too_many_dead = (n_part_tot - sum(beam%bunch%n_live) >= nint(lttp%dead_cutoff * n_part_tot))
    call mpi_Bcast(too_many_dead, 1, MPI_LOGICAL, master_rank$, MPI_COMM_WORLD, ierr)
    if (too_many_dead) then
      call ltt_print_mpi_info (lttp, ltt_com, 'PARTICLE LOSS GREATER THAN SET BY DEAD_CUTOFF. STOPPING NOW.', .true.)
      exit
    endif

    do i = 1, num_slaves
      ! Get data from slave.
      slave_rank = MPI_ANY_SOURCE
      call mpi_recv (nn, 1, MPI_INTEGER, slave_rank, base_tag$+ix_stage, MPI_COMM_WORLD, stat, ierr)
      if (ierr /= MPI_SUCCESS) call ltt_print_mpi_info (lttp, ltt_com, 'MPI ERROR!', .true.)
      slave_rank = stat(MPI_SOURCE)  ! Slave rank
      ip0 = ixp_slave(slave_rank-1)+1; ip1 = ixp_slave(slave_rank)
      n = ip1 + 1 - ip0
      dat_size = n * storage_size(beam%bunch(1)%particle(1)) / 8
      do nb = 1, size(beam%bunch)
        call mpi_recv(beam%bunch(nb)%particle(ip0:ip1), dat_size, MPI_BYTE, slave_rank, base_tag$+ix_stage, MPI_COMM_WORLD, stat, ierr)
        if (ierr /= MPI_SUCCESS) call ltt_print_mpi_info (lttp, ltt_com, 'MPI ERROR!', .true.)
      enddo
      call ltt_print_mpi_info (lttp, ltt_com, 'Master: Gathered data for stage ' // int_str(ix_stage) // ' from Slave: ' // int_str(slave_rank))
    enddo

    ! Track with particles assigned to master

    do nb = 1, size(beam%bunch)
      beam2%bunch(nb)%particle = beam%bunch(nb)%particle(1:ixp_slave(0))
    enddo

    nt0 = ix_stop_turn(ix_stage-1); nt1 = ix_stop_turn(ix_stage)
    call ltt_print_mpi_info (lttp, ltt_com, 'Master: Track beam for stage: ' // int_str(ix_stage) // ' (End turn: ' // int_str(nt1) // ')', .true.)
    call ltt_run_beam_mode(lttp, ltt_com, nt0, nt1, beam2)

    do nb = 1, size(beam%bunch)
      beam%bunch(nb)%particle(1:ixp_slave(0)) = beam2%bunch(nb)%particle
    enddo

    ! Write results

    call ltt_print_mpi_info (lttp, ltt_com, 'Master: Writing data for stage ' // int_str(ix_stage), .true.)
    i_turn = ix_stop_turn(ix_stage)
    call ltt_write_particle_data (lttp, ltt_com, i_turn, beam)
    call ltt_write_averages_data (lttp, i_turn, beam)
    call ltt_write_custom (lttp, ltt_com, i_turn, beam = beam)
    call ltt_write_beam_binary_file(lttp, ltt_com, i_turn, beam)
  enddo

  ! And end

  call ltt_print_mpi_info (lttp, ltt_com, 'Master: All done!', .true.)
  call ltt_print_mpi_info (lttp, ltt_com, 'Master: Final Ramper Ran State: ' // int_str(ltt_com%ramper_ran_state%ix))
  call run_timer ('ABS', now_time)
  call ltt_write_line('# tracking_time = ' // real_str((now_time-ltt_com%time_start)/60, 4, 2), lttp, 0)
  call mpi_finalize(ierr)

!------------------------------------------------------------------------------------------
else  ! BEAM slave

  call ltt_print_mpi_info (lttp, ltt_com, 'Slave Starting...', .true.)
  call ltt_print_mpi_info (lttp, ltt_com, 'Slave: Initial Ramper Ran State: ' // int_str(ltt_com%ramper_ran_state%ix))
  ip0 = ixp_slave(ltt_com%mpi_rank-1)+1; ip1 = ixp_slave(ltt_com%mpi_rank)
  n = ip1 + 1 - ip0

  ! Init positions from master
  call ltt_print_mpi_info (lttp, ltt_com, 'Slave: Waiting for init position info for ' // int_str(n) // ' particles.')
  call reallocate_beam(beam2, max(1, ltt_com%beam_init%n_bunch), n)
  dat_size = n * storage_size(beam2%bunch(1)%particle(1)) / 8
  do nb = 1, size(beam2%bunch)
    call mpi_recv (beam2%bunch(nb)%particle, dat_size, MPI_BYTE, master_rank$, base_tag$, MPI_COMM_WORLD, stat, ierr)
    if (ierr /= MPI_SUCCESS) call ltt_print_mpi_info (lttp, ltt_com, 'MPI ERROR!', .true.)
  enddo

  do ix_stage = 1, ubound(ix_stop_turn, 1)
    call mpi_barrier(MPI_COMM_WORLD, ierr)
    call mpi_Bcast(too_many_dead, 1, MPI_LOGICAL, master_rank$, MPI_COMM_WORLD, ierr)
    if (too_many_dead) exit

    nt0 = ix_stop_turn(ix_stage-1); nt1 = ix_stop_turn(ix_stage)
    call ltt_print_mpi_info (lttp, ltt_com, 'Slave: Track beam for stage: ' // int_str(ix_stage) // ' (End turn: ' // int_str(nt1) // ')', .true.)
    call ltt_run_beam_mode(lttp, ltt_com, nt0, nt1, beam2)  ! Beam tracking

    ! Send data to master
    call ltt_print_mpi_info (lttp, ltt_com, 'Slave: Send particle data for stage: ' // int_str(ix_stage))
    call mpi_send (ltt_com%mpi_rank, 1, MPI_INTEGER, master_rank$, base_tag$+ix_stage, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) call ltt_print_mpi_info (lttp, ltt_com, 'MPI ERROR!', .true.)
    do nb = 1, size(beam2%bunch)
      call mpi_send (beam2%bunch(nb)%particle, dat_size, MPI_BYTE, master_rank$, base_tag$+ix_stage, MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) call ltt_print_mpi_info (lttp, ltt_com, 'MPI ERROR!', .true.)
    enddo
  enddo

  call ltt_print_mpi_info (lttp, ltt_com, 'Slave: All done!', .true.)
  call ltt_print_mpi_info (lttp, ltt_com, 'Slave: Final Ramper Ran State: ' // int_str(ltt_com%ramper_ran_state%ix))
  call mpi_finalize(ierr)
endif

end program

