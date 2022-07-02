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

real(rp) particles_per_thread, now_time

integer num_slaves, slave_rank, stat(MPI_STATUS_SIZE)
integer i, n, nn, nb, ib, ix, ierr, rc, leng, bd_size, storage_size, dat_size, ix_stage
integer ip0, ip1, nt0, nt1, mpi_n_proc, iu, ios, n_particle, i_turn
integer, allocatable :: ix_stop_turn(:), ixp_slave(:)

integer, parameter :: base_tag$  = 1000

logical am_i_done, err_flag, ok
logical, allocatable :: slave_is_done(:), stop_here(:)

character(200) pwd
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

call ltt_read_params(lttp, ltt_com)
call run_timer ('ABS', ltt_com%time_start)

if (lttp%simulation_mode /= 'BEAM' .and. ltt_com%mpi_rank /= master_rank$) then
  call mpi_finalize(ierr)
  stop
endif

! Only the master should create a map file if a file is to be created.

call ltt_init_params(lttp, ltt_com)

if (ltt_com%mpi_rank == master_rank$) then
  call ltt_print_mpi_info (lttp, ltt_com, 'Master: Init tracking', .true.)
  call ltt_init_tracking (lttp, ltt_com)
  call mpi_Bcast (0, 1, MPI_INTEGER, master_rank$, MPI_COMM_WORLD, ierr)
  if (ierr /= MPI_SUCCESS) call ltt_print_mpi_info (lttp, ltt_com, 'MPI ERROR!', .true.)
  call ltt_print_inital_info (lttp, ltt_com)

else
  call mpi_Bcast (ix, 1, MPI_INTEGER, master_rank$, MPI_COMM_WORLD, ierr)
  if (ierr /= MPI_SUCCESS) call ltt_print_mpi_info (lttp, ltt_com, 'MPI ERROR!', .true.)
  call ltt_init_tracking (lttp, ltt_com)
endif

! Calculation start.

select case (lttp%simulation_mode)
case ('CHECK');  call ltt_run_check_mode(lttp, ltt_com)  ! A single turn tracking check
case ('SINGLE'); call ltt_run_single_mode(lttp, ltt_com) ! Single particle tracking
case ('STAT');   call ltt_run_stat_mode(lttp, ltt_com)              ! Lattice statistics (radiation integrals, etc.).
case ('BEAM');   ! Handled below. Only the BEAM simulation mode uses mpi.
case default;    print *, 'BAD SIMULATION_MODE: ' // lttp%simulation_mode
end select

! Not using mpi if there is only one thread.

if (.not. ltt_com%using_mpi) then
  print '(a, i0)', 'Number of threads is one! (Need to use mpirun or mpiexec if on a single machine.)'
  call ltt_init_beam_distribution(lttp, ltt_com, beam)
  call ltt_run_beam_mode(lttp, ltt_com, 0, lttp%n_turns, beam)
  stop
endif

!-----------------------------------------
! MPI BEAM simulation

! Init beam distribution in master

if (ltt_com%mpi_rank == master_rank$) then
  call ltt_init_beam_distribution(lttp, ltt_com, beam)
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

particles_per_thread = real(n_particle, rp) / real((num_slaves), rp)
allocate (ixp_slave(0:num_slaves), stop_here(lttp%n_turns))

ixp_slave(0) = 0
do i = 1, num_slaves
  ixp_slave(i) = ixp_slave(i-1) + nint(particles_per_thread)
  if (lttp%debug .and. ltt_com%mpi_rank == master_rank$) print *, 'ixp:', i, ixp_slave(i)
enddo
ixp_slave(num_slaves) = n_particle
if (lttp%debug .and. ltt_com%mpi_rank == master_rank$) print *, 'ixp:', num_slaves, ixp_slave(num_slaves)

!

stop_here = .false.

select case (lttp%averages_output_every_n_turns)
case (-1)  ! Nothing to do
case (0);  stop_here(lttp%n_turns) = .true.
case default
  do i = 1, lttp%n_turns/lttp%averages_output_every_n_turns + 1
    n = i * lttp%averages_output_every_n_turns 
    if (n > lttp%n_turns) exit
    stop_here(n) = .true.
  enddo
end select

select case (lttp%particle_output_every_n_turns)
case (-1)  ! Nothing to do
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
ix_stop_turn(0) = 0
n = 0
do i = 1, lttp%n_turns
  if (.not. stop_here(i)) cycle
  n = n + 1
  ix_stop_turn(n) = i
  if (lttp%debug .and. ltt_com%mpi_rank == master_rank$) print *, 'ix_stop:', n, ix_stop_turn(n)
enddo


!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
! Master:

if (ltt_com%mpi_rank == master_rank$) then
  allocate (slave_is_done(num_slaves))
  slave_is_done = .false.

  print '(a, i0)', 'Number of processes (including Master): ', mpi_n_proc

  ! Init positions to slaves

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

  call ltt_write_particle_data (lttp, ltt_com, 0, beam)
  call ltt_write_averages_data (lttp, 0, beam)
  call ltt_write_custom (lttp, ltt_com, 0, beam = beam)

  ! Loop over tracking states

  do ix_stage = 1, ubound(ix_stop_turn, 1)
    do i = 1, num_slaves
      ! Get data from slave.
      slave_rank = MPI_ANY_SOURCE

      call mpi_recv (i, 1, MPI_INTEGER, slave_rank, base_tag$+ix_stage, MPI_COMM_WORLD, stat, ierr)
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

    ! Write results

    call ltt_print_mpi_info (lttp, ltt_com, 'Master: Writing data for stage ' // int_str(ix_stage), .true.)
    i_turn = ix_stop_turn(ix_stage)
    call ltt_write_particle_data (lttp, ltt_com, i_turn, beam)
    call ltt_write_averages_data (lttp, i_turn, beam)
    call ltt_write_custom (lttp, ltt_com, i_turn, beam = beam)
  enddo

  if (lttp%beam_binary_output_file /= '') then
    call write_beam_file (lttp%beam_binary_output_file, beam)
  endif

  ! And end

  call ltt_print_mpi_info (lttp, ltt_com, 'Master: All done!', .true.)
  call run_timer ('ABS', now_time)
  call ltt_write_line('# tracking_time = ' // real_str((now_time-ltt_com%time_start)/60, 4, 2), lttp, 0)
  call mpi_finalize(ierr)

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
else  ! Is a slave

  call ltt_print_mpi_info (lttp, ltt_com, 'Slave Starting...', .true.)
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
    nt0 = ix_stop_turn(ix_stage-1); nt1 = ix_stop_turn(ix_stage)
    call ltt_print_mpi_info (lttp, ltt_com, 'Slave: Track beam for stage: ' // int_str(ix_stage), .true.)
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
  call mpi_finalize(ierr)
endif

end program

