module tao_mpi_mod

use tao_struct

implicit none

include 'mpif.h'

!
		
contains


!-------------------------------------------------------------------------
!+
! Subroutine tao_mpi_initialize() 
!                                          
! Initializes MPI processes.
!
! Output:
!    s%mpi    -- tao_mpi_struct: superuniverse MPI information
!
! Note: s%mpi%on will be set to False if only one task is found. 
!
! 
!-
subroutine tao_mpi_initialize() 
implicit none

integer :: mpierr, errorcode, tag, mpistatus(MPI_STATUS_SIZE)
integer numtasks, rank, hostname_length
character(MPI_MAX_PROCESSOR_NAME) hostname

character(20) :: r_name = 'tao_mpi_initialize'

!

call mpi_init(mpierr)                             ! Introduce yourself to the MPI daemon
if (mpierr /= MPI_SUCCESS) then
  print *,'Error starting MPI program. Terminating.'
  call mpi_abort(MPI_COMM_WORLD, errorcode, mpierr)
end if

call mpi_comm_rank(MPI_COMM_WORLD, rank, mpierr)  ! Get your rank number, store in myrank.  Master is rank 0.
call mpi_comm_size(MPI_COMM_WORLD, numtasks, mpierr)
call mpi_get_processor_name(hostname, hostname_length, mpierr)

!CALL mympi_register_derived_types()              ! Tell the daemon about derived types that will be 
												  ! shared among the cluster

! Set superuniverse mpi info
s%mpi%rank = rank
s%mpi%host_name = hostname(1:hostname_length)
s%mpi%max_rank = numtasks - 1
if(rank == 0) then
  s%mpi%master = .true.
else
  s%mpi%master = .false.
endif

! If there is only one task, turn MPI off
if (numtasks == 1) then
  s%mpi%on = .false.
else 
  s%mpi%on = .true. 
endif


end subroutine tao_mpi_initialize

!-------------------------------------------------------------------------
!+
! Subroutine tao_mpi_test_send_receive()
!   Simple program to test MPI send and receive 
!-
subroutine tao_mpi_test_send_receive()
implicit none

integer :: i, rank, mpierr, mpistatus(MPI_STATUS_SIZE), numtasks, tag
real(rp) :: dummy

call mpi_comm_rank(MPI_COMM_WORLD, rank, mpierr) 
call mpi_comm_size(MPI_COMM_WORLD, numtasks, mpierr)

! User tag. 
tag = 666 

if (rank == 0) then
  ! Master, receive from slaves
  do i = 1, numtasks - 1
    call MPI_RECV(dummy, 1, MPI_DOUBLE_PRECISION, i, tag,  MPI_COMM_WORLD, mpistatus, mpierr)  
    print *, '========received from' , i
  enddo
else
  ! Slave, send to master
  call MPI_SEND(dummy, 1, MPI_DOUBLE_PRECISION, 0, tag, MPI_COMM_WORLD, mpierr)
  print *, '======== sent from' , rank
endif


end subroutine


!-------------------------------------------------------------------------
!+
! Subroutine tao_broadcast_opt_vars_mpi(rank)
!
! Synchronizes the opt vars from process: rank  
!
! Input: 
!   rank, optional -- Integer: rank to broadcast from
!                              Default: 0
!
!-
subroutine tao_broadcast_opt_vars_mpi(rank)

use tao_var_mod, only: tao_get_opt_vars, tao_set_opt_vars

implicit none

real(rp), allocatable :: var_vec(:)

integer, optional :: rank
integer ::  master, mpierr, errorcode

!

if (present(rank)) then
  master = rank
else
  master = 0
endif

call tao_get_opt_vars (var_vec)

call MPI_BCAST (var_vec, size(var_vec), MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, mpierr)
if (mpierr /= MPI_SUCCESS) then
  print *,'Error broadcasting var_vec. Terminating.'
  call mpi_abort(MPI_COMM_WORLD, errorcode, mpierr)
end if

call tao_set_opt_vars (var_vec, s%global%optimizer_var_limit_warn)

end subroutine


!-------------------------------------------------------------------------
!+ subroutine tao_broadcast_chars_mpi(chars, rank)
!
! Input: 
!   rank, optional -- Integer: rank to broadcast from
!                              Default: 0
!
!-
subroutine tao_broadcast_chars_mpi(chars, rank)

implicit none

character(*) :: chars
integer, optional :: rank
integer ::  master, mpierr, errorcode

if (present(rank)) then
  master = rank
else
  master = 0
endif

call MPI_BCAST (chars, len(chars), MPI_CHARACTER, master, MPI_COMM_WORLD, mpierr)


end subroutine

!-------------------------------------------------------------------------
!+
! Subroutine tao_mpi_finalize()
!   Wrapper for mpi_finalize
!-
subroutine tao_mpi_finalize()
implicit none  
integer :: mpierr
  call mpi_finalize(mpierr)
end subroutine tao_mpi_finalize




end module
