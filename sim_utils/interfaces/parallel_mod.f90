!+
!   This module should be used for keeping track of the
!   rank of a process when doing parallel processing
!   a call to MPI_COMM_RANK(MPI_COMM_WORLD,global_rank,integer)
!   (where integer is any integer variable) will insure that only
!   the process designated as rank 0 prints to the screen and 
!   allows global_rank to be used to distinguish between processes
!   to allow differing processes to handle different portions of
!   the program
!-

module parallel_mod
  integer :: global_rank = 0
end module parallel_mod
