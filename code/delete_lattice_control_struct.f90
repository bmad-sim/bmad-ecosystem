!+
! Subroutine delete_lattice_control_struct (ring, ix_con)
! 
! Subroutine to delete a control structure from the lattice.
!
! Modules to use:
!   use bmad
!
! Input:
!   ix_con -- Integer: Index of ring%control_(:) control structure
!                   to be eliminated.
!
! Output:
!   ring -- Ring_struct: Ring with control structure fixed.
!-

#include "CESR_platform.inc"

subroutine delete_lattice_control_struct (ring, ix_con)

  use bmad_struct
  use bmad_interface, except => delete_lattice_control_struct

  implicit none

  type (ring_struct), target :: ring
  type (ele_struct), pointer :: lord, slave

  integer ix_con, i, ix_slave, ix_lord, ix_ic, n_ic, n_con

! 

  n_con = ring%n_control_max
  if (ix_con < 1 .or. ix_con > n_con) then
    print *, 'ERROR IN DELETE_LATTICE_CONTROL_STRUCT: CONTROL INDEX OUT OF RANGE:', ix_con
    call err_exit
  endif

  ix_slave = ring%control_(ix_con)%ix_slave
  ix_lord = ring%control_(ix_con)%ix_lord

! fix slave

  slave => ring%ele_(ix_slave)
  do i = slave%ic1_lord, slave%ic2_lord
    if (ring%ic_(i) == ix_con) then
      ix_ic = i
      exit
    endif
  enddo

  slave%n_lord = slave%n_lord - 1
  if (slave%n_lord == 0) then
    slave%control_type = free$
    slave%ic1_lord = 0
    slave%ic2_lord = -1
  else
    slave%ic2_lord = slave%ic2_lord - 1
  endif

! fix lord

  lord => ring%ele_(ix_lord)
  lord%n_slave = lord%n_slave - 1
  if (lord%n_slave == 0) then
    lord%ix1_slave = 0
    lord%ix2_slave = -1
  else
    lord%ix2_slave = lord%ix2_slave - 1
  endif

! fix ring%ic_

  n_ic = ring%n_ic_max
  ring%ic_(ix_ic+1:n_ic) = ring%ic_(ix_ic:n_ic-1)
  ring%n_ic_max = ring%n_ic_max - 1
  where (ring%ic_ > ix_con) ring%ic_ = ring%ic_ - 1


! fix ring%control_

  ring%control_(ix_con+1:n_con) = ring%control_(ix_con:n_con-1)
  ring%n_control_max = ring%n_control_max - 1

! fix other lord and slave elements

  where (ring%ele_%ix1_slave > ix_con) ring%ele_%ix1_slave = &
                                           ring%ele_%ix1_slave - 1
  where (ring%ele_%ix2_slave > ix_con) ring%ele_%ix2_slave = &
                                            ring%ele_%ix2_slave - 1
                                        
  where (ring%ele_%ic1_lord > ix_ic) ring%ele_%ic1_lord = &
                                            ring%ele_%ic1_lord - 1
  where (ring%ele_%ic2_lord > ix_ic) ring%ele_%ic2_lord = &
                                            ring%ele_%ic2_lord - 1


end subroutine
