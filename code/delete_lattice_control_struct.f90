!+
! Subroutine delete_lattice_control_struct (lat, ix_con)
! 
! Subroutine to delete a control structure from the lattice.
!
! Modules to use:
!   use bmad
!
! Input:
!   ix_con -- Integer: Index of lat%control(:) control structure
!                   to be eliminated.
!
! Output:
!   lat -- lat_struct: Lat with control structure fixed.
!-

#include "CESR_platform.inc"

subroutine delete_lattice_control_struct (lat, ix_con)

  use bmad_struct
  use bmad_interface, except_dummy => delete_lattice_control_struct

  implicit none

  type (lat_struct), target :: lat
  type (ele_struct), pointer :: lord, slave

  integer ix_con, i, ix_slave, ix_lord, ix_ic, n_ic, n_con

! 

  n_con = lat%n_control_max
  if (ix_con < 1 .or. ix_con > n_con) then
    print *, 'ERROR IN DELETE_LATTICE_CONTROL_STRUCT: CONTROL INDEX OUT OF RANGE:', ix_con
    call err_exit
  endif

  ix_slave = lat%control(ix_con)%ix_slave
  ix_lord = lat%control(ix_con)%ix_lord

! fix slave

  slave => lat%ele(ix_slave)
  do i = slave%ic1_lord, slave%ic2_lord
    if (lat%ic(i) == ix_con) then
      ix_ic = i
      exit
    endif
  enddo

  slave%n_lord = slave%n_lord - 1
  if (slave%n_lord == 0) then
    slave%lord_status  = not_a_lord$
    slave%slave_status = free$
    slave%ic1_lord = 0
    slave%ic2_lord = -1
  else
    slave%ic2_lord = slave%ic2_lord - 1
  endif

! fix lord

  lord => lat%ele(ix_lord)
  lord%n_slave = lord%n_slave - 1
  if (lord%n_slave == 0) then
    lord%ix1_slave = 0
    lord%ix2_slave = -1
  else
    lord%ix2_slave = lord%ix2_slave - 1
  endif

! fix lat%ic

  n_ic = lat%n_ic_max
  lat%ic(ix_ic+1:n_ic) = lat%ic(ix_ic:n_ic-1)
  lat%n_ic_max = lat%n_ic_max - 1
  where (lat%ic > ix_con) lat%ic = lat%ic - 1


! fix lat%control

  lat%control(ix_con+1:n_con) = lat%control(ix_con:n_con-1)
  lat%n_control_max = lat%n_control_max - 1

! fix other lord and slave elements

  where (lat%ele%ix1_slave > ix_con) lat%ele%ix1_slave = &
                                           lat%ele%ix1_slave - 1
  where (lat%ele%ix2_slave > ix_con) lat%ele%ix2_slave = &
                                            lat%ele%ix2_slave - 1
                                        
  where (lat%ele%ic1_lord > ix_ic) lat%ele%ic1_lord = &
                                            lat%ele%ic1_lord - 1
  where (lat%ele%ic2_lord > ix_ic) lat%ele%ic2_lord = &
                                            lat%ele%ic2_lord - 1


end subroutine
