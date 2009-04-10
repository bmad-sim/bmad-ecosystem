!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine init_wall (wall)

  use synrad_struct
  use synrad_interface

  implicit none

  type (wall_struct) wall

  integer i

!

  wall%seg(:)%power%power_tot = 0
  wall%seg(:)%power%power_per_len = 0
  wall%seg(:)%power%power_per_area = 0
  wall%seg(:)%power%ix_ele_source = 0
  wall%seg(:)%power%s_source = 0
  wall%seg(:)%power%n_source = 0

  do i = 1, size(wall%seg)
    nullify (wall%seg(i)%power%sources)
  enddo

end subroutine
