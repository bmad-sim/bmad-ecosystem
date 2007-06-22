!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine init_wall (wall)

  use sr_struct
  use sr_interface

  implicit none

  type (wall_struct) wall

  integer i

!

  wall%seg(:)%sr_power%power = 0
  wall%seg(:)%sr_power%power_per_len = 0
  wall%seg(:)%sr_power%power_per_area = 0
  wall%seg(:)%sr_power%ix_ele_source = 0
  wall%seg(:)%sr_power%s_source = 0
  wall%seg(:)%sr_power%n_source = 0

  do i = 1, size(wall%seg)
    nullify (wall%seg(i)%sr_power%rays)
  enddo

end subroutine
