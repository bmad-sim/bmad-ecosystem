!+
! Subroutine type_real_8_taylors (y, switch_z)
!
! Subroutine to type out the taylor series from a real_8 array.
!
! Modules needed:
!   use accelerator
!
! Input
!   y(6)     -- Real_8: 6 taylor series: (x, P_x, y, P_y, P_z, -z)
!   switch_z -- Logical, optional: If True then switch from PTC coordinate
!                       convention to BMAD's. Default is True.
!-

subroutine type_real_8_taylors (y, switch_z)

  use accelerator

  implicit none

  type (real_8) y(:)
  type (taylor_struct) b_taylor(6)

  logical, optional :: switch_z

!

  call init_taylor (b_taylor)
  call real_8_to_taylor (y, b_taylor, switch_z)
  call type_taylors (b_taylor)
  call kill_taylor (b_taylor)

end subroutine
