!+
! Subroutine display_size_and_resolution (ix_screen, x_size, y_size, x_res, y_res)
!
! Routine to return the size and resolution of a display screen.
! Note: All output numbers will be set to zero if there is a problem obtaining the display information.
! 
! Input:
!    ix_screen -- Screen index. Generally this should be 0.
!
! Output:
!    x_size  -- real(rp): Horizontal screen size in mm.
!    y_size  -- real(rp): Vertical screen size in mm.
!    x_res   -- real(rp): Horizontal resolution in pixels / mm.
!    y_res   -- real(rp): Vertical resolution in pixels / mm.
!-


subroutine display_size_and_resolution (ix_screen, x_size, y_size, x_res, y_res)

use precision_def

implicit none

interface
  subroutine display_size_and_res(ix_screen, x_size, y_size, x_res, y_res) bind (c)
    use, intrinsic :: iso_c_binding
    real(c_double) x_size, y_size, x_res, y_res
    integer(c_int), value :: ix_screen
  end subroutine
end interface

real(rp) x_size, y_size, x_res, y_res
integer ix_screen

!

call display_size_and_res(ix_screen, x_size, y_size, x_res, y_res)

end subroutine
