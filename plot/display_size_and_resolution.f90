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

use sim_utils, dummy => display_size_and_resolution

implicit none

interface
  subroutine display_size_and_res(ix_screen, x_size, y_size, x_res, y_res) bind (c)
    use, intrinsic :: iso_c_binding
    real(c_double) x_size, y_size, x_res, y_res
    integer(c_int), value :: ix_screen
  end subroutine
end interface

real(rp) x_size, y_size, x_res, y_res
real(rp) :: dpi
integer ix_screen, s, l, ios
character(*), parameter :: r_name = 'display_size_and_resolution'

!

call display_size_and_res(ix_screen, x_size, y_size, x_res, y_res)

call get_environment_variable('ACC_DPI_RESOLUTION', length=l, status=s)
if (s == 0) then
  block
    character(l) :: dpi_str
    call get_environment_variable('ACC_DPI_RESOLUTION', value=dpi_str)
    read (dpi_str, *, iostat = ios) dpi
    if (ios == 0) then
      x_res = dpi / 25.4_rp
      y_res = dpi / 25.4_rp
    else
      call out_io (s_error$, r_name, 'ACC_DPI_RESOLUTION ENVIRONMENT VARIABLE IS NOT A REAL NUMBER: ' // dpi_str)
    endif
  end block
end if

end subroutine
