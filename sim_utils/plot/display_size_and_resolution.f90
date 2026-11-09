!+
! Subroutine display_size_and_resolution (ix_screen, x_size, y_size, x_res, y_res, x_pix, y_pix)
!
! Routine to return the size and resolution of a display screen.
! Note: x_size and y_size will be set to zero if there is a problem obtaining the display information.
! If the display does not report a plausible physical size (for example, under WSL), the resolution
! is estimated by assuming a screen height of 11 inches. If the size in pixels is not known either,
! the resolution is set to the standard 96 pixels / inch.
! The environment variable ACC_DPI_RESOLUTION, if set, overrides the resolution (in pixels / inch).
! 
! Input:
!    ix_screen -- Screen index. Generally this should be 0.
!
! Output:
!    x_size  -- real(rp): Horizontal screen size in mm.
!    y_size  -- real(rp): Vertical screen size in mm.
!    x_res   -- real(rp): Horizontal resolution in pixels / mm.
!    y_res   -- real(rp): Vertical resolution in pixels / mm.
!    x_pix   -- integer, optional: Horizontal screen size in pixels. Zero if not known.
!    y_pix   -- integer, optional: Vertical screen size in pixels. Zero if not known.
!-


subroutine display_size_and_resolution (ix_screen, x_size, y_size, x_res, y_res, x_pix, y_pix)

use sim_utils, dummy => display_size_and_resolution
use, intrinsic :: iso_c_binding, only: c_int

implicit none

interface
  subroutine display_size_and_res(ix_screen, x_size, y_size, x_res, y_res, x_pix, y_pix) bind (c)
    use, intrinsic :: iso_c_binding
    real(c_double) x_size, y_size, x_res, y_res
    integer(c_int) x_pix, y_pix
    integer(c_int), value :: ix_screen
  end subroutine
end interface

real(rp) x_size, y_size, x_res, y_res
real(rp) :: dpi
real(rp), parameter :: assumed_screen_height = 11  ! inches. Used if the physical size is not known.
integer, optional :: x_pix, y_pix
integer ix_screen, s, l, ios
integer(c_int) nx, ny
character(*), parameter :: r_name = 'display_size_and_resolution'

!

call display_size_and_res(ix_screen, x_size, y_size, x_res, y_res, nx, ny)
if (present(x_pix)) x_pix = nx
if (present(y_pix)) y_pix = ny

! Some displays report a physical size of zero or nonsense. Estimate the resolution from the
! screen height in pixels. Assuming a fixed height in inches does better than assuming a fixed
! resolution since high resolution screens are common.

if (x_res < 50/25.4_rp .or. x_res > 400/25.4_rp .or. y_res < 50/25.4_rp .or. y_res > 400/25.4_rp) then
  if (ny > 0) then
    dpi = min(400.0_rp, max(72.0_rp, ny / assumed_screen_height))
  else
    dpi = 96
  endif
  x_res = dpi / 25.4_rp
  y_res = x_res
endif

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
