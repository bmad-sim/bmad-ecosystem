!+
! Subroutine tao_output_cmd (s, what)
!
! 
! Input:
!   s        -- tao_super_universe_struct
!
!  Output:
!   s        -- tao_super_universe_struct
!-

subroutine tao_output_cmd (s, what)

use tao_mod
use quick_plot

type (tao_super_universe_struct) s

character(*) what
character(20) :: r_name = 'tao_output_cmd'

integer ix

!

call string_trim (what, what, ix)
call string_trim (what(1:ix), what, ix)
select case (what)

! hard

case ('hard')
  call qp_open_page ('PS')
  call tao_plot_out (s%plot_page)   ! Update the plotting window
  call qp_close_page
  call qp_select_page (s%plot_page%id_window)  ! Back to X-windows
  call tao_plot_out (s%plot_page)   ! Update the plotting window
  call out_io (s_info$, r_name, 'Postscript file created: quick_plot.ps')

! error

case default

  call out_io (s_error$, r_name, 'UNKNOWN "WHAT": ' // what)

end select

end subroutine 
