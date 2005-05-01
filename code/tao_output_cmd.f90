!+
! Subroutine tao_output_cmd (what, who)
!
! 
! Input:
!
!  Output:
!-

subroutine tao_output_cmd (what, who)

use tao_mod
use quick_plot
use tao_single_mod
use tao_plot_mod

implicit none

character(*) what, who
character(16) action
character(20) :: r_name = 'tao_output_cmd'

integer ix

!

call string_trim (what, action, ix)

select case (action)

! hard

case ('hard')
  call qp_open_page ('PS')
  call tao_plot_out ()   ! Update the plotting window
  call qp_close_page
  call qp_select_page (s%plot_page%id_window)  ! Back to X-windows
  call tao_plot_out ()   ! Update the plotting window

  if (s%global%print_command == ' ') then
    call out_io (s_fatal$, r_name, &
        'P%PRINT_COMMAND NEEDS TO BE SET TO SEND THE PS FILE TO THE PRINTER!')
    return
  endif

  call system (trim(s%global%print_command) // ' quick_plot.ps')
  call out_io (s_blank$, r_name, 'Printing with command: ' // &
                                              s%global%print_command)

case ('gif')
  call qp_open_page ('GIF', x_len = s%plot_page%size(1), &
           y_len = s%plot_page%size(2), units = 'POINTS', plot_file = who)
  call tao_plot_out ()   ! Update the plotting window
  call qp_close_page
  call qp_select_page (s%plot_page%id_window)  ! Back to X-windows
  call tao_plot_out ()   ! Update the plotting window
  call out_io (s_info$, r_name, 'GIF file created: quick_plot.gif')

! ps

case ('ps')
  call qp_open_page ('PS', plot_file = who)
  call tao_plot_out ()   ! Update the plotting window
  call qp_close_page
  call qp_select_page (s%plot_page%id_window)  ! Back to X-windows
  call tao_plot_out ()   ! Update the plotting window

! variables

case ('var')
  call tao_var_write (s%global%var_out_file, .true.)

! error

case default

  call out_io (s_error$, r_name, 'UNKNOWN "WHAT": ' // what)

end select

end subroutine 
