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
use tao_top10_mod
use quick_plot
use tao_plot_mod
use io_mod

implicit none

character(*) what, who
character(16) action
character(20) :: r_name = 'tao_output_cmd'
character(100) file_name

integer i, ix

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
  call tao_var_write (s%global%var_out_file)

case ('lattice')
  file_name = who
  if (file_name == ' ') file_name = 'lat_universe_#.bmad'
  do i = 1, size(s%u)
    ix = index(file_name, '#')
    if (ix /= 0) write (file_name, '(a, i0, a)') file_name(1:ix-1), i, trim(file_name(ix+1:))
    call write_bmad_lattice_file (file_name, s%u(i)%model%lat)
    call out_io (s_info$, r_name, 'Writen: ' // file_name)
  enddo

case ('digested')
  file_name = who
  if (file_name == ' ') file_name = 'digested_lat_universe_#.bmad'
  do i = 1, size(s%u)
    ix = index(file_name, '#')
    if (ix /= 0) write (file_name, '(a, i0, a)') file_name(1:ix-1), i, trim(file_name(ix+1:))
    call write_digested_bmad_file (file_name, s%u(i)%model%lat)
    call out_io (s_info$, r_name, 'Writen: ' // file_name)
  enddo

! error

case default

  call out_io (s_error$, r_name, 'UNKNOWN "WHAT": ' // what)

end select

end subroutine 
