!+
! Subroutine tao_output_cmd (what, file_name)
!
! 
! Input:
!
!  Output:
!-

subroutine tao_output_cmd (what, file_name)

use tao_mod
use tao_top10_mod
use quick_plot
use tao_plot_mod
use io_mod

implicit none

character(*) what, file_name
character(16) action
character(20) :: r_name = 'tao_output_cmd'
character(100) file_name2

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
           y_len = s%plot_page%size(2), units = 'POINTS', plot_file = file_name)
  call tao_plot_out ()   ! Update the plotting window
  call qp_close_page
  call qp_select_page (s%plot_page%id_window)  ! Back to X-windows
  call tao_plot_out ()   ! Update the plotting window
  call out_io (s_info$, r_name, 'GIF file created: quick_plot.gif')

! ps

case ('ps')
  call qp_open_page ('PS', plot_file = file_name)
  call tao_plot_out ()   ! Update the plotting window
  call qp_close_page
  call qp_select_page (s%plot_page%id_window)  ! Back to X-windows
  call tao_plot_out ()   ! Update the plotting window

! variables

case ('var')
  if (file_name == ' ') then
    call tao_var_write (s%global%var_out_file)
  else
    call tao_var_write (file_name)
  endif

case ('lattice')
  file_name2 = file_name
  if (file_name2 == ' ') file_name2 = 'lat_universe_#.bmad'
  do i = 1, size(s%u)
    ix = index(file_name2, '#')
    if (ix /= 0) write (file_name2, '(a, i0, a)') file_name2(1:ix-1), i, trim(file_name2(ix+1:))
    call write_bmad_lattice_file (file_name2, s%u(i)%model%lat)
    call out_io (s_info$, r_name, 'Writen: ' // file_name2)
  enddo

case ('digested')
  file_name2 = file_name
  if (file_name2 == ' ') file_name2 = 'digested_lat_universe_#.bmad'
  do i = 1, size(s%u)
    ix = index(file_name2, '#')
    if (ix /= 0) write (file_name2, '(a, i0, a)') file_name2(1:ix-1), i, trim(file_name2(ix+1:))
    call write_digested_bmad_file (file_name2, s%u(i)%model%lat)
    call out_io (s_info$, r_name, 'Writen: ' // file_name2)
  enddo

! error

case default

  call out_io (s_error$, r_name, 'UNKNOWN "WHAT": ' // what)

end select

end subroutine 
