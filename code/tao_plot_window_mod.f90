module tao_plot_window_mod

use tao_mod
use quick_plot

logical :: has_been_created = .false.

contains
!------------------------------------------------------------------------------
!+
! Subroutine tao_create_plot_window ()
!
! Subroutine to create the plot window.
! This soubroutine knows not to create a second window if one already exists.
!-

subroutine tao_create_plot_window ()

implicit none

character(24) :: r_name = 'tao_create_plot_window'

  ! Open a plot window

  if (.not. s%global%plot_on) then
    call out_io (s_info$, r_name, 'global%plot_on = F so no plot window created.')
    return
  endif

  if (has_been_created) return

  has_been_created = .true.

  call qp_open_page (s%plotting%plot_display_type, s%plotting%id_window, &
                                s%plotting%size(1), s%plotting%size(2), 'POINTS')
  call qp_set_layout (page_border = s%plotting%border)

  ! Set default font size

  call qp_set_text_attrib ("TEXT",  height = s%plotting%text_height) 
  call qp_set_text_attrib ("MAIN_TITLE",  height = s%plotting%text_height) 
  call qp_set_text_attrib ("GRAPH_TITLE",  height = s%plotting%text_height) 
  call qp_set_text_attrib ("AXIS_LABEL",  height = s%plotting%text_height) 
  call qp_set_text_attrib ("AXIS_NUMBERS",  height = s%plotting%text_height) 
  call qp_set_text_attrib ("LEGEND",  height = s%plotting%text_height) 

end subroutine tao_create_plot_window

!contains
!------------------------------------------------------------------------------
!+
! Subroutine tao_destroy_plot_window
!
! Closes an already created plot window created by tao_create_plot_window
!
!-

subroutine tao_destroy_plot_window ()

implicit none

character(24) :: r_name = 'tao_destroy_plot_window'

  if (.not. has_been_created) return

  call qp_close_page
  has_been_created = .false.

end subroutine tao_destroy_plot_window

end module tao_plot_window_mod
