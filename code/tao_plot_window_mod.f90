module tao_plot_window_mod

use tao_interface
use quick_plot

logical :: has_been_created = .false.

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
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

  call qp_open_page (s%plot_page%plot_display_type, s%plot_page%id_window, &
                                s%plot_page%size(1), s%plot_page%size(2), 'POINTS')
  call qp_set_layout (page_border = s%plot_page%border)

  ! Set default font size

  call qp_set_text_attrib ("TEXT",  height = s%plot_page%text_height) 
  call qp_set_text_attrib ("MAIN_TITLE",  height = s%plot_page%text_height) 
  call qp_set_text_attrib ("GRAPH_TITLE",  height = s%plot_page%text_height) 
  call qp_set_text_attrib ("AXIS_LABEL",  height = s%plot_page%text_height) 
  call qp_set_text_attrib ("AXIS_NUMBERS",  height = s%plot_page%text_height) 
  call qp_set_text_attrib ("LEGEND",  height = s%plot_page%text_height) 

end subroutine tao_create_plot_window

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
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
  call qp_end
  has_been_created = .false.

end subroutine tao_destroy_plot_window

end module tao_plot_window_mod
