!+
! Subroutine tao_create_plot_window ()
!
! Subroutine to create the plot window.
! This soubroutine knows not to create a second window if one already exists.
!-

subroutine tao_create_plot_window ()

use tao_mod
use quick_plot

implicit none

logical :: has_been_created = .false.
character(24) :: r_name = 'tao_create_plot_window'

! Open a plot window

if (.not. s%global%plot_on) then
  call out_io (s_info$, r_name, 'global%plot_on = F so no plot window created.')
  return
endif

if (has_been_created) return

has_been_created = .true.

call qp_open_page ('X', s%plot_page%id_window, s%plot_page%size(1), s%plot_page%size(2), 'POINTS')
call qp_set_layout (page_border = s%plot_page%border)

! Set default font size

call qp_set_text_attrib ("TEXT",  height = s%plot_page%text_height) 
call qp_set_text_attrib ("MAIN_TITLE",  height = s%plot_page%text_height) 
call qp_set_text_attrib ("GRAPH_TITLE",  height = s%plot_page%text_height) 
call qp_set_text_attrib ("AXIS_LABEL",  height = s%plot_page%text_height) 
call qp_set_text_attrib ("AXIS_NUMBERS",  height = s%plot_page%text_height) 
call qp_set_text_attrib ("LEGEND",  height = s%plot_page%text_height) 

end subroutine
