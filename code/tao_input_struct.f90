!+
! module tao_input_struct
!
! Module to define the structures needed for the namelist input.
!-

module tao_input_struct

use tao_mod

!-------------------------------------------------------------
! data input structures

type tao_d2_data_input
  character(16) name           ! name of data
end type

type tao_d1_data_input
  character(16) name           ! type of data
end type

type tao_data_input
  character(16) :: name
  character(16) :: data_type
  character(16) :: ele_name
  character(16) :: ele2_name
  character(16) :: merit_type
  real(rp) :: meas_value
  real(rp) :: weight
  logical :: good_data
end type

!-------------------------------------------------------------
! variable input structures

type tao_v1_var_input
  character(16) name           ! name of variable
end type

type tao_var_input
  character(16) :: name
  character(16) :: ele_name
  character(16) attribute       ! attribute to vary
  character(16) universe
  real(rp) :: weight
  real(rp) :: step
  real(rp) low_lim
  real(rp) high_lim
  character(16) :: merit_type
end type

!-------------------------------------------------------------
! plot input structures

type tao_place_input
  character(16) region
  character(16) plot
end type

type tao_plot_page_input
  real(rp) size(2)
  real(rp) text_height
  type (qp_rect_struct) border
end type

type tao_curve_input
  character(16) data_source
  character(16) data_type
  real(rp) units_factor
  integer symbol_every
  integer ix_universe
  logical draw_line
  logical use_y2
  type (qp_line_struct) line
  type (qp_symbol_struct) symbol
end type

type tao_graph_input
  character(16) name
  character(16) type
  character(80) title
  integer box(4)
  integer ix_universe
  integer n_curve
  type (qp_rect_struct) margin
  type (qp_axis_struct) y
  type (qp_axis_struct) y2
end type 

type tao_plot_input
  character(16) name
  type (tao_plot_who_struct) who(10)
  type (qp_axis_struct) x
  character(16) x_axis_type
  logical convert  
  integer n_graph
  logical independent_graphs
end type

!-------------------------------------------------------------
! other structures

type tao_design_lat_input
  character(200) file
  character(16) :: parser = 'bmad'
end type

type tao_key_input
  character(16) ele_name
  character(16) attrib_name
  real(rp) delta
  character(16) universe
  real(rp) small_step
  real(rp) low_lim
  real(rp) high_lim
  real(rp) weight
  logical good_opt
end type

end module
