!+
! module tao_input_struct
!
! Module to define the structures needed for the namelist input.
!-

module tao_input_struct

use tao_struct
use bmad_struct
use bmad_interface

!-------------------------------------------------------------
! data input structures

type tao_d2_data_input
  character(40) name           ! name of data
end type

type tao_d1_data_input
  character(40) name           ! type of data
end type

type tao_data_input
  character(100) :: data_type
  character(40) :: ele0_name
  character(40) :: ele_name
  character(20) :: merit_type
  real(rp) :: meas
  real(rp) :: weight
  logical :: good_user
  character(20) data_source
  integer ix_bunch
  character(8) relative
  character(40) :: name
end type

!-------------------------------------------------------------
! variable input structures

type tao_v1_var_input
  character(40) name           ! name of variable
end type

type tao_var_input
  character(40) :: ele_name
  character(40) attribute       ! attribute to vary
  character(16) universe
  real(rp) :: weight
  real(rp) :: step
  real(rp) low_lim
  real(rp) high_lim
  character(40) :: merit_type
  character(8) :: good_user
  character(40) :: name
  character(8) key_bound
  real(rp) key_delta
end type

!-------------------------------------------------------------
! plot input structures

type tao_region_input
  character(40) name             ! Eg: 'top', 'bottom'.
  real(rp) location(4)           ! location on page.
end type

type tao_place_input
  character(40) region
  character(40) plot
end type

type tao_curve_input
  character(40) name
  character(40) data_source
  character(100) data_type_x
  character(100) data_type
  character(100) data_index
  real(rp) x_axis_scale_factor
  real(rp) y_axis_scale_factor
  integer symbol_every
  integer ix_universe
  logical draw_line
  logical draw_symbols
  logical draw_symbol_index
  logical use_y2
  logical draw_interpolated_curve
  logical smooth_line_calc
  type (qp_line_struct) line
  type (qp_symbol_struct) symbol
  character(40) ele_ref_name
  integer ix_ele_ref
  integer ix_bunch
end type

type tao_graph_input
  character(40) name
  character(40) type
  character(80) title
  character(60) component
  integer box(4)
  integer ix_universe
  integer n_curve
  type (qp_point_struct) legend_origin  ! For backwards compatibility
  type (qp_point_struct) text_legend_origin
  type (qp_point_struct) curve_legend_origin
  type (tao_data_var_component_struct) who(n_who_maxx)
  type (qp_rect_struct) margin
  type (qp_axis_struct) x
  type (qp_axis_struct) y
  type (qp_axis_struct) y2
  logical clip
  logical draw_axes
  logical correct_xy_distortion
  logical draw_curve_legend     ! For identifying curves. 
end type 

type tao_plot_input
  character(40) name
  type (qp_axis_struct) x
  character(16) x_axis_type
  integer n_graph
  logical independent_graphs
  logical autoscale_gang_x      ! scale cmd scales graphs independently?
  logical autoscale_gang_y      ! scale cmd scales graphs independently?
end type

!-------------------------------------------------------------
! other structures

type tao_design_lat_input
  character(100) file
  character(100) file2
  character(16) parser
end type

type tao_connected_uni_input
  integer from_universe
  character(40) at_element !connected at end of element
  integer at_ele_index ! connected at end of element
  real(rp) :: at_s ! connected at position s
  logical match_to_design
end type

type tao_key_input
  character(40) ele_name
  character(40) attrib_name
  real(rp) delta
  character(16) universe
  real(rp) small_step
  real(rp) low_lim
  real(rp) high_lim
  real(rp) weight
  logical good_opt
  character(40) merit_type
end type

end module
