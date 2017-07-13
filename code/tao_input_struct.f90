!+
! module tao_input_struct
!
! Module to define the structures needed for the namelist input.
!-

module tao_input_struct

use tao_struct
use bmad_struct
use bmad_interface

integer, private, save :: dummy = 0 ! So ranlib will not complain about no symbols

!-------------------------------------------------------------
! data input structures

type tao_d2_data_input
  character(40) :: name = ''      ! name of data
end type

type tao_d1_data_input
  character(40) :: name = ''      ! type of data
end type

type tao_datum_input
  character(200) :: data_type = ''
  character(40) :: ele_ref_name = ''
  character(40) :: ele_start_name = ''
  character(40) :: ele_name = ''
  character(20) :: merit_type = ''
  real(rp) :: meas = real_garbage$  ! used to tag when %meas_value is set in file
  real(rp) :: weight = 0
  logical :: good_user = .true.
  character(20) :: data_source = ''
  character(20) :: eval_point = 'end'
  real(rp) :: s_offset = 0
  integer :: ix_bunch = 0
  real(rp) :: invalid_value = 0
end type

!-------------------------------------------------------------
! variable input structures

type tao_v1_var_input
  character(40) :: name           ! name of variable
end type

type tao_var_input
  character(40) :: ele_name = ''
  character(40) :: attribute = ''  ! attribute to vary
  character(16) :: universe = ''
  real(rp) :: weight = 0
  real(rp) :: step = 0
  real(rp) :: low_lim = -1e30
  real(rp) :: high_lim = 1e30
  character(40) :: merit_type = ''
  logical :: good_user = .true.
  logical :: key_bound = .false.
  real(rp) :: key_delta = 0
end type

!-------------------------------------------------------------
! plot input structures

type tao_region_input
  character(40) :: name = ''        ! Eg: 'top', 'bottom'.
  real(rp) :: location(4) = 0       ! location on page.
end type

type tao_place_input
  character(40) :: region = ''
  character(40) :: plot = ''
end type

type tao_curve_input
  character(40) :: name = ''
  character(40) :: data_source = 'lat'
  character(100) :: data_type_x = ''
  character(100) :: data_type_z = ''
  character(100) :: data_type = ''
  character(100) :: data_index = ''
  character(40) :: legend_text = ''
  character(40) :: units = ''
  character(60) :: component = ''
  real(rp) :: y_axis_scale_factor = 1
  real(rp) :: z_color0 = 0, z_color1 = 0
  integer :: symbol_every = 1
  integer :: ix_universe = -1
  logical :: draw_line = .true.
  logical :: draw_symbols = .true.
  logical :: draw_symbol_index = .false.
  logical :: use_y2 = .false.
  logical :: use_z_color = .false.
  logical :: autoscale_z_color = .true.
  logical :: draw_interpolated_curve = .true.
  logical :: smooth_line_calc = .true.
  character(40) :: ele_ref_name = ''
  integer :: ix_branch = 0
  integer :: ix_ele_ref = -1
  integer :: ix_bunch = 0
  type (qp_line_struct) :: line = qp_line_struct()
  type (qp_symbol_struct) :: symbol = qp_symbol_struct()
  type (tao_histogram_struct) :: hist = tao_histogram_struct()
end type

type tao_graph_input
  character(40) :: name = ''
  character(40) :: type = 'data' 
  character(80) :: title = ''
  character(60) :: component = ''
  character(2) :: floor_plan_view = 'zx'
  character(16) :: floor_plan_orbit_color = 'RED'
  integer :: box(4) = [1, 1, 1, 1]
  integer :: ix_universe = -1
  integer :: ix_branch = 0
  integer :: n_curve = 0
  real(rp) :: x_axis_scale_factor = 1
  real(rp) :: symbol_size_scale = 0
  real(rp) :: floor_plan_rotation = 0    ! Rotation of floor plan plot: 1.0 -> 360^deg 
  real(rp) :: floor_plan_orbit_scale = 0 ! Scale factor for drawing orbits. 0 -> Do not draw.
  logical :: clip = .true.
  logical :: floor_plan_size_is_absolute = .false.
  logical :: floor_plan_draw_only_first_pass = .false.
  logical :: correct_xy_distortion = .false.
  logical :: draw_axes = .true.
  logical :: draw_grid = .true.
  logical :: draw_curve_legend = .true.
  logical :: draw_only_good_user_data_or_vars = .true.
  type (qp_point_struct) :: text_legend_origin = qp_point_struct(5.0_rp, 0.0_rp, 'POINTS/GRAPH/RT')
  type (qp_point_struct) :: curve_legend_origin = qp_point_struct(5.0_rp, -2.0_rp, 'POINTS/GRAPH/LT')
  type (qp_rect_struct) :: margin = qp_rect_struct(0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, '%GRAPH')
  type (qp_rect_struct) :: scale_margin = qp_rect_struct(0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, '%GRAPH')
  type (qp_axis_struct) :: x = qp_axis_struct()
  type (qp_axis_struct) :: y = qp_axis_struct()
  type (qp_axis_struct) :: y2 = qp_axis_struct()
end type 

type tao_plot_input
  character(40) :: name
  character(60) :: description
  character(16) :: x_axis_type
  integer :: n_graph
  logical :: autoscale_gang_x      ! scale cmd scales graphs independently?
  logical :: autoscale_gang_y      ! scale cmd scales graphs independently?
  logical :: autoscale_x 
  logical :: autoscale_y 
  type (qp_axis_struct) :: x
end type

!-------------------------------------------------------------
! other structures

type tao_design_lat_input
  character(200) :: file = ''
  character(200) :: file2 = ''
  character(16) :: language = ''
  character(40) :: use_line = ''
  logical :: one_turn_map_calc = .false.
  logical :: dynamic_aperture_calc = .false.
end type

type tao_key_input
  character(40) :: ele_name
  character(40) :: attrib_name
  real(rp) :: delta
  character(16) :: universe
  real(rp) :: small_step
  real(rp) :: low_lim
  real(rp) :: high_lim
  real(rp) :: weight
  logical :: good_opt
  character(40) :: merit_type
end type

type tao_plot_page_input
  character(8) :: plot_display_type = ''  
  character(80) :: ps_scale             ! scaling when creating PS files.
  real(rp) size(2)                      ! width and height of window in pixels.
  real(rp) :: text_height = 12              ! In points. Scales the height of all text
  real(rp) :: main_title_text_scale  = 1.3  ! Relative to text_height
  real(rp) :: graph_title_text_scale = 1.1  ! Relative to text_height
  real(rp) :: axis_number_text_scale = 0.9  ! Relative to text_height
  real(rp) :: axis_label_text_scale  = 1.0  ! Relative to text_height
  real(rp) :: legend_text_scale      = 0.7  ! Relative to text_height
  real(rp) :: key_table_text_scale   = 0.9  ! Relative to text_height
  real(rp) :: floor_plan_shape_scale = 1.0
  real(rp) :: lat_layout_shape_scale = 1.0
  real(rp) :: curve_legend_line_len  = 50   ! Points
  real(rp) :: curve_legend_text_offset = 10 ! Points
  integer :: n_curve_pts = 401           ! Number of points for plotting a smooth curve
  type (tao_title_struct) :: title(2)       ! Titles at top of page.
  type (qp_rect_struct) :: border           ! Border around plots edge of page.
end type

type tao_dynamic_aperture_input
  real(rp) :: min_angle = 0
  real(rp) :: max_angle = pi
  integer :: n_angle   = 9
  integer :: n_turn = 100         ! Number of turns a particle must survive
  real(rp) :: x_init = 1e-3_rp    ! Initial x coordinate to start with for theta_xy = 0.
  real(rp) :: y_init = 1e-3_rp    ! Initial y coordinate to start with for theta_xy = pi/2.
  real(rp) :: accuracy = 1e-5_rp  ! Resolution of bracketed aperture
  real(rp) :: pz(100) = real_garbage$  
end type

contains

!------------------------------------------------------------------------------

subroutine tao_set_plotting (plot_page, plot_input, use_cmd_line_geom, reverse)

implicit none

type (tao_plot_page_input) plot_page
type (tao_plot_page_struct) plot_input

integer ix
logical use_cmd_line_geom
logical, optional :: reverse
character(40) str
character(16), parameter :: r_name = 'tao_set_plotting'

! 

if (logic_option(.false., reverse)) then
  if (plot_input%plot_display_type /= '') plot_page%plot_display_type = plot_input%plot_display_type
  plot_page%ps_scale                     = plot_input%ps_scale
  plot_page%size                         = plot_input%size
  plot_page%text_height                  = plot_input%text_height
  plot_page%main_title_text_scale        = plot_input%main_title_text_scale
  plot_page%graph_title_text_scale       = plot_input%graph_title_text_scale
  plot_page%axis_number_text_scale       = plot_input%axis_number_text_scale
  plot_page%axis_label_text_scale        = plot_input%axis_label_text_scale
  plot_page%floor_plan_shape_scale       = plot_input%floor_plan_shape_scale
  plot_page%lat_layout_shape_scale       = plot_input%lat_layout_shape_scale
  plot_page%legend_text_scale            = plot_input%legend_text_scale
  plot_page%key_table_text_scale         = plot_input%key_table_text_scale
  plot_page%curve_legend_line_len        = plot_input%curve_legend_line_len
  plot_page%curve_legend_text_offset     = plot_input%curve_legend_text_offset
  plot_page%n_curve_pts                  = plot_input%n_curve_pts
  plot_page%title                        = plot_input%title 
  plot_page%border                       = plot_input%border
endif

! 

if (plot_page%plot_display_type /= '') plot_input%plot_display_type = plot_page%plot_display_type
plot_input%ps_scale                     = plot_page%ps_scale
plot_input%size                         = plot_page%size
plot_input%text_height                  = plot_page%text_height
plot_input%main_title_text_scale        = plot_page%main_title_text_scale
plot_input%graph_title_text_scale       = plot_page%graph_title_text_scale
plot_input%axis_number_text_scale       = plot_page%axis_number_text_scale
plot_input%axis_label_text_scale        = plot_page%axis_label_text_scale
plot_input%floor_plan_shape_scale       = plot_page%floor_plan_shape_scale
plot_input%lat_layout_shape_scale       = plot_page%lat_layout_shape_scale
plot_input%legend_text_scale            = plot_page%legend_text_scale
plot_input%key_table_text_scale         = plot_page%key_table_text_scale
plot_input%curve_legend_line_len        = plot_page%curve_legend_line_len
plot_input%curve_legend_text_offset     = plot_page%curve_legend_text_offset
plot_input%n_curve_pts                  = plot_page%n_curve_pts
plot_input%title                        = plot_page%title 
plot_input%border                       = plot_page%border

! Plot window geometry specified on cmd line?

if (use_cmd_line_geom .and. s%com%plot_geometry /= '') then
   str = s%com%plot_geometry
   ix = index(str, 'x')
   if (ix == 0) then
     call out_io (s_error$, r_name, 'Malformed -geometry argument. No "x" present: ' // str)
   else
     if (.not. is_integer(str(1:ix-1)) .or. .not. is_integer(str(ix+1:))) then
       call out_io (s_error$, r_name, 'Malformed -geometry argument: ' // str)
     else
       read (str(:ix-1), *) plot_input%size(1)
       read (str(ix+1:), *) plot_input%size(2)
     endif
   endif
 endif
 
end subroutine tao_set_plotting

!------------------------------------------------------------------------------

end module
