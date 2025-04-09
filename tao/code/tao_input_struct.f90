!+
! module tao_input_struct
!
! Module to define the structures needed for the namelist input.
!-

module tao_input_struct

use tao_struct
use bmad_interface

integer, private, save :: dummy = 0 ! So ranlib will not complain about no symbols

!-------------------------------------------------------------
! data input structures

type tao_d2_data_input
  character(40) :: name = ''         ! name of data
end type

type tao_d1_data_input
  character(40) :: name = ''         ! type of data
end type

type tao_datum_input
  character(600) :: data_type = ''   ! long due to expressions.
  character(40) :: ele_ref_name = ''
  character(40) :: ele_start_name = ''
  character(40) :: ele_name = ''
  character(20) :: merit_type = ''
  real(rp) :: meas = real_garbage$   ! used to tag when %meas_value is set in file
  real(rp) :: weight = real_garbage$
  logical :: good_user = .true.
  logical :: good_opt = .true.
  character(20) :: data_source = ''
  character(20) :: eval_point = 'end'
  real(rp) :: s_offset = 0
  integer :: ix_bunch = 0
  type (spin_axis_struct) :: spin_axis = spin_axis_struct()
  real(rp) :: invalid_value = 0
  real(rp) :: error_rms = 0
end type

!-------------------------------------------------------------
! variable input structures

type tao_v1_var_input
  character(40) :: name            ! name of variable
end type

type tao_var_input
  character(40) :: ele_name = ''
  character(40) :: attribute = ''  ! attribute to vary
  character(16) :: universe = ''
  real(rp) :: weight = real_garbage$
  real(rp) :: step = 0
  real(rp) :: low_lim = -1e30
  real(rp) :: high_lim = 1e30
  character(40) :: merit_type = ''
  logical(4) :: good_user = .true.
  logical(4) :: key_bound = .false.
  real(rp) :: key_delta = 0
  real(rp) :: meas = real_garbage$
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
  character(40) :: data_source = ''
  character(100) :: data_type_x = ''
  character(100) :: data_type_z = ''                       ! Deprecated. Use c%z_color%...
  character(200) :: data_type = ''
  character(100) :: data_index = ''
  character(40) :: legend_text = ''
  character(40) :: units = ''                              ! Unused.
  character(60) :: component = ''
  real(rp) :: y_axis_scale_factor = 1
  real(rp) :: z_color0 = invalid$, z_color1 = invalid$     ! Deprecated. Use c%z_color%...
  integer :: symbol_every = 1
  integer :: ix_universe = -1
  integer :: n_turn = -1
  logical :: draw_line = .true.
  logical :: draw_symbols = .true.
  logical :: draw_symbol_index = .false.
  logical :: draw_error_bars = .false.
  logical :: use_y2 = .false.
  logical :: use_z_color = .false.                ! Deprecated. Use c%z_color%...
  logical :: autoscale_z_color = .true.           ! Deprecated. Use c%z_color%...
  logical :: smooth_line_calc = .true.
  character(40) :: ele_ref_name = ''
  integer :: ix_branch = 0
  integer :: ix_bunch = 0
  type (qp_line_struct) :: line = qp_line_struct()
  type (qp_symbol_struct) :: symbol = qp_symbol_struct()
  type (tao_histogram_struct) :: hist = tao_histogram_struct()
  type (tao_curve_orbit_struct) :: orbit = tao_curve_orbit_struct()
  type (tao_curve_color_struct) :: z_color = tao_curve_color_struct()
end type

type tao_graph_input
  character(40) :: name = ''
  character(40) :: type = 'data' 
  character(80) :: title = ''
  character(60) :: component = ''
  character(100) :: text_legend(10) = ''
  character(2) :: floor_plan_view = ''                 ! deprecated. Use g%floor_plan%...
  character(16) :: floor_plan_orbit_color = ''         ! deprecated. Use g%floor_plan%...
  integer :: box(4) = [1, 1, 1, 1]
  integer :: ix_universe = -1
  integer :: ix_branch = 0
  integer :: n_curve = -1
  real(rp) :: x_axis_scale_factor = 1
  real(rp) :: symbol_size_scale = 0
  real(rp) :: floor_plan_rotation = real_garbage$      ! deprecated. Use g%floor_plan%...
  real(rp) :: floor_plan_orbit_scale = -1              ! deprecated. Use g%floor_plan%...
  logical :: floor_plan_flip_label_side = .false.      ! deprecated. Use g%floor_plan%...
  logical :: floor_plan_size_is_absolute = .false.     ! deprecated. Use g%floor_plan%...
  logical :: floor_plan_draw_only_first_pass = .false. ! deprecated. Use g%floor_plan%...
  logical :: correct_xy_distortion = .true.            ! deprecated. Use g%floor_plan%...
  logical :: clip = .true.
  logical :: draw_title = .true.
  logical :: draw_axes = .true.
  logical :: draw_grid = .true.
  logical :: draw_curve_legend = .true.
  logical :: draw_only_good_user_data_or_vars = .true.
  logical :: allow_wrap_around = .true.
  type (tao_floor_plan_struct) :: floor_plan = tao_floor_plan_struct()
  type (qp_point_struct) :: text_legend_origin = qp_point_struct(5.0_rp, 0.0_rp, 'POINTS/GRAPH/RT')
  type (qp_point_struct) :: curve_legend_origin = qp_point_struct(5.0_rp, -2.0_rp, 'POINTS/GRAPH/LT')
  type (qp_rect_struct) :: margin = qp_rect_struct(0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, '%GRAPH')
  type (qp_rect_struct) :: scale_margin = qp_rect_struct(0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, '%GRAPH')
  type (qp_axis_struct) :: x = qp_axis_struct()
  type (qp_axis_struct) :: y = qp_axis_struct()
  type (qp_axis_struct) :: x2 = qp_axis_struct()
  type (qp_axis_struct) :: y2 = qp_axis_struct()
  type (qp_legend_struct) :: curve_legend = qp_legend_struct(1.0_rp, 30.0_rp, 6.0_rp, .true., .true., .true.)
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
  integer :: n_curve_pts 
  type (qp_axis_struct) :: x
end type

!-------------------------------------------------------------
! other structures

type tao_design_lat_input
  character(400) :: file = ''
  character(400) :: file2 = ''
  character(16) :: language = ''
  character(40) :: use_line = ''
  logical :: one_turn_map_calc = .false.
  logical :: dynamic_aperture_calc = .false.
  logical :: reverse_lattice = .false.
  character(40) :: start_branch_at = ''
  character(80) :: slice_lattice = ''
  character(40) :: use_element_range(2) = ''
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
  type (tao_title_struct) :: title = tao_title_struct()         ! Title  at top of page.
  type (tao_title_struct) :: subtitle = tao_title_struct()      ! Subtitle at top of page.
  type (qp_rect_struct) :: border = qp_rect_struct()            ! Border around plots edge of page.
  character(8) :: plot_display_type = ''  
  real(rp) :: size(2) = 0.0                 ! width and height of window in pixels.
  real(rp) :: text_height = 12              ! In points. Scales the height of all text
  real(rp) :: main_title_text_scale  = 1.3  ! Relative to text_height
  real(rp) :: graph_title_text_scale = 1.1  ! Relative to text_height
  real(rp) :: axis_number_text_scale = 0.9  ! Relative to text_height
  real(rp) :: axis_label_text_scale  = 1.0  ! Relative to text_height
  real(rp) :: legend_text_scale      = 0.9  ! Relative to text_height
  real(rp) :: key_table_text_scale   = 0.9  ! Relative to text_height
  real(rp) :: floor_plan_shape_scale = 1.0
  real(rp) :: floor_plan_text_scale  = 1.0  ! Scale used = floor_plan_text_scale * legend_text_scale
  real(rp) :: lat_layout_shape_scale = 1.0
  real(rp) :: lat_layout_text_scale  = 1.0  ! Scale used = lat_layout_text_scale * legend_text_scale
  real(rp) :: curve_legend_line_len  = real_garbage$      ! OLD STYLE. Points.
  real(rp) :: curve_legend_text_offset = real_garbage$    ! OLD STYLE. Points.
  integer :: n_curve_pts = n_curve_pts_init$   ! Number of points for plotting a smooth curve
  logical :: delete_overlapping_plots = .true. ! Delete overlapping plots when a plot is placed?
  logical :: draw_graph_title_suffix = .true.
end type

type tao_ele_shape_input
  character(60) :: ele_id = ''       ! element "key::name" to match to.
  character(40) :: shape = ''        ! Shape to draw
  character(16) :: color = 'black'   ! Color of shape
  real(rp) :: size = 0               ! plot vertical height 
  character(16) :: label = 'name'    ! Can be: 'name', 's', 'none' 
  logical :: draw = .true.           ! Draw the shape?
  logical :: multi = .false.         ! Can be part of a multi-shape.
  integer :: line_width = 1          ! Width of lines used to draw the shape.
  real(rp) :: offset = 0             ! Vertical offset. 
end type

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine tao_set_plotting (plot_input, plot_page, use_cmd_line_geom, reverse)

implicit none

type (tao_plot_page_input) plot_input
type (tao_plot_page_struct) plot_page

integer ix
logical use_cmd_line_geom
logical, optional :: reverse
character(40) str
character(16), parameter :: r_name = 'tao_set_plotting'

! 

if (logic_option(.false., reverse)) then
  if (plot_page%plot_display_type /= '') plot_input%plot_display_type = plot_page%plot_display_type
  plot_input%size                         = plot_page%size
  plot_input%text_height                  = plot_page%text_height
  plot_input%main_title_text_scale        = plot_page%main_title_text_scale
  plot_input%graph_title_text_scale       = plot_page%graph_title_text_scale
  plot_input%axis_number_text_scale       = plot_page%axis_number_text_scale
  plot_input%axis_label_text_scale        = plot_page%axis_label_text_scale
  plot_input%floor_plan_shape_scale       = plot_page%floor_plan_shape_scale
  plot_input%floor_plan_text_scale        = plot_page%floor_plan_text_scale
  plot_input%lat_layout_shape_scale       = plot_page%lat_layout_shape_scale
  plot_input%lat_layout_text_scale        = plot_page%lat_layout_text_scale
  plot_input%legend_text_scale            = plot_page%legend_text_scale
  plot_input%key_table_text_scale         = plot_page%key_table_text_scale
  plot_input%n_curve_pts                  = plot_page%n_curve_pts
  plot_input%title                        = plot_page%title 
  plot_input%subtitle                     = plot_page%subtitle 
  plot_input%border                       = plot_page%border
  plot_input%delete_overlapping_plots     = plot_page%delete_overlapping_plots
  plot_input%draw_graph_title_suffix      = plot_page%draw_graph_title_suffix
endif

! 

if (plot_input%plot_display_type /= '') plot_page%plot_display_type = plot_input%plot_display_type
plot_page%size                         = plot_input%size
plot_page%text_height                  = plot_input%text_height
plot_page%main_title_text_scale        = plot_input%main_title_text_scale
plot_page%graph_title_text_scale       = plot_input%graph_title_text_scale
plot_page%axis_number_text_scale       = plot_input%axis_number_text_scale
plot_page%axis_label_text_scale        = plot_input%axis_label_text_scale
plot_page%floor_plan_shape_scale       = plot_input%floor_plan_shape_scale
plot_page%floor_plan_text_scale        = plot_input%floor_plan_text_scale
plot_page%lat_layout_shape_scale       = plot_input%lat_layout_shape_scale
plot_page%lat_layout_text_scale        = plot_input%lat_layout_text_scale
plot_page%legend_text_scale            = plot_input%legend_text_scale
plot_page%key_table_text_scale         = plot_input%key_table_text_scale
plot_page%n_curve_pts                  = plot_input%n_curve_pts
plot_page%title                        = plot_input%title 
plot_page%subtitle                     = plot_input%subtitle 
plot_page%border                       = plot_input%border
plot_page%delete_overlapping_plots     = plot_input%delete_overlapping_plots
plot_page%draw_graph_title_suffix      = plot_input%draw_graph_title_suffix

if (plot_input%curve_legend_line_len /= real_garbage$) then
  call out_io(s_error$, r_name, '"plot_page%curve_legend_line_len" has been replaced by "default_graph%curve_legend%line_len".', &
                                 'Please modify the input file. This setting will be ignored.')
endif

if (plot_input%curve_legend_text_offset /= real_garbage$) then
  call out_io(s_error$, r_name, '"plot_page%curve_legend_text_offset" has been replaced by "default_graph%curve_legend%text_offset".', &
                                 'Please modify the input file. This setting will be ignored.')
endif

! Plot window geometry specified on cmd line?

if (use_cmd_line_geom .and. s%init%geometry_arg /= '') then
   str = s%init%geometry_arg
   ix = index(str, 'x')
   if (ix == 0) then
     call out_io (s_error$, r_name, 'Malformed -geometry argument. No "x" present: ' // str)
   else
     if (.not. is_integer(str(1:ix-1)) .or. .not. is_integer(str(ix+1:))) then
       call out_io (s_error$, r_name, 'Malformed -geometry argument: ' // str)
     else
       read (str(:ix-1), *) plot_page%size(1)
       read (str(ix+1:), *) plot_page%size(2)
     endif
   endif
 endif
 
end subroutine tao_set_plotting

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

elemental function tao_ele_shape_struct_to_input (shape_struct) result (shape_input)

type (tao_ele_shape_struct), intent(in) :: shape_struct
type (tao_ele_shape_input) shape_input

!

shape_input%ele_id     = shape_struct%ele_id
shape_input%shape      = shape_struct%shape
shape_input%color      = shape_struct%color
shape_input%size       = shape_struct%size
shape_input%label      = shape_struct%label
shape_input%draw       = shape_struct%draw
shape_input%multi      = shape_struct%multi
shape_input%line_width = shape_struct%line_width
shape_input%offset     = shape_struct%offset

end function tao_ele_shape_struct_to_input

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

function tao_ele_shape_input_to_struct (shape_input, namelist_name) result (shape_struct)

type (tao_ele_shape_input), target :: shape_input
type (tao_ele_shape_struct), target :: shape_struct
type (tao_ele_shape_input), pointer :: s_in
type (tao_ele_shape_struct), pointer :: s_st

integer ix

character(*), optional :: namelist_name
character(40) shape, prefix
character(100) err_line
character(*), parameter :: r_name = 'tao_ele_shape_input_to_struct'

!

if (present(namelist_name)) then
  err_line = '  IN NAMELIST: ' // namelist_name
else
  err_line = ''
endif

! Transfer from shape_input to shape_struct

s_in => shape_input
s_st => shape_struct

s_st%ele_id     = downcase(s_in%ele_id)
s_st%shape      = downcase(s_in%shape)
s_st%color      = downcase(s_in%color)
s_st%label      = downcase(s_in%label)
s_st%size       = s_in%size
s_st%draw       = s_in%draw
s_st%multi      = s_in%multi
s_st%line_width = s_in%line_width
s_st%offset     = s_in%offset

if (s_st%ele_id == 'wall::building') s_st%ele_id = 'building_wall::*'    ! Convert old style to new

! Instances where we must preserve ele_id case.

if (s_st%ele_id(1:6) == 'data::' .or. s_st%ele_id(1:5) == 'var::' .or. &
    s_st%ele_id(1:5) == 'lat::' .or. s_st%ele_id(1:15) == 'building_wall::') then
  ix = index(s_st%ele_id, '::')
  s_st%ele_id = s_st%ele_id(1:ix+1) // s_in%ele_id(ix+2:)
endif

! Convert from old shape names to new names

if (s_st%ele_id(1:15) == 'building_wall::' .and. s_st%shape == '-') s_st%shape = 'solid_line' ! Convert old style to new

if (s_st%label == '') s_st%label = 'name'
if (index('false', trim(s_st%label)) == 1) s_st%label = 'none'
if (index('true', trim(s_st%label)) == 1) s_st%label = 'name'

if (s_st%shape(1:4) == 'var_')            then;   s_st%shape = s_st%shape(1:3) // ':' // s_st%shape(5:)
elseif (s_st%shape(1:5) == 'vvar_')       then;   s_st%shape = s_st%shape(1:4) // ':' // s_st%shape(6:)
elseif (s_st%shape(1:9) == 'asym_var_')   then;   s_st%shape = s_st%shape(1:8) // ':' // s_st%shape(10:)
elseif (s_st%shape(1:10) == 'asym_vvar_') then;   s_st%shape = s_st%shape(1:9) // ':' // s_st%shape(11:)
elseif (s_st%shape == '-')                then;   s_st%shape = 'SOLID_LINE'
endif

! Error checks

ix = index(s_st%shape, ':')
if (ix == 0) then
  prefix = ''
  shape = s_st%shape
else
  prefix = s_st%shape(1:ix-1)
  shape = s_st%shape(ix+1:)
endif

if (prefix == 'pattern') then
  if (all(shape /= s%plot_page%pattern(:)%name)) then
    call out_io (s_error$, r_name, 'ELE_SHAPE: ', trim(s_st%shape) // ' DOES NOT HAVE AN ASSOCIATED PATTERN.', err_line)
                                   
  endif
  return
endif

select case (prefix)
case ('vvar', 'asym_vvar')
  if (s_st%ele_id(1:6) /= 'data::' .and. s_st%ele_id(1:5) /= 'var::') then
    call out_io (s_error$, r_name, 'ELE_SHAPE WITH ', trim(s_st%shape) // &
                                   '  MUST BE ASSOCIATED WITH A DATUM OR VARIABLE! NOT: ' // s_st%ele_id, err_line)
  endif
case ('', 'var', 'asym_var')
case default
  call out_io (s_fatal$, r_name, 'UNKNOWN ELE_SHAPE PREFIX: ' // s_st%shape, err_line)
end select

select case (shape)
case ('box', 'xbox', 'bow_tie', 'rbow_tie', 'circle', 'diamond', 'x', 'r_triangle', 'l_triangle', 'u_triangle', 'd_triangle')
case ('solid_line', 'dashed_line', 'dash_dot_line', 'dotted_line')
  if (s_st%ele_id(1:15) /= 'building_wall::') then
    call out_io (s_fatal$, r_name, 'SHAPE "' // trim(shape) // '" MAY ONLY BE USED FOR building_wall ELEMENTS.', err_line)
  endif
case default
  call out_io (s_fatal$, r_name, 'UNKNOWN SHAPE: ' // s_st%shape, err_line)
end select

end function tao_ele_shape_input_to_struct 

end module
