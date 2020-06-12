!+
! Module tao_struct
!
! Module defining the basic tao structures
!
! If any pointers or allocatables are added remember to add a corresponding
! statment to tao_init\deallocate_everything.
!-

module tao_struct

use, intrinsic :: iso_c_binding
use equal_mod
use quick_plot_struct
use dynamic_aperture_mod
use complex_taylor_mod
use srdt_mod

integer, parameter :: model$ = 1, base$ = 2, design$ = 3
integer, parameter :: ix_common_uni$ = 0
integer, parameter :: apparent_emit$ = 1, projected_emit$ = 2

character(2), parameter :: tao_floor_plan_view_name(6) = [character(2):: 'xy', 'xz', 'yx', 'yz', 'zx', 'zy']
character(8), parameter :: tao_lat_type_name(3) = [character(8):: 'model ', 'base  ', 'design']
character(8), parameter :: tao_data_source_name(5) = [character(8):: 'lat', 'beam', 'data', 'var', 'multi_turn_orbit']
character(20), parameter :: tao_graph_type_name(6) = [character(20):: 'data', 'floor_plan', 'lat_layout', &
                                                       'phase_space', 'histogram', 'dynamic_aperture']
character(16), parameter :: tao_x_axis_type_name(10) = [character(16):: 'index', 'lat', 'var', &
                                   'ele_index', 's', 'none', 'phase_space', 'histogram', 'data', 'floor']
character(12), parameter :: tao_data_type_z_name(14) = [character(12):: 'x', 'px', 'y', 'py', 'z', 'pz', 'time', &
                  'intensity', 'intensity_x', 'intensity_y', 'phase_x', 'phase_y', 'Ja', 'energy']
character(8), parameter :: tao_var_merit_type_name(2) = [character(8):: 'target ', 'limit']
character(8), parameter :: tao_data_merit_type_name(7) = [character(8):: 'target', 'min', 'max', &
                                             'abs_min', 'abs_max', 'max-min', 'average']
character(12), parameter :: tao_optimizer_name(6) = [character(12):: 'de', 'lm', 'lmdif', 'custom', 'svd', 'geodesic_lm']
character(24), parameter :: tao_shape_shape_name(13) = [character(24):: 'Box', 'Var:Box', 'VVar:Box', 'Asym_Var:Box', &
                                    'Asym_VVar:Box', 'Xbox', 'Diamond', 'Bow_Tie', 'RBow_Tie', 'Circle', 'X', &
                                    'Pattern:<pattern-name>', '-']

character(8), parameter :: tao_shape_label_name(3) = [character(8):: 'name', 's', 'none']
character(24), parameter :: tao_wave_data_name(27) = [character(40):: 'orbit.x', 'orbit.y', 'beta.a', 'beta.b', &
    'eta.x', 'eta.y', 'phase.a', 'phase.b', 'cbar.12', 'cbar.11', 'cbar.22', &
    'ping_a.amp_x', 'ping_a.phase_x', 'ping_a.amp_y', 'ping_a.phase_y', &
    'ping_a.amp_sin_y', 'ping_a.amp_cos_y', 'ping_a.amp_sin_rel_y', 'ping_a.amp_cos_rel_y', &
    'ping_b.amp_y', 'ping_b.phase_y', 'ping_b.amp_x', 'ping_b.phase_x', &
    'ping_b.amp_sin_x', 'ping_b.amp_cos_x', 'ping_b.amp_sin_rel_x', 'ping_b.amp_cos_rel_x']

integer, parameter :: n_char_show = 1000

logical, save, target :: forever_true$ = .true.  ! Used for pointer init.

interface assignment (=)
  module procedure tao_lattice_equal_tao_lattice
  module procedure tao_lattice_branch_equal_tao_lattice_branch
  module procedure tao_lattice_branches_equal_tao_lattice_branches
end interface

!---------

type tao_cmd_history_struct  ! record the command history
  character(:), allocatable :: cmd     ! the command
  integer :: ix = 0      ! command index (1st command has ix = 1, etc.)
  logical cmd_file       ! Did command come from a command file
end type

!-----------------------------------------------------------------------
! A tao_real_pointer_struct is just a pointer to a real number.
! This is used to construct arrays of reals.

type tao_real_pointer_struct
  real(rp), pointer :: r => null()
  logical, pointer :: good_value => null()
  logical, pointer :: good_user => null()
end type

type tao_logical_array_struct
  logical, pointer :: l => null()
end type

type tao_integer_array_struct
  integer, pointer :: i => null()
end type

type tao_string_array_struct
  character(:), pointer :: s => null()
end type

! Note: Expressions may not have an associated lattice element or an associated longituinal position.

type tao_expression_info_struct
  logical :: good = .true.                    ! Expression is valid.
  type (ele_struct), pointer :: ele => null() ! Associated ele if it exists
  real(rp) :: s = real_garbage$               ! Longitudinal position of expression.
end type

type tao_eval_stack1_struct
  integer :: type = 0
  character(60) :: name = ''
  real(rp) :: scale = 1    ! Scale factor for ping data 
  real(rp), allocatable :: value(:)
  type (tao_expression_info_struct), allocatable :: info(:)
  type (tao_real_pointer_struct), allocatable :: value_ptr(:)  ! Used to point to data, lattice parameters, etc
end type

!----------------------------------------------------------------------

type tao_ele_pointer_struct
  type (ele_pointer_struct), allocatable :: eles(:)
  integer :: n_loc = 0
end type


type tao_ele_shape_struct    ! for the element layout plot
  character(60) :: ele_id = ''       ! element "key::name" to match to.
  character(40) :: shape = ''        ! Shape to draw
  character(16) :: color = 'black'   ! Color of shape
  real(rp) :: size = 0               ! plot vertical height 
  character(16) :: label = 'name'    ! Can be: 'name', 's', 'none' 
  logical :: draw = .true.           ! Draw the shape?
  logical :: multi = .false.         ! Can be part of a multi-shape.
  integer :: line_width = 1          ! Width of lines used to draw the shape.
  integer :: ix_key = 0              ! Extracted from ele_id. 0 => all classes (quadrupole, etc.)
  character(40) :: name_ele = ''     ! Name of element.
  type (tao_ele_pointer_struct), allocatable :: uni(:)
end type

type tao_drawing_struct
  type (tao_ele_shape_struct), allocatable :: ele_shape(:)
end type

type tao_shape_pattern_point_struct
  real(rp) :: s = real_garbage$, x = real_garbage$, radius = 0
end type

type tao_shape_pattern_struct
  character(40) :: name = ''
  type (qp_line_struct) :: line = qp_line_struct(1, 'Not_Set', 'solid')  ! Line color and pattern set by shape using this pattern.
  type (tao_shape_pattern_point_struct), allocatable :: pt(:)
end type

!-----------------------------------------------------------------------
! Plot structures.

type tao_title_struct
  character(100) :: string = ''      ! title character string.
  real(rp) :: x = 0.5, y = 0.97      ! x, y rwt lower left corner
  character(16) :: units = '%PAGE'   ! %BOX, POINTS, etc...
  character(2) :: justify = 'CC'     ! Left, Center, or Right justification.
  logical :: draw_it = .true.        ! draw the title?
end type

type tao_data_var_component_struct    ! Components to plot
  character(16) :: name = ''        ! Eg: 'meas', 'ref', 'model', etc.
  real(rp) :: sign = 1              ! +1 or -1
end type

type tao_histogram_struct
  logical :: density_normalized = .false.
  logical :: weight_by_charge = .true.
  real(rp) :: minimum = 0, maximum = 0
  real(rp) :: width = 0, center = 0
  integer :: number = 0
end type

! A curve is defined by a set of (x,y) points and the axis parameters.
! for example the horizontal orbit is one curve.

type tao_curve_struct
  character(40) :: name = ''             ! Name identifying the curve.
  character(40) :: data_source  = ''     ! 'lat', 'beam', 'data' (deprecated: 'dat'), 'var', 'multi_turn_orbit'
  character(100) :: data_index  = ''     ! Used for calculating %ix_symb(:).
  character(100) :: data_type_x = ''     ! Used for data slices and phase space plots.
  character(100) :: data_type_z = ''     ! Used for data phase space plots.
  character(200) :: data_type   = ''     ! 'orbit.x', etc.
  character(40) :: ele_ref_name = ''     ! Reference element.
  character(40) :: legend_text = ''      ! String to draw in a curve legend. 
  character(40) :: message_text = ''     ! Informational message to draw with graph.
  character(40) :: units = ''            ! Data units.
  character(60) :: component = ''        ! Who to plot. Eg: 'meas - design'
  character(80) :: why_invalid = '???'   ! Informative string to print.
  type (tao_graph_struct), pointer :: g  ! pointer to parent graph
  type (tao_histogram_struct) hist
  real(rp), allocatable :: x_line(:)     ! Coords for drawing a curve
  real(rp), allocatable :: y_line(:) 
  real(rp), allocatable :: y2_line(:)    ! Second array needed for beam chamber curve. 
  integer, allocatable :: ix_line(:)     ! Used by wave and aperture curves.
  real(rp), allocatable :: x_symb(:)     ! Coords for drawing the symbols
  real(rp), allocatable :: y_symb(:) 
  real(rp), allocatable :: z_symb(:)     ! Symbol color
  real(rp), allocatable :: err_symb(:)   ! Error bars
  real(rp), allocatable :: symb_size(:)  ! Symbol size. Used with symbol_size_scale. 
  integer, allocatable :: ix_symb(:)     ! Corresponding index in d1_data%d(:) array.
  real(rp) :: y_axis_scale_factor = 1    ! y-axis conversion from internal to plotting units.
  real(rp) :: s = 0                      ! longitudinal position
  real(rp) :: z_color0 = 0, z_color1 = 0 ! Min and max values for mapping z-axis to color.
  type (qp_line_struct) line             ! Line attributes
  type (qp_symbol_struct) symbol         ! Symbol attributes
  integer :: ix_universe = -1            ! Universe where data is. -1 => use s%com%default_universe
  integer :: symbol_every = 1            ! Symbol every how many points.
  integer :: ix_branch = 0
  integer :: ix_ele_ref = -1             ! Index in lattice of reference element.
  integer :: ix_ele_ref_track = -1       ! = ix_ele_ref except for super_lord elements.
  integer :: ix_bunch = 0                ! Bunch to plot.
  logical :: use_y2 = .false.            ! Use y2 axis?
  logical :: draw_line = .true.          ! Draw a line through the data points?
  logical :: draw_symbols = .true.       ! Draw a symbol at the data points?
  logical :: draw_symbol_index = .false. ! Draw the symbol index number curve%ix_symb?
  logical :: draw_error_bars = .false.   ! Draw error bars based upon data%error_rms if drawing data?
  logical :: smooth_line_calc = .true.   ! Calculate data between element edge points?
  logical :: use_z_color = .false.       ! For phase space plots.
  logical :: autoscale_z_color = .true.  ! Set %z_color0, %z_color1 automatically to the limits of %data_type_z
  logical :: valid = .false.             ! valid data? 
end type

! This is used with floor_plan drawings.

type tao_floor_plan_struct
  character(2) :: view = 'zx'                ! or 'xz'.
  real(rp) :: rotation = 0                   ! Rotation of floor plan plot: 1.0 -> 360^deg
  logical :: correct_distortion = .true.     ! T -> Shrink one axis so x-scale = y-scale.
  logical :: flip_label_side = .false.       ! Draw element label on other side of element?
  logical :: size_is_absolute = .false.      ! Are shape sizes in meters or window pixels?
  logical :: draw_only_first_pass = .false.  ! Draw only first pass with multipass elements?
  logical :: draw_building_wall = .true.     ! Draw the building wall?
  real(rp) :: orbit_scale = 0                ! Scale factor for drawing orbits. 0 -> Do not draw.
  character(16) :: orbit_color = 'red'
  character(16) :: orbit_pattern = 'solid'
  integer :: orbit_width = 1
end type

! A graph is a collection of overlayed curves with associated graph title, etc.
! For example a graph could contain just the horizontal orbit or could
! contain both overlayed horizontal and vertical orbits.

type tao_graph_struct
  character(40) :: name = ''          ! Name identifying the graph
  character(40) :: type = ''          ! 'data', 'lat_layout', 'phase_space', 'histogram', 'dynamic_aperture'
  character(100) :: title = ''
  character(100) :: title_suffix = ''
  character(100) :: text_legend(10) = ''            ! Array for holding descriptive info.
  character(100) :: text_legend_out(10) = ''        ! Array for holding descriptive info.
  character(60) :: component = ''                   ! Who to plot. Eg: 'meas - design'
  character(80) :: why_invalid = '???'              ! Informative string to print.
  type (tao_curve_struct), allocatable :: curve(:)
  type (tao_plot_struct), pointer :: p ! pointer to parent plot
  type (tao_floor_plan_struct) :: floor_plan = tao_floor_plan_struct()
  type (qp_point_struct) text_legend_origin
  type (qp_point_struct) curve_legend_origin
  type (qp_axis_struct) x                           ! X-axis parameters.
  type (qp_axis_struct) y                           ! Y-axis attributes.
  type (qp_axis_struct) x2                          ! X2-axis attributes (Only used for floor_plan).
  type (qp_axis_struct) y2                          ! Y2-axis attributes.
  type (qp_rect_struct) margin                      ! Margin around the graph.
  type (qp_rect_struct) scale_margin                ! Margin for scaling
  real(rp) :: x_axis_scale_factor = 1               ! x-axis conversion from internal to plotting units.
  real(rp) :: symbol_size_scale = 0                 ! Symbol size scale factor for phase_space plots.
  integer box(4)                                    ! Defines which box the plot is put in.
  integer :: ix_branch = 0                          ! Branch in lattice.
  integer :: ix_universe = -1                       ! Used for lat_layout plots.
  logical :: clip = .false.                         ! Clip plot at graph boundary.
  logical :: y2_mirrors_y = .true.                  ! Y2-axis same as Y-axis?
  logical :: limited = .false.                      ! True if at least one data point past graph bounds.
  logical :: draw_axes = .true.                     ! Draw axes, labels, etc?
  logical :: draw_curve_legend = .true.             ! Legend for displaying curve info.
  logical :: draw_grid = .true.                     ! Draw a grid?
  logical :: draw_title = .true.
  logical :: draw_only_good_user_data_or_vars = .true.
  logical :: allow_wrap_around = .true.             ! "Wrap" curves to extend past lattice boundaries?
  logical :: is_valid = .false.                     ! EG: Bad x_axis_type.
end type

! A plot is collection of graphs.
! For example a plot could contain three graphs. One for Cbar11, 
! One for Cbar12, and one for Cbar22.


type tao_plot_struct
  character(40) :: name = ' '                      ! Identifying name
  character(100) :: description = ''               ! Descriptive string.
  type (tao_graph_struct), allocatable :: graph(:) ! individual graphs of a plot
  type (qp_axis_struct) :: x = qp_axis_struct()    ! X-axis parameters.
  type (tao_plot_region_struct), pointer :: r => null() ! pointer to parent.
  integer :: ix_plot = -1                          ! Index in s%plot_page%template(:) or %region(:) arrays.
  integer :: n_curve_pts = -1                      ! Overrides s%plot_page%n_curve_pts.
  character(8) :: type = 'normal'                  ! or 'wave'
  character(16) :: x_axis_type = ''                ! 'index', 'ele_index', 's', 'none', 'floor', 'phase_space', etc.
  logical :: autoscale_x = .false.                 ! Horizontal autoscale.
  logical :: autoscale_y = .false.                 ! Vertical autoscale.
  logical :: autoscale_gang_x = .true.             ! scale cmd scales graphs together?
  logical :: autoscale_gang_y = .true.             ! scale cmd scales graphs together?
  logical :: list_with_show_plot_command = .true.  ! False used for default plots to shorten the output of "show plot"
  logical :: phantom = .false.                     ! True used to insert info into output of "show plot"
  logical :: default_plot = .false.                ! One of Tao's default plots? 
end type

! A region defines a plot and where to position the plot on the plot page
! %location = (x1, x2, y1, y2) gives the plotting region in percent of the 
!   part of the page inside the page_border with respect to the lower left corner.
! Eg: %location = (0.0, 1.0, 0.5, 1.0) gives the top half of the page inside the border.

type tao_plot_region_struct
  character(40) :: name = ''     ! Region name. Eg: 'r13', etc.
  type (tao_plot_struct) plot    ! Plot associated with this region
  real(rp) location(4)           ! [x1, x2, y1, y2] location on page.
  logical :: visible = .false.   ! To draw or not to draw.
  logical :: list_with_show_plot_command = .true.  ! False used for default plots to shorten the output of "show plot"
end type

! The tao_plot_page_struct defines the whole plotting window. 
! Note that the qp_com structure of quick_plot also is used to hold 
! plot page info.

type tao_plot_page_struct
  type (tao_title_struct) title(2)          ! Titles at top of page.
  type (qp_rect_struct) border              ! Border around plots edge of page.
  type (tao_drawing_struct) :: floor_plan
  type (tao_drawing_struct) :: lat_layout
  type (tao_shape_pattern_struct), allocatable :: pattern(:)
  type (tao_plot_struct), allocatable :: template(:)  ! Templates for the plots.
  type (tao_plot_region_struct), allocatable :: region(:)
  character(8) :: plot_display_type = 'X'   ! 'X' or 'TK'
  character(80) ps_scale                    ! scaling when creating PS files.
  real(rp) size(2)                          ! width and height of plot window in pixels.
  real(rp) :: text_height = 12              ! In points. Scales the height of all text
  real(rp) :: main_title_text_scale  = 1.3  ! Relative to text_height
  real(rp) :: graph_title_text_scale = 1.1  ! Relative to text_height
  real(rp) :: axis_number_text_scale = 0.9  ! Relative to text_height
  real(rp) :: axis_label_text_scale  = 1.0  ! Relative to text_height
  real(rp) :: legend_text_scale      = 0.9  ! Relative to text_height
  real(rp) :: key_table_text_scale   = 0.9  ! Relative to text_height
  real(rp) :: curve_legend_line_len  = 30   ! Points
  real(rp) :: curve_legend_text_offset = 6 ! Points
  real(rp) :: floor_plan_shape_scale = 1.0
  real(rp) :: lat_layout_shape_scale = 1.0
  integer :: n_curve_pts = 401              ! Default number of points for plotting a smooth curve.
  integer :: id_window = -1                 ! X window id number.
  logical :: delete_overlapping_plots = .true. ! Delete overlapping plots when a plot is placed?
end type

! Arrays of structures

type tao_plot_array_struct
  type (tao_plot_struct), pointer :: p
end type

type tao_graph_array_struct
  type (tao_graph_struct), pointer :: g
end type

type tao_curve_array_struct
  type (tao_curve_struct), pointer :: c
end type

!-----------------------------------------------------------------------
! The data_struct defines the fundamental data structure representing 
! one datum point.
! The universe_struct will hold an array of data_struct structures: u%data(:).
!
! %exists       -- The datum can exist. Non-existent datums can serve 
!                    as place holders in the u%data array.
! %good_model   -- A valid model value was computed. For example, good_model would 
!                    be False for some orbit data if the particle being tracked is lost.
! %good_base    -- A valid base value was computed. 
! %good_design  -- A valid design value was computed. 
! %good_meas    -- Set by the routine that reads in a data set. Good_meas may be 
!                    false, say, if a detector amplifyer is overloaded.
! %good_ref     -- Like good_meas this is set for a reference data set.
! %good_user    -- What the user has selected using the use, veto, and restore 
!                    commands.
! %good_opt     -- Not modified by Tao. Setting is reserved to be done by extension code.
! %good_plot    -- Not modified by Tao. Setting is reserved to be done by extension code.
! %useit_plot   -- Datum is valid for plotting:
!                    = %exists & %good_plot (w/o measured & reference data)
!                    = %exists & %good_plot & %good_user & %good_meas (w/ meas data)
!                    = %exists & %good_plot & %good_user & %good_ref (w/ reference data)
!                    = %exists & %good_plot & %good_user & %good_meas & %good_ref 
!                                                          (w/ measured & reference data)
! %useit_opt    -- Datum is possibly valid for optimization (minimizing the merit function):
!                    = %exists & %good_user & %good_opt (w/o reference data)
!                    = %exists & %good_user & %good_opt & %good_ref (w/ ref data)
! A datum is used in the optimization if both %useit_opt & %good_meas are true.

type tao_data_struct
  character(40) :: ele_name = ''           ! Name of the lattice element where datum is evaluated.
  character(40) :: ele_start_name = ''     ! Name of starting lattice element when there is a range 
  character(40) :: ele_ref_name = ''       ! Name of reference lattice element
  character(300) :: data_type = ''         ! Type of data: 'orbit.x', etc.
  character(40) :: merit_type = ''         ! Type of constraint: 'target', 'max', 'min', etc.
  character(40) :: id = ''                 ! Used by Tao extension code. Not used by Tao directly.
  character(20) :: data_source = ''        ! 'lat', 'beam', 'data' or 'var'. Last two used for expressions.
  integer :: ix_uni = -1                   ! Universe index of datum.
  integer :: ix_bunch = 0                  ! Bunch number to get the data from.
  integer :: ix_branch = 0                 ! Index of the associated lattice branch.
  integer :: ix_ele = -1                   ! Index of the lattice element corresponding to ele_name
  integer :: ix_ele_start = -1             ! Index of lattice elment when there is a range 
  integer :: ix_ele_ref = -1               ! Index of lattice elment when there is a reference.
  integer :: ix_ele_merit = -1             ! Index of lattice elment where merit is evaluated.
  integer :: ix_d1 = -1                    ! Index number in u%d2_data(i)%d1_data(j)%d(:) array.
  integer :: ix_data = -1                  ! Index of this datum in the u%data(:) array of data_structs.
  integer :: ix_dModel = -1                ! Row number in the dModel_dVar derivative matrix.
  integer :: eval_point = anchor_end$      ! Where to evaluate the data relative to the lattice element.
  real(rp) :: meas_value = 0               ! Measured datum value. 
  real(rp) :: ref_value = 0                ! Measured datum value from the reference data set.
  real(rp) :: model_value = 0              ! Datum value as calculated from the model.
  real(rp) :: design_value = 0             ! What the datum value is in the design lattice.
  real(rp) :: old_value = 0                ! The model_value at some previous time.
  real(rp) :: base_value = 0               ! The value as calculated from the base model.
  real(rp) :: error_rms = 0                ! RMS of data errors.
  real(rp) :: delta_merit = 0              ! Diff used to calculate the merit function term 
  real(rp) :: weight = 0                   ! Weight for the merit function term.
  real(rp) :: invalid_value = 0            ! Value used in merit calc if good_model = F (or possibly good_design & good_base).
  real(rp) :: merit = 0                    ! Merit function term value: weight * delta^2
  real(rp) :: s = real_garbage$            ! longitudinal position of ele.
  real(rp) :: s_offset = 0                 ! Offset of the evaluation point.
  type (spin_axis_struct) :: spin_axis = spin_axis_struct()  ! Input axes for spin g-matrix calculations.
  logical :: err_message_printed = .false. ! Used to prevent zillions of error messages being generated
  logical :: exists = .false.              ! See above
  logical :: good_model = .false.          ! See above
  logical :: good_base = .false.           ! See above
  logical :: good_design = .false.         ! See above
  logical :: good_meas = .false.           ! See above
  logical :: good_ref = .false.            ! See above
  logical :: good_user = .true.            ! See above
  logical :: good_opt = .true.             ! See above
  logical :: good_plot = .true.            ! See above
  logical :: useit_plot = .false.          ! See above
  logical :: useit_opt = .false.           ! See above
  type (tao_d1_data_struct), pointer :: d1 => null() ! Pointer to the parent d1_data_struct 
  type (tao_eval_stack1_struct), allocatable :: stack(:)
end type tao_data_struct

! A d1_data_struct represents, say, all the horizontal orbit data.
! The d1_data_struct has a pointer to the appropriate section in 
!   the u%data array. 

type tao_d1_data_struct
  character(40) name          ! Eg: 'x', etc.
  type (tao_d2_data_struct), pointer :: d2 => null() ! ptr to parent d2_data
  type (tao_data_struct), pointer :: d(:) => null()
                              ! Pointer to the appropriate section in u%data
end type

! A d2_data_struct represents all of a type of data. Eg: All orbit data.
! The d2_data_struct has pointers to the approprite d1_data_structs
! %ix_data and %ix_ref are used if the external data files are 
!   sequentially numbered.

type tao_d2_data_struct
  character(40) :: name = ''             ! Name to be used with commands.
  character(200) :: data_file_name = ''  ! Data file name .
  character(200) :: ref_file_name = ''   ! Reference file name.
  character(24) :: data_date = ''        ! Data measurement date.
  character(24) :: ref_date = ''         ! Reference data measurement date.
  character(80) :: descrip(10) = ''      ! Array for descriptive information.
  type (tao_d1_data_struct), allocatable :: d1(:) ! Points to children 
  integer ix_universe                    ! Index of universe this is in.
  integer ix_d2_data                     ! Index in u%d2_data(:) array.
  integer ix_ref                         ! Index of the reference data set. 
  logical :: data_read_in = .false.      ! A data set has been read in?
  logical :: ref_read_in = .false.       ! A reference data set has been read in?
end type

! A tao_data_array_struct is just a pointer to a tao_data_struct.
! This is used to construct arrays of tao_data_structs.

type tao_data_array_struct
  type (tao_data_struct), pointer :: d
end type

type tao_d1_data_array_struct
  type (tao_d1_data_struct), pointer :: d1
end type

type tao_d2_data_array_struct
  type (tao_d2_data_struct), pointer :: d2
end type

!-----------------------------------------------------------------------
! tao_var_slave_struct is for defining an array of pointers to variables
! in the tao_var_struct

type tao_var_slave_struct
  integer :: ix_uni = 1            ! universe index.
  integer :: ix_branch = 0
  integer :: ix_ele = -1           ! Index of element in the u%lattice%ele(:) array.
  real(rp), pointer :: model_value => null() ! Pointer to the variable in the model lat.
  real(rp), pointer :: base_value => null()  ! Pointer to the variable in the base lat.
end type  

! The var_struct defined the fundamental variable structure.
! The super_universe_struct will hold an array of var_structs: s%var(:).
!
! %exists     -- The variable exists. Non-existent variables can serve as place
!                  holders in the s%var array.
! %good_var   -- The variable can be varied. Used by the lm optimizer to
!                  veto variables that do not change the merit function.
! %good_user  -- What the user has selected using the use, veto, and restore 
!                  commands.
! %good_opt   -- Not modified by Tao. Setting is reserved to be done by extension code.
! %good_plot  -- Not modified by Tao. Setting is reserved to be done by extension code.
! %useit_opt  -- Variable is to be used for optimizing:
!                  %useit_opt = %exists & %good_user & %good_opt & %good_var
! %useit_plot -- If True variable is used in plotting variable values.
!                  %useit_plot = %exists & %good_plot & %good_user
!
! With common_lattice = True => var%slave(:)%model_value will point to the working universe.

type tao_var_struct
  character(40) :: ele_name = ''    ! Associated lattice element name.
  character(40) :: attrib_name = '' ! Name of the attribute to vary.
  character(40) :: id = ''          ! Used by Tao extension code. Not used by Tao directly.
  type (tao_var_slave_struct), allocatable :: slave(:)
  type (tao_var_slave_struct) :: common_slave
  integer :: ix_v1 = 0              ! Index of this var in the s%v1_var(i)%v(:) array.
  integer :: ix_var = 0             ! Index number of this var in the s%var(:) array.
  integer :: ix_dvar = -1           ! Column in the dData_dVar derivative matrix.
  integer :: ix_attrib = 0          ! Index in ele%value(:) array if appropriate.
  integer :: ix_key_table = 0       ! Has a key binding?
  real(rp), pointer :: model_value => null()     ! Model value.
  real(rp), pointer :: base_value => null()      ! Base value.
  real(rp) :: design_value = 0      ! Design value from the design lattice.
  real(rp) :: scratch_value = 0     ! Scratch space used by Tao.
  real(rp) :: old_value = 0         ! Scratch space used by Tao.
  real(rp) :: meas_value = 0        ! The value when the data measurement was taken.
  real(rp) :: ref_value = 0         ! Value when the reference measurement was taken.
  real(rp) :: correction_value = 0  ! Value determined by a fit to correct the lattice.
  real(rp) :: high_lim = -1d30      ! High limit for the model_value.
  real(rp) :: low_lim = 1d30        ! Low limit for the model_value.
  real(rp) :: step = 0              ! Sets what is a small step for varying this var.
  real(rp) :: weight = 0            ! Weight for the merit function term.
  real(rp) :: delta_merit = 0       ! Diff used to calculate the merit function term.
  real(rp) :: merit = 0             ! merit_term = weight * delta^2.
  real(rp) :: dMerit_dVar = 0       ! Merit derivative.     
  real(rp) :: key_val0 = 0          ! Key base value
  real(rp) :: key_delta = 0         ! Change in value when a key is pressed.
  real(rp) :: s = 0                 ! longitudinal position of ele.
  real(rp) :: extend_val            ! For extension code. Not used by Tao.
  character(40) :: merit_type = ''  ! 'target' or 'limit'
  logical :: exists = .false.       ! See above
  logical :: good_var = .false.     ! See above
  logical :: good_user = .true.     ! See above
  logical :: good_opt = .false.     ! See above
  logical :: good_plot = .false.    ! See above
  logical :: useit_opt = .false.    ! See above
  logical :: useit_plot = .false.   ! See above
  logical :: key_bound = .false.    ! Variable bound to keyboard key?
  type (tao_v1_var_struct), pointer :: v1 => null() ! Pointer to the parent.
end type tao_var_struct  

! A v1_var_struct represents, say, all the quadrupole power supplies.
! The v1_var_struct has a pointer to a section in the s%var array. 

type tao_v1_var_struct
  character(40) :: name = ''       ! V1 variable name. Eg: 'quad_k1'.
  integer ix_v1_var                ! Index to s%v1_var(:) array
  type (tao_var_struct), pointer :: v(:) => null() 
                                   ! Pointer to the appropriate section in s%var.
end type

! A tao_var_array_struct is just a pointer to a tao_var_struct.
! This is used to construct arrays of tao_var_structs.

type tao_var_array_struct
  type (tao_var_struct), pointer :: v
end type

type tao_v1_var_array_struct
  type (tao_v1_var_struct), pointer :: v1
end type


!------------------------------------------------------------------------
! Building wall structure

type tao_building_wall_orientation_struct
  real(rp) :: theta = 0
  real(rp) :: x_offset = 0
  real(rp) :: z_offset = 0
end type

type tao_building_wall_point_struct
  real(rp) :: z = 0, x = 0                  ! Global floor position
  real(rp) :: radius = 0                    ! Arc radius. +r -> CW rotation, same as bends. 
  real(rp) :: z_center = 0, x_center = 0    ! Arc center.
end type

type tao_building_wall_section_struct
  character(40) :: name = ''
  character(16) :: constraint = ''   ! "left_side" or "right_side" constraint.
  type (tao_building_wall_point_struct), allocatable :: point(:)
end type

type tao_building_wall_struct
  type (tao_building_wall_orientation_struct) :: orientation = tao_building_wall_orientation_struct()
  type (tao_building_wall_section_struct), allocatable :: section(:)
end type

!------------------------------------------------------------------------
! global parameters that the user has direct access to.
! Also see: tao_common_struct.
! If this structure is changed, change tao_set_global_cmd.

type tao_global_struct
  real(rp) :: lm_opt_deriv_reinit = -1   ! Reinit derivative matrix cutoff
  real(rp) :: de_lm_step_ratio = 1       ! Scaling for step sizes between DE and LM optimizers.
  real(rp) :: de_var_to_population_factor = 5.0_rp ! DE population = max(n_var*factor, 20)
  real(rp) :: lmdif_eps = 1e-12          ! tollerance for lmdif optimizer.
  real(rp) :: lmdif_negligible_merit = 1d-30
  real(rp) :: svd_cutoff = 1e-5          ! SVD singular value cutoff.
  real(rp) :: unstable_penalty = 1e-3    ! Used in unstable_ring datum merit calculation.
  real(rp) :: merit_stop_value = -1      ! Merit value below which an optimizer will stop.
  real(rp) :: dmerit_stop_value = 0      ! Fractional Merit change below which an optimizer will stop.
  real(rp) :: random_sigma_cutoff = -1   ! cut-off in sigmas.
  real(rp) :: delta_e_chrom = 0          ! delta E used from chrom calc.
  integer :: n_opti_cycles = 20          ! number of optimization cycles
  integer :: n_opti_loops = 1            ! number of optimization loops
  integer :: phase_units = radians$      ! Phase units on output.
  integer :: bunch_to_plot = 1           ! Which bunch to plot
  integer :: random_seed = 0             ! Use system clock by default
  integer :: n_top10_merit = 10          ! Number of top merit constraints to print.
  integer :: srdt_gen_n_slices = 10      ! number times to slice elements for summation RDT calculation
  integer :: srdt_sxt_n_slices = 20      ! number times to slice sextupoles for summation RDT calculation
  logical :: srdt_use_cache = .true.     ! create cache for SRDT calculations.  Can use lots of memory if srdt_*_n_slices large.
  character(12) :: quiet = 'off'                    ! "all", or "output". Print I/O when running a command file?
  character(16) :: random_engine = 'pseudo'         ! Non-beam random number engine
  character(16) :: random_gauss_converter = 'exact' ! Non-beam
  character(16) :: track_type    = 'single'         ! or 'beam'  
  character(40) :: prompt_string = 'Tao'
  character(16) :: prompt_color = 'DEFAULT'         ! See read_a_line routine for possible settings.
  character(16) :: optimizer     = 'lm'             ! optimizer to use.
  character(40) :: print_command = 'lpr'
  character(80) :: var_out_file  = 'var#.out'
  character(100) :: history_file = '~/.history_tao'
  logical :: beam_timer_on = .false.              ! For timing the beam tracking calculation.
  logical :: box_plots = .false.                  ! For debugging plot layout issues.
  logical :: command_file_print_on = .true.       ! Depracated. No longer used. 
  logical :: concatenate_maps = .false.           ! False => tracking using DA. 
  logical :: debug_on = .false.                   ! For debugging.
  logical :: derivative_recalc = .true.           ! Recalc before each optimizer run?
  logical :: derivative_uses_design = .false.     ! Derivative calc uses design lattice instead of model?
  logical :: disable_smooth_line_calc             ! Global disable of the smooth line calculation.
  logical :: draw_curve_off_scale_warn = .true.   ! Display warning on graphs?
  logical :: external_plotting = .false.          ! Used with matplotlib and gui.
  logical :: label_lattice_elements = .true.      ! For lat_layout plots
  logical :: label_keys = .true.                  ! For lat_layout plots
  logical :: lattice_calc_on = .true.             ! Turn on/off beam and single particle calculations.
  logical :: only_limit_opt_vars = .false.        ! Only apply limits to variables used in optimization.
  logical :: opt_with_ref = .false.               ! use reference data in optimization?
  logical :: opt_with_base = .false.              ! use base data in optimization?
  logical :: optimizer_allow_user_abort = .true.  ! See Tao manual for more details.
  logical :: optimizer_var_limit_warn = .true.    ! Warn when vars reach a limit with optimization.
  logical :: orm_analysis = .false.               ! ORM using MDSA? 
  logical :: plot_on = .true.                     ! Do plotting?
  logical :: rad_int_calc_on = .true.             ! Radiation integrals calculation on/off.
  logical :: rf_on = .false.                      ! RFcavities on or off? Does not affect lcavities.
  logical :: single_step = .false.                ! For debugging and demonstrations: Single step through a command file?
  logical :: stop_on_error = .true.               ! For debugging: False prevents tao from exiting on an error.
  logical :: svd_retreat_on_merit_increase = .true.
  logical :: var_limits_on = .true.               ! Respect the variable limits?
  logical :: wait_for_CR_in_single_mode = .false. ! For use with a python GUI.
end type

!

type tao_alias_struct
  character(40) :: name = ''
  character(200) :: expanded_str = ''
end type

type tao_command_file_struct
  character(200) :: full_name = ''
  character(200) :: dir = './'
  integer :: ix_unit
  character(40) cmd_arg(9)          ! Command file arguments.
  logical :: paused = .false.       ! Is the command file paused?
  integer :: n_line = 0             ! Current line number
end type

! tao_common_struct is for the global parameters that the user should not have direct access to.
! Also see tao_global_struct.

type tao_common_struct
  type (tao_alias_struct) alias(200)
  type (tao_alias_struct) key(100)
  type (tao_universe_struct), pointer :: u_working          ! Index of working universe.
  type (tao_command_file_struct), allocatable :: cmd_file(:)
  type (named_number_struct), allocatable :: symbolic_num(:)    ! Named numbers
  type (tao_plot_region_struct), allocatable :: plot_place_buffer(:)  ! Used when %external_plotting is on.
  real(rp), allocatable :: covar(:,:), alpha(:,:)
  real(rp) :: dummy_target = 0           ! Dummy varaible
  integer ix_ref_taylor, ix_ele_taylor   ! Taylor map end points
  integer :: n_alias = 0
  integer :: cmd_file_level = 0          ! For nested command files. 0 -> no command file.
  integer :: ix_key_bank = 0             ! For single mode.
  integer :: n_universes = 1   
  integer :: default_universe = 1        ! Default universe to work with.
  integer :: default_branch = 0          ! Default lattice branch to work with.
  integer :: ix_history = 0 ! present index to command history array
  integer :: n_history      ! present history index
  logical :: initialized = .false.              ! Does tao_init() need to be called?
  logical :: cmd_file_paused
  logical :: use_cmd_here  = .false.            ! Used for the cmd history stack
  logical :: cmd_from_cmd_file = .false.        ! was command from a command file?
  logical :: use_saved_beam_in_tracking = .false.
  logical :: single_mode = .false.
  logical :: combine_consecutive_elements_of_like_name
  logical :: common_lattice = .false.      
  logical :: have_tracked_beam = .false.      ! Used to catch error when beam plotting without having tracked a beam.
  logical :: init_plot_needed      = .true.   ! reinitialize plotting?
  logical :: init_beam             = .true.   ! Used by custom programs to control Tao init
  logical :: init_var              = .true.   ! Used by custom programs to control Tao init
  logical :: init_read_lat_info    = .true.   ! Used by custom programs to control Tao init
  logical :: parse_cmd_args        = .true.   ! Used by custom programs to control Tao init
  logical :: optimizer_running     = .false. 
  logical :: have_datums_using_expressions = .false.
  logical :: print_to_terminal = .true.        ! Print command prompt to the terminal? For use with GUIs.
  logical :: lattice_calc_done = .false.       ! Used by GUI for deciding when to refresh.
  character(100) :: cmd                        ! Used for the cmd history
  character(16) :: init_name = 'Tao'           ! label for initialization          
  character(200) :: hook_init_file = ''        ! 
  character(200) :: hook_lat_file = ''             ! To be set by tao_hook_parse_command_args
  character(200) :: hook_beam_file = ''            ! To be set by tao_hook_parse_command_args
  character(200) :: hook_data_file = ''            ! To be set by tao_hook_parse_command_args
  character(200) :: hook_plot_file = ''            ! To be set by tao_hook_parse_command_args
  character(200) :: hook_startup_file = ''         ! To be set by tao_hook_parse_command_args
  character(200) :: hook_var_file = ''             ! To be set by tao_hook_parse_command_args
  character(200) :: hook_building_wall_file = ''   ! To be set by tao_hook_parse_command_args
  character(200) :: saved_cmd_line = ''            ! Saved part of command line when there are mulitple commands on a line
  character(200) :: init_file_arg_path = ''        ! Path part of init_tao_file
  character(200) :: lattice_file_arg = ''          ! -lattice_file        command line argument.
  character(200) :: hook_init_file_arg = ''        ! -hook_init_file      command line argument
  character(200) :: init_file_arg = ''             ! -init_file           command line argument.
  character(200) :: beam_file_arg = ''             ! -beam_file           command line argument.
  character(200) :: beam_track_data_file_arg = ''  ! -beam_track_data_file       command line argument.
  character(200) :: beam_init_position_file_arg = '' ! -beam_init_position_file command line argument.
  character(200) :: data_file_arg = ''             ! -data_file           command line argument.
  character(200) :: plot_file_arg = ''             ! -plot_file           command line argument.
  character(200) :: startup_file_arg = ''          ! -startup_file        command line argument.
  character(200) :: var_file_arg = ''              ! -var_file            command line argument.
  character(200) :: building_wall_file_arg = ''    ! -building_wtrack_data_file  command line argument.
  character(16) :: geometry_arg = ''               ! -geometry            command line argument.
  character(80) :: slice_lattice_arg = ''          ! -slice_lattice       command line argument.
  character(12) :: log_startup_arg = ''            ! -log_startup         command line argument
  character(12) :: no_stopping_arg = ''            ! -no_stopping         command line argument
  character(12) :: noplot_arg = ''                 ! -noplot              command line argument
  character(12) :: no_rad_int_arg = ''             ! -no_rad_int          command line argument
  character(12) :: debug_arg = ''                  ! -debug               command line argument
  character(12) :: disable_smooth_line_calc_arg = '' ! -disable_smooth_line_calc
  character(12) :: rf_on_arg = ''                  ! -rf_on               command line argument
  character(12) :: prompt_color_arg = ''           ! -prompt_color        command line argument
  character(12) :: quiet_arg = ''                  ! -quiet               command line argument
  character(12) :: noinit_arg = ''                 ! -noinit              command line argument
  character(80) :: single_mode_buffer = ''
  character(100) :: unique_name_suffix
  character(16) :: valid_plot_who(10)            ! model, base, ref etc...
end type

!-----------------------------------------------------------------------
! scratch space

type tao_beam_shake_struct
  real(rp) cbar(2,2)
  real(rp) k_11a, k_12a, k_12b, k_22b
  real(rp) amp_a, amp_b, amp_na, amp_nb
  real(rp) :: one = 1.0
  logical :: coupling_calc_done = .false.
  logical :: amp_calc_done = .false.
end type

type tao_spin_map_struct
  type (taylor_struct) :: orbit_taylor(6) = taylor_struct()  ! Not yet used.
  type (taylor_struct) :: spin_taylor(0:3) = taylor_struct() ! Not yet used.
  type (spin_axis_struct) :: axis_dat = spin_axis_struct()   ! Axes from data_struct
  type (spin_axis_struct) :: axis0 = spin_axis_struct()      ! Initial axes.
  type (spin_axis_struct) :: axis1 = spin_axis_struct()      ! Final axes.
  integer :: ix_ele = 0, ix_ref = 0, ix_uni = 0
  integer :: ix_branch = 0
  real(rp) :: mat8(8,8) = 0
end type

type tao_spin_polarization_struct
  logical valid_value
  character(60) why_invalid
  real(rp) :: pol_limit = real_garbage$  ! Equalibrium Polarization calculated via the Derbenev-Kondratenko-Mane formula.
  real(rp) :: pol_rate = real_garbage$   ! Polarization rate (1/sec).
  real(rp) :: depol_rate = real_garbage$ ! Depolarization rate (1/sec).
end type

type tao_scratch_space_struct
  type (tao_beam_shake_struct), allocatable :: cc(:)
  type (ele_pointer_struct), allocatable :: eles(:)
  type (tao_d1_data_array_struct), allocatable :: d1_array(:)
  type (tao_v1_var_array_struct), allocatable :: v1_array(:)
  type (tao_eval_stack1_struct), allocatable :: stack(:)
  type (tao_var_array_struct), allocatable :: var_array(:)
  type (all_pointer_struct), allocatable :: attribs(:)
  type (tao_data_var_component_struct), allocatable :: comp(:)
  type (tao_expression_info_struct), allocatable :: info(:)
  type (tao_expression_info_struct), allocatable :: info_x(:), info_y(:), info_ix(:)
  type (tao_spin_map_struct), allocatable :: spin_map(:)
  type (tao_spin_polarization_struct) spin
  logical, allocatable :: picked(:)
  logical, allocatable :: this_u(:)
  real(rp), allocatable :: axis1(:), axis2(:), axis3(:)
  real(rp), allocatable :: x(:), y(:), err(:)
  real(rp), allocatable :: y_value(:)
  complex(rp), allocatable :: srdt_cache(:,:,:)
end type

type (tao_scratch_space_struct), save, target :: scratch

!-----------------------------------------------------------------------

type tao_lat_mode_struct
  real(rp) chrom
  real(rp) growth_rate
end type

type tao_sigma_mat_struct
  real(rp) sigma(6,6)
end type

! For caching lattice calculations associated with plotting.

type tao_plot_cache_struct
  type (ele_struct) ele
  type (coord_struct) orbit
end type

! The %bunch_params(:) array has a 1-to-1 correspondence with the lattice elements.
! NOTE: If tao_lattice_branch_struct is modified then the routine
! tao_lattice_branch_equal_tao_lattice_branch must be modified as well.

type tao_lattice_branch_struct
  type (summation_rdt_struct) srdt
  type (bunch_params_struct), allocatable :: bunch_params(:)
  type (tao_sigma_mat_struct), allocatable :: linear(:) ! Sigma matrix derived from linear lattice.
  type (coord_struct), allocatable :: orbit(:)
  type (coord_struct) orb0                     ! For saving beginning orbit
  type (lat_struct) :: high_E_lat, low_E_lat  ! For chrom calc.
  type (tao_plot_cache_struct), allocatable :: plot_cache(:)  ! Plotting data cache
  integer track_state
  logical has_open_match_element
  logical :: plot_cache_valid = .false.        ! Valid plotting data cache?
  real(rp) :: cache_x_min = 0, cache_x_max = 0
  integer :: cache_n_pts = 0
  type (normal_modes_struct) modes             ! Synchrotron integrals stuff
  type (tao_lat_mode_struct) a, b
  integer ix_rad_int_cache                     ! Radiation integrals cache index.
  type (normal_modes_struct) modes_rf_on       ! Synchrotron integrals stuff
end type

! Structure to hold a single lat_struct (model, base, or design) in
! a universe along with stuff like radiation integrals, etc.

type tao_lattice_struct
  character(8) :: name                         ! "model", "base", or "design".
  type (lat_struct) lat                        ! lattice structures
  type (tao_universe_struct), pointer :: u => null()  ! Parent universe
  type (rad_int_all_ele_struct) rad_int
  type (tao_lattice_branch_struct), allocatable :: tao_branch(:)
end type

! Universe wide structure for information that does not fit anywhere else.

!-----------------------------------------------------------------------
! tao_element_struct is for saving per-element information.

type tao_element_struct
  type (beam_struct) beam         ! Beam distribution at element.
  logical save_beam_internally    ! Save beam here? Beam also saved at fork elements and at track ends.
  logical save_beam_to_file       ! Save beam to a file? Beam also saved at fork elements and at track ends.
end type

! Information for a particular lattice branch of a particular universe.

type tao_universe_branch_struct
  type (tao_element_struct), allocatable :: ele(:) ! Per element information
  character(40) :: particle_track_start = '', particle_track_end = ''
  integer :: ix_particle_track_start = 0
  integer :: ix_particle_track_end = -1
end type

! Beam information for a particular universe 

type tao_beam_struct
  type (beam_init_struct) :: beam_init            ! Beam distrubution at beginning of lattice
  type (beam_struct) beam_at_start                ! Initial
  integer :: ix_track_start = 0                   ! Element start index of tracking
  integer :: ix_track_end = -1                    ! Element end index of tracking
  logical :: init_starting_distribution = .false. ! Init beam
  character(40) :: track_start = '', track_end = ''
  character(200) :: track_data_file = ''          ! track data from previous simulation for reanalysis.
  character(160) :: saved_at = ''
  character(160) :: dump_file = ''
  character(160) :: dump_at = ''
end type

! Logicals that determine what calculations need to be done.
! Keep data and plotting separate since when optimizing will only do a calc if the data needs it

type tao_universe_calc_struct
  logical :: rad_int_for_data = .false.           ! Do the radiation integrals need to be computed for
  integer :: srdt_for_data = 0                    ! 0 = false, 1 = 1st order, 2 = 1st & 2nd order
  logical :: rad_int_for_plotting = .false.       !   data or plotting?
  logical :: chrom_for_data = .false.             ! Does the chromaticity need to be computed for
  logical :: chrom_for_plotting = .false.         !   data or plotting? 
  logical :: beam_sigma_for_data = .false.        ! Do the beam sigmas need to be computed for
  logical :: beam_sigma_for_plotting = .false.    !   data or plotting? 
  logical :: dynamic_aperture = .false.           ! Do the dynamic_aperture calc?
  logical :: one_turn_map = .false.               ! Compute the one turn map?
  logical :: lattice = .true.                     ! Used to indicate which lattices need tracking done.
  logical :: twiss = .true.                       ! calc linear transfer matrix?
  logical :: track = .true.                       ! tracking needs to be done?
  logical :: spin_matrices = .false.              ! Calculate G and D spin matrices?
end type

!-----------------------------------------------------------------------
! MPI information structure

type tao_mpi_struct
  logical :: on = .false.           ! Is MPI on?
  logical :: master = .true.        ! Is this the master task? If yes, rank == 0
  integer :: rank = 0               ! Rank of task (rank is 0, 1, 2, ... n_tasks-1 ) 
  integer :: max_rank = 0           ! Maximum rank, should be n_tasks-1
  character(160) :: host_name  =''  ! Name of the host machine
end type

!-----------------------------------------------------------------------
! tao_dynamic_aperture_struct

type tao_dynamic_aperture_struct
  type (aperture_param_struct) param
  type (aperture_scan_struct), allocatable :: scan(:) ! One scan for each pz.
  real(rp), allocatable :: pz(:)
  integer :: ix_branch = 0
end type

!-----------------------------------------------------------------------
! Wave analysis structures

type tao_wave_kick_pt_struct
  real(rp) :: phi_s, phi_r, phi, amp
  real(rp) :: s                 ! s-position of kick
  integer :: ix_dat_before_kick ! Index of datum in data array just before the kick.
  type (ele_struct), pointer :: ele             ! lattice element at position of kick.
end type  

type tao_wave_struct     ! Struct for wave analysis
  character(40) data_type
  real(rp) rms_rel_a, rms_rel_b, rms_rel_as, rms_rel_bs, rms_rel_ar, rms_rel_br
  real(rp) rms_rel_k, rms_rel_ks, rms_rel_kr 
  real(rp) rms_phi, rms_phi_s, rms_phi_r
  real(rp) amp_ba_s, amp_ba_r, chi_a, chi_c, chi_ba
  real(rp) amp_a(2), amp_b(2), amp_ba(2)
  real(rp) coef_a(4), coef_b(4), coef_ba(4)
  integer n_func   ! Number of functions used in the fit.
  integer :: ix_a1 = -1, ix_a2 = -1, ix_b1 = -1, ix_b2 = -1
  integer i_a1, i_a2, i_b1, i_b2, n_a, n_b
  integer i_curve_wrap_pt      ! Index of last point before wrap in curve array. 
  integer, allocatable :: ix_data(:) ! Translates from plot point to datum index
  integer n_kick
  type (tao_wave_kick_pt_struct), allocatable :: kick(:)
  type (tao_graph_struct) :: base_graph   ! Graph before curves extended to 1.5 periods.
  type (tao_plot_region_struct), pointer :: region    ! Where the wave plot is
  type (tao_d1_data_struct), pointer :: d1_dat        ! D1 data for analysis
end type

!-----------------------------------------------------------------------
! For scaling ping amplitude

type tao_ping_scale_struct
  real(rp) :: a_mode_meas = 1
  real(rp) :: a_mode_ref = 1
  real(rp) :: b_mode_meas = 1
  real(rp) :: b_mode_ref = 1
end type

!-----------------------------------------------------------------------
! A universe is a snapshot of a machine

type tao_universe_struct
  type (tao_universe_struct), pointer :: common => null()
  type (tao_lattice_struct), pointer :: model, design, base
  type (tao_beam_struct) beam
  type (tao_dynamic_aperture_struct) :: dynamic_aperture
  type (tao_universe_branch_struct), pointer :: uni_branch(:) ! Per element information
  type (tao_d2_data_struct), allocatable :: d2_data(:)   ! The data types 
  type (tao_data_struct), allocatable :: data(:)         ! Array of all data.
  type (tao_ping_scale_struct) ping_scale
  type (lat_struct) scratch_lat                          ! Scratch area.
  type (tao_universe_calc_struct) calc                   ! What needs to be calculated?
  real(rp), allocatable :: dModel_dVar(:,:)              ! Derivative matrix.
  integer :: ix_uni = -1                    ! Universe index.
  integer :: n_d2_data_used = -1            ! Number of used %d2_data(:) components.
  integer :: n_data_used = -1               ! Number of used %data(:) components.
  logical :: reverse_tracking = .false.     ! Reverse tracking direction?
  logical :: is_on = .true.                 ! universe turned on
  logical :: design_same_as_previous = .false.  ! Design lat same as the previous uni?
  logical :: picked_uni = .false.           ! Scratch logical.
end type

!-----------------------------------------------------------------------
! The super_universe is the structure that holds an array of universes.
! Essentially this holds all the information known to the program.

type tao_super_universe_struct
  type (tao_global_struct) global                          ! User accessible global variables.
  type (tao_common_struct) :: com
  type (tao_plot_page_struct) :: plot_page                 ! Defines the plot window.
  type (tao_v1_var_struct), allocatable :: v1_var(:)       ! The variable types
  type (tao_var_struct), allocatable :: var(:)             ! array of all variables.
  type (tao_universe_struct), allocatable :: u(:)          ! array of universes.
  type (tao_mpi_struct) mpi
  integer, allocatable :: key(:)
  type (tao_building_wall_struct) :: building_wall
  type (tao_wave_struct) :: wave 
  integer :: n_var_used = 0
  integer :: n_v1_var_used = 0
  type (tao_cmd_history_struct) :: history(1000) ! command history
end type

type (tao_super_universe_struct), save, target :: s

!-----------------------------------------------------------------------
contains

subroutine tao_deallocate_plot_cache(plot_cache)

implicit none

type (tao_plot_cache_struct), allocatable :: plot_cache(:)
integer i

!

if (.not. allocated(plot_cache)) return

do i = 1, size(plot_cache)
  call deallocate_ele_pointers(plot_cache(i)%ele)
enddo

deallocate(plot_cache)

end subroutine tao_deallocate_plot_cache

!-----------------------------------------------------------------------
! contains

subroutine tao_lattice_branches_equal_tao_lattice_branches (tlb1, tlb2)

implicit none

type (tao_lattice_branch_struct), intent(inout) :: tlb1(:)
type (tao_lattice_branch_struct), intent(in) :: tlb2(:)
integer i

!

do i = 1, size(tlb1)
  tlb1(i) = tlb2(i)
enddo

end subroutine tao_lattice_branches_equal_tao_lattice_branches

!-----------------------------------------------------------------------
! contains

subroutine tao_lattice_branch_equal_tao_lattice_branch (tlb1, tlb2)

implicit none

type (tao_lattice_branch_struct), intent(inout) :: tlb1
type (tao_lattice_branch_struct), intent(in) :: tlb2
integer n1, n2, i

!

if (allocated(tlb2%bunch_params)) then
  tlb1%bunch_params           = tlb2%bunch_params
else
  if (allocated(tlb1%bunch_params)) deallocate(tlb1%bunch_params)
endif

if (allocated(tlb2%linear)) then
  tlb1%linear                 = tlb2%linear
else
  if (allocated(tlb1%linear)) deallocate(tlb1%linear)
endif

if (allocated(tlb2%orbit)) then
  tlb1%orbit                  = tlb2%orbit
else
  if (allocated(tlb1%orbit)) deallocate(tlb1%orbit)
endif

tlb1%srdt                   = tlb2%srdt
tlb1%orb0                   = tlb2%orb0
tlb1%high_E_lat             = tlb2%high_E_lat
tlb1%low_E_lat              = tlb2%low_E_lat
tlb1%track_state            = tlb2%track_state
tlb1%has_open_match_element = tlb2%has_open_match_element
tlb1%plot_cache_valid       = tlb2%plot_cache_valid
tlb1%cache_x_min            = tlb2%cache_x_min
tlb1%cache_x_max            = tlb2%cache_x_max
tlb1%cache_n_pts            = tlb2%cache_n_pts
tlb1%modes                  = tlb2%modes
tlb1%a                      = tlb2%a
tlb1%b                      = tlb2%b
tlb1%ix_rad_int_cache       = tlb2%ix_rad_int_cache
tlb1%modes_rf_on            = tlb2%modes_rf_on


n1 = -1
if (allocated(tlb1%plot_cache)) n1 = size(tlb1%plot_cache)

n2 = -1
if (allocated(tlb2%plot_cache)) n2 = size(tlb2%plot_cache)

if (n1 > -1 .and. n1 /= n2) call tao_deallocate_plot_cache(tlb1%plot_cache)

if (n2 > -1) then
  if (.not. allocated(tlb1%plot_cache)) allocate (tlb1%plot_cache(n2))
  do i = 1, n2
    tlb1%plot_cache(i)%ele   = tlb2%plot_cache(i)%ele
    tlb1%plot_cache(i)%orbit = tlb2%plot_cache(i)%orbit
  enddo
endif

end subroutine tao_lattice_branch_equal_tao_lattice_branch 

!-----------------------------------------------------------------------
! contains

subroutine tao_lattice_equal_tao_lattice (lat1, lat2)

implicit none

type (tao_lattice_struct), intent(inout) :: lat1
type (tao_lattice_struct), intent(in) :: lat2
integer i

!

lat1%lat            = lat2%lat
lat1%rad_int        = lat2%rad_int
lat1%tao_branch     = lat2%tao_branch

end subroutine tao_lattice_equal_tao_lattice

end module
