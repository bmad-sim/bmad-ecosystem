!+
! Module tao_struct
!
! Module defining the basic tao structures
!
! If any pointers or allocatables are added remember to add a corresponding
! statment to tao_init\deallocate_everything.
!-

module tao_struct

use bmad_struct, only: rp, lat_struct, coord_struct, radians$, ele_struct, normal_modes_struct
use equal_mod
use quick_plot, only: qp_line_struct, qp_symbol_struct, qp_axis_struct, &
                      qp_rect_struct, qp_point_struct
use beam_def_struct, only: beam_init_struct, beam_struct, bunch_params_struct
use tao_parameters
use rad_int_common, only: rad_int_common_struct

integer, parameter :: model$ = 1, base$ = 2, design$ = 3
integer, parameter :: ix_common_uni$ = 0

logical, save, target :: forever_true$ = .true.

interface assignment (=)
  module procedure tao_lat_equal_tao_lat
end interface

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
  character(40), pointer :: s => null()
end type

type tao_eval_stack1_struct
  integer type
  character(60) :: name = ''  
  real(rp), allocatable :: value(:)
  logical, allocatable :: good(:)
  type (tao_real_pointer_struct), allocatable :: value_ptr(:)
end type

!----------------------------------------------------------------------

type tao_ele_shape_struct    ! for the element layout plot
  character(60) ele_name     ! element name
  character(16) shape        ! plot shape
  character(16) color        ! plot color
  real(rp) dy_pix            ! plot vertical height 
  character(16) label_type 
end type

type tao_wave_kick_pt_struct
  real(rp) :: phi_s, phi_r, phi, amp
  integer :: ix_dat
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
  integer i_wrap_pt      ! Index of last point before wrap in curve array. 
  integer, allocatable :: ix_data(:) ! Translates from plot point to datum index
  integer n_kick
  type (tao_wave_kick_pt_struct), allocatable :: kick(:)
end type

!-----------------------------------------------------------------------
! Plot structures.

type tao_title_struct
  character(100) string      ! title character string.
  real(rp) x, y              ! x, y rwt lower left corner
  character(16) units        ! %BOX, POINTS, etc...
  character(2) justify       ! Left, Center, or Right justification.
  logical draw_it            ! draw the title
end type

type tao_data_var_component_struct    ! Components to plot
  character(16) name         ! Eg: 'meas', 'ref', 'model', etc.
  real(rp) sign              ! +1 or -1
end type

! A curve is defined by a set of (x,y) points and the axis parameters.
! for example the horizontal orbit is one curve.

type tao_curve_struct
  character(40) :: name = ''               ! Name identifying the curve.
  character(40) :: data_source  = ''       ! 'lat', 'dat', 'var', etc.
  character(100) :: data_index  = ''       ! Used for calculating %ix_symb(:).
  character(100) :: data_type_x = ''       ! Used for data slices and phase space plots.
  character(100) :: data_type   = ''       ! 'orbit.x', etc.
  character(40) :: ele_ref_name = ''       ! Reference element.
  character(40) :: legend_text = ''        ! String to print in a curve legend. 
  type (tao_graph_struct), pointer :: g    ! pointer to parent graph 
  real(rp), allocatable :: x_line(:)       ! Coords for drawing a curve
  real(rp), allocatable :: y_line(:) 
  real(rp), allocatable :: x_symb(:)       ! Coords for drawing the symbols
  real(rp), allocatable :: y_symb(:) 
  integer, allocatable :: ix_symb(:)       ! Corresponding index in d1_data%d(:) array.
  real(rp) :: y_axis_scale_factor = 1  ! y-axis conversion from internal to plotting units.
  type (qp_line_struct) line       ! Line attributes
  type (qp_symbol_struct) symbol   ! Symbol attributes
  integer :: ix_universe = -1      ! Universe where data is. -1 => use s%global%u_view
  integer :: symbol_every = 1      ! Symbol every how many points.
  integer :: ix_branch = 0
  integer :: ix_ele_ref = -1       ! Index in lattice of reference element.
  integer :: ix_ele_ref_track = -1 ! = ix_ele_ref except for super_lord elements.
  integer :: ix_bunch = 0          ! Bunch to plot.
  logical :: use_y2 = .false.      ! Use y2 axis?
  logical :: draw_line = .true.    ! Draw a line through the data points?
  logical :: draw_symbols = .true. ! Draw a symbol at the data points?
  logical :: draw_symbol_index = .false. ! Draw the symbol index number curve%ix_symb?
  logical :: smooth_line_calc = .true.   ! Calculate data between element edge points?
end type

! A graph is a collection of overlayed curves with associated graph title, etc.
! For example a graph could contain just the horizontal orbit or could
! contain both overlayed horizontal and vertical orbits.

type tao_graph_struct
  character(40) name            ! Name identifying the graph
  character(40) type            ! 'data', 'lat_layout', 'key_table', 'phase_space'
  character(100) title
  character(100) title_suffix 
  character(100) text_legend(n_legend_maxx) ! Array for holding descriptive info.
  character(60) component       ! Who to plot. Eg: 'meas - design'
  character(80) why_invalid     ! Informative string to print.
  type (tao_curve_struct), allocatable :: curve(:)
  type (tao_plot_struct), pointer :: p ! pointer to parent plot
  type (qp_point_struct) text_legend_origin
  type (qp_point_struct) curve_legend_origin
  type (qp_axis_struct) x       ! X-axis parameters.
  type (qp_axis_struct) y       ! Y-axis attributes.
  type (qp_axis_struct) y2      ! Y-axis attributes.
  type (qp_rect_struct) margin  ! Margin around the graph.
  real(rp) :: x_axis_scale_factor = 1  ! x-axis conversion from internal to plotting units.
  integer box(4)                ! Defines which box the plot is put in.
  integer :: ix_branch = 0
  integer :: ix_universe = -1   ! Used for lat_layout plots.
  logical clip                  ! Clip plot at graph boundary.
  logical valid                 ! valid if all curve y_dat computed OK.
  logical y2_mirrors_y          ! Y2-axis same as Y-axis?
  logical limited               ! True if at least one data point past graph bounds.
  logical draw_axes             ! Draw axes, labels, etc?
  logical correct_xy_distortion ! T -> Shrink floor plan along one axis to give both axes the same scale.
  logical draw_curve_legend     ! Legend for displaying curve info.
  logical :: visible = .true.   ! To draw or not to draw. 
end type

! A plot is collection of graphs.
! For example a plot could contain three graphs. One for Cbar11, 
! One for Cbar12, and one for Cbar22.

type tao_plot_struct
  character(40) :: name = ' '                 ! Identifying name
  type (tao_graph_struct), allocatable :: graph(:)
                                              ! individual graphs of a plot
  type (qp_axis_struct) x                     ! X-axis parameters.
  type (tao_plot_region_struct), pointer :: r ! pointer to parent.
  character(16) x_axis_type                   ! 'index', 'ele_index', 's', 'none', 
                                              !         'floor', or 'phase_space'
  logical :: autoscale_gang_x = .true.        ! scale cmd scales graphs together?
  logical :: autoscale_gang_y = .true.        ! scale cmd scales graphs together?
end type

! A region defines a plot and where to position the plot on the plot page
! %location = (x1, x2, y1, y2) gives the plotting region in percent of the 
!   part of the page inside the page_border with respect to the lower left corner.
! Eg: %location = (0.0, 1.0, 0.5, 1.0) gives the top half of the page inside the border.

type tao_plot_region_struct
  character(40) :: name = ''     ! Eg: 'top', 'bottom'.
  type (tao_plot_struct) plot    ! Plot associated with this region
  real(rp) location(4)           ! location on page.
  logical visible                ! To draw or not to draw.
end type

! The plot_page defines the whole plotting window. 
! The plot_page contains a collection of regions along with 
! other info (border margins etc.).
! Note that the qp_com structure of quick_plot also is used to hold 
! plot page info.

type tao_plot_page_struct               
  character(80) ps_scale             ! scaling when creating PS files.
  real(rp) :: shape_height_max = 40  ! maximum half height for drawing elements.
  real(rp) size(2)                   ! width and height of window in pixels.
  real(rp) :: text_height = 12              ! In points. Scales the height of all text
  real(rp) :: main_title_text_scale  = 1.3  ! Relative to text_height
  real(rp) :: graph_title_text_scale = 1.1  ! Relative to text_height
  real(rp) :: axis_number_text_scale = 0.9  ! Relative to text_height
  real(rp) :: axis_label_text_scale  = 1.0  ! Relative to text_height
  real(rp) :: legend_text_scale      = 0.8  ! Relative to text_height
  real(rp) :: key_table_text_scale   = 0.9  ! Relative to text_height
  real(rp) :: curve_legend_line_len  = 50   ! Points
  real(rp) :: curve_legend_text_offset = 10 ! Points
  integer :: n_curve_pts = 401       ! Number of points for plotting a smooth curve
  integer id_window                  ! X window id number.
  type (tao_title_struct) title(2)   ! Titles at top of page.
  type (qp_rect_struct) border       ! Border around plots edge of page.
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
! %good_opt     -- Not modified by Tao proper and reserved for use by extension code.
! %good_plot    -- Not modified by Tao proper and reserved for use by extension code.
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
  character(40) ele_name        ! Name of the lattice element where datum is evaluated.
  character(40) ele_start_name  ! Name of starting lattice element when there is a range 
  character(40) ele_ref_name    ! Name of reference lattice element
  character(100) data_type   ! Type of data: 'orbit.x', etc.
  character(40) merit_type   ! Type of constraint: 'target', 'max', 'min', etc.
  character(20) data_source  ! 'lat', or 'beam'
  integer ix_bunch           ! Bunch number to get the data from.
  integer ix_branch          ! Index of the lattice branch of the element
  integer ix_ele             ! Index of the lattice element corresponding to ele_name
  integer ix_ele_start       ! Index of lattice elment when there is a range 
  integer ix_ele_ref         ! Index of lattice elment when there is a reference.
  integer ix_ele_merit       ! Index of lattice elment where merit is evaluated.
  integer ix_d1              ! Index number in u%d2_data(i)%d1_data(j)%d(:) array.
  integer ix_data            ! Index of this datum in the u%data(:) array of data_structs.
  integer ix_dModel          ! Row number in the dModel_dVar derivative matrix.
  real(rp) meas_value        ! Measured datum value. 
  real(rp) ref_value         ! Measured datum value from the reference data set.
  real(rp) model_value       ! Datum value as calculated from the model.
  real(rp) design_value      ! What the datum value is in the design lattice.
  real(rp) old_value         ! The model_value at some previous time.
  real(rp) base_value        ! The value as calculated from the base model.
  real(rp) fit_value         ! The value as calculated from a fitting procedure.
  real(rp) delta_merit       ! Diff used to calculate the merit function term 
  real(rp) weight            ! Weight for the merit function term
  real(rp) invalid_value     ! Value used in merit calc if good_model = False.
  real(rp) merit             ! Merit function term value: weight * delta^2
  real(rp) s                 ! longitudinal position of ele.
  logical exists             ! See above
  logical good_model         ! See above
  logical good_base          ! See above
  logical good_design        ! See above
  logical good_meas          ! See above
  logical good_ref           ! See above
  logical good_user          ! See above
  logical good_opt           ! See above
  logical good_plot          ! See above
  logical useit_plot         ! See above
  logical useit_opt          ! See above
  type (tao_d1_data_struct), pointer :: d1 => null() 
                             ! Pointer to the parent d1_data_struct 
  type (tao_eval_stack1_struct), allocatable :: stack(:)
end type tao_data_struct

! A d1_data_struct represents, say, all the horizontal orbit data.
! The d1_data_struct has a pointer to the appropriate section in 
!   the u%data array. 

type tao_d1_data_struct
  character(40) name        ! Eg: 'x', etc.
  integer ix_data           ! index of the 0th element in u%data.
  type (tao_d2_data_struct), pointer :: d2 => null() ! ptr to parent d2_data
  type (tao_data_struct), pointer :: d(:) => null()  
                            ! Pointer to the appropriate section in u%data
end type

! A d2_data_struct represents all of a type of data. Eg: All orbit data.
! The d2_data_struct has pointers to the approprite d1_data_structs
! %ix_data and %ix_ref are used if the external data files are 
!   sequentially numbered.

type tao_d2_data_struct
  character(40) name              ! Name to be used with commands.
  character(200) data_file_name   ! Data file name .
  character(200) ref_file_name    ! Reference file name.
  character(20) data_date         ! Data measurement date.
  character(20) ref_date          ! Reference data measurement date.
  character(80) descrip(n_descrip_maxx) ! Array for descriptive information.
  type (tao_d1_data_struct), allocatable :: d1(:) ! Points to children 
  integer ix_uni                  ! Index of universe this is in.
  integer ix_data                 ! Index of the data set.
  integer ix_ref                  ! Index of the reference data set. 
  logical data_read_in            ! A data set has been read in?
  logical ref_read_in             ! A reference data set has been read in?
end type

! A tao_data_array_struct is just a pointer to a tao_data_struct.
! This is used to construct arrays of tao_data_structs.

type tao_data_array_struct
  type (tao_data_struct), pointer :: d
end type

type tao_d1_data_array_struct
  type (tao_d1_data_struct), pointer :: d1
end type

!-----------------------------------------------------------------------
! tao_this_var_struct is for defining an array of pointers to variables
! in the tao_var_struct

type tao_this_var_struct
  integer ix_uni            ! universe index.
  integer :: ix_branch = 0
  integer ix_ele            ! Index of element in the u%lattice%ele(:) array.
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
! %good_opt   -- Not modified by Tao proper and reserved for use by extension code.
! %good_plot  -- Not modified by Tao proper and reserved for use by extension code.
! %useit_opt  -- Variable is to be used for optimizing:
!                  = %exists & %good_user & %good_opt & %good_var
! %useit_plot -- Variable value to be plotted:
!                  = %exists & %good_plot & %good_user
!
! With common_lattice = True => var%this(:)%model_value will point to the working universe.

type tao_var_struct
  character(40) ele_name    ! Associated lattice element name.
  character(40) attrib_name ! Name of the attribute to vary.
  type (tao_this_var_struct), allocatable :: this(:)
  type (tao_this_var_struct) :: common
  integer ix_v1             ! Index of this var in the s%v1_var(i)%v(:) array.
  integer ix_var            ! Index number of this var in the s%var(:) array.
  integer ix_dvar           ! Column in the dData_dVar derivative matrix.
  integer ix_attrib         ! Index in ele%value(:) array if appropriate.
  integer ix_key_table      ! Has a key binding?
  real(rp), pointer :: model_value      ! Model value.
  real(rp), pointer :: base_value       ! Base value.
  real(rp) design_value     ! Design value from the design lattice.
  real(rp) old_value        ! The model_value at some previous time.
  real(rp) meas_value       ! The value when the data measurement was taken.
  real(rp) ref_value        ! Value when the reference measurement was taken.
  real(rp) correction_value ! Value determined by a fit to correct the lattice.
  real(rp) high_lim         ! High limit for the model_value.
  real(rp) low_lim          ! Low limit for the model_value.
  real(rp) step             ! Sets what is a small step for varying this var.
  real(rp) weight           ! Weight for the merit function term.
  real(rp) delta_merit      ! Diff used to calculate the merit function term.
  real(rp) merit            ! merit_term = weight * delta^2.
  real(rp) dMerit_dVar      ! Merit derivative.     
  real(rp) key_val0         ! Key base value
  real(rp) key_delta        ! Change in value when a key is pressed.
  real(rp) s                ! longitudinal position of ele.
  character(40) merit_type  ! 'target' or 'limit'
  logical exists            ! See above
  logical good_var          ! See above
  logical good_user         ! See above
  logical good_opt          ! See above
  logical good_plot         ! See above
  logical useit_opt         ! See above
  logical useit_plot        ! See above
  type (tao_v1_var_struct), pointer :: v1 => null() ! Pointer to the parent.
end type tao_var_struct  

! A v1_var_struct represents, say, all the quadrupole power supplies.
! The v1_var_struct has a pointer to a section in the s%var array. 

type tao_v1_var_struct
  character(40) :: name = ' '  ! Eg: 'quad_k1'
  integer ix_var0              ! Index of the 0th element in s%var
  integer ix_v1                ! Index to s%v1_var(:) array
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
! Tunnel wall structure

integer, parameter :: arc$ = 1, point$ = 2
integer, parameter :: left_side$ = 1, right_side$ = 2

type tao_wall_point_struct
  integer type              ! arc$ or point$
  real(rp) z, x                 ! floor position
  real(rp) r, theta1, theta2    ! For arcs
end type

type tao_wall_struct
  integer side    ! left_side$ or right_side$
  type (tao_wall_point_struct), allocatable :: point(:)
end type

!------------------------------------------------------------------------
! global parameters that the user has direct access to.
! Also see: tao_common_struct.

type tao_global_struct
  real(rp) :: y_axis_plot_dmin = 1e-4    ! Minimum y_max-y_min allowed for a graph.
  real(rp) :: lm_opt_deriv_reinit = -1   ! Reinit derivative matrix cutoff
  real(rp) :: de_lm_step_ratio = 1       ! Scaling for step sizes between DE and LM optimizers.
  real(rp) :: lmdif_eps = 1e-12          ! tollerance for lmdif optimizer.
  real(rp) :: unstable_penalty = 1e-3    ! Used in unstable_ring datum merit calculation.
  real(rp) :: merit_finish = 1           ! Merit value below which an optimizer will stop.
  real(rp) :: floor_plan_rotation = 0    ! Rotation of floor plan plot: 1.0 -> 360^deg 
  integer :: u_view = 1                  ! Which universe we are viewing.
  integer :: n_opti_cycles = 20          ! number of optimization cycles
  integer :: n_opti_loops = 1            ! number of optimization loops
  integer :: phase_units = radians$      ! Phase units on output.
  integer :: bunch_to_plot = 1           ! Which bunch to plot
  integer :: n_curve_pts = -1            ! For backward compatability
  integer :: random_seed = 0             ! Use system clock by default
  integer :: n_top10 = 10                ! Number of top constraints to print.
  real(rp) :: random_sigma_cutoff = 4    ! cut-off in sigmas.
  character(16) :: random_engine = 'pseudo'         ! Non-beam random number engine
  character(16) :: random_gauss_converter = 'exact' ! Non-beam
  character(16) :: track_type    = 'single'         ! or 'beam'  
  character(16) :: prompt_string = 'Tao'
  character(16) :: optimizer     = 'de'             ! optimizer to use.
  character(16) :: default_key_merit_type 
  character(40) :: print_command = 'lpr'
  character(80) :: var_out_file  = 'var#.out'
  logical :: var_limits_on = .true.             ! Respect the variable limits?
  logical :: auto_scale = .false.               ! Automatically scale and x-scale the plots?
  logical :: opt_with_ref = .false.             ! use reference data in optimization?
  logical :: opt_with_base = .false.            ! use base data in optimization?
  logical :: label_lattice_elements = .true.    ! For lat_layout plots
  logical :: label_keys = .true.                ! For lat_layout plots
  logical :: derivative_recalc = .true.         ! Recalc before each optimizer run?
  logical :: init_plot_needed = .true.          ! reinitialize plotting?
  logical :: orm_analysis = .false.             ! orm using mdsa? 
  logical :: plot_on = .true.                   ! Do plotting?
  logical :: lattice_calc_on = .true.           ! Turn on/off calculations.
  logical :: command_file_print_on = .true.     ! print to terminal when using a cmd file?
  logical :: box_plots = .false.                ! For debugging plot layout issues.
  logical :: beam_timer_on = .false.            ! For timing the beam tracking calculation.
  logical :: optimizer_var_limit_warn = .true.  ! Warn when vars reach a limit with optimization.
end type

!

type tao_alias_struct
  character(40) :: name
  character(100) :: string
end type

type tao_command_file_struct
  character(200) name
  integer :: ix_unit
  character(40) cmd_arg(9)          ! Command file arguments.
  logical :: paused = .false.       ! Is the command file paused?
  integer :: n_line = 0             ! Current line number
end type

! tao_common_struct is for those global parameters that the user 
! should not have direct access to.
! Also see tao_global_struct.

type tao_common_struct
  type (tao_alias_struct) alias(100)
  type (tao_ele_shape_struct), allocatable :: ele_shape_floor_plan(:)
  type (tao_ele_shape_struct), allocatable :: ele_shape_lat_layout(:)
  type (tao_universe_struct), pointer :: u_working          ! Index of working universe.
  type (tao_command_file_struct), allocatable :: cmd_file(:)
  real(rp), allocatable :: covar(:,:), alpha(:,:)
  real(rp) :: dummy_target = 0           ! Dummy varaible
  integer ix_ref_taylor, ix_ele_taylor         ! Taylor map end points
  integer :: n_alias = 0
  integer :: cmd_file_level = 0          ! For nested command files. 0 -> no command file.
  integer :: ix_key_bank = 0             ! For single mode.
  integer :: n_universes = 1   
  logical :: cmd_file_paused
  logical :: use_cmd_here  = .false.                   ! Used for the cmd history stack
  logical :: multi_commands_here = .false.
  logical :: cmd_from_cmd_file = .false.               ! was command from a command file?
  logical :: use_saved_beam_in_tracking = .false.
  logical :: single_mode = .false.
  logical :: lattice_recalc = .true.         ! recalculate the lattice?
  logical :: combine_consecutive_elements_of_like_name
  logical :: common_lattice = .false.      
  character(100) :: cmd                                ! Used for the cmd history
  character(16) :: init_name = 'Tao'                   ! label for initialization
  character(200) :: init_lat_file = ''                 ! '-lat' argument.
  character(100) :: init_tao_file                      ! '-init' argument.
  character(100) :: default_init_tao_file = 'tao.init'          
  character(100) :: beam_all_file = ''  ! Command line input beam data file.
  character(100) :: beam0_file    = ''  ! Command line input beam data file.
  character(16) :: aperture_limit_on
  character(40) :: unique_name_suffix
  character(16) :: valid_plot_who(10)          ! model, base, ref etc...
end type

!------------------------------------------------------------------------
! for connected universes

type tao_connected_uni_struct
  logical connected       ! This universe is injected from another
  logical match_to_design ! match the design lattices
  integer from_uni        ! The universe whose beam injects into this universe
  integer from_uni_ix_ele ! element index where the connection occurs
  real(rp) from_uni_s     ! s position in from_uni where the connection occurs
  type (ele_struct) :: match_ele ! element used to match universes
  type (beam_struct) injecting_beam ! used for beam injection
end type

!-----------------------------------------------------------------------
! This says which datumns to evaluate at this ele
! The point of this is to ave time by not looping through all the data at every
! elements finding which datums need to be evaluated. Instead, do the searching
! beforehand and just keep a log of where to evaluate.

  type tao_ix_data_struct
    ! list of all datums evaluated at this ele
    integer, allocatable :: ix_datum(:)
  endtype

!-----------------------------------------------------------------------

type tao_lat_mode_struct
  real(rp) chrom
  real(rp) growth_rate
end type

! The %bunch_params(:) array has a 1-to-1 correspondence with the lattice elements.
! The %bunch_params2(:) array, if used, is for drawing smooth data lines and has 
! a lot more elements than the %bunch_params(:) array

type tao_lattice_branch_struct
  type (bunch_params_struct), allocatable :: bunch_params(:)
  type (coord_struct), allocatable :: orbit(:)
end type

type tao_lattice_struct
  type (lat_struct) lat                           ! lattice structures
  type (tao_lattice_branch_struct), allocatable :: lat_branch(:)
  type (bunch_params_struct), allocatable :: bunch_params2(:)
  type (normal_modes_struct) modes                ! Synchrotron integrals stuff
  type (rad_int_common_struct) rad_int
  type (tao_lat_mode_struct) a, b
  integer n_bunch_params2                          ! bunch_params2 array size.
end type

!-----------------------------------------------------------------------
! tao_element_struct is for saving per-element information.

type tao_element_struct
  type (beam_struct) beam         ! Beam distribution at element.
  logical save_beam               ! Save beam here?
  integer n_particle_lost_here    ! How many particles are lost here.
  integer ixx                     ! Scratch variable
end type

type tao_universe_branch_struct
  type (tao_element_struct), allocatable :: ele(:) ! Per element information
  type (beam_struct) beam0                         ! Beam at the beginning of lattice
  type (beam_init_struct) :: beam_init             ! Beam distrubution
                                                   !  at beginning of lattice
  integer ix_track_start                 ! Element start index of tracking
  integer ix_track_end                   ! Element end index of tracking
  logical :: init_beam0 = .false.        ! Init beam
  character(80) :: beam_all_file = ''  ! Input beam data file for entire lattice.
  character(80) :: beam0_file    = ''  ! Input beam data file at the start of the lattice.
end type

!-----------------------------------------------------------------------
! A universe is a snapshot of a machine

type tao_universe_struct
  type (tao_universe_struct), pointer :: common => null()
  type (tao_lattice_struct), pointer :: model, design, base
  type (tao_universe_branch_struct), pointer :: uni_branch(:) ! Per element information
  type (beam_struct) current_beam                  ! Beam at the current position
  type (tao_connected_uni_struct)   :: connect     ! Connection data put in 'to' uni.
  type (tao_d2_data_struct), allocatable :: d2_data(:)   ! The data types 
  type (tao_data_struct), allocatable :: data(:)         ! Array of all data.
  type (coord_struct) model_orb0                         ! For saving beginning orbit
  type (tao_ix_data_struct), allocatable :: ix_data(:)   ! which data to evaluate at this ele
  real(rp), allocatable :: dModel_dVar(:,:)              ! Derivative matrix.
  character(100) beam_saved_at
  integer ix_uni                         ! Universe index.
  integer n_d2_data_used
  integer n_data_used
  integer ix_rad_int_cache
  logical do_synch_rad_int_calc
  logical do_chrom_calc
  logical is_on                          ! universe turned on
  logical calc_beam_emittance            ! for a lat calculate emittance
  logical universe_recalc                ! Allows for fine control of lattice calculations
  logical :: mat6_recalc_on = .true.     ! calc linear transfer matrix
  logical picked_uni                     ! Scratch logical.
end type

! The super_universe is the structure that holds an array of universes.
! Essentially this and tao_com hold all the information known to the program.

type tao_super_universe_struct
  type (tao_global_struct) global                          ! global variables.
  type (tao_plot_struct), allocatable :: template_plot(:)  ! Templates for the plots.
  type (tao_plot_page_struct) :: plot_page                 ! Defines the plot window.
  type (tao_plot_region_struct), allocatable :: plot_region(:)
  type (tao_v1_var_struct), allocatable :: v1_var(:)       ! The variable types
  type (tao_var_struct), allocatable :: var(:)             ! array of all variables.
  type (tao_universe_struct), allocatable :: u(:)          ! array of universes.
  integer, allocatable :: key(:)
  type (tao_wall_struct), allocatable :: wall(:)
  type (tao_wave_struct) :: wave 
  integer n_var_used
  integer n_v1_var_used
end type

!-----------------------------------------------------------------------
! The grand global scheme

type (tao_super_universe_struct), save, target :: s
type (tao_common_struct), save, target :: tao_com

!-----------------------------------------------------------------------
contains

subroutine tao_lat_equal_tao_lat (lat1, lat2)

implicit none

type (tao_lattice_struct), intent(inout) :: lat1
type (tao_lattice_struct), intent(in) :: lat2
integer ix2

!

lat1%lat          = lat2%lat
lat1%lat_branch   = lat2%lat_branch
lat1%modes        = lat2%modes
lat1%rad_int      = lat2%rad_int
lat1%a            = lat2%a
lat1%b            = lat2%b

if (allocated(lat2%bunch_params2)) then
  ix2 = size(lat2%bunch_params2)
  if (allocated(lat1%bunch_params2)) then
    if (size(lat1%bunch_params2) /= ix2) deallocate(lat1%bunch_params2)
  endif
  if (.not. allocated(lat1%bunch_params2)) allocate(lat1%bunch_params2(ix2))
else
  if (allocated(lat1%bunch_params2)) deallocate(lat1%bunch_params2)
endif

lat1%n_bunch_params2 = lat2%n_bunch_params2

end subroutine

end module
