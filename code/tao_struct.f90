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
use macroparticle_mod, only: macro_init_struct, macro_beam_struct
use macro_utils_mod, only: macro_bunch_params_struct
use beam_def_struct, only: beam_init_struct, beam_struct, bunch_params_struct
use tao_parameters
use rad_int_common, only: rad_int_common_struct

interface assignment (=)
  module procedure tao_lat_equal_tao_lat
end interface


!-----------------------------------------------------------------------
! misc.

  integer, parameter :: n_key_maxx = 200

!----------------------------------------------------------------------

type tao_ele_shape_struct    ! for the element layout plot
  character(40) key_name     ! Element key name
  character(40) ele_name     ! element name
  character(16) shape        ! plot shape
  character(16) color        ! plot color
  real(rp) dy_pix            ! plot vertical height 
  Logical :: draw_name = .true.
  integer key                ! Element key index to match to
end type

type tao_keyboard_struct
  real(rp) val0                            ! Base value
  real(rp) delta                           ! Change in value
  integer ix_var                           ! Index to variable array.
end type

!-----------------------------------------------------------------------
! Plot structures.

type tao_title_struct
  character(100) string      ! title character string.
  character(2) justify       ! Left, Center, or Right justification.
  real(rp) x, y              ! x, y rwt lower left corner
  character(16) units        ! %BOX, POINTS, etc...
  logical draw_it            ! draw the title
end type

type tao_plot_who_struct     ! Who to plot
  character(16) name         ! Eg: 'meas', 'ref', 'model', etc.
  integer sign               ! +1 or -1
end type

! A curve is defined by a set of (x,y) points and the axis parameters.
! for example the horizontal orbit is one curve.

type tao_curve_struct
  character(40) :: name                    ! Name identifying the curve.
  character(40) :: data_source             ! "calculation", "data_array", or "var_array"
  character(40) :: data_type = ' '         ! "orbit.x", etc.
  character(40) :: ele_ref_name            ! Reference element.
  type (tao_graph_struct), pointer :: g    ! pointer to parent graph 
  real(rp), allocatable :: x_line(:)       ! coords for drawing a curve
  real(rp), allocatable :: y_line(:) 
  real(rp), allocatable :: x_symb(:)       ! coords for drawing the symbols
  real(rp), allocatable :: y_symb(:) 
  integer, allocatable :: ix_symb(:)       ! corresponding index in d1_data%d(:) array.
  real(rp) x_axis_scale_factor ! x-axis conversion from internal to plotting units.
  real(rp) y_axis_scale_factor ! y-axis conversion from internal to plotting units.
  type (qp_line_struct) line   ! Line attributes
  type (qp_symbol_struct) symbol ! Symbol attributes
  integer ix_universe          ! universe to take the data from. 0 => use s%global%u_view
  integer symbol_every         ! symbol every how many points.
  integer ix_ele_ref           ! Index in lattice of reference element.
  logical use_y2               ! Use y2 axis?
  logical draw_line            ! draw a line through the data points?
  logical draw_symbols         ! draw a line through the data points?
  logical convert              ! Eg: covert coupling to cbar?
  logical draw_interpolated_curve
end type

! A graph is a collection of overlayed curves with associated graph title, etc.
! For example a graph could contain just the horizontal orbit or could
! contain both overlayed horizontal and vertical orbits.

type tao_graph_struct
  character(40) name           ! Name identifying the graph
  character(40) type           ! "data", "lat_layout", "key_table", "phase_space"
  character(80) title
  character(80) title_suffix 
  character(80) legend(n_legend_maxx) ! Array for holding descriptive info.
  type (qp_point_struct) legend_origin
  type (tao_plot_who_struct) who(n_who_maxx)  ! Who to plot. Eg: Data - Design
  type (qp_axis_struct) y      ! Y-axis attributes.
  type (qp_axis_struct) y2     ! Y-axis attributes.
  type (qp_rect_struct) margin ! margin around the graph.
  type (tao_curve_struct), allocatable :: curve(:)
  type (tao_plot_struct), pointer :: p ! pointer to parent plot
  real(rp) x_min, x_max         ! min and max of floor_plan drawing.
  real(rp) y_min, y_max         ! min and max of floor_plan drawing.
  logical clip                 ! clip plot at graph boundary.
  integer box(4)               ! Defines which box the plot is put in.
  integer ix_universe          ! Used for lat_layout plots.
  logical valid                ! valid if all curve y_dat computed OK.
  logical y2_mirrors_y         ! Y2-axis same as Y-axis?
  logical limited              ! True if at least one data point past graph bounds.
end type

! A plot is collection of graphs.
! For example a plot could contain three graphs. One for Cbar11, 
! One for Cbar12, and one for Cbar22.

type tao_plot_struct
  character(40) :: name = ' '   ! Identifying name
  type (tao_graph_struct), allocatable :: graph(:)
                                ! individual graphs of a plot
  type (qp_axis_struct) x       ! X-axis parameters.
  type (tao_plot_region_struct), pointer :: r ! pointer to parent.
  real(rp) x_divisions          ! Nominal number of x-axis divisions.
  character(16) x_axis_type     ! 'index', 'ele_index', 's', 'none', or 'phase_space'
  logical independent_graphs    ! scale cmd scales graphs independently?
end type

! A region defines a plot and where to position the plot on the plot page
! %location = (x1, x2, y1, y2) gives the plotting region in percent of the 
!   part of the page inside the page_border with respect to the lower left corner.
! Eg: %location = (0.0, 1.0, 0.5, 1.0) gives the top half of the page inside the border.

type tao_plot_region_struct
  character(40) name             ! Eg: 'top', 'bottom'.
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
  type (tao_title_struct) title(2)   ! Titles at top of page.
  type (tao_plot_region_struct), pointer :: region(:) => null() 
  type (qp_rect_struct) border       ! Border around plots edge of page.
  type (tao_ele_shape_struct) ele_shape(20)
  character(80) ps_scale             ! scaling when creating PS files.
  real(rp) size(2)                   ! width and height of window in pixels.
  real(rp) text_height
  real(rp) :: text_scale = 1
  integer id_window                  ! X window id number.
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
! %exists     -- The datum can exist. Non-existant datums can serve 
!                  as place holders in the u%data array.
! %good_meas  -- Set by the routine that reads in a data set. Good_meas may be 
!                  false, say, if a detector amplifyer is overloaded.
! %good_ref   -- Like good_meas this is set for a reference data set.
! %good_user  -- What the user has selected using the use, veto, and restore 
!                  commands.
! %good_opt   -- Convenient way to veto data to use with optimization without 
!                  touching the other logicals.
! %good_plot  -- Conveninet way to veto data to plot without 
!                  touching the other logicals.
! %useit_plot -- Datum is valid for plotting:
!                  = %exists & %good_plot (w/o measured & reference data)
!                  = %exists & %good_plot & %good_user & %good_meas (w/ meas data)
!                  = %exists & %good_plot & %good_user & %good_ref (w/ reference data)
!                  = %exists & %good_plot & %good_user & %good_meas & %good_ref 
!                                                        (w/ measured & reference data)
! %useit_opt  -- Datum is valid for optimizing (minimizing the merit function):
!                  = %exists & %good_meas & %good_user & %good_opt (w/o reference data)
!                  = %exists & %good_meas & %good_user & %good_opt & %good_ref (otherwise)

type tao_data_struct
  character(40) name        ! Datum name. Eg: "X Orbit @ Det 10"
  character(40) ele_name    ! Name of the element in the Lattice corresponding to the datum.
  character(40) ele0_name   ! Name lattice element when there is a range 
  character(40) data_type   ! Type of data: "orbit.x", etc.
  character(16) merit_type  ! Type of constraint: 'target', 'max', 'min', etc.
  integer ix_ele            ! Index of the element in the lattice element array.
  integer ix_ele0           ! Index of lattice elment when there is a range or reference.
  integer ix_ele_merit      ! Index of lattice elment where merit is evaluated.
  integer ix_d1             ! Index number in u%d2_data(i)%d1_data(j)%d(:) array.
  integer ix_data           ! Index of this datum in the u%data(:) array of data_structs.
  integer ix_dModel         ! Row number in the dModel_dVar derivative matrix.
  real(rp) meas_value       ! Measured datum value. 
  real(rp) ref_value        ! Measured datum value from the reference data set.
  real(rp) model_value      ! Datum value as calculated from the model.
  real(rp) design_value     ! What the datum value is in the design lattice.
  real(rp) old_value        ! The model_value at some previous time.
  real(rp) base_value       ! The value as calculated from the base model.
  real(rp) fit_value        ! The value as calculated from a fitting procedure.
  real(rp) delta_merit      ! Diff used to calculate the merit function term 
  real(rp) weight           ! Weight for the merit function term
  real(rp) merit            ! Merit function term value: weight * delta^2
  real(rp) conversion_factor ! Typically used to convert coupling to cbar
  real(rp) s                ! longitudinal position of ele.
  logical exists            ! See above
  logical good_meas         ! See above
  logical good_ref          ! See above
  logical good_user         ! See above
  logical good_opt          ! See above
  logical good_plot         ! See above
  logical useit_plot        ! See above
  logical useit_opt         ! See above
  type (tao_d1_data_struct), pointer :: d1 => null() 
                            ! Pointer to the parent d1_data_struct 
end type tao_data_struct

! A d1_data_struct represents, say, all the horizontal orbit data.
! The d1_data_struct has a pointer to the appropriate section in 
!   the u%data array. 

type tao_d1_data_struct
  character(40) name        ! Eg: "X", etc.
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
  type (tao_d1_data_struct), pointer :: d1(:) => null() ! Points to children 
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

! A tao_real_array_struct is just a pointer to a real number.
! This is used to construct arrays of reals.

type tao_real_array_struct
  real(rp), pointer :: r
end type

type tao_logical_array_struct
  logical, pointer :: l
end type

!-----------------------------------------------------------------------
! The var_struct defined the fundamental variable structure.
! The super_universe_struct will hold an array of var_structs: s%var(:).
!
! %exists     -- The variable exists. Non-existant variables can serve as place
!                  holders in the s%var array.
! %good_var   -- The variable can be varied. Eg: Permanent magnet quads are 
!                  generally considered not to be variables.
! %good_user  -- What the user has selected using the use, veto, and restore 
!                  commands.
! %good_plot  -- Conveninet way to veto variables to plot without 
!                  touching the other logicals.
! %useit_opt  -- Variable is to be used for optimizing:
!                  = %exists & %good_user & %good_opt
! %useit_plot -- Variable value to be plotted:
!                  = %exists & %good_plot

type tao_this_var_struct
  integer ix_uni            ! universe index.
  integer ix_ele            ! Index of element in the u%lattice%ele(:) array.
  real(rp), pointer :: model_ptr => null() ! Pointer to the model value.
  real(rp), pointer :: base_ptr => null()  ! Pointer to the base value.
end type  

type tao_var_struct
  character(40) name        ! Variable name.
  character(40) alias       ! Short alias name.
  character(40) ele_name    ! Associated lattice element name.
  character(40) attrib_name ! Name of the attribute to vary.
  type (tao_this_var_struct), allocatable :: this(:)
  integer ix_v1             ! Index of this var in the s%v1_var(i)%v(:) array.
  integer ix_var            ! Index number of this var in the s%var(:) array.
  integer ix_dvar           ! Column in the dData_dVar derivative matrix.
  real(rp) model_value      ! Model value.
  real(rp) base_value       ! Base value.
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
  real(rp) conversion_factor! Not currently used for anything
  real(rp) s                ! longitudinal position of ele.
  character(16) merit_type  ! 'target' or 'limit'
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
  character(40) :: name = ' '  ! Eg: "quad_k1"
  integer ix_var0              ! Index of the 0th element in s%var
  type (tao_var_struct), pointer :: v(:) => null() 
                               ! Pointer to the appropriate section in s%var.
end type

! A tao_var_array_struct is just a pointer to a tao_var_struct.
! This is used to construct arrays of tao_var_structs.

type tao_var_array_struct
  type (tao_var_struct), pointer :: v
end type


!------------------------------------------------------------------------
! global switches

type tao_global_struct
  real(rp) :: y_axis_plot_dmin = 1e-4    ! Minimum y_max-y_min allowed for a graph.
  real(rp) :: lm_opt_deriv_reinit = -1   ! Reinit derivative matrix cutoff
  real(rp) :: de_lm_step_ratio = 1       ! Scaling for step sizes between DE and LM optimizers.
  real(rp) :: lmdif_eps = 1e-12          ! tollerance for lmdif optimizer.
  integer :: u_view = 1                  ! Which universe we are viewing.
  integer :: n_opti_cycles = 20          ! number of optimization cycles
  integer :: ix_key_bank = 0             ! For single mode.
  integer :: n_key_table_max = 0         ! Maximum key table index.
  integer :: n_lat_layout_label_rows = 1 ! How many rows with a lat_layout
  integer :: phase_units = radians$      ! Phase units on output.
  integer :: bunch_to_plot = 1           ! Which bunch to plot
  integer :: n_curve_pts = 401           ! Number of points for plotting a smooth curve
  integer :: random_seed = 0             ! use system clock by default
  integer :: n_write_file = 0            ! used for indexing 'show write' files
  character(16) :: track_type    = 'single'    ! or 'beam' or 'macro' 
  character(16) :: prompt_string = 'Tao'
  character(16) :: optimizer     = 'de'        ! optimizer to use.
  character(16) :: default_key_merit_type 
  character(80) :: write_file    = 'tao_show.dat'
  character(16) :: valid_plot_who(10)          ! model, base, ref etc...
  character(40) :: print_command = 'awprint'
  character(80) :: init_file     = 'tao.init'  ! used with 'reinitialize' command
  character(80) :: beam_file     = ''          ! 
  character(80) :: var_out_file  = 'var#.out'
  logical :: var_limits_on = .true.      ! Respect the variable limits?
  logical :: plot_on = .true.            ! Do plotting?
  logical :: auto_scale = .false.        ! Automatically scale and x-scale the plots?
  logical :: opt_with_ref = .false.      ! use reference data in optimization?
  logical :: opt_with_base = .false.     ! use base data in optimization?
  logical :: single_mode = .false.
  logical :: optimizer_running 
  logical :: init_opt_wrapper = .true.
  logical :: label_lattice_elements = .true. ! For lat_layout plots
  logical :: label_keys = .true.             ! For lat_layout plots
  logical :: derivative_recalc = .true.      ! Recalc before each optimizer run?
  logical :: lattice_recalc = .true.         ! recalculate the lattice?
  logical :: init_plot_needed = .true.       ! reinitialize plotting?
  logical :: matrix_recalc_on = .true.       ! calc linear transfer matrix
  logical :: save_beam_everywhere = .false.  ! Save the beam info at all elements?
  logical :: use_saved_beam_in_tracking = .false.
end type

!------------------------------------------------------------------------
! for coupling universes

type tao_coupled_uni_struct
  logical coupled   ! This universe us coupled from another
  logical match_to_design ! match the design lattices
  logical use_coupling_ele ! to use the coupling_ele
  integer from_uni   ! The universe whose beam injects into this universe
  integer from_uni_ix_ele ! element index where coupling occurs
  real(rp) from_uni_s ! s position in from_uni where coupling occurs
  type (ele_struct) :: coupling_ele ! element used to match universes
  type (beam_struct) injecting_beam ! used for beam injection
  type (macro_beam_struct) injecting_macro_beam ! used for macroparticle injection
end type

!-----------------------------------------------------------------------
! This says which datumns to evaluate at this ele
! The point of this is to ave time by not looping through all the data at every
! elements finding which datums need to be evaluated. Instead, do the searching
! beforehand and just keep a log of where to evaluate.

  type tao_ix_data_struct
    ! list of all datums evaluated at this ele
    integer, pointer :: ix_datum(:) => null()
  endtype

!-----------------------------------------------------------------------
! Macroparticle beam structures

type tao_macro_beam_struct
  type (macro_beam_struct) beam             ! macroparticle beam
  type (macro_init_struct) macro_init ! macro distribution at beginning of lat
  type (macro_bunch_params_struct) params ! macro bunch parameters for viewed bunch
  logical calc_emittance     ! for a lat calculate emittance
  integer, pointer :: ix_lost(:,:,:) => null()
                                      ! ^ if .ne. -1 then this macro lost at this ele
                                      ! ix_lost(bunch,slice,macro)
end type

!-----------------------------------------------------------------------
! The %bunch_params(:) array has a 1-to-1 correspondence with the lattice elements.
! The %bunch_params2(:) array, if used, is for drawing smooth data lines and has 
! a lot more elements than the %bunch_params(:) array

type tao_lat_mode_struct
  real(rp) chrom
  real(rp) growth_rate
end type

type tao_lattice_struct
  type (lat_struct) lat                           ! lattice structures
  type (coord_struct), allocatable :: orb(:)
  type (normal_modes_struct) modes                ! Synchrotron integrals stuff
  type (rad_int_common_struct) rad_int
  type (tao_lat_mode_struct) a, b
  type (bunch_params_struct), allocatable :: bunch_params(:)
  type (bunch_params_struct), allocatable :: bunch_params2(:)
  integer n_bunch_params2                          ! bunch_params2 array size.
end type

!-----------------------------------------------------------------------
! A universe is a snapshot of a machine

type tao_universe_struct
  type (tao_lattice_struct) model, design, base
  type (beam_struct), allocatable :: beam_at_element(:) ! beam at element.
  type (beam_struct) current_beam                  ! beam at the current position
  type (beam_init_struct) :: beam_init             ! beam distrubution
                                                   !  at beginning of lattice
  type (tao_macro_beam_struct) macro_beam          ! macroparticle beam 
  type (tao_coupled_uni_struct)   :: coupling      !used for coupled lattices
  type (tao_d2_data_struct), pointer :: d2_data(:) => null()  ! The data types 
  type (tao_data_struct), pointer :: data(:) => null()        ! array of all data.
  type (tao_ix_data_struct), pointer :: ix_data(:) ! which data to evaluate at this ele
  real(rp), pointer :: dModel_dVar(:,:) => null()             ! Derivative matrix.
  integer ix_uni                                   ! Universe index.
  integer n_d2_data_used
  integer n_data_used
  integer ix_rad_int_cache
  integer ixx                                      ! scratch variable
  logical do_synch_rad_int_calc
  logical do_chrom_calc
  logical is_on                                    ! universe turned on
  logical calc_beam_emittance                      ! for a lat calculate emittance
end type

! The super_universe is the structure that holds an array of universes.
! Essentially this holds all the information known to the program.

type tao_super_universe_struct
  type (tao_global_struct) global                          ! global variables.
  type (tao_plot_struct) :: template_plot(n_template_maxx) ! Templates for the plots.
  type (tao_plot_page_struct) :: plot_page                 ! Defines the plot window.
  type (tao_v1_var_struct), pointer :: v1_var(:) => null() ! The variable types
  type (tao_var_struct), pointer :: var(:) => null()       ! array of all variables.
  type (tao_universe_struct), pointer :: u(:) => null()    ! array of universes.
  type (tao_keyboard_struct), pointer :: key(:) => null()
  integer n_var_used
  integer n_v1_var_used
end type

!-----------------------------------------------------------------------
! Define common variables

type tao_alias_struct
  character(40) :: name
  character(100) :: string
end type

type tao_common_struct
  type (tao_alias_struct) alias(100)
  logical opti_init        ! init needed?
  logical opti_at_limit    ! Variable at limit?
  character(40) cmd_arg(9) ! Command file arguments.
  character(100) cmd       ! Used for the cmd history
  character(16) :: init_name = "Tao"  !label for initialization
  integer :: n_alias = 0
  integer :: cmd_file_level = 0 ! for nested command files
              ! unit numbers for a command files. 0 -> no command file.
  integer, pointer :: lun_command_file(:) => null() 
  logical :: use_cmd_here  = .false. ! Used for the cmd history stack
  logical cmd_from_cmd_file ! was command from a command file?
end type


type (tao_super_universe_struct), save, target :: s
type (tao_common_struct), save :: tao_com

!-----------------------------------------------------------------------
contains

subroutine tao_lat_equal_tao_lat (lat1, lat2)

  implicit none

  type (tao_lattice_struct), intent(inout) :: lat1
  type (tao_lattice_struct), intent(in) :: lat2

!

  lat1%lat   = lat2%lat
  lat1%orb   = lat2%orb
  lat1%modes = lat2%modes
  lat1%a     = lat2%a
  lat1%b     = lat2%b

end subroutine

end module
