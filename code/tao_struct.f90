!+
! Module tao_struct
!
! Module defining the basic tao structures
!-

module tao_struct

use bmad_struct, only: rp, ring_struct, coord_struct, radians$
use quick_plot, only: qp_line_struct, qp_symbol_struct, qp_axis_struct, qp_rect_struct
use tao_parameters
use tao_hook_mod

!-----------------------------------------------------------------------
! misc.

  integer, parameter :: n_key_maxx = 200

!----------------------------------------------------------------------

type tao_ele_shape_struct    ! for the element layout plot
  character(16) key_name     ! Element key name
  character(16) ele_name     ! element name
  character(16) shape        ! plot shape
  character(16) color        ! plot color
  integer dy_pix             ! plot vertical height 
  Logical :: plot_name = .true.
  integer key                ! Element key index to match to
end type

type tao_keyboard_struct
  real(rp) val0                            ! Base value
  real(rp) delta                           ! Change in value
  real(rp) :: normalizer = 1.0
  integer ix_var                           ! Index to variable array.
end type

type tao_io_struct
  character(40) routine_name
  character(40) severity
  character(100), pointer :: lines(:) => null()
end type

!-----------------------------------------------------------------------
! Plot structures.

type tao_title_struct
  character(100) string      ! title character string.
  character(1) justify       ! Left, Center, or Right justification.
end type

type tao_plot_who_struct     ! Who to plot
  character(16) name         ! Eg: 'meas', 'ref', 'model', etc.
  integer sign               ! +1 or -1
end type

! A curve is defined by a set of (x,y) points and the axis parameters.
! for example the horizontal orbit is one curve.

type tao_curve_struct
  character(16) :: data_source             ! "lat_layout", "data_array"
  character(16) :: data_type = ' '         ! "orbit:x", etc.
  real(rp), pointer :: x_line(:) => null() ! coords for drawing a curve
  real(rp), pointer :: y_line(:) => null()
  real(rp), pointer :: x_symb(:) => null() ! coords for drawing the symbols
  real(rp), pointer :: y_symb(:) => null()
  integer, pointer :: ix_symb(:) => null() ! corresponding index in d1_data%d(:) array.
  real(rp) units_factor        ! conversion from internal to plotting units.
  type (qp_line_struct) line   ! Line attributes
  type (qp_symbol_struct) symbol ! Symbol attributes
  type (tao_curve_hook) hook   ! Custom stuff. Defined in tao_hook.f90.
  integer ix_universe          ! universe to take the data from. 0 => use s%global%u_view
  integer symbol_every         ! symbol every how many points.
  logical use_y2               ! Use y2 axis?
  logical draw_line            ! draw a line through the data points?
  logical limited              ! True if at least one data point past limit.
  logical convert              ! Eg: covert coupling to cbar?
end type

! A graph is a collection of overlayed curves with associated graph title, etc.
! For example a graph could contain just the horizontal orbit or could
! contain both overlayed horizontal and vertical orbits.

type tao_graph_struct
  character(16) name           ! Name identifying the graph
  character(16) type           ! "data", "lat_layout", "key_table"
  character(80) title
  character(80) title_suffix 
  character(80) legend(n_legend_maxx) ! Array for holding descriptive information.
  type (qp_axis_struct) y      ! Y-axis attributes.
  type (qp_axis_struct) y2     ! Y-axis attributes.
  type (qp_rect_struct) margin ! margin around the graph.
  type (tao_curve_struct), pointer :: curve(:) => null()
  type (tao_graph_hook) hook   ! Custom stuff. Defined in tao_hook.f90.
  logical clip                 ! clip plot at graph boundary.
  integer box(4)               ! Defines which box the plot is put in.
  integer ix_universe          ! Used for lat_layout plots.
end type

! A region defines where to position a plot on the plot page
! %location = (x1, x2, y1, y2) gives the plotting region in percent of the 
!   part of the page inside the page_border with respect to the lower left corner.
! Eg: %location = (0.0, 1.0, 0.5, 1.0) gives the top half of the page inside the border.

type tao_plot_region_struct
  character(16) name        ! Eg: 'top', 'bottom'.
  real(rp) location(4)        
end type

! A plot is collection of graphs.
! For example a plot could contain three graphs. One for Cbar11, 
! One for Cbar12, and one for Cbar22.

type tao_plot_struct
  character(16) :: name = ' '           ! Identifying name
  type (tao_plot_region_struct) region  ! Where on the plot page to put the plot
  type (tao_plot_who_struct) who(10)    ! Who to plot. Eg: Data - Design
  type (tao_graph_struct), pointer :: graph(:) => null() 
                                  ! individual graphs of a plot
  type (tao_plot_hook) hook       ! Custom stuff. Defined in tao_hook.f90
  type (qp_axis_struct) x         ! X-axis parameters.
  real(rp) x_divisions            ! Nominal number of x-axis divisions.
  character(16) x_axis_type       ! 'index', 's'
  logical independent_graphs      ! Graph y-axis scales independent when using the scale cmd?
  logical visible                 ! To draw or not to draw.
  logical valid                   ! valid if all curve y_dat computed OK.
end type

! The plot_page defines the whole plotting window. 
! The plot_page contains a collection of plots along with 
! other info (border margins etc.).
! Note that the qp_com structure of quick_plot also is used to hold 
! plot page info.

type tao_plot_page_struct               
  type (tao_title_struct) title(2)            ! Titles at top of page.
  type (tao_plot_struct), pointer :: plot(:) => null() ! Individual plots on a page. 
  type (tao_plot_page_hook) hook     ! Custom stuff. Defined in tao_hook.f90
  type (qp_rect_struct) border       ! Border around plots edge of page.
  type (tao_ele_shape_struct) ele_shape(20)
  character(80) ps_scale             ! scaling when creating PS files.
  real(rp) size(2)                   ! width and height of window in pixels.
  integer id_window                  ! X window id number.
end type

!-----------------------------------------------------------------------
! The data_struct defines the fundamental data structure representing 
! one datum point.
! The universe_struct will hold an array of data_struct structures: u%data(:).
!
! %exists     -- The datum can exist. Non-existant datums can serve 
!                  as place holders in the u%data array.
! %good_data  -- Set by the routine that reads in a data set. Good_data may be 
!                  false, say, if a detector amplifyer is overloaded.
! %good_ref   -- Like good_data this is set for a reference data set.
! %good_user  -- What the user has selected using the use, veto, and restore 
!                  commands.
! %good_opt   -- Convenient way to veto data to use with optimization without 
!                  touching the other logicals.
! %good_plot  -- Conveninet way to veto data to plot without 
!                  touching the other logicals.
! %useit_plot -- Datum is valid for plotting:
!                  = %exists & %good_dat & %good_user & %good_plot (w/o reference data)
!                  = %exists & %good_dat & %good_user & %good_plot & %good_ref (otherwise)
! %useit_opt  -- Datum is valid for optimizing (minimizing the merit function):
!                  = %exists & %good_dat & %good_user & %good_opt (w/o reference data)
!                  = %exists & %good_dat & %good_user & %good_opt & %good_ref (otherwise)

type tao_data_struct
  character(16) name        ! Datum name. Eg: "X Orbit @ Det 10"
  character(16) ele_name    ! Name of the element in the Lattice corresponding to the datum.
  character(16) ele2_name   ! Name lattice element when there is a range 
  character(16) data_type   ! Type of data: "orbit:x", etc.
  character(16) merit_type  ! Type of constraint: 'target', 'max', 'min', etc.
  integer ix_ele            ! Index of the element in the lattice element array.
  integer ix_ele2           ! Index of lattice elment when there is a range.
  integer ix_ele_merit      ! Index of lattice elment where merit is evaluated.
  integer ix_d1             ! Index number of this datum.
  integer ix_data           ! Index of this datum in the u%data(:) array of data_structs.
  integer ix_dModel         ! Row number in the dModel_dVar derivative matrix.
  real(rp) meas_value       ! Measured datum value. 
  real(rp) ref_value        ! Measured datum value from the reference data set.
  real(rp) model_value      ! Datum value as calculated from the model.
  real(rp) design_value     ! What the datum value is in the design lattice.
  real(rp) old_value        ! The model_value at some previous time.
  real(rp) base_value       ! The value as calculated from the base model.
  real(rp) fit_value        ! The value as calculated from a fitting procedure.
  real(rp) delta            ! Diff used to calculate the merit function term 
  real(rp) weight           ! Weight for the merit function term
  real(rp) merit            ! Merit function term value: weight * delta^2
  real(rp) conversion_factor ! Typically used to convert coupling to cbar
  real(rp) s                ! longitudinal position of ele.
  logical exists            ! See above
  logical good_data         ! See above
  logical good_ref          ! See above
  logical good_user         ! See above
  logical good_opt          ! See above
  logical good_plot         ! See above
  logical useit_plot        ! See above
  logical useit_opt         ! See above
  type (tao_data_hook) hook ! Custom stuff. Defined in tao_hook.f90
  type (tao_d1_data_struct), pointer :: d1 => null() 
                            ! Pointer to the parent d1_data_struct 
end type tao_data_struct

! A d1_data_struct represents, say, all the horizontal orbit data.
! The d1_data_struct has a pointer to the appropriate section in 
!   the u%data array. 

type tao_d1_data_struct
  character(16) name        ! Eg: "X", etc.
  integer ix_data           ! index of the 0th element in u%data.
  type (tao_d1_data_hook) hook  ! Custom stuff. Defined in tao_hook.f90
  type (tao_d2_data_struct), pointer :: d2 => null() ! ptr to parent d2_data
  type (tao_data_struct), pointer :: d(:) => null()  
                            ! Pointer to the appropriate section in u%data
end type

! A d2_data_struct represents all of a type of data. Eg: All orbit data.
! The d2_data_struct has pointers to the approprite d1_data_structs
! %ix_data and %ix_ref are used if the external data files are 
!   sequentially numbered.

type tao_d2_data_struct
  character(16) name              ! Name to be used with commands.
  character(200) data_file_name   ! Data file name .
  character(200) ref_file_name    ! Reference file name.
  character(20) data_date         ! Data measurement date.
  character(20) ref_date          ! Reference data measurement date.
  character(80) Descrip(n_descrip_maxx) ! Array for descriptive information.
  type (tao_d2_data_hook) hook    ! Custom stuff. Defined in tao_hook.f90
  type (tao_d1_data_struct), pointer :: d1(:) => null() ! Points to children 
  integer ix_data                 ! Index of the data set.
  integer ix_ref                  ! Index of the reference data set. 
  logical data_read_in            ! A data set has been read in?
  logical ref_read_in             ! A reference data set has been read in?
end type

!-----------------------------------------------------------------------
! The var_struct defined the fundamental variable structure.
! The universe_struct will hold an array of var_struct structures: u%var(:).
!
! %exists     -- The variable exists. Non-existant variables can serve as place
!                  holders in the u%var array.
! %good_var   -- The variable can be varied. Eg: Permanent magnet quads are 
!                  generally considered not to be variables.
! %good_user  -- What the user has selected using the use, veto, and restore 
!                  commands.
!                  touching the other logicals.
! %good_plot  -- Conveninet way to veto variables to plot without 
!                  touching the other logicals.
! %useit_opt  -- Variable is to be used for optimizing:
!                  = %exists & %good_user & %good_opt
! %useit_plot -- Variable value to be plotted:
!                  = %exists & %good_plot

type tao_this_var_struct
  integer ix_uni            ! universe index.
  integer ix_ele            ! Index of element in the u%lattice%ele_(:) array.
  real(rp), pointer :: model_ptr => null() ! Pointer to the model value.
  real(rp), pointer :: base_ptr => null()  ! Pointer to the base value.
end type  

type tao_var_struct
  character(32) name        ! Variable name.
  character(16) alias       ! Short alias name.
  character(16) ele_name    ! Associated lattice element name.
  character(16) attrib_name ! Name of the attribute to vary.
  type (tao_this_var_struct), pointer :: this(:) => null()
  integer ix_v1             ! Index of this var in the v1_var_struct%v array.
  integer ix_var            ! Index number of this var in the u%var array.
  integer ix_dvar           ! Column in the dData_dVar derivative matrix.
  real(rp) model_value      ! Model value.
  real(rp) base_value      ! Model value.
  real(rp) design_value     ! Design value from the design lattice.
  real(rp) old_value        ! The model_value at some previous time.
  real(rp) meas_value       ! The value when the data measurement was taken.
  real(rp) ref_value        ! Value when the reference measurement was taken.
  real(rp) correction_value ! Value determined by a fit to correct the lattice.
  real(rp) high_lim         ! High limit for the model_value.
  real(rp) low_lim          ! Low limit for the model_value.
  real(rp) step             ! Sets what is a small step for varying this var.
  real(rp) weight           ! Weight for the merit function term.
  real(rp) delta            ! Diff used to calculate the merit function term.
  real(rp) merit            ! merit_term = weight * delta^2.
  real(rp) dMerit_dVar      ! Merit derivative.     
  character(16) merit_type  ! 'target' or 'limit'
  logical exists            ! See above
  logical good_var          ! See above
  logical good_user         ! See above
  logical good_opt          ! See above
  logical useit_opt         ! See above
  logical useit_plot        ! See above
  type (tao_var_hook) hook  ! Custom stuff. Defined in tao_hook.f90
  type (tao_v1_var_struct), pointer :: v1 => null() ! Pointer to the parent.
end type tao_var_struct  

! A v1_var_struct represents, say, all the quadrupole power supplies.
! The v1_var_struct has a pointer to a section in the u%var array. 

type tao_v1_var_struct
  character(16) :: name = ' '  ! Eg: "quad_k1"
  integer ix_var0              ! Index of the 0th element in u%var
  type (tao_v1_var_hook) hook  ! Custom stuff. Defined in tao_hook.f90
  type (tao_var_struct), pointer :: v(:) => null() 
                               ! Pointer to the appropriate section in u%var.
end type

!------------------------------------------------------------------------
! global switches

type tao_global_struct
  real(rp) :: y_axis_plot_dmin = 1e-4 
                                     ! Minimum y_max-y_min allowed for a graph.
  integer :: u_view = 1              ! Which universe we are viewing.
  integer :: n_opti_cycles = 20      ! number of optimization cycles
  integer :: lun_command_file = 0    ! unit number for a command file.
                                     !  0 -> no command file.
  integer :: ix_key_bank = 0         ! For single mode.
  integer :: phase_units = radians$  ! Phase units on output.
  character(16) :: prompt_string = 'Tao'
  character(16) :: optimizer = 'de'  ! optimizer to use.
  type (tao_global_hook) hook        ! Custom stuff. Defined in tao_hook.f90
  logical :: var_limits_on = .false. ! Respect the variable limits?
  logical :: plot_on = .true.        ! Do plotting?
  logical :: opt_with_ref = .false.  ! use reference data in optimization?
  logical :: opt_with_base = .false. ! use base data in optimization?
  logical :: single_mode = .false.
  logical :: optimizer_running 
  logical :: init_opt_wrapper = .true.
  logical :: label_lattice_elements = .true. ! For lat_layout plots
  logical :: label_keys = .true.             ! For lat_layout plots
  logical :: derivative_recalc = .true.      ! Recalc before each optimizer run?
  character(16) :: valid_plot_who(10) 
  character(40) :: print_command = 'awprint'
  character(80) :: var_out_file = 'var#.out'
  character(80) :: opt_var_out_file = 'opt_var#.out'
end type

!-----------------------------------------------------------------------
! A universe is a snapshot of a machine

type tao_universe_struct
  type (ring_struct) model, design, base           ! lattice structures
  type (coord_struct), allocatable :: model_orb(:), design_orb(:), base_orb(:)
  type (tao_d2_data_struct), pointer :: d2_data(:) => null()  ! The data types 
  type (tao_data_struct), pointer :: data(:) => null()        ! array of all data.
  real(rp), pointer :: dModel_dVar(:,:)                       ! Derivative matrix.
  integer n_d2_data_used
  integer n_data_used
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
  type (tao_keyboard_struct), pointer :: key(:)
  integer n_var_used
  integer n_v1_var_used
end type

!-----------------------------------------------------------------------
! Define common variables

type tao_alias_struct
  character(16) :: name
  character(100) :: string
end type

type tao_common_struct
  type (tao_alias_struct) alias(100)
  logical opti_init        ! init needed?
  logical opti_at_limit    ! Variable at limit?
  character(40) cmd_arg(9) ! Command file arguments.
  character(100) cmd
  integer :: n_alias = 0
  logical :: use_cmd_here  = .false. ! Used for the cmd history stack
end type


type (tao_super_universe_struct), save :: s
type (tao_common_struct), save :: tao_com

end module
