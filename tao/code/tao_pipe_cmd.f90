!+
! Subroutine tao_pipe_cmd (input_str)
!
! Print information in a form easily parsed by a scripting program like python.
! This routine formally was tao_python_cmd.
!
! Output will be printed to the terminal or written to a file depending upon the switches embedded
! in the input_str string argument. See the routine "end_stuff" below. For a few commands (for
! example, the "pipe lat_list" command), the output can be stored on the tao_c_interface_com%c_integer (for 
! integer output) or tao_c_interface_com%c_real (for real output) arrays for faster processing.
!
! Note: The syntax for "parameter list form" is:
!   {component_name};{type};{can_vary};{component_value(s)}
!
! {type} is the type of the parameter and is one of:
!   INT         ! Integer number
!   INT_ARR     ! Integer array.
!   REAL        ! Real number
!   REAL_ARR    ! Real array
!   COMPLEX     ! Complex number (Re;Im)
!   LOGIC       ! Logical: "T" or "F".
!   INUM        ! Integer whose allowed values can be obtained using the "pipe inum" command.
!   ENUM        ! String whose allowed values can be obtained using the "pipe enum" command.
!   ENUM_ARR    ! Array of enums.
!   FILE        ! Name of file.
!   CRYSTAL     ! Crystal name string. EG: "Si(111)"
!   DAT_TYPE    ! Data type string. EG: "orbit.x"
!   DAT_TYPE_Z  ! Data type string if plot%x_axis_type = 'data'. Otherwise is a data_type_z enum.
!   SPECIES     ! Species name string. EG: "H2SO4++"
!   ELE_PARAM   ! Lattice element parameter string. EG "K1"
!   STR         ! String that does not fall into one of the above string categories.
!   STR_ARR     ! String array
!   STRUCT      ! Structure. In this case {component_value} is of the form:
!                   {name1};{type1};{value1};{name2};{type2};{value2};...
!   COMPONENT   ! For curve component parameters.
!
! {can_vay} indicates if the component can be varied. It is one of:
!   T         ! Can vary
!   F         ! Cannot vary
!   I         ! Ignore (Do not display)
!
! If the {component_name} has a "^" symbol in it: The component is an enum or inum. Example: "graph^type"
! In this case, use the entire string when using "pipe enum" but suppress everything before the "^"
! when displaying the compoent.
!
! Input:
!   input_str  -- Character(*): What to show.
!-

subroutine tao_pipe_cmd (input_str)

use tao_interface, dummy => tao_pipe_cmd
use location_encode_mod, only: location_encode
use twiss_and_track_mod, only: twiss_and_track_at_s
use wall3d_mod, only: calc_wall_radius, wall3d_d_radius
use tao_command_mod, only: tao_next_switch, tao_cmd_split, tao_next_word
use tao_init_data_mod, only: tao_point_d1_to_data
use tao_init_variables_mod, only: tao_point_v1_to_var, tao_var_stuffit2
use tao_c_interface_mod, only: tao_c_interface_com, re_allocate_c_double
use tao_plot_mod, only: tao_set_floor_plan_axis_label
use tao_dmerit_mod, only: tao_dmodel_dvar_calc
use tao_input_struct, only: tao_ele_shape_input, tao_ele_shape_input_to_struct
use opti_de_mod, only: opti_de_param
use rad_6d_mod, only: emit_6d

implicit none

type (tao_universe_struct), pointer :: u
type (tao_d2_data_struct), pointer :: d2_ptr
type (tao_d1_data_struct), pointer :: d1_ptr
type (tao_d2_data_struct), allocatable :: d2_temp(:)
type (tao_d1_data_struct), allocatable :: d1_temp(:)
type (tao_data_struct), pointer :: data, d_ptr
type (tao_data_struct), allocatable :: d_temp(:)
type (tao_data_struct), target :: datum
type (tao_v1_var_array_struct), allocatable, target :: v1_array(:)
type (tao_v1_var_struct), pointer :: v1_ptr
type (tao_var_struct), pointer :: v_ptr, var
type (tao_var_array_struct), allocatable, target :: v_array(:)
type (tao_v1_var_struct), allocatable :: v1_temp(:)
type (tao_var_struct), allocatable :: v_temp(:)
type (tao_plot_array_struct), allocatable :: plots(:)
type (tao_graph_array_struct), allocatable :: graphs(:)
type (tao_curve_array_struct), allocatable :: curves(:)
type (tao_plot_region_struct), pointer :: pr
type (tao_plot_struct), allocatable :: plot_temp(:)
type (tao_plot_struct), pointer :: p
type (tao_graph_struct), pointer :: g
type (tao_graph_struct) :: graph
type (tao_graph_struct), allocatable :: graph_temp(:)
type (tao_curve_struct), allocatable :: curve_temp(:)
type (tao_curve_struct), pointer :: c
type (tao_lattice_struct), pointer :: tao_lat
type (tao_plot_region_struct), pointer :: region
type (tao_d2_data_array_struct), allocatable :: d2_array(:)
type (tao_d1_data_array_struct), allocatable :: d1_array(:)
type (tao_data_array_struct), allocatable :: d_array(:)
type (beam_struct), pointer :: beam
type (beam_init_struct), pointer :: beam_init
type (lat_struct), pointer :: lat
type (bunch_struct), pointer :: bunch
type (ele_struct), pointer :: ele, ele0, ele1, ele2, lord, slave
type (ele_struct), target :: this_ele
type (coord_struct), pointer :: orbit
type (coord_struct), target :: orb, orb_start, orb_end, orb_here
type (bunch_params_struct), pointer :: bunch_params
type (bunch_params_struct), pointer :: bunch_p
type (bunch_track_struct), pointer :: bunch_params_comb(:)
type (bunch_track_struct), pointer :: comb1
type (ele_pointer_struct), allocatable :: eles(:), eles2(:)
type (branch_struct), pointer :: branch, branch2
type (tao_model_branch_struct), pointer :: model_branch
type (random_state_struct) ran_state
type (ele_attribute_struct) attrib
type (ac_kicker_struct), pointer :: ac
type (cartesian_map_struct), pointer :: ct_map
type (cartesian_map_term1_struct), pointer :: ctt
type (cylindrical_map_struct), pointer :: cy_map
type (cylindrical_map_term1_struct), pointer :: cyt
type (em_field_struct) :: field
type (taylor_struct) taylor(6)
type (taylor_term_struct), pointer :: tt
type (floor_position_struct) floor, floor1, floor2, end1, end2, f_orb
type (tao_floor_plan_struct), pointer :: fp
type (wake_struct), pointer :: wake
type (wake_sr_mode_struct), pointer :: wsr
type (wake_lr_mode_struct), pointer :: lr_mode
type (wall3d_struct), pointer :: wall3d
type (wall3d_section_struct), pointer :: sec
type (gen_grad_map_struct), pointer :: gg_map
type (gen_grad1_struct), pointer :: gg
type (twiss_struct), pointer :: twiss_arr(:)
type (em_taylor_term_struct), pointer :: em_tt
type (grid_field_struct), pointer :: g_field
type (grid_field_pt1_struct), pointer :: g_pt
type (tao_drawing_struct), pointer :: drawing
type (tao_shape_pattern_struct), pointer :: pattern
type (tao_shape_pattern_struct), allocatable :: pat_temp(:)
type (tao_shape_pattern_point_struct), allocatable :: pat_pt_temp(:)
type (tao_ele_shape_struct), pointer :: shapes(:)
type (tao_ele_shape_struct), allocatable :: shapes_temp(:)
type (tao_ele_shape_struct), pointer :: shape
type (tao_ele_shape_input) shape_input
type (photon_element_struct), pointer :: ph
type (qp_axis_struct) x_ax, y_ax
type (tao_building_wall_section_struct), pointer :: bws
type (tao_building_wall_section_struct), allocatable :: bws_temp(:)
type (tao_building_wall_point_struct), pointer :: bwp(:)
type (tao_building_wall_point_struct), allocatable :: bwp_temp(:)
type (tao_dynamic_aperture_struct), pointer :: da
type (tao_wave_kick_pt_struct), pointer :: wk
type (tao_model_element_struct), pointer :: tao_ele
type (tao_lattice_branch_struct), pointer :: tao_branch
type (all_pointer_struct) a_ptr
type (control_var1_struct), pointer :: cvar

real(rp) z, s_pos, value, values(40), y1, y2, v_old(3), r_vec(3), dr_vec(3), w_old(3,3), v_vec(3), dv_vec(3)
real(rp) length, angle, cos_t, sin_t, cos_a, sin_a, ang, s_here, z1, z2, rdummy, time1, gamma
real(rp) x_bend(0:400), y_bend(0:400), dx_bend(0:400), dy_bend(0:400), dx_orbit(0:400), dy_orbit(0:400)
real(rp) a(0:n_pole_maxx), b(0:n_pole_maxx), a2(0:n_pole_maxx), b2(0:n_pole_maxx)
real(rp) knl(0:n_pole_maxx), tn(0:n_pole_maxx)
real(rp) mat6(6,6), vec0(6), array(7), perp(3), origin(3), r_wall
real(rp), allocatable :: real_arr(:), value_arr(:)

type (tao_spin_map_struct), pointer :: sm
real(rp) n0(3), l0(3), m0(3), qs, q, xi_sum, xi_diff
complex(rp) eval(6), evec(6,6), n_eigen(6,3)

integer :: i, j, k, ib, id, iv, iv0, ie, ip, is, iu, nn, md, ct, nl2, n, ix, ix2, iu_write, data_type
integer :: ix_ele, ix_ele1, ix_ele2, ix_branch, ix_bunch, ix_d2, n_who, ix_pole_max, attrib_type, loc
integer :: ios, n_loc, ix_line, n_d1, ix_min(20), ix_max(20), n_delta, why_not_free, ix_uni, ix_shape_min
integer line_width, n_bend, ic, num_ele, n_arr, n_add, n1, n2, i0, i1, i2, n_order
integer, allocatable :: index_arr(:), int_arr(:)
integer, target :: nl
integer, pointer :: nl_ptr

logical :: err, print_flag, opened, doprint, free, matched, track_only, use_real_array_buffer, can_vary
logical first_time, found_one, calc_ok, no_slaves, index_order, ok, not
logical, allocatable :: picked(:), logic_arr(:)

character(*) input_str
character(len(input_str)) line
character(n_char_show), allocatable, target :: li(:)
character(n_char_show), pointer :: li_ptr(:)
character(n_char_show) li2
character(300), allocatable :: name_arr(:)
character(200) file_name, all_who, tail_str
character(40) imt, jmt, rmt, lmt, amt, amt2, iamt, vamt, rmt2, ramt, cmt, label_name
character(40) who, max_loc, ele_name, name1(40), name2(40), a_name, name, attrib_name, command
character(40), allocatable :: str_arr(:)
character(40), allocatable :: name_list(:)
character(40) cmd, which, v_str, head, tail
character(40) switch, color, shape_shape
character(1) :: mode(3) = ['a', 'b', 'c']
character(*), parameter :: r_name = 'tao_pipe_cmd'

!

line = input_str
doprint = .true.
opened = .false.
tao_c_interface_com%n_real = 0
tao_c_interface_com%n_int = 0

do
  call tao_next_switch (line, [character(8):: '-append ', '-write', '-noprint'], .false., switch, err)
  if (err) return
  if (switch == '') exit

  select case (switch)
  case ('-noprint')
    doprint = .false.

  case ('-append', '-write')
    call tao_next_word(line, file_name)
    iu_write = lunget()

    if (switch == '-append') then
      open (iu_write, file = file_name, position = 'APPEND', status = 'UNKNOWN', recl = 500)
    else
      open (iu_write, file = file_name, status = 'REPLACE', recl = 500)
    endif

    opened = .true.
  end select
enddo

call string_trim(line, line, ix)
cmd = line(1:ix)
call string_trim(line(ix+1:), line, ix_line)

! Needed:
!   EM field
!   HOM
!   x_axis_type (variable parameter)

call match_word (cmd, [character(40) :: &
          'beam', 'beam_init', 'branch1', 'bunch_comb', 'bunch_params', 'bunch1', 'bmad_com',&
          'building_wall_list', 'building_wall_graph', 'building_wall_point', 'building_wall_section', &
          'constraints', 'da_params', 'da_aperture', &
          'data', 'data_d2_create', 'data_d2_destroy', 'data_d_array', 'data_d1_array', &
          'data_d2', 'data_d2_array', 'data_set_design_value', 'data_parameter', &
          'datum_create', 'datum_has_ele', 'derivative', &
          'ele:ac_kicker', 'ele:cartesian_map', 'ele:chamber_wall', 'ele:control_var', &
          'ele:cylindrical_map', 'ele:elec_multipoles', 'ele:floor', 'ele:gen_attribs', 'ele:gen_grad_map', &
          'ele:grid_field', 'ele:head', 'ele:lord_slave', 'ele:mat6', 'ele:methods', &
          'ele:multipoles', 'ele:orbit', 'ele:param', 'ele:photon', 'ele:spin_taylor', 'ele:taylor', & 
          'ele:twiss', 'ele:wake', 'ele:wall3d', 'em_field', 'enum', 'evaluate', 'floor_plan', 'floor_orbit', &
          'global', 'global:opti_de', 'global:optimization', 'global:ran_state', 'help', 'inum', &
          'lat_branch_list', 'lat_calc_done', 'lat_ele_list', 'lat_header', 'lat_list', 'lat_param_units', &
          'matrix', 'merit', 'orbit_at_s', 'place_buffer', &
          'plot_curve', 'plot_curve_manage', 'plot_graph', 'plot_graph_manage', 'plot_histogram', &
          'plot_lat_layout', 'plot_line', 'plot_list', &
          'plot_symbol', 'plot_template_manage', 'plot_transfer', 'plot1', &
          'ptc_com', 'ring_general', &
          'shape_list', 'shape_manage', 'shape_pattern_list', 'shape_pattern_manage', &
          'shape_pattern_point_manage', 'shape_set', 'show', 'space_charge_com', 'species_to_int', 'species_to_str', &
          'spin_invariant', 'spin_polarization', 'spin_resonance', 'super_universe', &
          'taylor_map', 'twiss_at_s', 'universe', &
          'var_v1_create', 'var_v1_destroy', 'var_create', 'var_general', 'var_v1_array', 'var_v_array', 'var', &
          'wall3d_radius', 'wave'], ix, matched_name = command)

if (ix == 0) then
  call out_io (s_error$, r_name, 'pipe what? ' // quote(cmd) // ' is not recognized.')
  return
endif

if (ix < 0) then
  call out_io (s_error$, r_name, 'pipe what? ' // quote(cmd) // ' is an ambiguous command.')
  return
endif

amt  = '(100a)'
amt2 = '(a, l1, 10a)'
imt  = '(a, 100(i0, a))'
jmt  = '(i0, a, i0)'
rmt  = '(a, 100(es22.14, a))'
ramt = '(a, 100(a, es22.14))'
rmt2 = '(a, l1, a, 100(es22.14, a))'
lmt  = '(a, 100(l1, a))'
vamt = '(a, i0, 3a)'

nl = 0
call re_allocate_lines (li, 200)

li_ptr => li   ! To get around ifort bug
nl_ptr => nl   ! To get around ifort bug

select case (command)

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% beam
!
! Output beam parameters that are not in the beam_init structure.
!
! Notes
! -----
! Command syntax:
!   pipe beam {ix_uni}@{ix_branch}
!
! Where:
!   {ix_uni} is a universe index. Defaults to s%global%default_universe.
!   {ix_branch} is a lattice branch index. Defaults to s%global%default_branch.
!
! Note: To set beam_init parameters use the "set beam" command.
!
! Parameters
! ----------
! ix_uni : optional
! ix_branch : ""
!
! Returns
! -------
! string_list 
!
! Examples
! -------- 
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/csr_beam_tracking/tao.init
!  args:
!    ix_uni: 1
!    ix_branch: 0

case ('beam')

  u => point_to_uni(line, .true., err); if (err) return
  ix_branch = parse_branch(line, u, .false., err); if (err) return

  nl=incr(nl); write (li(nl), lmt) 'always_reinit;LOGIC;T;',           u%beam%always_reinit
  nl=incr(nl); write (li(nl), lmt) 'track_beam_in_universe;LOGIC;F;',  u%beam%track_beam_in_universe
  nl=incr(nl); write (li(nl), amt) 'saved_at;STR;T;',                  trim(u%beam%saved_at)
  nl=incr(nl); write (li(nl), amt) 'dump_at;STR;T;',                   trim(u%beam%dump_at)
  nl=incr(nl); write (li(nl), amt) 'dump_file;STR;T;',                 trim(u%beam%dump_file)
  nl=incr(nl); write (li(nl), amt) 'track_start;STR;T;',               trim(u%model_branch(ix_branch)%beam%track_start)
  nl=incr(nl); write (li(nl), amt) 'track_end;STR;T;',                 trim(u%model_branch(ix_branch)%beam%track_end)
  nl=incr(nl); write (li(nl), rmt) 'comb_ds_save;REAL;T;',             u%model%tao_branch(ix_branch)%comb_ds_save
  if (allocated(u%model%tao_branch(ix_branch)%bunch_params_comb)) then
    nl=incr(nl); write (li(nl), rmt) 'ds_save;REAL;F;',                  u%model%tao_branch(ix_branch)%bunch_params_comb(1)%ds_save
  else
    nl=incr(nl); write (li(nl), rmt) 'ds_save;REAL;F;',                  -1.0_rp
  endif

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% beam_init
!
! Output beam_init parameters.
!
! Notes
! -----
! Command syntax:
!   pipe beam_init {ix_uni}@{ix_branch}
!
! Where:
!   {ix_uni} is a universe index. Defaults to s%global%default_universe.
!   {ix_branch} is a lattice branch index. Defaults to s%global%default_branch.
!
! Note: To set beam_init parameters use the "set beam_init" command
!
! Parameters
! ----------
! ix_uni : optional
! ix_branch : ""
!
! Returns
! -------
! string_list 
!
! Examples
! -------- 
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/csr_beam_tracking/tao.init
!  args:
!    ix_uni: 1
!    ix_branch: 0

case ('beam_init')

  u => point_to_uni(line, .true., err); if (err) return
  ix_branch = parse_branch(line, u, .false., err); if (err) return
  beam_init => u%model_branch(ix_branch)%beam%beam_init

  nl=incr(nl); write (li(nl), amt) 'distribution_type;STR_ARR;T',              (';', trim(beam_init%distribution_type(k)), k = 1, 3)
  nl=incr(nl); write (li(nl), amt) 'position_file;FILE;T;',                    trim(beam_init%position_file)
  nl=incr(nl); write (li(nl), rmt) 'sig_z_jitter;REAL;T;',                     beam_init%sig_z_jitter
  nl=incr(nl); write (li(nl), rmt) 'sig_pz_jitter;REAL;T;',                    beam_init%sig_pz_jitter
  nl=incr(nl); write (li(nl), amt) 'center_jitter;REAL_ARR;T',                 (';', re_str(beam_init%center_jitter(k), 8), k = 1, 6)
  nl=incr(nl); write (li(nl), amt) 'emit_jitter;REAL_ARR;T',                   (';', re_str(beam_init%emit_jitter(k), 8), k = 1, 2)
  nl=incr(nl); write (li(nl), imt) 'n_particle;INT;T;',                        beam_init%n_particle
  nl=incr(nl); write (li(nl), lmt) 'renorm_center;LOGIC;T;',                   beam_init%renorm_center
  nl=incr(nl); write (li(nl), lmt) 'renorm_sigma;LOGIC;T;',                    beam_init%renorm_sigma
  nl=incr(nl); write (li(nl), amt) 'random_engine;ENUM;T;',                    trim(beam_init%random_engine)
  nl=incr(nl); write (li(nl), amt) 'random_gauss_converter;ENUM;T;',           trim(beam_init%random_gauss_converter)
  nl=incr(nl); write (li(nl), rmt) 'random_sigma_cutoff;REAL;T;',              beam_init%random_sigma_cutoff
  nl=incr(nl); write (li(nl), rmt) 'a_norm_emit;REAL;T;',                      beam_init%a_norm_emit
  nl=incr(nl); write (li(nl), rmt) 'b_norm_emit;REAL;T;',                      beam_init%b_norm_emit
  nl=incr(nl); write (li(nl), rmt) 'a_emit;REAL;T;',                           beam_init%a_emit
  nl=incr(nl); write (li(nl), rmt) 'b_emit;REAL;T;',                           beam_init%b_emit
  nl=incr(nl); write (li(nl), rmt) 'dpz_dz;REAL;T;',                           beam_init%dPz_dz
  nl=incr(nl); write (li(nl), rmt) 'dt_bunch;REAL;T;',                         beam_init%dt_bunch
  nl=incr(nl); write (li(nl), rmt) 't_offset;REAL;T;',                         beam_init%t_offset
  nl=incr(nl); write (li(nl), amt) 'center;REAL_ARR;T',                        (';', re_str(beam_init%center(k), 8), k = 1, 6)
  nl=incr(nl); write (li(nl), amt) 'spin;REAL_ARR;T',                          (';', re_str(beam_init%spin(k), 10), k = 1, 3)
  nl=incr(nl); write (li(nl), rmt) 'sig_z;REAL;T;',                            beam_init%sig_z
  nl=incr(nl); write (li(nl), rmt) 'sig_pz;REAL;T;',                           beam_init%sig_pz
  nl=incr(nl); write (li(nl), rmt) 'bunch_charge;REAL;T;',                     beam_init%bunch_charge
  nl=incr(nl); write (li(nl), imt) 'n_bunch;INT;T;',                           beam_init%n_bunch
  nl=incr(nl); write (li(nl), imt) 'ix_turn;INT;T;',                           beam_init%ix_turn
  nl=incr(nl); write (li(nl), amt) 'species;SPECIES;T;',                       trim(beam_init%species)
  nl=incr(nl); write (li(nl), lmt) 'full_6d_coupling_calc;LOGIC;T;',           beam_init%full_6D_coupling_calc
  nl=incr(nl); write (li(nl), lmt) 'use_particle_start;LOGIC;T;',              beam_init%use_particle_start
  nl=incr(nl); write (li(nl), lmt) 'use_t_coords;LOGIC;T;',                    beam_init%use_t_coords
  nl=incr(nl); write (li(nl), lmt) 'use_z_as_t;LOGIC;T;',                      beam_init%use_z_as_t

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% bmad_com
!
! Output bmad_com structure components.
!
! Notes
! -----
! Command syntax:
!   pipe bmad_com
! 
! Returns
! -------
! string_list 
!
! Examples
! -------- 
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:

case ('bmad_com')

  nl=incr(nl); write (li(nl), rmt) 'max_aperture_limit;REAL;T;',                 bmad_com%max_aperture_limit
  nl=incr(nl); write (li(nl), amt) 'd_orb;REAL_ARR;T',                           (';', re_str(bmad_com%d_orb(k), 6), k = 1, 6)
  nl=incr(nl); write (li(nl), rmt) 'default_ds_step;REAL;T;',                    bmad_com%default_ds_step
  nl=incr(nl); write (li(nl), rmt) 'significant_length;REAL;T;',                 bmad_com%significant_length
  nl=incr(nl); write (li(nl), rmt) 'rel_tol_tracking;REAL;T;',                   bmad_com%rel_tol_tracking
  nl=incr(nl); write (li(nl), rmt) 'abs_tol_tracking;REAL;T;',                   bmad_com%abs_tol_tracking
  nl=incr(nl); write (li(nl), rmt) 'rel_tol_adaptive_tracking;REAL;T;',          bmad_com%rel_tol_adaptive_tracking
  nl=incr(nl); write (li(nl), rmt) 'abs_tol_adaptive_tracking;REAL;T;',          bmad_com%abs_tol_adaptive_tracking
  nl=incr(nl); write (li(nl), rmt) 'init_ds_adaptive_tracking;REAL;T;',          bmad_com%init_ds_adaptive_tracking
  nl=incr(nl); write (li(nl), rmt) 'min_ds_adaptive_tracking;REAL;T;',           bmad_com%min_ds_adaptive_tracking
  nl=incr(nl); write (li(nl), rmt) 'fatal_ds_adaptive_tracking;REAL;T;',         bmad_com%fatal_ds_adaptive_tracking
  nl=incr(nl); write (li(nl), rmt) 'autoscale_amp_abs_tol;REAL;T;',              bmad_com%autoscale_amp_abs_tol
  nl=incr(nl); write (li(nl), rmt) 'autoscale_amp_rel_tol;REAL;T;',              bmad_com%autoscale_amp_rel_tol
  nl=incr(nl); write (li(nl), rmt) 'autoscale_phase_tol;REAL;T;',                bmad_com%autoscale_phase_tol
  nl=incr(nl); write (li(nl), rmt) 'electric_dipole_moment;REAL;T;',             bmad_com%electric_dipole_moment
  nl=incr(nl); write (li(nl), rmt) 'synch_rad_scale;REAL;T;',                    bmad_com%synch_rad_scale
  nl=incr(nl); write (li(nl), rmt) 'sad_eps_scale;REAL;T;',                      bmad_com%sad_eps_scale
  nl=incr(nl); write (li(nl), rmt) 'sad_amp_max;REAL;T;',                        bmad_com%sad_amp_max
  nl=incr(nl); write (li(nl), imt) 'sad_n_div_max;INT;T;',                       bmad_com%sad_n_div_max
  nl=incr(nl); write (li(nl), imt) 'taylor_order;INT;T;',                        bmad_com%taylor_order
  nl=incr(nl); write (li(nl), imt) 'runge_kutta_order;INT;T;',                   bmad_com%runge_kutta_order
  nl=incr(nl); write (li(nl), imt) 'default_integ_order;INT;T;',                 bmad_com%default_integ_order
  nl=incr(nl); write (li(nl), imt) 'max_num_runge_kutta_step;INT;T;',            bmad_com%max_num_runge_kutta_step
  nl=incr(nl); write (li(nl), lmt) 'rf_phase_below_transition_ref;LOGIC;T;',     bmad_com%rf_phase_below_transition_ref
  nl=incr(nl); write (li(nl), lmt) 'sr_wakes_on;LOGIC;T;',                       bmad_com%sr_wakes_on
  nl=incr(nl); write (li(nl), lmt) 'lr_wakes_on;LOGIC;T;',                       bmad_com%lr_wakes_on
  nl=incr(nl); write (li(nl), lmt) 'auto_bookkeeper;LOGIC;T;',                   bmad_com%auto_bookkeeper
  nl=incr(nl); write (li(nl), lmt) 'csr_and_space_charge_on;LOGIC;T;',           bmad_com%csr_and_space_charge_on
  nl=incr(nl); write (li(nl), lmt) 'spin_tracking_on;LOGIC;T;',                  bmad_com%spin_tracking_on
  nl=incr(nl); write (li(nl), lmt) 'spin_sokolov_ternov_flipping_on;LOGIC;T;',   bmad_com%spin_sokolov_ternov_flipping_on
  nl=incr(nl); write (li(nl), lmt) 'radiation_damping_on;LOGIC;T;',              bmad_com%radiation_damping_on
  nl=incr(nl); write (li(nl), lmt) 'radiation_fluctuations_on;LOGIC;T;',         bmad_com%radiation_fluctuations_on
  nl=incr(nl); write (li(nl), lmt) 'conserve_taylor_maps;LOGIC;T;',              bmad_com%conserve_taylor_maps
  nl=incr(nl); write (li(nl), lmt) 'absolute_time_tracking;LOGIC;T;',            bmad_com%absolute_time_tracking
  nl=incr(nl); write (li(nl), lmt) 'convert_to_kinetic_momentum;LOGIC;T;',       bmad_com%convert_to_kinetic_momentum
  nl=incr(nl); write (li(nl), lmt) 'aperture_limit_on;LOGIC;T;',                 bmad_com%aperture_limit_on
  nl=incr(nl); write (li(nl), lmt) 'debug;LOGIC;T;',                             bmad_com%debug

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% branch1
!
! Output lattice branch information for a particular lattice branch.
!
! Notes
! -----
! Command syntax:
!   pipe branch1 {ix_uni}@{ix_branch}
!
! Where:
!   {ix_uni} is a universe index. Defaults to s%global%default_universe.
!   {ix_branch} is a lattice branch index. Defaults to s%global%default_branch.
!
! Parameters
! ----------
! ix_uni : ""
! ix_branch : ""
! 
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!    ix_uni: 1
!    ix_branch: 0

case ('branch1')

  u => point_to_uni(line, .true., err); if (err) return
  ix_branch = parse_branch(line, u, .false., err); if (err) return
  branch => u%model%lat%branch(ix_branch)

  nl=incr(nl); write (li(nl), amt) 'name;STR;F;',                               trim(branch%name)
  nl=incr(nl); write (li(nl), imt) 'ix_branch;INT;F;',                          branch%ix_branch
  nl=incr(nl); write (li(nl), imt) 'ix_from_branch;INT;F;',                     branch%ix_from_branch
  nl=incr(nl); write (li(nl), imt) 'ix_from_ele;INT;F;',                        branch%ix_from_ele

  nl=incr(nl); write (li(nl), rmt) 'param.n_part;REAL;F;',                      branch%param%n_part
  nl=incr(nl); write (li(nl), rmt) 'param.total_length;REAL;F;',                branch%param%total_length
  nl=incr(nl); write (li(nl), rmt) 'param.unstable_factor;REAL;F;',             branch%param%unstable_factor
  nl=incr(nl); write (li(nl), amt) 'param.particle;SPECIES;T;',                 trim(species_name(branch%param%particle))
  nl=incr(nl); write (li(nl), amt) 'param.default_tracking_species;SPECIES;T;', trim(species_name(branch%param%default_tracking_species))
  nl=incr(nl); write (li(nl), amt) 'param.geometry;ENUM;T;',                    trim(geometry_name(branch%param%geometry))
  nl=incr(nl); write (li(nl), lmt) 'param.stable;LOGIC;F;',                     branch%param%stable

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% bunch_comb
!
! Outputs bunch parameters at a comb point. 
! Also see the "write bunch_comb" and "show bunch -comb" commands.
!
! Notes
! -----
! Command syntax:
!   pipe bunch_comb {flags} {who} {ix_uni}@{ix_branch} {ix_bunch}
!
! Where:
!   {flags} are optional switches:
!       -array_out : If present, the output will be available in the 
!              tao_c_interface_com%c_real array.
!   {ix_uni} is a universe index. Defaults to s%global%default_universe.
!   {ix_branch} is a branch index. Defaults to s%global%default_branch.
!   {ix_bunch} is the bunch index. Defaults to 1.
!   {who} is one of:
!       x, px, y, py, z, pz, t, s, spin.x, spin.y, spin.z, p0c, beta     -- centroid 
!       x.Q, y.Q, z.Q, a.Q, b.Q, c.Q where Q is one of: beta, alpha, gamma, phi, 
!                                       eta, etap, sigma, sigma_p, emit, norm_emit
!     sigma.IJ where I, J in range [1,6]
!     rel_min.I, rel_max.I where I in range [1,6]
!     charge_live, n_particle_live, n_particle_lost_in_ele, ix_ele
!
!   Note: If ix_uni or ix_branch is present, "@" must be present.
!
! Example:
!   pipe bunch_comb py 2@1 1
!
! Parameters
! ----------
! who
! ix_uni : optional
! ix_branch : optional
! ix_bunch : default=1
! flags : default=-array_out
!
! Returns
! -------
! string_list
!   if '-array_out' not in flags
! real_array
!   if '-array_out' in flags
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/csr_beam_tracking/tao.init
!  args:
!    who: x.beta

case ('bunch_comb')

  use_real_array_buffer = .false.

  if (index('-array_out', line(1:ix_line)) == 1) then
    call string_trim(line(ix_line+1:), line, ix_line)
    use_real_array_buffer = .true.
  endif

  u => s%u(s%global%default_universe)
  ix_branch = s%global%default_branch
  ix_bunch = 1

  call string_trim (line, line, ix)
  which = line(:ix)
  call string_trim(line(ix+1:), line, ix)

  if (ix > 0) then
    if (index(line(1:ix), '@') /= 0) then
      head = line(1:ix)
      call string_trim(line(ix+1:), line, ix)
      u => point_to_uni(head, .true., err); if (err) return
      ix_branch = parse_branch(head, u, .false., err); if (err) return
    endif
  endif

  branch => u%model%lat%branch(ix_branch)

  if (.not. allocated(u%model%tao_branch(ix_branch)%bunch_params_comb)) then
    call invalid ('COMB ARRAY NOT ALLOCATED. PROBABLY CAUSED BY NO BUNCH TRACKING.')
    return
  endif
  bunch_params_comb => u%model%tao_branch(ix_branch)%bunch_params_comb

  if (ix > 0) then
    ix_bunch = parse_int(line, err, 1, size(bunch_params_comb), 1)
  endif

  comb1 => bunch_params_comb(ix_bunch)

  if (comb1%ds_save < 0) then
    call invalid ('COMB_DS_SAVE NOT POSITIVE.')
    return
  endif

  n = comb1%n_pt
  if (n < 0) then
    call invalid ('COMB POINTS NOT CALCULATED.')
    return
  endif

  !

  ix = index(which, '.')
  if (ix == 0) then
    head = which
  else
    head = which(1:ix)
    tail = which(ix+1:)
  endif

  select case (head)
  case ('x', 'px', 'y', 'py', 'z', 'pz')
    call match_word(which, coord_name, ix, .true., .true.)
    call real_array_out(comb1%pt%centroid%vec(ix), use_real_array_buffer, 0, n)

  case ('spin.')
    call match_word(tail, ['x', 'y', 'z'], ix, .true., .true.)
    if (ix < 0) then
      call invalid ('"WHO" NOT RECOGNIZED: ' // which)
      return
    endif
    call real_array_out(comb1%pt%centroid%spin(ix), use_real_array_buffer, 0, n)

  case ('x.', 'y.', 'z.', 'a.', 'b.', 'c.')
    select case (head)
    case ('x.');  twiss_arr => comb1%pt%x
    case ('y.');  twiss_arr => comb1%pt%y
    case ('z.');  twiss_arr => comb1%pt%z
    case ('a.');  twiss_arr => comb1%pt%a
    case ('b.');  twiss_arr => comb1%pt%b
    case ('c.');  twiss_arr => comb1%pt%c
    end select

    select case (tail)
    case ('beta');      call real_array_out(twiss_arr%beta, use_real_array_buffer, 0, n)
    case ('alpha');     call real_array_out(twiss_arr%alpha, use_real_array_buffer, 0, n)
    case ('gamma');     call real_array_out(twiss_arr%gamma, use_real_array_buffer, 0, n)
    case ('phi');       call real_array_out(twiss_arr%phi, use_real_array_buffer, 0, n)
    case ('eta');       call real_array_out(twiss_arr%eta, use_real_array_buffer, 0, n)
    case ('etap');      call real_array_out(twiss_arr%etap, use_real_array_buffer, 0, n)
    case ('sigma');     call real_array_out(twiss_arr%sigma, use_real_array_buffer, 0, n)
    case ('sigma_p');   call real_array_out(twiss_arr%sigma_p, use_real_array_buffer, 0, n)
    case ('emit');      call real_array_out(twiss_arr%emit, use_real_array_buffer, 0, n)
    case ('norm_emit'); call real_array_out(twiss_arr%norm_emit, use_real_array_buffer, 0, n)
    case default
      call invalid ('Bad {who}: ' // which)
      return
    end select

  case ('sigma.')
    i = parse_int(tail(1:1), err, 1, 6);   if (err) return
    j = parse_int(tail(2:2), err, 1, 6);   if (err) return
    call real_array_out (comb1%pt%sigma(i,j), use_real_array_buffer, 0, n)

  case ('rel_min.')
    i = parse_int(tail(1:1), err, 1, 6);   if (err) return
    call real_array_out (comb1%pt%rel_min(i), use_real_array_buffer, 0, n)

  case ('rel_max.')
    i = parse_int(tail(1:1), err, 1, 6);   if (err) return
    call real_array_out (comb1%pt%rel_max(i), use_real_array_buffer, 0, n)

  case ('s');                       call real_array_out(comb1%pt%centroid%s, use_real_array_buffer, 0, n)
  case ('t');                       call real_array_out(real(comb1%pt%centroid%t, rp), use_real_array_buffer, 0, n)
  case ('p0c');                     call real_array_out(comb1%pt%centroid%p0c, use_real_array_buffer, 0, n)
  case ('beta');                    call real_array_out(comb1%pt%centroid%beta, use_real_array_buffer, 0, n)
  case ('charge_live');             call real_array_out(comb1%pt%charge_live, use_real_array_buffer, 0, n)
  case ('n_particle_live');         call real_array_out(1.0_rp*comb1%pt%n_particle_live, use_real_array_buffer, 0, n)
  case ('n_particle_lost_in_ele');  call real_array_out(1.0_rp*comb1%pt%n_particle_lost_in_ele, use_real_array_buffer, 0, n)
  case ('ix_ele');                  call real_array_out(1.0_rp*comb1%pt%ix_ele, use_real_array_buffer, 0, n)
  case default
    call invalid ('Bad {who}: ' // which)
    return
  end select

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% bunch_params
!
! Outputs bunch parameters at the exit end of a given lattice element.
!
! Notes
! -----
! Command syntax:
!   pipe bunch_params {ele_id}|{which}
!
! Where:
!   {ele_id} is an element name or index.
!   {which} is one of: "model", "base" or "design"
!
! Example:
!   pipe bunch_params end|model  ! parameters at model lattice element named "end".
!
! Parameters
! ----------
! ele_id
! which : default=model
!
! Returns
! -------
! string_list 
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/csr_beam_tracking/tao.init
!  args:
!    ele_id: end
!    which: model

case ('bunch_params')

  u => point_to_uni(line, .true., err); if (err) return
  tao_lat => point_to_tao_lat(line, u, err, which, tail_str); if (err) return
  ele => point_to_ele(line, tao_lat%lat, err); if (err) return

  bunch_params => tao_lat%tao_branch(ele%ix_branch)%bunch_params(ele%ix_ele)

  call bunch_params_out(bunch_params)

  beam => u%model_branch(ele%ix_branch)%ele(ele%ix_ele)%beam
  nl=incr(nl); write (li(nl), lmt) 'beam_saved;LOGIC;T;', allocated(beam%bunch)

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% bunch1
!
! Outputs Bunch parameters at the exit end of a given lattice element.
!
! Notes
! -----
! Command syntax:
!   pipe bunch1 {ele_id}|{which} {ix_bunch} {coordinate}
!
! Where:
!   {ele_id} is an element name or index.
!   {which} is one of: "model", "base" or "design"
!   {ix_bunch} is the bunch index.
!   {coordinate} is one of: x, px, y, py, z, pz, "s", "t", "charge", "p0c", 
!                                                                 "state", "ix_ele"
!
! For example, if {coordinate} = "px", the phase space px coordinate of each particle
! of the bunch is displayed. The "state" of a particle is an integer. 
! A value of 1 means alive and any other value means the particle has been lost.
!
! Parameters
! ----------
! ele_id
! coordinate
! which : default=model
! ix_bunch : default=1
!
! Returns
! -------
! real_array
!   if coordinate in ['x', 'px', 'y', 'py', 'z', 'pz', 's', 't', 'charge', 'p0c']
! integer_array
!   if coordinate in ['state', 'ix_ele']
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/csr_beam_tracking/tao.init
!  args:
!    ele_id: end
!    coordinate: x
!    which: model
!    ix_bunch: 1

case ('bunch1')

  u => point_to_uni(line, .true., err); if (err) return
  tao_lat => point_to_tao_lat(line, u, err, which, tail_str); if (err) return
  ele => point_to_ele(line, tao_lat%lat, err); if (err) return

  beam => u%model_branch(ele%ix_branch)%ele(ele%ix_ele)%beam
  if (.not. allocated(beam%bunch)) then
    call invalid ('BEAM NOT SAVED AT ELEMENT.')
    return
  endif

  ix_bunch = parse_int(tail_str, err, 1, size(beam%bunch)); if (err) return

  call coord_out(beam%bunch(ix_bunch), tail_str)

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% building_wall_list
!
! Output List of building wall sections or section points
!
! Notes
! -----
! Command syntax:
!   pipe building_wall_list {ix_section}
!
! Where:
!   {ix_section} is a building wall section index.
!
! If {ix_section} is not present, a list of building wall sections is given.
! If {ix_section} is present, a list of section points is given.
! 
! Parameters
! ----------
! ix_section : optional
!
! Returns
! -------
! string_list 
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_wall
!  args:
!    ix_section:
!
! Example: 2
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_wall
!  args:
!    ix_section: 1

case ('building_wall_list')

  if (line == '') then
    do ib = 1, size(s%building_wall%section)
      shape => tao_pointer_to_building_wall_shape (s%building_wall%section(ib)%name)
      if (associated(shape)) then
        nl=incr(nl); write (li(nl), '(i0, 9a, i0)') ib, ';', trim(s%building_wall%section(ib)%name),';', &
              trim(s%building_wall%section(ib)%constraint), ';', trim(shape%shape), ';', trim(shape%color), ';', shape%line_width
      else
        nl=incr(nl); write (li(nl), '(i0, 5a)') ib, ';', trim(s%building_wall%section(ib)%name),';', &
              trim(s%building_wall%section(ib)%constraint), ';;;'
      endif
    enddo

  else
    ib = parse_int (line, err, 1, size(s%building_wall%section)); if (err) return
    if (allocated(s%building_wall%section(ib)%point)) then
      bwp => s%building_wall%section(ib)%point
      do ip = 1, size(bwp)
        nl=incr(nl); write (li(nl), '(i0, 10a)') ip, ';', re_str(bwp(ip)%z, 6), ';', re_str(bwp(ip)%x, 6), ';',  re_str(bwp(ip)%radius, 6), ';',  &
                                                                           re_str(bwp(ip)%z_center, 6), ';',  re_str(bwp(ip)%x_center, 6)
      enddo
    endif
  endif

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% building_wall_graph
!
! Output (x, y) points for drawing the building wall for a particular graph.
!
! Notes
! -----
! Command syntax:
!   pipe building_wall_graph {graph}
!
! Where:
!   {graph} is a plot region graph name.
!
! Note: The graph defines the coordinate system for the (x, y) points.
!
! Parameters
! ----------
! graph
!
! Returns
! -------
! string_list 
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_wall
!  args:
!    graph: floor_plan.g

case ('building_wall_graph')

  call tao_find_plots (err, line(1:ix_line), 'REGION', graph = graphs, only_visible = .false.)
  call string_trim(line(ix_line+1:), line, ix_line)

  if (err .or. size(graphs) /= 1) then
    call invalid ('Bad graph name')
    return
  endif

  g => graphs(1)%g
  u => tao_pointer_to_universe(g%ix_universe)
  lat => u%model%lat

  if (.not. allocated(s%building_wall%section)) then
    call invalid ('No building wall defined')
    return
  endif

  do ib = 1, size(s%building_wall%section)
    bwp => s%building_wall%section(ib)%point
    do j = 1, size(bwp)
      call tao_floor_to_screen (g, [bwp(j)%x, 0.0_rp, bwp(j)%z], end1%r(1), end1%r(2))
      nl=incr(nl); write (li(nl), '(2(i0,a), 3(es14.6, a))') ib, ';', j, ';', end1%r(1), ';', end1%r(2), ';', bwp(j)%radius
    enddo
  enddo

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% building_wall_point
!
! add or delete a building wall point
!
! Notes
! -----
! Command syntax:
!   pipe building_wall_point {ix_section}^^{ix_point}^^{z}^^{x}^^{radius}^^
!                                                               {z_center}^^{x_center}
!
! Where:
!   {ix_section}    -- Section index.
!   {ix_point}      -- Point index. Points of higher indexes will be moved up 
!                        if adding a point and down if deleting.
!   {z}, etc...     -- See tao_building_wall_point_struct components.
!                   -- If {z} is set to "delete" then delete the point.
! 
! Parameters
! ----------
! ix_section
! ix_point
! z
! x
! radius
! z_center
! x_center
!
! Returns
! -------
! None 
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_wall
!  args:
!    ix_section: 1
!    ix_point: 1
!    z: 0
!    x: 0
!    radius: 0
!    z_center: 0
!    x_center: 0

case ('building_wall_point')

  call split_this_line (line, name1, 7, err); if (err) return

  is = parse_int(name1(1), err, 1, size(s%building_wall%section));  if (err) return
  bws => s%building_wall%section(is)
  if (allocated(bws%point)) then
    n = size(bws%point)
  else
    n = 0
  endif

  select case (name1(3))
  case ('delete')
    ip = parse_int(name1(2), err, 1, n)
    call move_alloc(bws%point, bwp_temp)
    allocate (bws%point(n-1))
    bws%point(1:ip-1) = bwp_temp(1:ip-1)
    bws%point(ip:) = bwp_temp(ip+1:)

  case default
    ip = parse_int(name1(2), err, 1, n+1)
    if (allocated(bws%point)) then
      call move_alloc(bws%point, bwp_temp)
      allocate (bws%point(n+1))
      bws%point(1:ip-1) = bwp_temp(1:ip-1)
      bws%point(ip+1:) = bwp_temp(ip:)
    else
      allocate (bws%point(n+1))  ! n = 0 here and ip = 1
    endif

    bws%point(ip)%z        = parse_real(name1(3), err);  if (err) return
    bws%point(ip)%x        = parse_real(name1(4), err);  if (err) return
    bws%point(ip)%radius   = parse_real(name1(5), err);  if (err) return
    bws%point(ip)%z_center = parse_real(name1(6), err);  if (err) return
    bws%point(ip)%x_center = parse_real(name1(7), err);  if (err) return
  end select

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% building_wall_section
!
! Add or delete a building wall section
!
! Notes
! -----
! Command syntax:
!   pipe building_wall_section {ix_section}^^{sec_name}^^{sec_constraint}
!
! Where:
!   {ix_section}      -- Section index. Sections with higher indexes will be
!                          moved up if adding a section and down if deleting.
!   {sec_name}        -- Section name.
!   {sec_constraint}  -- A section constraint name or "delete". Must be one of:
!       delete          -- Delete section. Anything else will add the section.
!       none
!       left_side
!       right_side
! 
! Parameters
! ----------
! ix_section
! sec_name
! sec_constraint
! 
! Returns
! -------
! None
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!    ix_section: 1
!    sec_name: test
!    sec_constraint: none


case ('building_wall_section')

  n = size(s%building_wall%section)
  call split_this_line (line, name1, 3, err); if (err) return

  select case (name1(3))
  case ('delete')
    is = parse_int(name1(1), err, 1, n);  if (err) return
    call move_alloc(s%building_wall%section, bws_temp)
    allocate (s%building_wall%section(n-1))
    s%building_wall%section(1:is-1) = bws_temp(1:is-1)
    s%building_wall%section(is:) = bws_temp(is+1:)

  case default
    is = parse_int(name1(1), err, 1, n+1);  if (err) return
    call move_alloc(s%building_wall%section, bws_temp)
    allocate (s%building_wall%section(n+1))
    s%building_wall%section(1:is-1) = bws_temp(1:is-1)
    s%building_wall%section(is+1:) = bws_temp(is:)

    bws => s%building_wall%section(is)
    bws%name       = name1(2)
    bws%constraint = name1(3)
  end select

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% constraints
!
! Output optimization data and variable parameters that contribute to the merit function.
!
! Notes
! -----
! Command syntax:
!   pipe constraints {who}
!
! Where:
!   {who} is one of: "data" or "var"
!
! Data constraints output is:
!   data name
!   constraint type
!   evaluation element name
!   start element name
!   end/reference element name
!   measured value
!   ref value (only relavent if global%opt_with_ref = T)
!   model value
!   base value (only relavent if global%opt_with_base = T)
!   weight
!   merit value
!   location where merit is evaluated (if there is a range)
! Var constraints output is:
!   var name
!   Associated varible attribute
!   meas value
!   ref value (only relavent if global%opt_with_ref = T)
!   model value
!   base value (only relavent if global%opt_with_base = T)
!   weight
!   merit value
!   dmerit/dvar
!
! Parameters
! ----------
! who
!
! Returns
! -------
! string_list 
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_optics_matching
!  args:
!    who: data
!
! Example: 2
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!    who:var

case ('constraints')

  select case (line)
  case ('data')
    do i = lbound(s%u, 1), ubound(s%u, 1)
      do j = 1, size(s%u(i)%data)
        data => s%u(i)%data(j)
        if (.not. data%useit_opt) cycle

        ie = data%ix_ele_merit
        a_name = ''
        if (ie >= 0) then
          branch => s%u(i)%model%lat%branch(data%ix_branch)
          a_name = branch%ele(ie)%name
        endif

        nl=incr(nl); write (li(nl), '(10a, 6(es22.14, a), 2a)') trim(tao_datum_name(data)), ';', &
            trim(tao_constraint_type_name(data)), ';', &
            trim(data%ele_name), ';', trim(data%ele_start_name), ';', trim(data%ele_ref_name), ';', &
            data%meas_value, ';', data%ref_value, ';', data%model_value, ';', data%base_value, ';', &
            data%weight, ';', data%merit, ';', a_name
      enddo
    enddo

  case ('var')
    do i = 1, s%n_var_used
      var => s%var(i)
      if (.not. var%useit_opt) cycle
      nl=incr(nl); write (li(nl), '(4a, 7(es22.14, a))') trim(tao_var1_name(var)), ';', &
            trim(tao_var_attrib_name(var)), ';', &
            var%meas_value, ';', var%ref_value, ';', var%model_value, ';', var%base_value, ';', &
            var%weight, var%merit, var%dmerit_dvar
    enddo

  case default
    call invalid ('Bad {who}')
    return
  end select

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% da_aperture
!
! Output dynamic aperture data
!
! Notes
! -----
! Command syntax:
!   pipe da_aperture {ix_uni}
!
! Where:
!   {ix_uni} is a universe index. Defaults to s%global%default_universe.
!
! Parameters
! ----------
! ix_uni : optional
!
! Returns
! -------
! string_list
!
! Examples
! --------

case ('da_aperture')

  u => point_to_uni(line, .false., err); if (err) return
  da => u%dynamic_aperture

  if (.not. allocated(da%scan)) then
    call invalid ('Scan not done.')
    return
  endif

  do i = 1, size(da%scan)
    do j = 1, size(da%scan(i)%point)
      nl=incr(nl); write (li(nl), '(2(i0, a), 2(es14.6, a))') i, ';', j, ';', &
                                                       da%scan(i)%point(j)%x, ';', da%scan(i)%point(j)%y
    enddo
  enddo

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% da_params
!
! Output dynamic aperture input parameters
!
! Notes
! -----
! Command syntax:
!   pipe da_params {ix_uni}
!
! Where:
!   {ix_uni} is a universe index. Defaults to s%global%default_universe.
! 
! Parameters
! ----------
! ix_uni : optional
!
! Returns
! -------
! string_list
!
! Examples
! --------

case ('da_params')

  u => point_to_uni(line, .false., err); if (err) return
  da => u%dynamic_aperture

  if (.not. allocated(da%pz)) then
    call invalid('No pz points set.')
    return
  endif

  nl=incr(nl); write (li(nl), rmt)  'min_angle;REAL;T;',    da%param%min_angle
  nl=incr(nl); write (li(nl), rmt)  'max_angle;REAL;T;',    da%param%max_angle
  nl=incr(nl); write (li(nl), rmt)  'n_angle;REAL;T;',      da%param%n_angle
  nl=incr(nl); write (li(nl), rmt)  'n_turn;REAL;T;',       da%param%n_turn
  nl=incr(nl); write (li(nl), rmt)  'x_init;REAL;T;',       da%param%x_init
  nl=incr(nl); write (li(nl), rmt)  'y_init;REAL;T;',       da%param%y_init
  nl=incr(nl); write (li(nl), rmt)  'rel_accuracy;REAL;T;', da%param%rel_accuracy
  nl=incr(nl); write (li(nl), rmt)  'abs_accuracy;REAL;T;', da%param%abs_accuracy
  nl=incr(nl); write (li(nl), ramt) 'pz;REAL_ARR;T',        (';', da%pz(i), i = 1, size(da%pz))

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% data
!
! Output Individual datum parameters.
!
! Notes
! -----
! Command syntax:
!   pipe data {ix_uni}@{d2_name}.{d1_name}[{dat_index}]
!
! Where:
!   {ix_uni} is a universe index. Defaults to s%global%default_universe.
!   {d2_name} is the name of the d2_data structure the datum is in.
!   {d1_datum} is the name of the d1_data structure the datum is in.
!   {dat_index} is the index of the datum.
!
! Use the "pipe data-d1" command to get detailed info on a specific d1 array.
!
! Example:
!   pipe data 1@orbit.x[10]
! 
! Parameters
! ----------
! d2_name
! d1_name
! ix_uni : optional
! dat_index : default=1
!
! Returns
! -------
! string_list
!
! Examples
! --------
!
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_optics_matching
!  args:
!    ix_uni:
!    d2_name: twiss
!    d1_name: end 
!    dat_index: 1  
!
! Example: 2
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_optics_matching
!  args:
!    ix_uni: 1
!    d2_name: twiss
!    d1_name: end
!    dat_index: 1

case ('data')

  call tao_find_data (err, line, d_array = d_array)

  if (size(d_array) == 0) then
    call invalid ('Not a valid datum name.')
    return
  endif

  d_ptr => d_array(1)%d
  ix_uni = d_ptr%d1%d2%ix_universe

  nl=incr(nl); write (li(nl), amt) 'ele_name;STR;T;',                         trim(d_ptr%ele_name)
  nl=incr(nl); write (li(nl), amt) 'ele_start_name;STR;T;',                   trim(d_ptr%ele_start_name)
  nl=incr(nl); write (li(nl), amt) 'ele_ref_name;STR;T;',                     trim(d_ptr%ele_ref_name)
  nl=incr(nl); write (li(nl), amt) 'data_type;DAT_TYPE;T;',                   trim(d_ptr%data_type)
  nl=incr(nl); write (li(nl), amt) 'data^merit_type;ENUM;T;',                 trim(d_ptr%merit_type)
  nl=incr(nl); write (li(nl), amt) 'data_source;ENUM;T;',                     trim(d_ptr%data_source)
  nl=incr(nl); write (li(nl), amt) 'eval_point;ENUM;T;',                      trim(anchor_pt_name(d_ptr%eval_point))
  nl=incr(nl); write (li(nl), jmt) ix_uni, '^ix_bunch;INUM;T;',               d_ptr%ix_bunch
  nl=incr(nl); write (li(nl), jmt) ix_uni, '^ix_branch;INUM;T;',              d_ptr%ix_branch
  nl=incr(nl); write (li(nl), imt) 'ix_ele;INT;I;',                           d_ptr%ix_ele
  nl=incr(nl); write (li(nl), imt) 'ix_ele_start;INT;I;',                     d_ptr%ix_ele_start
  nl=incr(nl); write (li(nl), imt) 'ix_ele_ref;INT;I;',                       d_ptr%ix_ele_ref
  nl=incr(nl); write (li(nl), imt) 'ix_ele_merit;INT;F;',                     d_ptr%ix_ele_merit
  nl=incr(nl); write (li(nl), imt) 'ix_d1;INT;F;',                            d_ptr%ix_d1
  nl=incr(nl); write (li(nl), imt) 'ix_data;INT;F;',                          d_ptr%ix_data
  nl=incr(nl); write (li(nl), imt) 'ix_dmodel;INT;F;',                        d_ptr%ix_dModel
  nl=incr(nl); write (li(nl), rmt) 'meas_value;REAL;T;',                      d_ptr%meas_value
  nl=incr(nl); write (li(nl), rmt) 'ref_value;REAL;T;',                       d_ptr%ref_value
  nl=incr(nl); write (li(nl), rmt) 'model_value;REAL;F;',                     d_ptr%model_value
  nl=incr(nl); write (li(nl), rmt) 'design_value;REAL;F;',                    d_ptr%design_value
  nl=incr(nl); write (li(nl), rmt) 'old_value;REAL;F;',                       d_ptr%old_value
  nl=incr(nl); write (li(nl), rmt) 'base_value;REAL;F;',                      d_ptr%base_value
  nl=incr(nl); write (li(nl), rmt) 'error_rms;REAL;T;',                       d_ptr%error_rms
  nl=incr(nl); write (li(nl), rmt) 'delta_merit;REAL;F;',                     d_ptr%delta_merit
  nl=incr(nl); write (li(nl), rmt) 'weight;REAL;T;',                          d_ptr%weight
  nl=incr(nl); write (li(nl), rmt) 'invalid_value;REAL;T;',                   d_ptr%invalid_value
  nl=incr(nl); write (li(nl), rmt) 'merit;REAL;F;',                           d_ptr%merit
  nl=incr(nl); write (li(nl), rmt) 's;REAL;F;',                               d_ptr%s
  nl=incr(nl); write (li(nl), rmt) 's_offset;REAL;F;',                        d_ptr%s_offset
  nl=incr(nl); write (li(nl), lmt) 'exists;LOGIC;F;',                         d_ptr%exists
  nl=incr(nl); write (li(nl), lmt) 'good_model;LOGIC;F;',                     d_ptr%good_model
  nl=incr(nl); write (li(nl), lmt) 'good_base;LOGIC;F;',                      d_ptr%good_base
  nl=incr(nl); write (li(nl), lmt) 'good_design;LOGIC;F;',                    d_ptr%good_design
  nl=incr(nl); write (li(nl), lmt) 'good_meas;LOGIC;T;',                      d_ptr%good_meas
  nl=incr(nl); write (li(nl), lmt) 'good_ref;LOGIC;T;',                       d_ptr%good_ref
  nl=incr(nl); write (li(nl), lmt) 'good_user;LOGIC;T;',                      d_ptr%good_user
  nl=incr(nl); write (li(nl), lmt) 'good_opt;LOGIC;F;',                       d_ptr%good_opt
  nl=incr(nl); write (li(nl), lmt) 'good_plot;LOGIC;T;',                      d_ptr%good_plot
  nl=incr(nl); write (li(nl), lmt) 'useit_plot;LOGIC;F;',                     d_ptr%useit_plot
  nl=incr(nl); write (li(nl), lmt) 'useit_opt;LOGIC;F;',                      d_ptr%useit_opt

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% data_d_array
!
! Output list of datums for a given d1_data structure.
!
! Notes
! -----
! Command syntax:
!   pipe data_d_array {ix_uni}@{d2_name}.{d1_name}
!
! Where:
!   {ix_uni} is a universe index. Defaults to s%global%default_universe.
!   {d2_name} is the name of the containing d2_data structure.
!   {d1_name} is the name of the d1_data structure containing the array of datums.
!
! Example:
!   pipe data_d_array 1@orbit.x
! 
! Parameters
! ----------
! d2_name
! d1_name
! ix_uni : optional
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_optics_matching
!  args:
!    ix_uni: 1 
!    d2_name: twiss
!    d1_name: end


case ('data_d_array')


  call tao_find_data (err, line, d_array = d_array)

  if (.not. allocated(d_array)) then
    call invalid ('Not a valid d1_datum name.')
    return
  endif

  do i = 1, size(d_array)
    d_ptr => d_array(i)%d
    name = tao_constraint_type_name(d_ptr)
    nl=incr(nl); write(li(nl), '(i0, 11a, 3(es22.14, a), 3(l1, a), es22.14, a, l1)') d_ptr%ix_d1, ';', &
              trim(d_ptr%data_type), ';', trim(d_ptr%merit_type), ';', &
              trim(d_ptr%ele_ref_name), ';', trim(d_ptr%ele_start_name), ';', trim(d_ptr%ele_name), ';', &
              d_ptr%meas_value, ';', d_ptr%model_value, ';', d_ptr%design_value, ';', &
              d_ptr%useit_opt, ';', d_ptr%useit_plot, ';', d_ptr%good_user, ';', d_ptr%weight, ';', d_ptr%exists
  enddo


!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% data_d1_array
!
! Output list of d1 arrays for a given data_d2.
!
! Notes
! -----
! Command syntax:
!   pipe data_d1_array {d2_datum}
!
! {d2_datum} should be of the form
!   {ix_uni}@{d2_datum_name}
! 
! Parameters
! ----------
! d2_datum
! ix_uni : optional
!
! Returns
! -------
! string_list 
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_optics_matching
!  args:
!    ix_uni: 1 
!    d2_datum: twiss

case ('data_d1_array')

  call tao_find_data (err, line, d2_array = d2_array)

  if (err .or. .not. allocated(d2_array)) then
    call invalid ('Not a valid d2 data name')
    return
  endif

  d2_ptr => d2_array(1)%d2
  do i = lbound(d2_ptr%d1, 1), ubound(d2_ptr%d1, 1)
    d1_ptr => d2_ptr%d1(i)
    call location_encode (line, d1_ptr%d%useit_opt, d1_ptr%d%exists, lbound(d1_ptr%d, 1))
    nl=incr(nl); write (li(nl), '(a, i0, 5a, i0, a, i0)') 'd1[', i, '];STR2;F;', trim(d1_ptr%name), ';', trim(line), ';', &
                                                                                     lbound(d1_ptr%d, 1), ';', ubound(d1_ptr%d, 1)
  enddo

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% data_d2
!
! Output information on a d2_datum.
!
! Notes
! -----
! Command syntax:
!   pipe data_d2 {ix_uni}@{d2_name}
!
! Where:
!   {ix_uni} is a universe index. Defaults to s%global%default_universe.
!   {d2_name} is the name of the d2_data structure.
!
! Parameters
! ----------
! d2_name
! ix_uni : optional
!
! Returns
! -------
! string_list 
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_optics_matching
!  args:
!    ix_uni: 1
!    d2_name: twiss

case ('data_d2')

  call tao_find_data (err, line, d2_array = d2_array)

  if (err .or. .not. allocated(d2_array)) then
    call invalid ('Not a valid d2 data name')
    return
  endif

  d2_ptr => d2_array(1)%d2

  nl=incr(nl); write (li(nl), imt) 'n_d1;INT;F;',                             size(d2_ptr%d1)
  nl=incr(nl); write (li(nl), imt) 'ix_d2_data;INT;F;',                       d2_ptr%ix_d2_data
  nl=incr(nl); write (li(nl), amt) 'name;STR;T;',                             trim(d2_ptr%name)
  nl=incr(nl); write (li(nl), amt) 'data_file_name;FILE;F;',                  trim(d2_ptr%data_file_name)
  nl=incr(nl); write (li(nl), amt) 'ref_file_name;FILE;F;',                   trim(d2_ptr%ref_file_name)
  nl=incr(nl); write (li(nl), amt) 'data_date;STR;T;',                        trim(d2_ptr%data_date)
  nl=incr(nl); write (li(nl), amt) 'ref_date;STR;T;',                         trim(d2_ptr%ref_date)
  nl=incr(nl); write (li(nl), imt) 'ix_universe;INUM;T;',                     d2_ptr%ix_universe
  nl=incr(nl); write (li(nl), imt) 'ix_ref;INT;F;',                           d2_ptr%ix_ref
  nl=incr(nl); write (li(nl), lmt) 'data_read_in;LOGIC;F;',                   d2_ptr%data_read_in
  nl=incr(nl); write (li(nl), lmt) 'ref_read_in;LOGIC;F;',                    d2_ptr%ref_read_in

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% data_d2_array
!
! Output data d2 info for a given universe.
!
! Notes
! -----
! Command syntax:
!   pipe data_d2_array {ix_uni}
!
! Where:
!   {ix_uni} is a universe index. Defaults to s%global%default_universe.
!
! Example:
!   pipe data_d2_array 1
! 
! Parameters
! ----------
! ix_uni
! 
! Returns
! -------
! string_list  
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!    ix_uni : 1 

case ('data_d2_array')

  u => point_to_uni(line, .false., err); if (err) return

  do i = 1, u%n_d2_data_used
    d2_ptr => u%d2_data(i)
    if (d2_ptr%name == '') cycle
    nl=incr(nl); write (li(nl), '(a)') d2_ptr%name
  enddo

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% data_d2_create
!
! Create a d2 data structure along with associated d1 and data arrays.
!
! Notes
! -----
! Command syntax:
!   pipe data_d2_create {ix_uni}@{d2_name}^^{n_d1_data}^^{d_data_arrays_name_min_max}
!
! Where:
!   {ix_uni} is a universe index. Defaults to s%global%default_universe.
!   {d2_name} is the name of the d2_data structure to create.
!   {n_d1_data} is the number of associated d1 data structures.
!   {d_data_arrays_name_min_max} has the form
!     {name1}^^{lower_bound1}^^{upper_bound1}^^....
!                                            ^^{nameN}^^{lower_boundN}^^{upper_boundN}
!   where {name} is the data array name and 
!   {lower_bound} and {upper_bound} are the bounds of the array.
! 
! Example:
!   pipe data_d2_create 2@orbit^^2^^x^^0^^45^^y^^1^^47
! This example creates a d2 data structure called "orbit" with 
! two d1 structures called "x" and "y".
! 
! The "x" d1 structure has an associated data array with indexes in the range [0, 45].
! The "y" d1 structure has an associated data arrray with indexes in the range [1, 47].
! 
! Use the "set data" command to set created datum parameters.
!
! Note: When setting multiple data parameters, 
!       temporarily toggle s%global%lattice_calc_on to False
!   ("set global lattice_calc_on = F") to prevent Tao trying to 
!       evaluate the partially created datum and generating unwanted error messages.
! 
! Parameters
! ----------
! d2_name
! n_d1_data
! d_data_arrays_name_min_max
! ix_uni : optional
!    
! Returns
! -------
! None
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_optics_matching
!  args:
!    ix_uni: 1
!    d2_name: orbit
!    n_d1_data: 2 
!    d_data_arrays_name_min_max: x^^0^^45^^y^^1^^47

case ('data_d2_create')

  call split_this_line (line, name1, -1, err);  if (err) return

  if (.not. is_integer(name1(2))) then
    call invalid ('Number of d1 arrays missing or invalid')
    return
  endif
  read (name1(2), *) n_d1

  do i = 1, n_d1
    j = 3 * i
    if (.not. is_integer(name1(j+1)) .or. .not. is_integer(name1(j+2))) then
      call invalid('Malformed data parameters: ' // trim(name1(j)) // '^^' // trim(name1(j+1)) // '^^' // trim(name1(j+2)))
      return
    endif
    name2(i) = name1(j)
    read (name1(j+1), *) ix_min(i)
    read (name1(j+2), *) ix_max(i)
  enddo

  if (name1(j+3) /= '') then
    call invalid ('Extra stuff on line: ' // name1(j+3))
    return
  endif

  ! Now create the d2 structure

  name = name1(1)
  u => point_to_uni(name, .true., err); if (err) return

  call tao_find_data(err, name1(1), d2_array, print_err = .false.)
  if (size(d2_array) /= 0) then
    call destroy_this_data_d2 (name1(1))
    call out_io (s_warn$, r_name, '"pipe ' // trim(input_str) // '": Data with this name already exists.', &
                                   'This old data has been destroyed to make room for the new data.')
  endif


  if (allocated(u%d2_data)) then
    n2 = size(u%d2_data)
    if (u%n_d2_data_used + 1 > n2) then
      call move_alloc(u%d2_data, d2_temp)
      allocate (u%d2_data(n2+1))
      u%d2_data(1:n2) = d2_temp
    endif
  else
    allocate (u%d2_data(1))
  endif

  n_delta = sum(ix_max(1:n_d1)) - sum(ix_min(1:n_d1)) + n_d1

  if (allocated(u%d2_data)) then
    n = size(u%data)
    if (u%n_data_used + n_delta > n) then
      call move_alloc(u%data, d_temp)
      allocate (u%data(u%n_data_used + n_delta))
      u%data(1:u%n_data_used) = d_temp(1:u%n_data_used)
      do i = 1, size(u%data)
        u%data(i)%ix_data = i
        u%data(i)%ix_uni = u%ix_uni
      enddo
    endif
  else
    allocate (u%data(n_delta))
  endif

  i2 = 0   ! In case no d2 structures have yet been defined.

  do i = 1, u%n_d2_data_used
    n1 = size(u%d2_data(i)%d1)
    do j = 1, n1
      d1_ptr => u%d2_data(i)%d1(j)
      d1_ptr%d2 => u%d2_data(i)
      i1 = lbound(d1_ptr%d, 1)
      i1 = d1_ptr%d(i1)%ix_data
      i2 = ubound(d1_ptr%d, 1)
      i2 = d1_ptr%d(i2)%ix_data
      call tao_point_d1_to_data (d1_ptr, u%data(i1:i2), u%data(i1)%ix_d1)
    enddo
  enddo

  u%n_data_used = u%n_data_used + n_delta
  nn = u%n_d2_data_used + 1
  u%n_d2_data_used = nn
  u%d2_data(nn)%ix_d2_data = nn
  u%d2_data(nn)%name = name
  u%d2_data(nn)%ix_universe = u%ix_uni
  if (allocated(u%d2_data(nn)%d1)) deallocate(u%d2_data(nn)%d1) ! Can happen if data has been destroyed.
  allocate (u%d2_data(nn)%d1(n_d1))

  do j = 1, n_d1
    d1_ptr => u%d2_data(nn)%d1(j)
    d1_ptr%d2 => u%d2_data(nn)
    d1_ptr%name = name2(j)
    i1 = i2 + 1
    i2 = i2 + 1 + ix_max(j) - ix_min(j)
    call tao_point_d1_to_data (d1_ptr, u%data(i1:i2), ix_min(j))
    ! Load some sensible defaults so they don't need to be set manually
    u%data(i1:i2)%merit_type = 'target'
    u%data(i1:i2)%data_source = 'lat'
  enddo

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% data_d2_destroy
!
! Destroy a d2 data structure along with associated d1 and data arrays.
!
! Notes
! -----
! Command syntax:
!   pipe data_d2_destroy {ix_uni}@{d2_name}
!
! Where:
!   {ix_uni} is a universe index. Defaults to s%global%default_universe.
!   {d2_name} is the name of the d2_data structure to destroy.
!
! Example:
!   pipe data_d2_destroy 2@orbit
! This destroys the orbit d2_data structure in universe 2.
! 
! Parameters
! ----------
! d2_name
! ix_uni : optional
!    
! Returns
! -------
! None
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!    d2_name: orbit

case ('data_d2_destroy')

call destroy_this_data_d2(line)

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% data_parameter
!
! Output an array of values for a particular datum parameter for a given array of datums, 
!
! Notes
! -----
! Command syntax:
!   pipe data_parameter {data_array} {parameter}
!
! {parameter} may be any tao_data_struct parameter.
! Example:
!   pipe data_parameter orbit.x model_value
!
! Parameters
! ----------
! data_array
! parameter
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_optics_matching
!  args:
!    data_array: twiss.end 
!    parameter: model_value


case ('data_parameter')

  call split_this_line (line, name1, 2, err, n, .true.); if (err) return
  call tao_find_data (err, name1(1), d_array = d_array)

  n = size(d_array)
  if (n == 0) then
    call invalid ('Not a valid datum name.')
    return
  endif

  select case (name1(2))
  case ('data_type')
    allocate (name_arr(n))

  case ('ele_name', 'ele_start_name', 'ele_ref_name', 'merit_type', 'id', 'data_source')
    allocate (str_arr(n))

  case ('ix_uni', 'ix_bunch', 'ix_branch', 'ix_ele', 'ix_ele_start', 'ix_ele_ref', &
        'ix_ele_merit', 'ix_d1', 'ix_data', 'ix_dModel', 'eval_point')
    allocate (int_arr(n))

  case ('meas_value', 'ref_value', 'model_value', 'design_value', 'old_value', 'base_value', &
        'error_rms', 'delta_merit', 'weight', 'invalid_value', 'merit', 's', 's_offset')
    allocate (real_arr(n))

  case ('err_message_printed', 'exists', 'good_model', 'good_base', 'good_design', 'good_meas', &
        'good_ref', 'good_user', 'good_opt', 'good_plot', 'useit_plot', 'useit_opt')
    allocate (logic_arr(n))
  end select

  !

  do j = 1, n
    select case (name1(2))
    case ('data_type');            name_arr(j) = d_array(j)%d%data_type

    case ('ele_name');             str_arr(j) = d_array(j)%d%ele_name
    case ('ele_start_name');       str_arr(j) = d_array(j)%d%ele_start_name
    case ('ele_ref_name');         str_arr(j) = d_array(j)%d%ele_ref_name
    case ('merit_type');           str_arr(j) = d_array(j)%d%merit_type
    case ('id');                   str_arr(j) = d_array(j)%d%id
    case ('data_source');          str_arr(j) = d_array(j)%d%data_source

    case ('ix_uni');               int_arr(j) = d_array(j)%d%ix_uni
    case ('ix_bunch');             int_arr(j) = d_array(j)%d%ix_bunch
    case ('ix_branch');            int_arr(j) = d_array(j)%d%ix_branch
    case ('ix_ele');               int_arr(j) = d_array(j)%d%ix_ele
    case ('ix_ele_start');         int_arr(j) = d_array(j)%d%ix_ele_start
    case ('ix_ele_ref');           int_arr(j) = d_array(j)%d%ix_ele_ref
    case ('ix_ele_merit');         int_arr(j) = d_array(j)%d%ix_ele_merit
    case ('ix_d1');                int_arr(j) = d_array(j)%d%ix_d1
    case ('ix_data');              int_arr(j) = d_array(j)%d%ix_data
    case ('ix_dModel');            int_arr(j) = d_array(j)%d%ix_dModel
    case ('eval_point');           int_arr(j) = d_array(j)%d%eval_point

    case ('meas_value');           real_arr(j) = d_array(j)%d%meas_value
    case ('ref_value');            real_arr(j) = d_array(j)%d%ref_value
    case ('model_value');          real_arr(j) = d_array(j)%d%model_value
    case ('design_value');         real_arr(j) = d_array(j)%d%design_value
    case ('old_value');            real_arr(j) = d_array(j)%d%old_value
    case ('base_value');           real_arr(j) = d_array(j)%d%base_value
    case ('error_rms');            real_arr(j) = d_array(j)%d%error_rms
    case ('delta_merit');          real_arr(j) = d_array(j)%d%delta_merit
    case ('weight');               real_arr(j) = d_array(j)%d%weight
    case ('invalid_value');        real_arr(j) = d_array(j)%d%invalid_value
    case ('merit');                real_arr(j) = d_array(j)%d%merit
    case ('s');                    real_arr(j) = d_array(j)%d%s
    case ('s_offset');             real_arr(j) = d_array(j)%d%s_offset

    case ('err_message_printed');  logic_arr(j) = d_array(j)%d%err_message_printed
    case ('exists');               logic_arr(j) = d_array(j)%d%exists
    case ('good_model');           logic_arr(j) = d_array(j)%d%good_model
    case ('good_base');            logic_arr(j) = d_array(j)%d%good_base
    case ('good_design');          logic_arr(j) = d_array(j)%d%good_design
    case ('good_meas');            logic_arr(j) = d_array(j)%d%good_meas
    case ('good_ref');             logic_arr(j) = d_array(j)%d%good_ref
    case ('good_user');            logic_arr(j) = d_array(j)%d%good_user
    case ('good_opt');             logic_arr(j) = d_array(j)%d%good_opt
    case ('good_plot');            logic_arr(j) = d_array(j)%d%good_plot
    case ('useit_plot');           logic_arr(j) = d_array(j)%d%useit_plot
    case ('useit_opt');            logic_arr(j) = d_array(j)%d%useit_opt
    end select
  enddo

  if (allocated(name_arr)) then
    do j = 1, n
      nl=incr(nl); write(li(nl), '(i0, 2a)') j, ';', trim(name_arr(j))
    enddo

  elseif (allocated(str_arr)) then
    do j = 0, (n-1)/10
      iv0 = 10 * j + 1
      nl=incr(nl); write(li(nl), '(i0, 30a)') iv0, (';', trim(str_arr(i)), i = iv0, min(iv0+9, n))
    enddo

  elseif (allocated(int_arr)) then
    do j = 0, (n-1)/10
      iv0 = 10 * j + 1
      nl=incr(nl); write(li(nl), '(i0, 20(a, i0))') iv0, (';', int_arr(i), i = iv0, min(iv0+9, n))
    enddo

  elseif (allocated(real_arr)) then
    do j = 0, (n-1)/10
      iv0 = 10 * j + 1
      nl=incr(nl); write(li(nl), '(i0, 30a)') iv0, (';', re_str(real_arr(i), 10), i = iv0, min(iv0+9, n))
    enddo

  elseif (allocated(logic_arr)) then
    do j = 0, (n-1)/10
      iv0 = 10 * j + 1
      nl=incr(nl); write(li(nl), '(i0, 20(a, l1))') iv0, (';', logic_arr(i), i = iv0, min(iv0+9, n))
    enddo
  endif

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% data_set_design_value
!
! Set the design (and base & model) values for all datums.
!
! Notes
! -----
! Command syntax:
!   pipe data_set_design_value
!
! Example:
!   pipe data_set_design_value
! 
! Note: Use the "data_d2_create" and "datum_create" first to create datums.
! 
! Returns
! -------
! None
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_optics_matching
!  args:

case ('data_set_design_value')

  do i = lbound(s%u, 1), ubound(s%u, 1)
    u => s%u(i)
    u%scratch_lat = u%model%lat
    u%model%lat = u%design%lat
  enddo

  s%u%calc%lattice = .true.
  call tao_lattice_calc (calc_ok)

  do i = lbound(s%u, 1), ubound(s%u, 1)
    u => s%u(i)
    u%design%tao_branch = u%model%tao_branch
    u%data%design_value = u%data%model_value
    u%data%good_design  = u%data%good_model

    u%model%lat = u%base%lat
  enddo

  s%u%calc%lattice = .true.
  call tao_lattice_calc (calc_ok)

  do i = lbound(s%u, 1), ubound(s%u, 1)
    u => s%u(i)
    u%base%tao_branch = u%model%tao_branch
    u%data%base_value = u%data%model_value
    u%data%good_base  = u%data%good_model

    u%model%lat = u%scratch_lat
  enddo

  s%u%calc%lattice = .true.
  call tao_lattice_calc (calc_ok)

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% datum_create
!
! Create a datum.
!
! Notes
! -----
! Command syntax:
!   pipe datum_create {datum_name}^^{data_type}^^{ele_ref_name}^^{ele_start_name}^^
!                       {ele_name}^^{merit_type}^^{meas}^^{good_meas}^^{ref}^^
!                       {good_ref}^^{weight}^^{good_user}^^{data_source}^^
!                       {eval_point}^^{s_offset}^^{ix_bunch}^^{invalid_value}^^
!                       {spin_axis%n0(1)}^^{spin_axis%n0(2)}^^{spin_axis%n0(3)}^^
!                       {spin_axis%l(1)}^^{spin_axis%l(2)}^^{spin_axis%l(3)}
! 
! Note: The 3 values for spin_axis%n0, as a group, are optional. 
!       Also the 3 values for spin_axis%l are, as a group, optional.
! Note: Use the "pipe data_d2_create" command first to create a d2 structure 
!       with associated d1 arrays.
! Note: After creating all your datums, use the "pipe data_set_design_value" routine
!       to set the design (and model) values.
! 
! Parameters
! ----------
! datum_name          ! EG: orb.x[3]
! data_type           ! EG: orbit.x
! ele_ref_name : optional
! ele_start_name : optional
! ele_name : optional
! merit_type : optional
! meas : default=0
! good_meas : default=F
! ref : default=0
! good_ref : default=F
! weight : default=0
! good_user : default=T
! data_source : default=lat
! eval_point : default=END
! s_offset : default=0
! ix_bunch : default=0
! invalid_value : default=0
! spin_axis%n0(1) : optional
! spin_axis%n0(2) : optional
! spin_axis%n0(3) : optional
! spin_axis%l(1) : optional
! spin_axis%l(2) : optional
! spin_axis%l(3) : optional
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_optics_matching
!  args:
!    datum_name: twiss.end[6]
!    data_type: beta.y
!    ele_ref_name:
!    ele_start_name:
!    ele_name: P1
!    merit_type: target
!    meas: 0
!    good_meas: T
!    ref: 0
!    good_ref: T
!    weight: 0.3
!    good_user: T
!    data_source: lat
!    eval_point: END
!    s_offset: 0
!    ix_bunch: 1
!    invalid_value: 0


case ('datum_create')

  allocate (name_arr(23))
  call split_this_line (line, name_arr, -1, err, n_arr);  if (err) return
  if (n_arr /= 17 .and. n_arr /= 20 .and. n_arr /= 23) then
    call invalid('NUMBER OF COMPONENTS ON LINE NOT CORRECT.')
  endif

  call tao_find_data (err, name_arr(1), d_array = d_array);
  if (err .or. size(d_array) /= 1) then
    call invalid('BAD DATUM NAME')
    return
  endif

  d_ptr => d_array(1)%d
  u => s%u(d_ptr%ix_uni)

  d_ptr%data_type = name_arr(2)

  ele_name = upcase(name_arr(3))
  d_ptr%ele_start_name = ele_name
  if (ele_name /= '') then
    call lat_ele_locator (ele_name, u%model%lat, eles, n_loc, err)
    if (err .or. n_loc /= 1) then
      call invalid('UNIQUE LATTICE START ELEMENT NOT FOUND FOR: ' // ele_name)
      return
    endif
    d_ptr%ix_ele_start = eles(1)%ele%ix_ele
  endif

  ele_name = upcase(name_arr(4))
  d_ptr%ele_ref_name = ele_name
  if (ele_name /= '') then
    call lat_ele_locator (ele_name, u%model%lat, eles, n_loc, err)
    if (err .or. n_loc /= 1) then
      call invalid('UNIQUE LATTICE REF ELEMENT NOT FOUND FOR: ' // ele_name)
      return
    endif
    d_ptr%ix_ele_ref = eles(1)%ele%ix_ele
  endif

  ele_name = upcase(name_arr(5))
  d_ptr%ele_name = ele_name
  if (ele_name /= '') then
    call lat_ele_locator (ele_name, u%model%lat, eles, n_loc, err)
    if (err .or. n_loc == 0) then
      call invalid('LATTICE ELEMENT NOT FOUND FOR: ' // ele_name)
      return
    elseif (n_loc > 1) then
      call invalid('UNIQUE LATTICE ELEMENT NOT FOUND FOR: ' // ele_name)
      return
    endif
    d_ptr%ix_ele = eles(1)%ele%ix_ele
    d_ptr%ix_branch = eles(1)%ele%ix_branch
  endif

  d_ptr%merit_type  = name_arr(6)
  if (d_ptr%merit_type == '') d_ptr%merit_type = 'target'
  d_ptr%meas_value  = real_val(name_arr(7), 0.0_rp, err);      if (err) return
  d_ptr%good_meas   = logic_val(name_arr(8), .false., err);    if (err) return
  d_ptr%ref_value   = real_val(name_arr(9), 0.0_rp, err);      if (err) return
  d_ptr%good_ref    = logic_val(name_arr(10), .false., err);    if (err) return
  d_ptr%weight      = real_val(name_arr(11), 0.0_rp, err);      if (err) return
  d_ptr%good_user   = logic_val(name_arr(12), .true., err);    if (err) return
  d_ptr%data_source = str_val(name_arr(13), 'lat')
  if (name_arr(12) == '') then
    d_ptr%eval_point  = anchor_end$
  else
    call match_word (name_arr(14), anchor_pt_name(1:), d_ptr%eval_point)
  endif
  d_ptr%s_offset = real_val(name_arr(15), 0.0_rp, err);          if (err) return
  d_ptr%ix_bunch = int_val(name_arr(16), 1, err);                if (err) return
  d_ptr%invalid_value = real_val(name_arr(17), 0.0_rp, err);     if (err) return
  if (n_arr > 17) then
    d_ptr%spin_map%axis_input%n0(1) = real_val(name_arr(18), 0.0_rp, err);   if (err) return
    d_ptr%spin_map%axis_input%n0(2) = real_val(name_arr(19), 0.0_rp, err);   if (err) return
    d_ptr%spin_map%axis_input%n0(3) = real_val(name_arr(20), 0.0_rp, err);   if (err) return
  endif
  if (n_arr == 23) then
    d_ptr%spin_map%axis_input%l(1) = real_val(name_arr(21), 0.0_rp, err);   if (err) return
    d_ptr%spin_map%axis_input%l(2) = real_val(name_arr(22), 0.0_rp, err);   if (err) return
    d_ptr%spin_map%axis_input%l(3) = real_val(name_arr(23), 0.0_rp, err);   if (err) return
  endif

  d_ptr%exists = tao_data_sanity_check (d_ptr, d_ptr%exists, '')
  if (tao_chrom_calc_needed(d_ptr%data_type, d_ptr%data_source)) u%calc%chrom_for_data = .true.

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% datum_has_ele
!
! Output whether a datum type has an associated lattice element
!
! Notes
! -----
! Command syntax:
!   pipe datum_has_ele {datum_type}
! 
! Parameters
! ----------
! datum_type
!
! Returns
! -------
! string_list 
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_optics_matching
!  args:
!    datum_type: twiss.end 

case ('datum_has_ele')

  select case (tao_datum_has_associated_ele(line))
  case (no$);             nl=incr(nl); li(nl) = 'no'
  case (yes$);            nl=incr(nl); li(nl) = 'yes'
  case (maybe$);          nl=incr(nl); li(nl) = 'maybe'
  case (provisional$);    nl=incr(nl); li(nl) = 'provisional'
  end select

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% derivative
!
! Output optimization derivatives
!
! Notes
! -----
! Command syntax:
!   pipe derivative
!
! Note: To save time, this command will not recalculate derivatives. 
! Use the "derivative" command beforehand to recalcuate if needed.
! 
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_optics_matching
!  args:

case ('derivative')

  call tao_dmodel_dvar_calc(.false., err);  if (err) return

  do iu = lbound(s%u, 1), ubound(s%u, 1)
    u => s%u(iu)
    n2 = ubound(u%dModel_dVar, 2)
    do id = 1, ubound(u%dModel_dVar, 1)
      do iv = 0, (n2-1)/10
        iv0 = 10 * iv + 1
        nl=incr(nl); write(li(nl), '(3(i0, a), 20a)') iu, ';', id, ';', iv0, (';', re_str(u%dModel_dVar(id, j), 10), j = iv0, min(iv0+9, n2))
      enddo
    enddo
  enddo

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% ele:ac_kicker
!
! Output element ac_kicker parameters
!
! Notes
! -----
! Command syntax:
!   pipe ele:ac_kicker {ele_id}|{which}
!
! Where: 
!   {ele_id} is an element name or index.
!   {which} is one of: "model", "base" or "design"
!
! Example:
!   pipe ele:ac_kicker 3@1>>7|model
! This gives element number 7 in branch 1 of universe 3.
! 
! Parameters
! ----------
! ele_id
! which : default=model
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!   ele_id: 1@0>>1
!   which: model

case ('ele:ac_kicker')

  u => point_to_uni(line, .true., err); if (err) return
  tao_lat => point_to_tao_lat(line, u, err, which, tail_str); if (err) return
  ele => point_to_ele(line, tao_lat%lat, err); if (err) return

  if (.not. associated(ele%ac_kick)) return
  ac => ele%ac_kick

  if (allocated(ac%amp_vs_time)) then
    nl=incr(nl); write (li(nl), '(a)') 'has#amp_vs_time'
    do i = 1, size(ac%amp_vs_time)
      nl=incr(nl); write (li(nl), '(i0, 2(a, es22.14))') i, ';', ac%amp_vs_time(i)%amp, ';', ac%amp_vs_time(i)%time
    enddo

  else
    nl=incr(nl); write (li(nl), '(a)') 'has#frequencies'
    do i = 1, size(ac%frequency)
      nl=incr(nl); write (li(nl), '(i0, 3(a, es22.14))') i, ';', &
                      ac%frequency(i)%f, ';', ac%frequency(i)%amp, ';', ac%frequency(i)%phi
    enddo
  endif

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% ele:cartesian_map
!
! Output element cartesian_map parameters
!
! Notes
! -----
! Command syntax:
!   pipe ele:cartesian_map {ele_id}|{which} {index} {who}
!
! Where:
!   {ele_id} is an element name or index
!   {which} is one of: "model", "base" or "design"
!   {index} is the index number in the ele%cartesian_map(:) array
!   {who} is one of: "base", or "terms"
!
! Example:
!   pipe ele:cartesian_map 3@1>>7|model 2 base
! This gives element number 7 in branch 1 of universe 3.
! 
! Parameters
! ----------
! ele_id
! index 
! who
! which : default=model
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_em_field
!  args:
!   ele_id: 1@0>>1
!   which: model
!   index: 1
!   who: base

case ('ele:cartesian_map')

  u => point_to_uni(line, .true., err); if (err) return
  tao_lat => point_to_tao_lat(line, u, err, which, tail_str); if (err) return
  ele => point_to_ele(line, tao_lat%lat, err); if (err) return

  if (.not. associated(ele%cartesian_map)) then
    call invalid ('cartesian_map not allocated')
    return
  endif
  ix = parse_int (tail_str, err, 1, size(ele%cartesian_map));  if (err) return
  ct_map => ele%cartesian_map(ix)

  select case (tail_str)
  case ('base')
    nl=incr(nl); write (li(nl), amt) 'file;FILE;T;',                          trim(ct_map%ptr%file)
    nl=incr(nl); write (li(nl), rmt) 'field_scale;REAL;T;',                   ct_map%field_scale
    nl=incr(nl); write (li(nl), ramt) 'r0;REAL_ARR;T',                        (';', ct_map%r0(i), i = 1, 3)
    name = attribute_name(ele, ct_map%master_parameter)
    if (name(1:1) == '!') name = '<None>'
    nl=incr(nl); write (li(nl), amt) 'master_parameter;ELE_PARAM;T;',         trim(name)
    nl=incr(nl); write (li(nl), amt) 'ele_anchor_pt;ENUM;T;',                 trim(anchor_pt_name(ct_map%ele_anchor_pt))
    nl=incr(nl); write (li(nl), amt) 'nongrid^field_type;ENUM;T;',            trim(em_field_type_name(ct_map%field_type))

  case ('terms')
    do i = 1, size(ct_map%ptr%term)
      ctt => ct_map%ptr%term(i)
      nl=incr(nl); write (li(nl), '(i0, 7(a, es22.14), 4a)') i, ';', &
            ctt%coef, ';', ctt%kx, ';', ctt%ky, ';', ctt%kz, ';', ctt%x0, ';', ctt%y0, ';', ctt%phi_z, ';', &
            trim(cartesian_map_family_name(ctt%family)), ';', trim(cartesian_map_form_name(ctt%form))
    enddo

  case default
    call invalid ('{who} is not "base" or "terms"')
    return
  end select

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% ele:chamber_wall
!
! Output element beam chamber wall parameters
!
! Notes
! -----
! Command syntax:
!   pipe ele:chamber_wall {ele_id}|{which} {index} {who}
!
! Where:
!   {ele_id} is an element name or index.
!   {which} is one of: "model", "base" or "design"
!   {index} is index of the wall.
!   {who} is one of:
!     "x"       ! Return min/max in horizontal plane
!     "y"       ! Return min/max in vertical plane
! 
! Parameters
! ----------
! ele_id
! index
! who
! which : default=model
!
! Returns
! -------
! string_list 
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_wall3d
!  args:
!   ele_id: 1@0>>1
!   which: model
!   index: 1
!   who: x

case ('ele:chamber_wall')

  u => point_to_uni(line, .true., err); if (err) return
  tao_lat => point_to_tao_lat(line, u, err, which, tail_str); if (err) return
  ele => point_to_ele(line, tao_lat%lat, err); if (err) return

  if (.not. associated(ele%wall3d)) then
    call invalid ('No associated wall')
    return
  endif

  n = parse_int(tail_str, err, 1, size(ele%wall3d));        if (err) return
  wall3d => ele%wall3d(n)

  do i = 1, size(wall3d%section)
    select case (tail_str)
    case ('x')
      call calc_wall_radius (wall3d%section(i)%v,  1.0_rp, 0.0_rp, z1, rdummy)
      call calc_wall_radius (wall3d%section(i)%v, -1.0_rp, 0.0_rp, z2, rdummy)
    case ('y')
      call calc_wall_radius (wall3d%section(i)%v, 0.0_rp,  1.0_rp, z1, rdummy)
      call calc_wall_radius (wall3d%section(i)%v, 0.0_rp, -1.0_rp, z2, rdummy)
    case default
      call invalid ('{who} is not "x" or "y"')
      return
    end select

    nl=incr(nl); write (li(nl), '(i0, 3(a, es14.6))') i, ';', wall3d%section(i)%s, ';', z1, ';', -z2
  enddo

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% ele:control_var
!
! Output list of element control variables.
! Used for group, overlay and ramper type elements.
!
! Notes
! -----
! Command syntax:
!   pipe ele:control_var {ele_id}|{which}
!
! Where:
!   {ele_id} is an element name or index.
!   {which} is one of: "model", "base" or "design"
!
! Example:
!   pipe ele:control_var 3@1>>7|model
! This gives control info on element number 7 in branch 1 of universe 3.
! 
! Parameters
! ----------
! ele_id
! which : default=model
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!   ele_id: 1@0>>873
!   which: model

case ('ele:control_var')

  u => point_to_uni(line, .true., err); if (err) return
  tao_lat => point_to_tao_lat(line, u, err, which, tail_str); if (err) return
  ele => point_to_ele(line, tao_lat%lat, err); if (err) return

  if (.not. associated(ele%control)) then
    call invalid ('ele%control not allocated')
    return
  endif

  if (.not. allocated(ele%control%var)) then
    call invalid ('ele%control%var not allocated')
    return
  endif

  ! Group controller var has an old_value. Overlay and ramper vars do not.

  if (ele%key == group$) then
    do i = 1, size(ele%control%var)
      cvar => ele%control%var(i)
      nl=incr(nl); write (li(nl), '(i0, 2a, 2(a, es22.14))') i, ';', trim(cvar%name), ';', cvar%value, ';', cvar%old_value
    enddo
  else
    do i = 1, size(ele%control%var)
      cvar => ele%control%var(i)
      nl=incr(nl); write (li(nl), '(i0, 2a, 2(a, es22.14))') i, ';', trim(cvar%name), ';', cvar%value
    enddo
  endif

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% ele:cylindrical_map
!
! Output element cylindrical_map
!
! Notes
! -----
! Command syntax:
!   pipe ele:cylindrical_map {ele_id}|{which} {index} {who}
!
! Where 
!   {ele_id} is an element name or index.
!   {which} is one of: "model", "base" or "design"
!   {index} is the index number in the ele%cylindrical_map(:) array
!   {who} is one of: "base", or "terms"
!
! Example:
!   pipe ele:cylindrical_map 3@1>>7|model 2 base
! This gives map #2 of element number 7 in branch 1 of universe 3.
! 
! Parameters
! ----------
! ele_id
! index
! who
! which : default=model
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_em_field
!  args:
!   ele_id: 1@0>>5
!   which: model
!   index: 1
!   who: base

case ('ele:cylindrical_map')

  u => point_to_uni(line, .true., err); if (err) return
  tao_lat => point_to_tao_lat(line, u, err, which, tail_str); if (err) return
  ele => point_to_ele(line, tao_lat%lat, err); if (err) return

  if (.not. associated(ele%cylindrical_map)) then
    call invalid ('cylindrical_map not allocated')
    return
  endif
  ix = parse_int (tail_str, err, 1, size(ele%cylindrical_map));  if (err) return
  cy_map => ele%cylindrical_map(ix)

  select case (tail_str)
  case ('base')
    nl=incr(nl); write (li(nl), amt)  'file;FILE;T;',                          trim(cy_map%ptr%file)
    nl=incr(nl); write (li(nl), imt)  'm;INT;T;',                              cy_map%m
    nl=incr(nl); write (li(nl), imt)  'harmonic;INT;T;',                       cy_map%harmonic
    nl=incr(nl); write (li(nl), rmt)  'phi0_fieldmap;REAL;T;',                 cy_map%phi0_fieldmap
    nl=incr(nl); write (li(nl), rmt)  'theta0_azimuth;REAL;T;',                cy_map%theta0_azimuth
    nl=incr(nl); write (li(nl), rmt)  'field_scale;REAL;T;',                   cy_map%field_scale
    nl=incr(nl); write (li(nl), rmt)  'dz;REAL;T;',                            cy_map%dz
    nl=incr(nl); write (li(nl), ramt) 'r0;REAL_ARR;T',                         (';', cy_map%r0(i), i = 1, 3)
    name = attribute_name(ele, cy_map%master_parameter)
    if (name(1:1) == '!') name = '<None>'
    nl=incr(nl); write (li(nl), amt)  'master_parameter;ELE_PARAM;T;',         trim(name)
    nl=incr(nl); write (li(nl), amt)  'ele_anchor_pt;ENUM;T;',                 trim(anchor_pt_name(cy_map%ele_anchor_pt))
    nl=incr(nl); write (li(nl), imt)  'number_of_terms;INT;F;',                size(cy_map%ptr%term)

  case ('terms')
    do i = 1, size(cy_map%ptr%term)
      cyt => cy_map%ptr%term(i)
      nl=incr(nl); write (li(nl), '(i0, 2a)') i, cmplx_str(cyt%e_coef), cmplx_str(cyt%b_coef)
    enddo

  case default
    call invalid ('{who} is not "base" or "terms"')
    return
  end select

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% ele:elec_multipoles
!
! Output element electric multipoles
!
! Notes
! -----
! Command syntax:
!   pipe ele:elec_multipoles {ele_id}|{which}
!
! Where:
!   {ele_id} is an element name or index.
!   {which} is one of: "model", "base" or "design"
!
! Example:
!   pipe ele:elec_multipoles 3@1>>7|model
! This gives element number 7 in branch 1 of universe 3.
! 
! Parameters
! ----------
! ele_id
! which : default=model
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!   ele_id: 1@0>>1
!   which: model
 

case ('ele:elec_multipoles')

  u => point_to_uni(line, .true., err); if (err) return
  tao_lat => point_to_tao_lat(line, u, err, which, tail_str); if (err) return
  ele => point_to_ele(line, tao_lat%lat, err); if (err) return

  nl=incr(nl); write (li(nl), lmt) 'multipoles_on;LOGIC;T;', ele%multipoles_on
  if (attribute_index(ele, 'SCALE_MULTIPOLES') == scale_multipoles$) then
    nl=incr(nl); write (li(nl), lmt) 'scale_multipoles;LOGIC;T;', ele%scale_multipoles
  endif

  can_vary = (which == 'model')

  nl=incr(nl); li(nl) = 'An_elec;Bn_elec;An_elec (Scaled);Bn_elec (Scaled)'

  if (.not. associated(ele%a_pole_elec)) then
    call end_stuff(li, nl)
    return
  endif

  call multipole_ele_to_ab (ele, .false., ix_pole_max, a, b, electric$)

  do i = 0, n_pole_maxx
    if (ele%a_pole_elec(i) == 0 .and. ele%b_pole_elec(i) == 0) cycle
    nl=incr(nl); write (li(nl), '(i0, 4(a, es22.14))') i, ';', ele%a_pole_elec(i), ';', ele%b_pole_elec(i), ';', a(i), ';', b(i)
  enddo

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% ele:floor
!
! Output element floor coordinates. The output gives four lines. "Reference" is
! without element misalignments and "Actual" is with misalignments. The lines with "-W"
! give the W matrix. The exception is that if ele is a multipass_lord, there will be 4*N
! lines where N is the number of slaves.
!
! Notes
! -----
! Command syntax:
!   pipe ele:floor {ele_id}|{which} {where}
!
! Where:
!   {ele_id} is an element name or index.
!   {which} is one of: "model", "base" or "design"
!   {where} is an optional argument which, if present, is one of
!     beginning  ! Upstream end.
!     center     ! Middle of the element. This is the surface of element when used 
!                !  with photonic reflecting elements such as crystal and mirror elements.
!     end        ! Downstream end (default).
!
! Example:
!   pipe ele:floor 3@1>>7|model
! This gives element number 7 in branch 1 of universe 3.
! 
! Parameters
! ----------
! ele_id
! which : default=model
! where : default=end
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!   ele_id: 1@0>>1
!   which: model
!   where: 
!
! Example: 2
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!   ele_id: 1@0>>1
!   which: model
!   where: center

case ('ele:floor')

  u => point_to_uni(line, .true., err); if (err) return
  tao_lat => point_to_tao_lat(line, u, err, which, tail_str); if (err) return
  ele => point_to_ele(line, tao_lat%lat, err); if (err) return

  if (tail_str == '') tail_str = 'end'
  call match_word (tail_str, [character(12):: 'beginning', 'center', 'end'], loc)
  if (loc == 0) then
    call invalid ('BAD "where" SWITCH. SHOULD BE ONE OF "", "beginning", "center", or "end".')
    return
  endif

  can_vary = (ele%ix_ele == 0 .and. which == 'model')
  call write_this_ele_floor(ele, loc, can_vary, '')

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% ele:gen_attribs
!
! Output element general attributes
!
! Notes
! -----
! Command syntax:
!   pipe ele:gen_attribs {ele_id}|{which}
!
! Where: 
!   {ele_id} is an element name or index.
!   {which} is one of: "model", "base" or "design"
!
! Example:
!   pipe ele:gen_attribs 3@1>>7|model
! This gives element number 7 in branch 1 of universe 3.
! 
! Parameters
! ----------
! ele_id
! which : default=model
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!   ele_id: 1@0>>1
!   which: model

case ('ele:gen_attribs')

  u => point_to_uni(line, .true., err); if (err) return
  tao_lat => point_to_tao_lat(line, u, err, which, tail_str); if (err) return
  ele => point_to_ele(line, tao_lat%lat, err); if (err) return

  do i = 1, num_ele_attrib$
    attrib = attribute_info(ele, i)
    a_name = attrib%name
    if (a_name == null_name$) cycle
    if (attrib%state == private$) cycle

    free = attribute_free (ele, a_name, .false., why_not_free = why_not_free)
    if (.not. free .and. why_not_free == field_master_dependent$) free = .true.
    attrib_type = attribute_type(a_name)
    if (which /= 'model') free = .false.

    select case (attrib_type)
    case (is_logical$)
      nl=incr(nl); write (li(nl), '(2a, l1, a, l1)') trim(a_name), ';LOGIC;', free, ';', is_true(ele%value(i))
    case (is_integer$)
      nl=incr(nl); write (li(nl), '(2a, l1, a, i0)') trim(a_name), ';INT;', free, ';', nint(ele%value(i))
    case (is_real$)
      nl=incr(nl); write (li(nl), '(2a, l1, a, es22.14)') trim(a_name), ';REAL;', free, ';', ele%value(i)
      nl=incr(nl); write (li(nl), '(4a)') 'units#', trim(a_name), ';STR;F;', attrib%units
    case (is_switch$)
      name = switch_attrib_value_name (a_name, ele%value(i), ele)
      nl=incr(nl); write (li(nl), '(2a, l1, 2a)') trim(a_name), ';ENUM;', free, ';', trim(name)
    end select
  enddo

  nl=incr(nl); write (li(nl), amt) 'lord_status;ENUM;F;', trim(control_name(ele%lord_status))
  nl=incr(nl); write (li(nl), amt) 'slave_status;ENUM;F;', trim(control_name(ele%slave_status))
  

  if (attribute_name(ele, aperture_at$) == 'APERTURE_AT') then
    nl=incr(nl); write (li(nl), amt) 'aperture_at;ENUM;T;', trim(aperture_at_name(ele%aperture_at))
    nl=incr(nl); write (li(nl), lmt) 'offset_moves_aperture;LOGIC;T;',          ele%offset_moves_aperture
  endif

  if (attribute_name(ele, aperture_type$) == 'APERTURE_TYPE') then
    nl=incr(nl); write (li(nl), amt) 'aperture_type;ENUM;T;', trim(aperture_type_name(ele%aperture_type))
  endif

  if (attribute_index(ele, 'FIELD_MASTER') /= 0) then
    nl=incr(nl); write (li(nl), lmt) 'field_master;LOGIC;T;',                   ele%field_master
  endif

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% ele:gen_grad_map
!
! Output element gen_grad_map 
!
! Notes
! -----
! Command syntax:
!   pipe ele:gen_grad_map {ele_id}|{which} {index} {who}
!
! Where: 
!   {ele_id} is an element name or index.
!   {which} is one of: "model", "base" or "design"
!   {index} is the index number in the ele%gen_grad_map(:) array
!   {who} is one of: "base", or "derivs".
!
! Example:
!   pipe ele:gen_grad_map 3@1>>7|model 2 base
! This gives element number 7 in branch 1 of universe 3.
! 
! Parameters
! ----------
! ele_id
! index
! who
! which : default=model
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_em_field
!  args:
!   ele_id: 1@0>>9
!   which: model
!   index: 1
!   who: derivs

case ('ele:gen_grad_map')

  u => point_to_uni(line, .true., err); if (err) return
  tao_lat => point_to_tao_lat(line, u, err, which, tail_str); if (err) return
  ele => point_to_ele(line, tao_lat%lat, err); if (err) return

  if (.not. associated(ele%gen_grad_map)) then
    call invalid ('gen_grad_map not allocated')
    return
  endif
  ix = parse_int (tail_str, err, 1, size(ele%gen_grad_map));  if (err) return
  gg_map => ele%gen_grad_map(ix)

  select case (tail_str)
  case ('base')
    nl=incr(nl); write (li(nl), amt)  'file;FILE;T;',                          trim(gg_map%file)
    nl=incr(nl); write (li(nl), rmt)  'field_scale;REAL;T;',                   gg_map%field_scale
    nl=incr(nl); write (li(nl), ramt) 'r0;REAL_ARR;T',                         (';', gg_map%r0(i), i = 1, 3)
    nl=incr(nl); write (li(nl), rmt)  'dz;REAL;T;',                            gg_map%dz
    name = attribute_name(ele, gg_map%master_parameter)
    if (name(1:1) == '!') name = '<None>'
    nl=incr(nl); write (li(nl), amt) 'master_parameter;ELE_PARAM;T;',        trim(name)
    nl=incr(nl); write (li(nl), amt) 'ele_anchor_pt;ENUM;T;',                 trim(anchor_pt_name(gg_map%ele_anchor_pt))
    nl=incr(nl); write (li(nl), amt) 'nongrid^field_type;ENUM;T;',            trim(em_field_type_name(gg_map%field_type))
    nl=incr(nl); write (li(nl), lmt) 'curved_ref_frame;LOGIC;T;',             gg_map%curved_ref_frame
    nl=incr(nl); write (li(nl), imt) 'iz0;INT;F;',                            gg_map%iz0
    nl=incr(nl); write (li(nl), imt) 'iz1;INT;F;',                            gg_map%iz1
    nl=incr(nl); write (li(nl), imt) 'size_of_gg;INT;F;',                     size(gg_map%gg)

  case ('derivs')
    do i = 1, size(gg_map%gg)
      gg => gg_map%gg(i)
      do j = gg_map%iz0, gg_map%iz1
        do k = 0, ubound(gg%deriv,2)
          nl=incr(nl); write (li(nl), '(3(i0,a), es22.14, a, es22.14)') i, ';', j, ';', k, ';', &
                              j*gg_map%dz, ';', gg%deriv(j,k)
        enddo
      enddo
    enddo
  end select

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% ele:grid_field
!
! Output element grid_field
!
! Notes
! -----
! Command syntax:
!   pipe ele:grid_field {ele_id}|{which} {index} {who}
!
! Where:
!   {ele_id} is an element name or index.
!   {which} is one of: "model", "base" or "design"
!   {index} is the index number in the ele%grid_field(:) array.
!   {who} is one of: "base", or "points"
!
! Example:
!   pipe ele:grid_field 3@1>>7|model 2 base
! This gives grid #2 of element number 7 in branch 1 of universe 3.
! 
! Parameters
! ----------
! ele_id
! index 
! who 
! which : default=model
!    
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_grid
!  args:
!   ele_id: 1@0>>1
!   which: model
!   index: 1
!   who: base 

case ('ele:grid_field')

  u => point_to_uni(line, .true., err); if (err) return
  tao_lat => point_to_tao_lat(line, u, err, which, tail_str); if (err) return
  ele => point_to_ele(line, tao_lat%lat, err); if (err) return

  if (.not. associated(ele%grid_field)) then
    call invalid ('grid_field not allocated')
    return
  endif
  ix = parse_int (tail_str, err, 1, size(ele%grid_field));  if (err) return
  g_field => ele%grid_field(ix)

  select case (tail_str)
  case ('base')
    nl=incr(nl); write (li(nl), ramt) 'dr;REAL_ARR;T',                        (';', g_field%dr(i), i = 1, 3)
    nl=incr(nl); write (li(nl), ramt) 'r0;REAL_ARR;T',                        (';', g_field%r0(i), i = 1, 3)
    name = attribute_name(ele, g_field%master_parameter)
    if (name(1:1) == '!') name = '<None>'
    nl=incr(nl); write (li(nl), amt) 'master_parameter;ELE_PARAM;T;',         trim(name)
    nl=incr(nl); write (li(nl), amt) 'ele_anchor_pt;ENUM;T;',                 trim(anchor_pt_name(g_field%ele_anchor_pt))
    nl=incr(nl); write (li(nl), amt) 'field_type;ENUM;T;',                    trim(em_field_type_name(g_field%field_type))
    nl=incr(nl); write (li(nl), amt) 'grid_field^geometry;ENUM;T;',           trim(grid_field_geometry_name(g_field%geometry))
    nl=incr(nl); write (li(nl), imt) 'harmonic;INT;T;',                       g_field%harmonic
    nl=incr(nl); write (li(nl), rmt) 'phi0_fieldmap;REAL;T;',                 g_field%phi0_fieldmap
    nl=incr(nl); write (li(nl), rmt) 'field_scale;REAL;T;',                   g_field%field_scale
    nl=incr(nl); write (li(nl), imt) 'interpolation_order;INUM;T;',           g_field%interpolation_order
    nl=incr(nl); write (li(nl), lmt) 'curved_ref_frame;LOGIC;T;',             g_field%curved_ref_frame
    nl=incr(nl); write (li(nl), amt) 'file;FILE;T;',                          trim(g_field%ptr%file)

  case ('points')
    do i = lbound(g_field%ptr%pt, 1), ubound(g_field%ptr%pt, 1)
    do j = lbound(g_field%ptr%pt, 2), ubound(g_field%ptr%pt, 2)
    do k = lbound(g_field%ptr%pt, 3), ubound(g_field%ptr%pt, 3)
      g_pt => g_field%ptr%pt(i,j,k)
      select case (g_field%field_type)
      case (electric$)
        if (g_field%harmonic == 0) then
          nl=incr(nl); write (li(nl), '(2(i0, a), i0, 3a)') i, ';', j, ';', k, (real_part_str(g_pt%E(ix)), ix = 1, 3)
        else
          nl=incr(nl); write (li(nl), '(2(i0, a), i0, 3a)') i, ';', j, ';', k, (cmplx_str(g_pt%E(ix)), ix = 1, 3)
        endif
      case (magnetic$)
        if (g_field%harmonic == 0) then
          nl=incr(nl); write (li(nl), '(2(i0, a), i0, 3a)') i, ';', j, ';', k, (real_part_str(g_pt%B(ix)), ix = 1, 3)
        else
          nl=incr(nl); write (li(nl), '(2(i0, a), i0, 3a)') i, ';', j, ';', k, (cmplx_str(g_pt%B(ix)), ix = 1, 3)
        endif
      case (mixed$)
        if (g_field%harmonic == 0) then
          nl=incr(nl); write (li(nl), '(2(i0, a), i0, 6a)') i, ';', j, ';', k, &
                                                            (real_part_str(g_pt%B(ix)), ix = 1, 3), (real_part_str(g_pt%E(ix)), ix = 1, 3)
        else
          nl=incr(nl); write (li(nl), '(2(i0, a), i0, 6a)') i, ';', j, ';', k, &
                                                            (cmplx_str(g_pt%B(ix)), ix = 1, 3), (cmplx_str(g_pt%B(ix)), ix = 1, 3)
        endif
      end select
    enddo; enddo; enddo
  end select

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% ele:head
!
! Output "head" Element attributes
!
! Notes
! -----
! Command syntax:
!   pipe ele:head {ele_id}|{which}
!
! Where: 
!   {ele_id} is an element name or index.
!   {which} is one of: "model", "base" or "design"
!
! Example:
!   pipe ele:head 3@1>>7|model
! This gives element number 7 in branch 1 of universe 3.
! 
! Parameters
! ----------
! ele_id 
! which : default=model
!
! Returns
! -------
! string_list 
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!   ele_id: 1@0>>1
!   which: model

case ('ele:head')

  u => point_to_uni(line, .true., err); if (err) return
  tao_lat => point_to_tao_lat(line, u, err, which, tail_str); if (err) return
  ele => point_to_ele(line, tao_lat%lat, err); if (err) return

  can_vary = (ele%slave_status /= multipass_slave$ .and. ele%slave_status /= super_slave$ .and. ele%ix_ele /= 0)

  nl=incr(nl); write (li(nl), imt) 'universe;INT;F;',                 u%ix_uni
  nl=incr(nl); write (li(nl), jmt) u%ix_uni, '^ix_branch;INUM;F;',    ele%ix_branch
  nl=incr(nl); write (li(nl), imt) 'ix_ele;INT;I;',                   ele%ix_ele

  nl=incr(nl); write (li(nl), amt) 'key;ENUM;F;',                     trim(key_name(ele%key))
  nl=incr(nl); write (li(nl), amt) 'name;STR;F;',                     trim(ele%name)
  nl=incr(nl); write (li(nl), amt2) 'type;STR;', can_vary, ';',       ele%type
  nl=incr(nl); write (li(nl), amt2) 'alias;STR;', can_vary, ';',      ele%alias
  if (associated(ele%descrip)) then
    nl=incr(nl); write (li(nl), amt2) 'descrip;STR;', can_vary, ';',  ele%descrip
  else
    nl=incr(nl); write (li(nl), amt2) 'descrip;STR;', can_vary, ';',  ''
  endif
  nl=incr(nl); write (li(nl), '(2(a,l1))') 'is_on;LOGIC;', attribute_free(ele, 'is_on', .false.), ';', ele%is_on

  nl=incr(nl); write (li(nl), rmt) 's;REAL;F;',                     ele%s
  nl=incr(nl); write (li(nl), rmt) 's_start;REAL;F;',               ele%s_start
  nl=incr(nl); write (li(nl), rmt) 'ref_time;REAL;F;',              ele%ref_time

  nl=incr(nl); write (li(nl), lmt) 'has#methods;LOGIC;F;',          (ele%key /= overlay$ .and. ele%key /= group$ .and. ele%key /= girder$)
  nl=incr(nl); write (li(nl), lmt) 'has#ab_multipoles;LOGIC;F;',    (attribute_name(ele, a0$) == 'A0')
  nl=incr(nl); write (li(nl), lmt) 'has#kt_multipoles;LOGIC;F;',    (ele%key == multipole$)
  nl=incr(nl); write (li(nl), lmt) 'has#multipoles_elec;LOGIC;F;',  (attribute_name(ele, a0_elec$) == 'A0_ELEC')
  nl=incr(nl); write (li(nl), lmt) 'has#ac_kick;LOGIC;F;',          associated(ele%ac_kick)
  nl=incr(nl); write (li(nl), lmt) 'has#taylor;LOGIC;F;',           associated(ele%taylor(1)%term)
  nl=incr(nl); write (li(nl), lmt) 'has#spin_taylor;LOGIC;F;',      associated(ele%spin_taylor(1)%term)
  nl=incr(nl); write (li(nl), lmt) 'has#wake;LOGIC;F;',             associated(ele%wake)
  n = 0; if (associated(ele%cartesian_map)) n = size(ele%cartesian_map)
  nl=incr(nl); write (li(nl), imt) 'num#cartesian_map;INT;F;',    n
  n = 0; if (associated(ele%cylindrical_map)) n = size(ele%cylindrical_map)
  nl=incr(nl); write (li(nl), imt) 'num#cylindrical_map;INT;F;',  n
  n = 0; if (associated(ele%gen_grad_map)) n = size(ele%gen_grad_map)
  nl=incr(nl); write (li(nl), imt) 'num#gen_grad_map;INT;F;',     n
  n = 0; if (associated(ele%grid_field)) n = size(ele%grid_field)
  nl=incr(nl); write (li(nl), imt) 'num#grid_field;INT;F;',       n
  n = 0; if (associated(ele%wall3d)) n = size(ele%wall3d)
  nl=incr(nl); write (li(nl), imt) 'has#wall3d;INT;F;',           n
  nl=incr(nl); write (li(nl), lmt) 'has#control;LOGIC;F;',          associated(ele%control)
  nl=incr(nl); write (li(nl), lmt) 'has#twiss;LOGIC;F;',            (ele%a%beta /= 0)
  nl=incr(nl); write (li(nl), lmt) 'has#mat6;LOGIC;F;',             (attribute_name(ele, mat6_calc_method$) == 'MAT6_CALC_METHOD')
  nl=incr(nl); write (li(nl), lmt) 'has#floor;LOGIC;F;',            (ele%lord_status /= multipass_lord$)
  nl=incr(nl); write (li(nl), lmt) 'has#photon;LOGIC;F;',           associated(ele%photon)
  nl=incr(nl); write (li(nl), lmt) 'has#lord_slave;LOGIC;F;',       .true.

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% ele:lord_slave
!
! Output the lord/slave tree of an element.
!
! Notes
! -----
! Command syntax:
!   pipe ele:lord_slave {ele_id}|{which}
!
! Where: 
!   {ele_id} is an element name or index.
!   {which} is one of: "model", "base" or "design"
!
! Example:
!   pipe ele:lord_slave 3@1>>7|model
! This gives lord and slave info on element number 7 in branch 1 of universe 3.
! Note: The lord/slave info is independent of the setting of {which}.
! 
! The output is a number of lines.
! Each line gives information on an element (element index, etc.).
! Some lines begin with the word "Element". 
! After each "Element" line, there are a number of lines (possibly zero) 
! that begin with the word "Slave or "Lord".
! These "Slave" and "Lord" lines are the slaves and lords of the "Element" element.
! 
! Parameters
! ----------
! ele_id
! which : default=model
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!   ele_id: 1@0>>1
!   which: model
 

case ('ele:lord_slave')

  u => point_to_uni(line, .true., err); if (err) return
  tao_lat => point_to_tao_lat(line, u, err, which, tail_str); if (err) return
  ele => point_to_ele(line, tao_lat%lat, err); if (err) return

  call tao_control_tree_list(ele, eles)
  do i = size(eles), 1, -1  ! Show lords first
    ele => eles(i)%ele

    nl=incr(nl); write (li(nl), '(8a)') 'Element;', trim(ele_loc_name(ele, .true.)), ';', &
                                  trim(ele%name), ';', trim(key_name(ele%key)), ';', control_name(ele%lord_status)

    do j = 1, ele%n_lord
      lord => pointer_to_lord(ele, j)
      nl=incr(nl); write (li(nl), '(8a)') 'Lord;', trim(ele_loc_name(lord, .true.)), ';', &
                                  trim(lord%name), ';', trim(key_name(lord%key)), ';', control_name(lord%lord_status)
    enddo

    do j = 1, ele%n_slave+ele%n_slave_field
      slave => pointer_to_slave(ele, j)
      nl=incr(nl); write (li(nl), '(8a)') 'Slave;', trim(ele_loc_name(slave, .true.)), ';', &
                                  trim(slave%name), ';', trim(key_name(slave%key)), ';', control_name(slave%slave_status)
    enddo
  enddo

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% ele:mat6
!
! Output element mat6
!
! Notes
! -----
! Command syntax:
!   pipe ele:mat6 {ele_id}|{which} {who}
!
! Where: 
!   {ele_id} is an element name or index.
!   {which} is one of: "model", "base" or "design"
!   {who} is one of: "mat6", "vec0", or "err"
!
! Example:
!   pipe ele:mat6 3@1>>7|model mat6
! This gives element number 7 in branch 1 of universe 3.
! 
! Parameters
! ----------
! ele_id
! which : default=model
! who : default=mat6
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!   ele_id: 1@0>>1
!   which: model
!   who: mat6

case ('ele:mat6')

  u => point_to_uni(line, .true., err); if (err) return
  tao_lat => point_to_tao_lat(line, u, err, which, tail_str); if (err) return
  ele => point_to_ele(line, tao_lat%lat, err); if (err) return

  select case (tail_str)
  case ('mat6')
    do i = 1, 6
      nl=incr(nl); write (li(nl), '(i0, a, 6(a, es22.14))') i, ';REAL_ARR;F', (';', ele%mat6(i,j), j = 1, 6)
    enddo

  case ('vec0')
    nl=incr(nl); write (li(nl), ramt) 'vec0;REAL_ARR;F', (';', ele%vec0(i), i = 1, 6)

  case ('err')
    nl=incr(nl); write (li(nl), rmt) 'symplectic_error;REAL;F;', mat_symp_error(ele%mat6)

  case default
    call invalid ('Bad or missign {who} switch.')
    return
  end select

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% ele:methods
!
! Output element methods
!
! Notes
! -----
! Command syntax:
!   pipe ele:methods {ele_id}|{which}
!
! Where: 
!   {ele_id} is an element name or index.
!   {which} is one of: "model", "base" or "design"
!
! Example:
!   pipe ele:methods 3@1>>7|model
! This gives element number 7 in branch 1 of universe 3.
! 
! Parameters
! ----------
! ele_id
! which : default=model
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!   ele_id: 1@0>>1
!   which: model

case ('ele:methods')

  u => point_to_uni(line, .true., err); if (err) return
  tao_lat => point_to_tao_lat(line, u, err, which, tail_str); if (err) return
  ele => point_to_ele(line, tao_lat%lat, err); if (err) return

  if (attribute_name(ele, crystal_type$) == 'CRYSTAL_TYPE') then
    nl=incr(nl); write (li(nl), amt) 'crystal_type;STR;T;', trim(ele%component_name)
  endif

  if (attribute_name(ele, material_type$) == 'MATERIAL_TYPE') then
    nl=incr(nl); write (li(nl), amt) 'material_type;STR;T;', trim(ele%component_name)
  endif

  if (attribute_name(ele, origin_ele$) == 'ORIGIN_ELE') then
    nl=incr(nl); write (li(nl), amt) 'origin_ele;STR;T;', trim(ele%component_name)
  endif

  if (attribute_name(ele, physical_source$) == 'PHYSICAL_SOURCE') then
    nl=incr(nl); write (li(nl), amt) 'physical_source;STR;T;', trim(ele%component_name)
  endif

  if (attribute_name(ele, mat6_calc_method$) == 'MAT6_CALC_METHOD') then
    nl=incr(nl); write (li(nl), amt) 'mat6_calc_method;ENUM;T;', trim(mat6_calc_method_name(ele%mat6_calc_method))
  endif

  if (attribute_name(ele, tracking_method$) == 'TRACKING_METHOD') then
    nl=incr(nl); write (li(nl), amt) 'tracking_method;ENUM;T;', trim(tracking_method_name(ele%tracking_method))
  endif

  if (attribute_name(ele, spin_tracking_method$) == 'SPIN_TRACKING_METHOD') then
    nl=incr(nl); write (li(nl), amt) 'spin_tracking_method;ENUM;T;', trim(spin_tracking_method_name(ele%spin_tracking_method))
  endif

  if (attribute_name(ele, csr_method$) == 'CSR_METHOD') then
    nl=incr(nl); write (li(nl), amt) 'csr_method;ENUM;T;', trim(csr_method_name(ele%csr_method))
  endif

  if (attribute_name(ele, space_charge_method$) == 'SPACE_CHARGE_METHOD') then
    nl=incr(nl); write (li(nl), amt) 'space_charge_method;ENUM;T;', trim(space_charge_method_name(ele%space_charge_method))
  endif

  if (attribute_name(ele, ptc_integration_type$) == 'PTC_INTEGRATION_TYPE') then
    nl=incr(nl); write (li(nl), amt) 'ptc_integration_type;ENUM;T;', trim(ptc_integration_type_name(ele%ptc_integration_type))
  endif

  if (attribute_name(ele, field_calc$) == 'FIELD_CALC') then
    nl=incr(nl); write (li(nl), amt) 'field_calc;ENUM;T;', trim(field_calc_name(ele%field_calc))
  endif

  if (ele%key /= overlay$ .and. ele%key /= group$ .and. ele%key /= girder$) then
    nl=incr(nl); write (li(nl), imt) 'longitudinal_orientation;INT;F;',              ele%orientation
  endif

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% ele:multipoles
!
! Output element multipoles
!
! Notes
! -----
! Command syntax:
!   pipe ele:multipoles {ele_id}|{which}
!
! Where: 
!   {ele_id} is an element name or index.
!   {which} is one of: "model", "base" or "design"
!
! Example:
!   pipe ele:multipoles 3@1>>7|model
! This gives element number 7 in branch 1 of universe 3.
! 
! Parameters
! ----------
! ele_id
! which : default=model
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!   ele_id: 1@0>>1
!   which: model

case ('ele:multipoles')

  u => point_to_uni(line, .true., err); if (err) return
  tao_lat => point_to_tao_lat(line, u, err, which, tail_str); if (err) return
  ele => point_to_ele(line, tao_lat%lat, err); if (err) return

  nl=incr(nl); write (li(nl), lmt) 'multipoles_on;LOGIC;T;', ele%multipoles_on
  if (attribute_index(ele, 'SCALE_MULTIPOLES') == scale_multipoles$) then
    nl=incr(nl); write (li(nl), lmt) 'scale_multipoles;LOGIC;T;', ele%scale_multipoles
  endif

  if (ele%key == multipole$) then
    nl=incr(nl); li(nl) = 'KnL;Tn;KnL (w/Tilt);Tn (w/Tilt);An (equiv);Bn (equiv)'
  elseif (ele%key == ab_multipole$) then
    nl=incr(nl); li(nl) = 'An;Bn;An (w/Tilt);Bn (w/Tilt);KnL (equiv);Tn (equiv)'
  else
    nl=incr(nl); li(nl) = 'An;Bn;An (Scaled);Bn (Scaled);An (w/Tilt);Bn (w/Tilt);KnL (equiv);Tn (equiv)'
  endif

  if (.not. associated(ele%a_pole)) then
    call end_stuff(li, nl)
    return
  endif

  a = 0; b = 0; a2 = 0; b2 = 0; knl = 0; tn = 0
  if (ele%key == multipole$) then
    call multipole_ele_to_ab (ele, .false., ix_pole_max, a,  b)
    call multipole_ele_to_kt (ele, .true.,  ix_pole_max, knl, tn)
  else
    call multipole_ele_to_ab (ele, .false., ix_pole_max, a,  b)
    call multipole_ele_to_ab (ele, .true.,  ix_pole_max, a2, b2)
    call multipole_ele_to_kt (ele, .true.,  ix_pole_max, knl, tn)
  endif

  do i = 0, n_pole_maxx
    if (ele%a_pole(i) == 0 .and. ele%b_pole(i) == 0) cycle

    if (ele%key == multipole$) then
      nl=incr(nl); write (li(nl), '(i0, 6(a, es22.14))') i, ';', &
                      ele%a_pole(i), ';', ele%b_pole(i), ';', knl(i), ';', tn(i), ';', a(i), ';', b(i)

    elseif (ele%key == ab_multipole$) then
      nl=incr(nl); write (li(nl), '(i0, 6(a, es22.14))') i, ';', &
                      ele%a_pole(i), ';', ele%b_pole(i), ';', a2(i), ';', b2(i), ';', knl(i), ';', tn(i)

    else
      nl=incr(nl); write (li(nl), '(i0, 8(a, es22.14))') i, ';', &
                      ele%a_pole(i), ';', ele%b_pole(i), ';', a(i), ';', b(i), ';', a2(i), ';', b2(i), ';', knl(i), ';', tn(i)
    endif
  enddo

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% ele:orbit
!
! Output element orbit
!
! Notes
! -----
! Command syntax:
!   pipe ele:orbit {ele_id}|{which}
!
! Where: 
!   {ele_id} is an element name or index.
!   {which} is one of: "model", "base" or "design"
!
! Example:
!   pipe ele:orbit 3@1>>7|model
! This gives element number 7 in branch 1 of universe 3.
! 
! Parameters
! ----------
! ele_id
! which : default=model
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!   ele_id: 1@0>>1
!   which: model

case ('ele:orbit')

  u => point_to_uni(line, .true., err); if (err) return
  tao_lat => point_to_tao_lat(line, u, err, which, tail_str); if (err) return
  ele => point_to_ele(line, tao_lat%lat, err); if (err) return

  call orbit_out (tao_lat%tao_branch(ele%ix_branch)%orbit(ele%ix_ele))

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% ele:param
!
! Output lattice element parameter
!
! Notes
! -----
! Command syntax:
!   pipe ele:param {ele_id}|{which} {who}
!
! Where: 
!   {ele_id} is an element name or index.
!   {which} is one of: "model", "base" or "design"
!   {who} values are the same as {who} values for "pipe lat_list".
!         Note: Here {who} must be a single parameter and not a list.
!
! Example:
!   pipe ele:param 3@1>>7|model e_tot
! This gives E_tot of element number 7 in branch 1 of universe 3.
!
! Note: On output the {variable} component will always be "F" (since this 
! command cannot tell if a parameter is allowed to vary).
!
! Also see: "pipe lat_list".
!
! Parameters
! ----------
! ele_id
! who 
! which : default=model
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_photon
!  args:
!   ele_id: 1@0>>1
!   which: model
!   who: orbit.vec.1

case ('ele:param')

  u => point_to_uni(line, .true., err); if (err) return
  tao_lat => point_to_tao_lat(line, u, err, which, tail_str); if (err) return
  ele => point_to_ele(line, tao_lat%lat, err); if (err) return
  orbit => tao_lat%tao_branch(ele%ix_branch)%orbit(ele%ix_ele)

  data_type = is_real$

  select case (tail_str)
  case ('ele.mat6')
    n_add = 36
    do ix = 1, 6
      values(6*(ix-1)+1:6*ix) = ele%mat6(ix,:)
    enddo
  case ('ele.vec0')
    n_add = 6
    values(1:6) = ele%vec0
  case ('ele.c_mat')
    n_add = 4
    values(1:4) = [ele%c_mat(1,1), ele%c_mat(1,2), ele%c_mat(2,1), ele%c_mat(2,2)]
  case default
    n_add = 1
    values(1) = ele_param_value(tail_str, ele, orbit, data_type, err); if (err) return
  end select


  select case (data_type)
  case (is_real$)
    nl=incr(nl); write (li(nl), amt) trim(tail_str) // ';REAL;F',    (';', re_str(values(k), 8), k = 1, n_add)
  case (is_integer$)
    nl=incr(nl); write (li(nl), imt) trim(tail_str) // ';INT;F;',     nint(values(1))
  end select

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% ele:photon
!
! Output element photon parameters
!
! Notes
! -----
! Command syntax:
!   pipe ele:photon {ele_id}|{which} {who}
!
! Where: 
!   {ele_id} is an element name or index.
!   {which} is one of: "model", "base" or "design"
!   {who} is one of: "base", "material", or "curvature"
!
! Example:
!   pipe ele:photon 3@1>>7|model base
! This gives element number 7 in branch 1 of universe 3.
! 
! Parameters
! ----------
! ele_id
! who 
! which : default=model
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_photon
!  args:
!   ele_id: 1@0>>1
!   which: model
!   who: base
 

case ('ele:photon')

  u => point_to_uni(line, .true., err); if (err) return
  tao_lat => point_to_tao_lat(line, u, err, which, tail_str); if (err) return
  ele => point_to_ele(line, tao_lat%lat, err); if (err) return

  if (.not. associated(ele%photon)) then
    call invalid ('photon structure not allocated for element.')
    return
  endif

  ph => ele%photon
  select case (tail_str)
  case ('base')
    nl=incr(nl); write (li(nl), lmt) 'has#pixel;LOGIC;F;',  (allocated(ele%photon%pixel%pt))
    nl=incr(nl); write (li(nl), lmt) 'has#material;LOGIC;F;', &
                           (attribute_name(ele, material_type$) == 'MATERIAL_TYPE' .or. ele%key == crystal$)

  case ('material')
    if (ele%key == multilayer_mirror$) then
      nl=incr(nl); write (li(nl), amt) 'F0_m1;COMPLEX;F',          cmplx_str(ph%material%f0_m1)
      nl=incr(nl); write (li(nl), amt) 'F0_m2;COMPLEX;F',          cmplx_str(ph%material%f0_m2)
    else
      nl=incr(nl); write (li(nl), amt) 'F0_m2;COMPLEX;F',          cmplx_str(ph%material%f0_m2)
    endif
    nl=incr(nl); write (li(nl), amt) 'F_H;COMPLEX;F',              cmplx_str(ph%material%f_h)
    nl=incr(nl); write (li(nl), amt) 'F_Hbar;COMPLEX;F',           cmplx_str(ph%material%f_hbar)
    nl=incr(nl); write (li(nl), amt) 'Sqrt(F_H*F_Hbar);COMPLEX;F', cmplx_str(ph%material%f_hkl)

  case ('curvature')
    nl=incr(nl); write (li(nl), rmt) 'spherical_curvature;REAL;T;',      ph%curvature%spherical
    nl=incr(nl); write (li(nl), ramt) 'elliptical_curvature;REAL_ARR;T', (';', ph%curvature%elliptical(i), i = 1, 3)
    do i = 0, ubound(ph%curvature%xy, 2)
      nl=incr(nl); write (li(nl), ramt) 'xy(' // int_str(i) // ',:);REAL_ARR;T', (';', ph%curvature%xy(i,j), j = 0, ubound(ph%curvature%xy, 1))
    enddo
  end select

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% ele:spin_taylor
!
! Output element spin_taylor parameters
!
! Notes
! -----
! Command syntax:
!   pipe ele:spin_taylor {ele_id}|{which}
!
! Where: 
!   {ele_id} is an element name or index.
!   {which} is one of: "model", "base" or "design"
!
! Example:
!   pipe ele:spin_taylor 3@1>>7|model
! This gives element number 7 in branch 1 of universe 3.
! 
! Parameters
! ----------
! ele_id
! which : default=model
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_spin
!  args:
!   ele_id: 1@0>>2
!   which: model
case ('ele:spin_taylor')

  u => point_to_uni(line, .true., err); if (err) return
  tao_lat => point_to_tao_lat(line, u, err, which, tail_str); if (err) return
  ele => point_to_ele(line, tao_lat%lat, err); if (err) return

  if (.not. associated(ele%spin_taylor(1)%term)) then
    call invalid('Spin Taylor map not allocated')
    return
  endif

  do i = 0, 3
    do j = 1, size(ele%spin_taylor(i)%term)
      tt => ele%spin_taylor(i)%term(j)
      nl=incr(nl); write (li(nl), '(i0, a, es22.14, 6(a, i0))') i, ';term;', tt%coef, (';', tt%expn(k), k = 1, 6)
    enddo
  enddo

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% ele:taylor
!
! Output element taylor map 
!
! Notes
! -----
! Command syntax:
!   pipe ele:taylor {ele_id}|{which}
!
! Where: 
!   {ele_id} is an element name or index.
!   {which} is one of: "model", "base" or "design"
!
! Example:
!   pipe ele:taylor 3@1>>7|model
! This gives element number 7 in branch 1 of universe 3.
! 
! Parameters
! ----------
! ele_id
! which : default=model
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_taylor
!  args:
!   ele_id: 1@0>>34
!   which: model

case ('ele:taylor')

  u => point_to_uni(line, .true., err); if (err) return
  tao_lat => point_to_tao_lat(line, u, err, which, tail_str); if (err) return
  ele => point_to_ele(line, tao_lat%lat, err); if (err) return

  if (attribute_name(ele, taylor_map_includes_offsets$) == 'TAYLOR_MAP_INCLUDES_OFFSETS') then
    nl=incr(nl); write (li(nl), lmt) 'taylor_map_includes_offsets;LOGIC;T;',    ele%taylor_map_includes_offsets
  endif

  if (.not. associated(ele%taylor(1)%term)) then
    call invalid('Taylor map not allocated')
    return
  endif

  do i = 1, 6
    nl=incr(nl); write (li(nl), '(i0, a, es22.14)') i, ';ref;', ele%taylor(i)%ref
    do j = 1, size(ele%taylor(i)%term)
      tt => ele%taylor(i)%term(j)
      nl=incr(nl); write (li(nl), '(2(i0, a), es22.14, 6(a, i0))') i, ';', j, ';', tt%coef, (';', tt%expn(k), k = 1, 6)
    enddo
  enddo

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% ele:twiss
!
! Output element Twiss parameters
!
! Notes
! -----
! Command syntax:
!   pipe ele:twiss {ele_id}|{which}
!
! Where: 
!   {ele_id} is an element name or index.
!   {which} is one of: "model", "base" or "design"
!
! Example:
!   pipe ele:twiss 3@1>>7|model
! This gives element number 7 in branch 1 of universe 3.
! 
! Parameters
! ----------
! ele_id
! which : default=model
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!   ele_id: 1@0>>1
!   which: model

case ('ele:twiss')

  u => point_to_uni(line, .true., err); if (err) return
  tao_lat => point_to_tao_lat(line, u, err, which, tail_str); if (err) return
  ele => point_to_ele(line, tao_lat%lat, err); if (err) return

  if (ele%a%beta == 0) return
  free = attribute_free(ele, 'BETA_A', .false.) .and. (which == 'model')

  nl=incr(nl); write (li(nl), lmt) 'mode_flip;LOGIC;F;', ele%mode_flip

  call twiss_out (ele%a, '', 'a', can_vary = free)
  call twiss_out (ele%b, '', 'b', can_vary = free)
  call xy_disp_out (ele%x, 'x', can_vary = free)
  call xy_disp_out (ele%y, 'y', can_vary = free)

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% ele:wake
!
! Output element wake.
!
! Notes
! -----
! Command syntax:
!   pipe ele:wake {ele_id}|{which} {who}
!
! Where: 
!   {ele_id} is an element name or index.
!   {which} is one of: "model", "base" or "design"
!   {Who} is one of:
!       "sr_long"        "sr_long_table"
!       "sr_trans"       "sr_trans_table"
!       "lr_mode_table"  "base"
!
! Example:
!   pipe ele:wake 3@1>>7|model
! This gives element number 7 in branch 1 of universe 3.
! 
! Parameters
! ----------
! ele_id
! who
! which : default=model
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_wake
!  args:
!   ele_id: 1@0>>1
!   which: model
!   who: sr_long

case ('ele:wake')

  u => point_to_uni(line, .true., err); if (err) return
  tao_lat => point_to_tao_lat(line, u, err, which, tail_str); if (err) return
  ele => point_to_ele(line, tao_lat%lat, err); if (err) return

  if (.not. associated(ele%wake)) then
    call invalid ('No wake associated')
    return
  endif

  wake => ele%wake

  select case (tail_str)
  case ('base')
    nl=incr(nl); write (li(nl), rmt) 'sr%z_max;REAL;T;',              wake%sr%z_max
    nl=incr(nl); write (li(nl), rmt) 'sr%amp_scale;REAL;T;',          wake%sr%amp_scale
    nl=incr(nl); write (li(nl), rmt) 'sr%z_scale;REAL;T;',            wake%sr%z_scale
    nl=incr(nl); write (li(nl), lmt) 'sr%scale_with_length;LOGIC;T;', wake%sr%scale_with_length
    nl=incr(nl); write (li(nl), rmt) 'lr%freq_spread;REAL;T;',        wake%lr%freq_spread
    nl=incr(nl); write (li(nl), rmt) 'lr%amp_scale;REAL;T;',          wake%lr%amp_scale
    nl=incr(nl); write (li(nl), rmt) 'lr%time_scale;REAL;T;',         wake%lr%time_scale
    nl=incr(nl); write (li(nl), lmt) 'lr%self_wake_on;LOGIC;T;',      wake%lr%self_wake_on
    nl=incr(nl); write (li(nl), lmt) 'has#sr_long;LOGIC;F;',          allocated(wake%sr%long)
    nl=incr(nl); write (li(nl), lmt) 'has#sr_trans;LOGIC;F;',         allocated(wake%sr%trans)
    nl=incr(nl); write (li(nl), lmt) 'has#lr_mode;LOGIC;F;',          allocated(wake%lr%mode)

  case ('sr_long')
    nl=incr(nl); write (li(nl), rmt) 'z_ref;REAL;T;',   wake%sr%z_ref_long

  case ('sr_long_table')
    do i = 1, size(wake%sr%long)
      wsr => wake%sr%long(i)
      nl=incr(nl); write (li(nl), '(4(es16.8), 4a)') wsr%amp, ';', wsr%damp, ';', wsr%k, ';', wsr%phi, ';', &
          sr_longitudinal_position_dep_name(wsr%position_dependence)
    enddo

  case ('sr_trans')
    nl=incr(nl); write (li(nl), rmt) 'z_ref;REAL;T;',   wake%sr%z_ref_trans

  case ('sr_trans_table')
    do i = 1, size(wake%sr%trans)
      wsr => wake%sr%trans(i)
      nl=incr(nl); write (li(nl), '(4(es16.8), 4a)') wsr%amp, ';', wsr%damp, ';', wsr%k, ';', wsr%phi, ';', &
          sr_transverse_polarization_name(wsr%polarization), ';', sr_transverse_position_dep_name(wsr%position_dependence)
    enddo

  case ('lr_mode_table')
    do i = 1, size(wake%lr%mode)
      lr_mode => wake%lr%mode(i)
      v_str = 'none'
      if (lr_mode%polarized) write (v_str, '(f8.3)') lr_mode%angle
      nl=incr(nl); write (li(nl), '(4(es22.14, a), 2a)') &
                lr_mode%freq, ';', lr_mode%R_over_Q, ';', lr_mode%Q, ';', lr_mode%m, ';', v_str
    enddo

  case default
    call invalid ('Bad or missign {who} switch.')
    return
  end select

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% ele:wall3d
!
! Output element wall3d parameters.
!
! Notes
! -----
! Command syntax:
!   pipe ele:wall3d {ele_id}|{which} {index} {who}
!
! Where: 
!   {ele_id} is an element name or index.
!   {which} is one of: "model", "base" or "design"
!   {index} is the index number in the ele%wall3d(:) array 
!             The size of this array is obtained from "pipe ele:head".
!   {who} is one of: "base", or "table".
! Example:
!   pipe ele:wall3d 3@1>>7|model 2 base
! This gives element number 7 in branch 1 of universe 3.
! 
! Parameters
! ----------
! ele_id
! index 
! who 
! which : default=model
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_wall3d
!  args:
!   ele_id: 1@0>>1
!   which: model
!   index: 1
!   who: table

case ('ele:wall3d')

  u => point_to_uni(line, .true., err); if (err) return
  tao_lat => point_to_tao_lat(line, u, err, which, tail_str); if (err) return
  ele => point_to_ele(line, tao_lat%lat, err); if (err) return

  if (.not. associated(ele%wall3d)) then
    call invalid ('wall3d not allocated')
    return
  endif
  ix = parse_int (tail_str, err, 1, size(ele%wall3d));  if (err) return
  wall3d => ele%wall3d(ix)

  select case (tail_str)
  case ('base')
    nl=incr(nl); write (li(nl), amt) 'name;STR;T;',                    trim(wall3d%name)
    nl=incr(nl); write (li(nl), amt) 'ele_anchor_pt;ENUM;T;',          trim(anchor_pt_name(wall3d%ele_anchor_pt))
    select case (ele%key)
    case (capillary$)
    case (diffraction_plate$, mask$)
      nl=incr(nl); write (li(nl), rmt) 'thickness;REAL;T;',            wall3d%thickness
      nl=incr(nl); write (li(nl), amt) 'clear_material;SPECIES;T;',    trim(wall3d%clear_material)
      nl=incr(nl); write (li(nl), amt) 'opaque_material;SPECIES;T;',   trim(wall3d%opaque_material)
    case default
      nl=incr(nl); write (li(nl), lmt) 'superimpose;LOGIC;T;',         wall3d%superimpose
    end select

  case ('table')
    do i = 1, size(wall3d%section)
      sec => wall3d%section(i)
      nl=incr(nl); write (li(nl), imt)   'section;INT;F;',    i
      nl=incr(nl); write (li(nl), rmt)   's;REAL;T;',         sec%s
      nl=incr(nl); write (li(nl), ramt)  'r0;REAL_ARR;T',    (';', sec%r0(j), j = 1, size(sec%r0))
      if (ele%key /= capillary$) then
        nl=incr(nl); write (li(nl), amt) 'wall3d_section^type;ENUM;T;',    trim(wall3d_section_type_name(sec%type))
      endif
      nl=incr(nl); write (li(nl), imt)   'vertex;INT;F;',    i
      do j = 1, size(sec%v)
        nl=incr(nl); write (li(nl), '(i0, 5(a, es22.14))') j, ';', sec%v(j)%x, ';', sec%v(j)%y, ';', &
                                              sec%v(j)%radius_x, ';', sec%v(j)%radius_y, ';', sec%v(j)%tilt
      enddo
    enddo

  case default
    call invalid ('Bad or missign {who} switch.')
    return
  end select

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% evaluate
!
! Output the value of an expression. The result may be a vector.
!
! Notes
! -----
! Command syntax:
!   pipe evaluate {flags} {expression}
!
! Where:
!   Optional {flags} are:
!       -array_out : If present, the output will be available in the 
!                     tao_c_interface_com%c_real array.
!   {expression} is expression to be evaluated.
!
! Example:
!   pipe evaluate 3+data::cbar.11[1:10]|model
! 
! Parameters
! ----------
! expression :
! flags : default=-array_out
!   If -array_out, the output will be available in the tao_c_interface_com%c_real.
!
! Returns
! -------
! string_list
!   if '-array_out' not in flags
! real_array
!   if '-array_out' in flags
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!    expression: data::cbar.11[1:10]|model
 

case ('evaluate')

  use_real_array_buffer = .false.

  if (index('-array_out', line(1:ix_line)) == 1) then
    call string_trim(line(ix_line+1:), line, ix_line)
    use_real_array_buffer = .true.
  endif

  if (index(line, 'chrom') /= 0 .or. index(line, 'rad') /= 0) then
    if (index(line, 'chrom') /= 0) s%com%force_chrom_calc = .true.
    if (index(line, 'rad') /= 0) s%com%force_rad_int_calc = .true.
    s%u%calc%lattice = .true.
    call tao_lattice_calc(ok)
  endif

  call tao_evaluate_expression (line, 0, .false., value_arr, err)
  if (err) then
    call invalid ('Invalid expression')
    return
  endif

  call real_array_out (value_arr, use_real_array_buffer)

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% em_field
!
! Output EM field at a given point generated by a given element.
!
! Notes
! -----
! Command syntax:
!   pipe em_field {ele_id}|{which} {x} {y} {z} {t_or_z}
!
! Where:
!   {which} is one of: "model", "base" or "design"
!   {x}, {y}  -- Transverse coords.
!   {z}       -- Longitudinal coord with respect to entrance end of element.
!   {t_or_z}  -- time or phase space z depending if lattice is setup for 
!             --   absolute time tracking.
! 
! Parameters
! ----------
! ele_id
! x
! y
! z
! t_or_z
! which : default=model
!    
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!    ele_id: 1@0>>22
!    which: model
!    x: 0
!    y: 0
!    z: 0
!    t_or_z: 0

case ('em_field')

  u => point_to_uni(line, .true., err); if (err) return
  tao_lat => point_to_tao_lat(line, u, err, which, tail_str); if (err) return
  ele => point_to_ele(line, tao_lat%lat, err); if (err) return

  call init_coord (orb, ele = ele, element_end = downstream_end$)

  call split_this_line (tail_str, name1, 4, err, space_sep = .true.); if (err) return
  orb%vec(1)  = parse_real(name1(1), err);  if (err) return
  orb%vec(3)  = parse_real(name1(2), err);  if (err) return
  z           = parse_real(name1(3), err);  if (err) return
  if (bmad_com%absolute_time_tracking) then
    orb%t = parse_real(name1(4), err);  if (err) return
  else
    orb%vec(5) = parse_real(name1(4), err);  if (err) return
  endif

  call em_field_calc (ele, ele%branch%param, z, orb, .false., field, err_flag = err);  if (err) return
  nl=incr(nl); write (li(nl), '(6(es22.14, a))') (field%B(i), ';',  i = 1, 3), (field%E(i), ';',  i = 1, 2), field%E(3)

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% enum
!
! Output list of possible values for enumerated numbers.
!
! Notes
! -----
! Command syntax:
!   pipe enum {enum_name}
!
! Example:
!   pipe enum tracking_method
! 
! Parameters
! ----------
! enum_name
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!    enum_name: tracking_method

case ('enum')

  if (index(line, 'color') /= 0) then
    do i = lbound(qp_color_name, 1), ubound(qp_color_name, 1)
      nl=incr(nl); write(li(nl), '(i0, 2a)') i, ';', trim(qp_color_name(i))
    enddo
    call end_stuff(li, nl)
    return
  endif

  select case (line)
  case ('axis^type')
    nl=incr(nl); write(li(nl), '(a)') '1;LINEAR'
    nl=incr(nl); write(li(nl), '(a)') '2;LOG'

  case ('bounds')
    nl=incr(nl); write(li(nl), '(a)') '1;GENERAL'
    nl=incr(nl); write(li(nl), '(a)') '2;ZERO_AT_END'
    nl=incr(nl); write(li(nl), '(a)') '3;ZERO_SYMMETRIC'

  case ('building^constraint')
    nl=incr(nl); write(li(nl), '(a)') '1;none'
    nl=incr(nl); write(li(nl), '(a)') '2;left_side'
    nl=incr(nl); write(li(nl), '(a)') '3;right_side'

  case ('data^merit_type')
    do i = 1, size(tao_data_merit_type_name)
      nl=incr(nl); write(li(nl), '(i0, 2a)') i, ';', trim(tao_data_merit_type_name(i))
    enddo

  case ('data_source')
    do i = 1, size(tao_data_source_name)
      nl=incr(nl); write(li(nl), '(i0, 2a)') i, ';', trim(tao_data_source_name(i))
    enddo

  case ('distribution_type')
    do i = 1, size(beam_distribution_type_name)
      nl=incr(nl); write(li(nl), '(i0, 2a)') i, ';', trim(beam_distribution_type_name(i))
    enddo

  case ('floor_plan_view_name')
    do i = 1, size(tao_floor_plan_view_name)
      nl=incr(nl); write(li(nl), '(i0, 2a)') i, ';', trim(tao_floor_plan_view_name(i))
    enddo

  case ('graph^type')
    do i = 1, size(tao_graph_type_name)
      nl=incr(nl); write(li(nl), '(i0, 2a)') i, ';', trim(tao_graph_type_name(i))
    enddo

  case ('line^pattern', 'orbit_pattern')
    do i = 1, size(qp_line_pattern_name)
      nl=incr(nl); write(li(nl), '(i0, 2a)') i, ';', trim(qp_line_pattern_name(i))
    enddo

  case ('lord_status')
    nl=incr(nl); li(nl) = '4;Group_Lord'
    nl=incr(nl); li(nl) = '5;Super_Lord' 
    nl=incr(nl); li(nl) = '6;Overlay_Lord' 
    nl=incr(nl); li(nl) = '7;Girder_Lord' 
    nl=incr(nl); li(nl) = '8;Multipass_Lord'
    nl=incr(nl); li(nl) = '10;Not_a_Lord' 
    nl=incr(nl); li(nl) = '12;Control_Lord' 
    nl=incr(nl); li(nl) = '13;Ramper_Lord'

  case ('optimizer')
    do i = 1, size(tao_optimizer_name)
      nl=incr(nl); write(li(nl), '(i0, 2a)') i, ';', trim(tao_optimizer_name(i))
    enddo

  case ('orbit_lattice')
    nl=incr(nl); li(nl) = '1;model'
    nl=incr(nl); li(nl) = '2;design'
    nl=incr(nl); li(nl) = '3;base'

  case ('photon_type')
    do i = 1, size(photon_type_name)
      nl=incr(nl); write(li(nl), '(i0, 2a)') i, ';', trim(photon_type_name(i))
    enddo

  case ('plot^type')
    nl=incr(nl); li(nl) = '1;normal'
    nl=incr(nl); li(nl) = '2;wave'

  case ('random_engine')
    nl=incr(nl); li(nl) = '1;pseudo'
    nl=incr(nl); li(nl) = '2;quasi'

  case ('random_gauss_converter')
    nl=incr(nl); li(nl) = '1;exact'
    nl=incr(nl); li(nl) = '2;quick'

  case ('shape^label')
    do i = 1, size(tao_shape_label_name)
      nl=incr(nl); write(li(nl), '(i0, 2a)') i, ';', trim(tao_shape_label_name(i))
    enddo

  case ('shape^shape')
    do i = 1, size(tao_shape_shape_name)
      nl=incr(nl); write(li(nl), '(i0, 2a)') i, ';', trim(tao_shape_shape_name(i))
    enddo

  case ('slave_status')
    nl=incr(nl); li(nl) = '1;Minor_Slave'
    nl=incr(nl); li(nl) = '2;Super_Slave' 
    nl=incr(nl); li(nl) = '3;Free' 
    nl=incr(nl); li(nl) = '9;Multipass_Slave' 
    nl=incr(nl); li(nl) = '11;Slice_Slave' 

  case ('fill_pattern')
    do i = 1, size(qp_symbol_fill_pattern_name)
      nl=incr(nl); write(li(nl), '(i0, 2a)') i, ';', trim(qp_symbol_fill_pattern_name(i))
    enddo

  case ('symbol^type')
    do i = lbound(qp_symbol_type_name, 1), ubound(qp_symbol_type_name, 1)
      nl=incr(nl); write(li(nl), '(i0, 2a)') i, ';', trim(qp_symbol_type_name(i))
    enddo

  case ('track_type')
    nl=incr(nl); li(nl) = 'single'
    nl=incr(nl); li(nl) = 'beam'

  case ('var^merit_type')
    do i = 1, size(tao_var_merit_type_name)
      nl=incr(nl); write(li(nl), '(i0, 2a)') i, ';', trim(tao_var_merit_type_name(i))
    enddo

  case ('view')
    nl=incr(nl); li(nl) = 'zx'
    nl=incr(nl); li(nl) = 'xz'
    nl=incr(nl); li(nl) = 'xy'
    nl=incr(nl); li(nl) = 'yx'
    nl=incr(nl); li(nl) = 'zy'
    nl=incr(nl); li(nl) = 'yz'

  case ('wave_data_type')
    do i = 1, size(tao_wave_data_name)
      nl=incr(nl); write(li(nl), '(i0, 2a)') i, ';', trim(tao_wave_data_name(i))
    enddo

  case ('x_axis_type')
    do i = 1, size(tao_x_axis_type_name)
      nl=incr(nl); write(li(nl), '(i0, 2a)') i, ';', trim(tao_x_axis_type_name(i))
    enddo

  case ('data_type_z')
    do i = 1, size(tao_data_type_z_name)
      nl=incr(nl); write(li(nl), '(i0, 2a)') i, ';', trim(tao_data_type_z_name(i))
    enddo

  case default

    name = upcase(line)
    if (name == 'EVAL_POINT') name = 'ELE_ORIGIN'  ! Cheat since data%eval_point is not recognized by switch_attrib_value_name

    a_name = switch_attrib_value_name(name, 1.0_rp, this_ele, name_list = name_list)
    if (.not. allocated(name_list)) then
      call invalid ('Not a valid switch name.')
      return
    endif

    do i = lbound(name_list, 1), ubound(name_list, 1)
      if (index(name_list(i), '!') /= 0 .or. name_list(i) == '') cycle
      nl=incr(nl); write(li(nl), '(i0, 2a)') i, ';', trim(name_list(i))
    enddo
  end select

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% floor_plan
!
! Output (x,y) points and other information that can be used for drawing a floor_plan.
!
! Notes
! -----
! Command syntax:
!   pipe floor_plan {graph}
! 
! Parameters
! ----------
! graph
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_optics_matching
!  args:
!    graph: r13.g

case ('floor_plan')

  call tao_find_plots (err, line, 'BOTH', graph = graphs, blank_means_all = .true., only_visible = .false.)

  if (err .or. size(graphs) /= 1) then
    call invalid ('Bad graph name')
    return
  endif

  g => graphs(1)%g

  if (g%ix_universe == -2) then
    do iu = 1, size(s%u)
      call this_floor_plan(iu, g)
    enddo
  else
    iu = tao_universe_index(g%ix_universe)
    call this_floor_plan(iu, g)
  endif

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% floor_orbit
!
! Output (x, y) coordinates for drawing the particle orbit on a floor plan.
!
! Notes
! -----
! Command syntax:
!   pipe floor_orbit {graph}
! 
! Parameters
! ----------
! graph
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_floor_orbit
!  args:
!    graph: r33.g 

case ('floor_orbit')

  call tao_find_plots (err, line, 'REGION', graph = graphs, only_visible = .false.)

  if (err .or. size(graphs) /= 1) then
    call invalid ('Bad graph name')
    return
  endif

  g => graphs(1)%g
  u => tao_pointer_to_universe(g%ix_universe)
  lat => u%model%lat

  if (g%floor_plan%orbit_scale == 0) then
    call invalid ('graph%floor_plan%orbit_scale is zero!')
    return
  endif

  do ib = 0, ubound(lat%branch, 1)
    branch => lat%branch(ib)
    do ie = 1, branch%n_ele_track
      ele => branch%ele(ie)
      if (ele%value(l$) == 0 .and. ele%key /= patch$) cycle

      orb_start = u%model%tao_branch(ele%ix_branch)%orbit(ele%ix_ele-1)
      orb_end   = u%model%tao_branch(ele%ix_branch)%orbit(ele%ix_ele)

      floor%r = [0.0_rp, 0.0_rp, 0.0_rp]
      floor1 = coords_local_curvilinear_to_floor (floor, ele, .true.)

      floor%r = [0.0_rp, 0.0_rp, ele%value(l$)]
      floor2 = coords_local_curvilinear_to_floor (floor, ele, .true.)

      call tao_floor_to_screen_coords (g, floor1, end1)
      call tao_floor_to_screen_coords (g, floor2, end2)

      ! Bends can be tricky if they are not in the X-Z plane.
      ! Bends are parameterized by a set of points (x_bend, y_bend) along their
      ! centerline and a set of vectors (dx_bend, dy_bend) tangent to the centerline.

      if (ele%key == sbend$) then

        ! Start at entrance end (not upstream end)
        if (ele%orientation == 1) then
          floor = floor1
        else
          floor = floor2
        endif

        v_old = floor%r
        call floor_angles_to_w_mat (floor%theta, floor%phi, 0.0_rp, w_old)

        n_bend = min(abs(int(100 * ele%value(angle$))) + 1, ubound(x_bend, 1))
        if (n_bend < 1) return   ! A crazy angle can cause int(100*angle) to be negative !!
        ang    = ele%value(angle$) * ele%orientation
        length = ele%value(l$)     * ele%orientation

        do j = 0, n_bend
          angle = j * ang / n_bend
          cos_t = cos(ele%value(ref_tilt_tot$))
          sin_t = sin(ele%value(ref_tilt_tot$))
          cos_a = cos(angle)
          sin_a = sin(angle)
          if (ele%value(g$) == 0) then
            r_vec = length * j * [0, 0, 1]
          else
            r_vec = ele%value(rho$) * [cos_t * (cos_a - 1), sin_t * (cos_a - 1), sin_a]
          endif
          dr_vec = [-cos_t * sin_a, -sin_t * sin_a, cos_a]
          ! This keeps dr_vec pointing to the inside (important for the labels).
          if (cos_t < 0) dr_vec = -dr_vec
          v_vec = matmul (w_old, r_vec) + v_old
          dv_vec = matmul (w_old, dr_vec)
          call tao_floor_to_screen (g, v_vec, x_bend(j), y_bend(j))
          call tao_floor_to_screen (g, dv_vec, dx_bend(j), dy_bend(j))

          s_here = j * ele%value(l$) / n_bend
          call twiss_and_track_intra_ele (ele, ele%branch%param, 0.0_rp, s_here, &
                                                           .true., .true., orb_start, orb_here)
          f_orb%r(1:2) = g%floor_plan%orbit_scale * orb_here%vec(1:3:2)
          f_orb%r(3) = s_here
          f_orb = coords_local_curvilinear_to_floor (f_orb, ele, .false.)
          call tao_floor_to_screen (g, f_orb%r, dx_orbit(j), dy_orbit(j))
        enddo

        do ix = 0, 100
          i0 = 50*ix
          i1 = min(50*(ix+1), n_bend)
          nl=incr(nl); write (li(nl), '(2(i0, a), 1000(a, es14.6))') ib, ';', ie, ';x', (';', dx_orbit(i), i = i0, i1)
          nl=incr(nl); write (li(nl), '(2(i0, a), 1000(a, es14.6))') ib, ';', ie, ';y', (';', dy_orbit(i), i = i0, i1)
          if (i1 == n_bend) exit
        enddo

      elseif (ele%key == patch$) then
        ele0 => pointer_to_next_ele (ele, -1)
        floor%r(1:2) = g%floor_plan%orbit_scale * orb_start%vec(1:3:2)
        floor%r(3) = ele0%value(l$)
        floor1 = coords_local_curvilinear_to_floor (floor, ele0, .false.)
        call tao_floor_to_screen_coords (g, floor1, f_orb)
        dx_orbit(0) = f_orb%r(1)
        dy_orbit(0) = f_orb%r(2)

        floor%r(1:2) = g%floor_plan%orbit_scale * orb_end%vec(1:3:2)
        floor%r(3) = ele%value(l$)
        floor1 = coords_local_curvilinear_to_floor (floor, ele, .false.)
        call tao_floor_to_screen_coords (g, floor1, f_orb)
        dx_orbit(1) = f_orb%r(1)
        dy_orbit(1) = f_orb%r(2)

      else
        n = 10
        do ic = 0, n
          s_here = ic * ele%value(l$) / n
          call twiss_and_track_intra_ele (ele, ele%branch%param, 0.0_rp, s_here, &
                                                     .true., .true., orb_start, orb_here)
          floor%r(1:2) = g%floor_plan%orbit_scale * orb_here%vec(1:3:2)
          floor%r(3) = s_here
          floor1 = coords_local_curvilinear_to_floor (floor, ele, .false.)
          call tao_floor_to_screen_coords (g, floor1, f_orb)
          dx_orbit(ic) = f_orb%r(1)
          dy_orbit(ic) = f_orb%r(2)
        enddo

        nl=incr(nl); write (li(nl), '(2(i0, a), 1000(a, es14.6))') ib, ';', ie, ';x', (';', dx_orbit(i), i = 0, n)
        nl=incr(nl); write (li(nl), '(2(i0, a), 1000(a, es14.6))') ib, ';', ie, ';y', (';', dy_orbit(i), i = 0, n)
      endif
    enddo
  enddo

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% global
!
! Output global parameters.
!
! Notes
! -----
! Command syntax:
!   pipe global
!
! Output syntax is parameter list form. See documentation at the beginning of this file.
! 
! Note: The follow is intentionally left out:
!   optimizer_allow_user_abort
!   quiet
!   single_step
!   prompt_color
!   prompt_string
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:

case ('global')

  nl=incr(nl); write (li(nl), rmt) 'de_lm_step_ratio;REAL;T;',                s%global%de_lm_step_ratio
  nl=incr(nl); write (li(nl), rmt) 'de_var_to_population_factor;REAL;T;',     s%global%de_var_to_population_factor
  nl=incr(nl); write (li(nl), rmt) 'lm_opt_deriv_reinit;REAL;T;',             s%global%lm_opt_deriv_reinit
  nl=incr(nl); write (li(nl), rmt) 'lmdif_eps;REAL;T;',                       s%global%lmdif_eps
  nl=incr(nl); write (li(nl), rmt) 'lmdif_negligible_merit;REAL;T;',          s%global%lmdif_negligible_merit
  nl=incr(nl); write (li(nl), rmt) 'svd_cutoff;REAL;T;',                      s%global%svd_cutoff
  nl=incr(nl); write (li(nl), rmt) 'unstable_penalty;REAL;T;',                s%global%unstable_penalty
  nl=incr(nl); write (li(nl), rmt) 'merit_stop_value;REAL;T;',                s%global%merit_stop_value
  nl=incr(nl); write (li(nl), rmt) 'dmerit_stop_value;REAL;T;',               s%global%dmerit_stop_value
  nl=incr(nl); write (li(nl), rmt) 'random_sigma_cutoff;REAL;T;',             s%global%random_sigma_cutoff
  nl=incr(nl); write (li(nl), rmt) 'delta_e_chrom;REAL;T;',                   s%global%delta_e_chrom
  nl=incr(nl); write (li(nl), imt) 'n_opti_cycles;INT;T;',                    s%global%n_opti_cycles
  nl=incr(nl); write (li(nl), imt) 'n_opti_loops;INT;T;',                     s%global%n_opti_loops
  nl=incr(nl); write (li(nl), amt) 'phase_units;ENUM;T;',                     trim(angle_units_name(s%global%phase_units))
  nl=incr(nl); write (li(nl), imt) 'bunch_to_plot;INT;T;',                    s%global%bunch_to_plot
  nl=incr(nl); write (li(nl), imt) 'random_seed;INT;T;',                      s%global%random_seed
  nl=incr(nl); write (li(nl), imt) 'n_threads;INT;T;',                        s%global%n_threads
  nl=incr(nl); write (li(nl), imt) 'n_top10_merit;INT;T;',                    s%global%n_top10_merit
  nl=incr(nl); write (li(nl), imt) 'n_opti_loops;INT;T;',                     s%global%n_opti_loops
  nl=incr(nl); write (li(nl), imt) 'n_opti_cycles;INT;T;',                    s%global%n_opti_cycles
  nl=incr(nl); write (li(nl), imt) 'srdt_gen_n_slices;INT;T;',                s%global%srdt_gen_n_slices
  nl=incr(nl); write (li(nl), imt) 'srdt_sxt_n_slices;INT;T;',                s%global%srdt_sxt_n_slices
  nl=incr(nl); write (li(nl), lmt) 'srdt_use_cache;LOGIC;T;',                 s%global%srdt_use_cache
  nl=incr(nl); write (li(nl), amt) 'random_engine;STR;T;',                    trim(s%global%random_engine)
  nl=incr(nl); write (li(nl), amt) 'random_gauss_converter;STR;T;',           trim(s%global%random_gauss_converter)
  nl=incr(nl); write (li(nl), amt) 'track_type;ENUM;T;',                      trim(s%global%track_type)
  nl=incr(nl); write (li(nl), amt) 'optimizer;ENUM;T;',                       trim(s%global%optimizer)
  nl=incr(nl); write (li(nl), amt) 'print_command;STR;T;',                    trim(s%global%print_command)
  nl=incr(nl); write (li(nl), amt) 'var_out_file;FILE;T;',                    trim(s%global%var_out_file)

  nl=incr(nl); write (li(nl), lmt) 'external_plotting;LOGIC;I;',              s%global%external_plotting
  nl=incr(nl); write (li(nl), lmt) 'opt_with_ref;LOGIC;T;',                   s%global%opt_with_ref
  nl=incr(nl); write (li(nl), lmt) 'opt_with_base;LOGIC;T;',                  s%global%opt_with_base
  nl=incr(nl); write (li(nl), lmt) 'label_lattice_elements;LOGIC;T;',         s%global%label_lattice_elements
  nl=incr(nl); write (li(nl), lmt) 'label_keys;LOGIC;T;',                     s%global%label_keys
  nl=incr(nl); write (li(nl), lmt) 'concatenate_maps;LOGIC;T;',               s%global%concatenate_maps
  nl=incr(nl); write (li(nl), lmt) 'derivative_recalc;LOGIC;T;',              s%global%derivative_recalc
  nl=incr(nl); write (li(nl), lmt) 'derivative_uses_design;LOGIC;T;',         s%global%derivative_uses_design
  nl=incr(nl); write (li(nl), lmt) 'plot_on;LOGIC;T;',                        s%global%plot_on
  nl=incr(nl); write (li(nl), lmt) 'lattice_calc_on;LOGIC;T;',                s%global%lattice_calc_on
  nl=incr(nl); write (li(nl), lmt) 'svd_retreat_on_merit_increase;LOGIC;T;',  s%global%svd_retreat_on_merit_increase
  nl=incr(nl); write (li(nl), lmt) 'stop_on_error;LOGIC;T;',                  s%global%stop_on_error
  nl=incr(nl); write (li(nl), lmt) 'box_plots;LOGIC;T;',                      s%global%box_plots
  nl=incr(nl); write (li(nl), lmt) 'beam_timer_on;LOGIC;T;',                  s%global%beam_timer_on
  nl=incr(nl); write (li(nl), lmt) 'var_limits_on;LOGIC;T;',                  s%global%var_limits_on
  nl=incr(nl); write (li(nl), lmt) 'only_limit_opt_vars;LOGIC;T;',            s%global%only_limit_opt_vars
  nl=incr(nl); write (li(nl), lmt) 'opt_match_auto_recalc;LOGIC;T;',          s%global%opt_match_auto_recalc
  nl=incr(nl); write (li(nl), lmt) 'opti_write_var_file;LOGIC;T;',            s%global%opti_write_var_file
  nl=incr(nl); write (li(nl), lmt) 'optimizer_var_limit_warn;LOGIC;T;',       s%global%optimizer_var_limit_warn
  nl=incr(nl); write (li(nl), lmt) 'optimizer_allow_user_abort;LOGIC;T;',     s%global%optimizer_allow_user_abort
  nl=incr(nl); write (li(nl), lmt) 'rf_on;LOGIC;T;',                          s%global%rf_on
  nl=incr(nl); write (li(nl), lmt) 'symbol_import;LOGIC;T;',                  s%global%symbol_import
  nl=incr(nl); write (li(nl), lmt) 'draw_curve_off_scale_warn;LOGIC;T;',      s%global%draw_curve_off_scale_warn
  nl=incr(nl); write (li(nl), lmt) 'wait_for_cr_in_single_mode;LOGIC;T;',     s%global%wait_for_CR_in_single_mode
  nl=incr(nl); write (li(nl), lmt) 'disable_smooth_line_calc;LOGIC;T;',       s%global%disable_smooth_line_calc
  nl=incr(nl); write (li(nl), lmt) 'debug_on;LOGIC;T;',                       s%global%debug_on

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% global:optimization
!
! Output optimization parameters.
! Also see global:opti_de.
!
! Notes
! -----
! Command syntax:
!   pipe global:optimization
!
! Output syntax is parameter list form. See documentation at the beginning of this file.
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:

case ('global:optimization')

  nl=incr(nl); write (li(nl), rmt) 'de_lm_step_ratio;REAL;T;',                s%global%de_lm_step_ratio
  nl=incr(nl); write (li(nl), rmt) 'de_var_to_population_factor;REAL;T;',     s%global%de_var_to_population_factor
  nl=incr(nl); write (li(nl), rmt) 'lm_opt_deriv_reinit;REAL;T;',             s%global%lm_opt_deriv_reinit
  nl=incr(nl); write (li(nl), rmt) 'lmdif_eps;REAL;T;',                       s%global%lmdif_eps
  nl=incr(nl); write (li(nl), rmt) 'lmdif_negligible_merit;REAL;T;',          s%global%lmdif_negligible_merit
  nl=incr(nl); write (li(nl), rmt) 'svd_cutoff;REAL;T;',                      s%global%svd_cutoff
  nl=incr(nl); write (li(nl), rmt) 'unstable_penalty;REAL;T;',                s%global%unstable_penalty
  nl=incr(nl); write (li(nl), rmt) 'merit_stop_value;REAL;T;',                s%global%merit_stop_value
  nl=incr(nl); write (li(nl), rmt) 'dmerit_stop_value;REAL;T;',               s%global%dmerit_stop_value

  nl=incr(nl); write (li(nl), imt) 'n_top10_merit;INT;T;',                    s%global%n_top10_merit
  nl=incr(nl); write (li(nl), imt) 'n_opti_loops;INT;T;',                     s%global%n_opti_loops
  nl=incr(nl); write (li(nl), imt) 'n_opti_cycles;INT;T;',                    s%global%n_opti_cycles

  nl=incr(nl); write (li(nl), amt) 'optimizer;ENUM;T;',                       trim(s%global%optimizer)
  nl=incr(nl); write (li(nl), amt) 'var_out_file;FILE;T;',                    trim(s%global%var_out_file)
  nl=incr(nl); write (li(nl), lmt) 'opti_write_var_file;LOGIC;T;',            s%global%opti_write_var_file

  nl=incr(nl); write (li(nl), lmt) 'derivative_recalc;LOGIC;T;',              s%global%derivative_recalc
  nl=incr(nl); write (li(nl), lmt) 'derivative_uses_design;LOGIC;T;',         s%global%derivative_uses_design
  nl=incr(nl); write (li(nl), lmt) 'opt_with_ref;LOGIC;T;',                   s%global%opt_with_ref
  nl=incr(nl); write (li(nl), lmt) 'opt_with_base;LOGIC;T;',                  s%global%opt_with_base
  nl=incr(nl); write (li(nl), lmt) 'opt_match_auto_recalc;LOGIC;T;',          s%global%opt_match_auto_recalc
  nl=incr(nl); write (li(nl), lmt) 'optimizer_var_limit_warn;LOGIC;T;',       s%global%optimizer_var_limit_warn
  nl=incr(nl); write (li(nl), lmt) 'optimizer_allow_user_abort;LOGIC;T;',     s%global%optimizer_allow_user_abort
  nl=incr(nl); write (li(nl), lmt) 'svd_retreat_on_merit_increase;LOGIC;T;',  s%global%svd_retreat_on_merit_increase
  nl=incr(nl); write (li(nl), lmt) 'var_limits_on;LOGIC;T;',                  s%global%var_limits_on
  nl=incr(nl); write (li(nl), lmt) 'only_limit_opt_vars;LOGIC;T;',            s%global%only_limit_opt_vars


!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% global:opti_de
!
! Output DE optimization parameters.
!
! Notes
! -----
! Command syntax:
!   pipe global:opti_de
!
! Output syntax is parameter list form. See documentation at the beginning of this file.
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:

case ('global:opti_de')

  nl=incr(nl); write (li(nl), rmt) 'CR;REAL;T;',                              opti_de_param%CR
  nl=incr(nl); write (li(nl), rmt) 'F;REAL;T;',                               opti_de_param%F
  nl=incr(nl); write (li(nl), rmt) 'l_best;REAL;T;',                          opti_de_param%l_best
  nl=incr(nl); write (li(nl), lmt) 'binomial_cross;LOGIC;T;',                 opti_de_param%binomial_cross
  nl=incr(nl); write (li(nl), lmt) 'use_2nd_diff;LOGIC;T;',                   opti_de_param%use_2nd_diff
  nl=incr(nl); write (li(nl), lmt) 'randomize_F;LOGIC;T;',                    opti_de_param%randomize_F
  nl=incr(nl); write (li(nl), lmt) 'minimize_merit;LOGIC;T;',                 opti_de_param%minimize_merit

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% help
!
! Output list of "help xxx" topics
!
! Notes
! -----
! Command syntax:
!   pipe help
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:

case ('help')

  call tao_help ('help-list', '', li, n)

  nl2 = 0
  do i = 1, n
    if (li(i) == '') cycle
    call string_trim(li(i), line, ix)
    nl=incr(nl); name1(nl) = line(1:ix)
    call string_trim(line(ix+1:), line, ix)
    if (ix == 0) cycle
    nl2=nl2+1; name2(nl2) = line
  enddo

  li(1:nl) = name1(1:nl)
  li(nl+1:nl+nl2) = name2(1:nl2)
  nl = nl + nl2

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% inum
!
! Output list of possible values for an INUM parameter.
! For example, possible index numbers for the branches of a lattice.
!
! Notes
! -----
! Command syntax:
!   pipe inum {who}
! 
! Parameters
! ----------
! who
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!    who: ix_universe

case ('inum')

  ix = index(line, '^')
  if (ix /= 0) then
    head = line(:ix-1)
    line = line(ix+1:)
  endif

  select case (line)
  case ('ix_branch')
    u => point_to_uni(head, .false., err);  if (err) return
    do i = 0, size(u%design%lat%branch)
      nl=incr(nl); write (li(nl), '(i0)') i
    enddo

  case ('ix_universe')
    do i = 1, ubound(s%u, 1)
      nl=incr(nl); write (li(nl), '(i0)') i
    enddo

  case ('ix_bunch')
    u => point_to_uni(head, .false., err);  if (err) return
    ix_branch = parse_branch(line, u, .false., err); if (err) return
    do i = 0, u%model_branch(ix_branch)%beam%beam_init%n_bunch
      nl=incr(nl); write (li(nl), '(i0)') i
    enddo

  case ('interpolation_order')
    nl=incr(nl); write (li(nl), '(i0)') 1
    nl=incr(nl); write (li(nl), '(i0)') 3

  case('tick_side', 'number_side')
    nl=incr(nl); write (li(nl), '(i0)') -1
    nl=incr(nl); write (li(nl), '(i0)') +1

  case default
    call invalid ('Not a recognized inum')
  end select

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% lat_calc_done
!
! Output if a lattice recalculation has been proformed since the last 
!   time "pipe lat_calc_done" was called.
!
! Notes
! -----
! Command syntax:
!   pipe lat_calc_done
! 
! Parameters
! ----------
! branch_name
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!    branch_name: 1@0

case ('lat_calc_done')

  nl=incr(nl); write (li(nl), '(l1)') s%com%lattice_calc_done
  s%com%lattice_calc_done = .false.

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% lat_ele_list
!
! Output lattice element list.
!
! Notes
! -----
! Command syntax:
!   pipe lat_ele_list {branch_name}
!
! {branch_name} should have the form:
!   {ix_uni}@{ix_branch}
! 
! Parameters
! ----------
! branch_name : optional
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!    branch_name: 1@0

case ('lat_ele_list')

  u => point_to_uni(line, .true., err); if (err) return
  ix_branch = parse_branch(line, u, .false., err); if (err) return
  branch => u%model%lat%branch(ix_branch)

  do i = 0, branch%n_ele_max
    nl=incr(nl); write (li(nl), '(i0, 2a)') i, ';', branch%ele(i)%name
  enddo

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% lat_header
!
! Output lattice "header" info like the lattice and machine names.
!
! Notes
! -----
! Command syntax:
!   pipe lat_header {ix_uni}
! 
! Output syntax is parameter list form. See documentation at the beginning of this file.
! 
! Parameters
! ----------
! ix_uni : optional
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!    ix_uni: 1

case ('lat_header') 

  u => point_to_uni(line, .false., err); if (err) return
  lat => u%model%lat

  nl=incr(nl); write (li(nl), amt) 'use_name;STR;T;',                   trim(lat%use_name)
  nl=incr(nl); write (li(nl), amt) 'lattice;STR;T;',                    trim(lat%lattice)
  nl=incr(nl); write (li(nl), amt) 'machine;STR;T;',                    trim(lat%machine)
  nl=incr(nl); write (li(nl), amt) 'input_file_name;STR;T;',            trim(lat%input_file_name)
  nl=incr(nl); write (li(nl), amt) 'title;STR;T;',                      trim(lat%title)
  nl=incr(nl); write (li(nl), amt) 'photon_type;ENUM;T;',               trim(photon_type_name(lat%photon_type))

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% lat_branch_list
!
! Output lattice branch list
!
! Notes
! -----
! Command syntax:
!   pipe lat_branch_list {ix_uni}
! 
! Output syntax:
!   branch_index;branch_name;n_ele_track;n_ele_max
! 
! Parameters
! ----------
! ix_uni : optional
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!    ix_uni: 1

case ('lat_branch_list')  ! lat_general is deprecated.

  u => point_to_uni(line, .false., err); if (err) return
  lat => u%model%lat

  do i = 0, ubound(lat%branch, 1)
    branch => lat%branch(i)
    nl=incr(nl); write (li(nl), '(i0, 3a, 2(i0, a))') i, ';', trim(branch%name), ';', branch%n_ele_track, ';', branch%n_ele_max
  enddo

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% lat_list
!
! Output list of parameters at ends of lattice elements
!
! Notes
! -----
! Command syntax:
!   pipe lat_list {flags} {ix_uni}@{ix_branch}>>{elements}|{which} {who}
!
! Where:
!  Optional {flags} are:
!   -no_slaves   - If present, multipass_slave and super_slave elements will not 
!                -   be matched to.
!   -track_only  - If present, lord elements will not be matched to.
!   -index_order - If present, order elements by element index instead of the 
!                -   standard s-position.
!   -array_out   - If present, the output will be available in the 
!     tao_c_interface_com%c_real or tao_c_interface_com%c_integer arrays. 
!     See the code below for when %c_real vs %c_integer is used.
!     Note: Only a single {who} item permitted when -array_out is present.
!
!   {which} is one of: "model", "base" or "design"
! 
!   {who} is a comma deliminated list of:
!     orbit.floor.x, orbit.floor.y, orbit.floor.z    ! Floor coords at particle orbit.
!     orbit.spin.1, orbit.spin.2, orbit.spin.3,
!     orbit.vec.1, orbit.vec.2, orbit.vec.3, orbit.vec.4, orbit.vec.5, orbit.vec.6,
!     orbit.t, orbit.beta,
!     orbit.state,     ! Note: state is an integer. alive$ = 1, anything else is lost.
!     orbit.energy, orbit.pc,
!     ele.name, ele.key, ele.ix_ele, ele.ix_branch
!     ele.a.beta, ele.a.alpha, ele.a.eta, ele.a.etap, ele.a.gamma, ele.a.phi,
!     ele.b.beta, ele.b.alpha, ele.b.eta, ele.b.etap, ele.b.gamma, ele.b.phi,
!     ele.x.eta, ele.x.etap,
!     ele.y.eta, ele.y.etap,
!     ele.ref_time, ele.ref_time_start
!     ele.s, ele.l
!     ele.e_tot, ele.p0c
!     ele.mat6      ! Output: mat6(1,:), mat6(2,:), ... mat6(6,:)
!     ele.vec0      ! Output: vec0(1), ... vec0(6)
!     ele.c_mat     ! Output: c_mat11, c_mat12, c_mat21, c_mat22.
!     ele.gamma_c   ! Parameter associated with coupling c-matrix.
!     ele.XXX       ! Where XXX is a Bmad syntax element attribute. 
!                   !   EG: ele.beta_a, ele.k1, etc.
! 
!   {elements} is a string to match element names to.
!     Use "*" to match to all elements.
!
! Examples:
!   pipe lat_list -track 3@0>>Q*|base ele.s,orbit.vec.2
!   pipe lat_list 3@0>>Q*|base real:ele.s    
! 
! Also see: "pipe ele:param"
!
! Parameters
! ----------
! elements 
! who 
! ix_uni : optional
! ix_branch : optional
! which : default=model
! flags : optional, default=-array_out -track_only
!
! Returns
! -------
! string_list
!   if ('-array_out' not in flags) or (who in ['ele.name', 'ele.key'])
! integer_array
!    if '-array_out' in flags and who in ['orbit.state', 'ele.ix_ele']
! real_array
!    if ('-array_out' in flags) or ('real:' in who) 
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!    ix_uni: 1  
!    ix_branch: 0 
!    elements: Q* 
!    which: model
!    who: orbit.floor.x
!
! Example: 2
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!    ix_uni: 1  
!    ix_branch: 0 
!    elements: Q* 
!    which: design
!    who: ele.ix_ele

case ('lat_list')

  no_slaves = .false.
  track_only = .false.
  index_order = .false.
  use_real_array_buffer = .false.

  do
    if (ix_line == 0) then
      call invalid ('List of elements not present.')
      return
    endif

    if (index('-no_slaves', line(1:ix_line)) == 1) then
      call string_trim(line(ix_line+1:), line, ix_line)
      no_slaves = .true.
      cycle
    endif

    if (index('-track_only', line(1:ix_line)) == 1) then
      call string_trim(line(ix_line+1:), line, ix_line)
      track_only = .true.
      cycle
    endif

    if (index('-index_order', line(1:ix_line)) == 1) then
      call string_trim(line(ix_line+1:), line, ix_line)
      index_order = .true.
      cycle
    endif

    if (index('-array_out', line(1:ix_line)) == 1) then
      call string_trim(line(ix_line+1:), line, ix_line)
      use_real_array_buffer = .true.
      call re_allocate(real_arr, 1000)
      cycle
    endif

    exit
  enddo

  u => point_to_uni(line, .true., err); if (err) return
  tao_lat => point_to_tao_lat(line, u, err, tail_str = all_who); if (err) return

  if (all_who(1:5) == 'real:') then  ! Old style
    use_real_array_buffer = .true.
    call re_allocate(real_arr, 1000)
    all_who = all_who(6:)
  endif

  call upcase_string(line)
  lat => tao_lat%lat
  call lat_ele_locator (line, lat, eles, n_loc, err, order_by_index = index_order);  if (err) return

  n_who = 0
  do
    ix = index(all_who, ',')
    n_who = n_who + 1
    if (ix == 0) then
      name1(n_who) = all_who
      exit
    else
      name1(n_who) = all_who(1:ix-1)
      all_who = adjustl(all_who(ix+1:))
    endif
  enddo

  if (use_real_array_buffer .and. n_who /= 1) then
    call invalid ('Number of "who" must be 1 for real buffered output.')
  endif

  n_arr = 0

  do ie = 1, n_loc
    ele => eles(ie)%ele
    if (track_only .and. ele%ix_ele > lat%n_ele_track) cycle
    if (no_slaves .and. (ele%slave_status == super_slave$ .or. ele%slave_status == multipass_slave$)) cycle
    orbit => tao_lat%tao_branch(ele%ix_branch)%orbit(ele%ix_ele)

    do i = 1, n_who

      n_add = 1

      select case (name1(i))
      case ('ele.mat6')
        n_add = 36
        do ix = 1, 6
          values(6*(ix-1)+1:6*ix) = ele%mat6(ix,:)
        enddo
      case ('ele.vec0')
        n_add = 6
        values(1:6) = ele%vec0
      case ('ele.c_mat')
        n_add = 4
        values(1:4) = [ele%c_mat(1,1), ele%c_mat(1,2), ele%c_mat(2,1), ele%c_mat(2,2)]
      case ('ele.name')
        nl=incr(nl); li(nl) = ele%name
        cycle
      case ('ele.key')
        nl=incr(nl); li(nl) = key_name(ele%key)
        cycle
      case default
        values(1) = ele_param_value(name1(i), ele, orbit, data_type, err)
        if (err) return
      end select

      !

      if (use_real_array_buffer) then
        n_arr = n_arr
        if (n_arr+n_add > size(real_arr)) call re_allocate(real_arr, 2*(n_arr+n_add))
        real_arr(n_arr+1:n_arr+n_add) = values(1:n_add)
        n_arr = n_arr + n_add

      else
        do ix = 1, n_add
          if (data_type == is_integer$) then
            if (i == 1 .and. ix == 1) then
              nl=incr(nl); write (li(nl), '(i0)') nint(values(ix))
            else
              write (li(nl), '(2a, i0)') trim(li(nl)), ';', nint(values(ix))
            endif
          else
            if (i == 1 .and. ix == 1) then
              nl=incr(nl); write (li(nl), '(es22.14)') values(ix)
            else
              write (li(nl), '(2a, es22.14)') trim(li(nl)), ';', values(ix)
            endif
          endif
        enddo
      endif

    enddo ! i = 1, n_who

  enddo ! ie = 1, n_loc

  if (use_real_array_buffer) then
    if (data_type == is_integer$) then
      if (.not. allocated(tao_c_interface_com%c_integer)) allocate (tao_c_interface_com%c_integer(n_arr))
      if (size(tao_c_interface_com%c_integer) < n_arr) then
        deallocate (tao_c_interface_com%c_integer)
        allocate (tao_c_interface_com%c_integer(n_arr))
      endif

      tao_c_interface_com%n_int = n_arr
      tao_c_interface_com%c_integer(1:n_arr) = nint(real_arr(1:n_arr))

    else
      call re_allocate_c_double(tao_c_interface_com%c_real, n_arr, .false.)
      tao_c_interface_com%n_real = n_arr
      tao_c_interface_com%c_real(1:n_arr) = real_arr(1:n_arr)
    endif
  endif

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% lat_param_units
!
! Output units of a parameter associated with a lattice or lattice element.
!
! Notes
! -----
! Command syntax:
!   pipe lat_param_units {param_name}
! 
! Parameters
! ----------
! param_name
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!    param_name: L   

case ('lat_param_units')

  name = upcase(line)
  a_name = attribute_units(name)
  nl=incr(nl); write(li(nl), '(a)') a_name

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% matrix
!
! Output matrix value from the exit end of one element to the exit end of the other.
!
! Notes
! -----
! Command syntax:
!   pipe matrix {ele1_id} {ele2_id}
!
! Where:
!   {ele1_id} is the start element.
!   {ele2_id} is the end element.
! If {ele2_id} = {ele1_id}, the 1-turn transfer map is computed.
! Note: {ele2_id} should just be an element name or index without universe, 
!       branch, or model/base/design specification.
!
! Example:
!   pipe matrix 2@1>>q01w|design q02w
! 
! Parameters
! ----------
! ele1_id
! ele2_id
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!    ele1_id: 1@0>>q01w|design
!    ele2_id: q02w


case ('matrix')

  u => point_to_uni(line, .true., err); if (err) return
  tao_lat => point_to_tao_lat(line, u, err, which, tail_str); if (err) return
  ele => point_to_ele(line, tao_lat%lat, err); if (err) return

  call lat_ele_locator (tail_str, tao_lat%lat, eles, n_loc, err, ix_dflt_branch = ele%ix_branch)
  if (err .or. n_loc == 0) then
    call invalid ('Bad ele2_id: ' // line)
    return
  endif
  if (n_loc > 1) then
    call invalid ('More than one element matches name: ' // line)
    return
  endif

  call transfer_matrix_calc (tao_lat%lat, mat6, vec0, ele%ix_ele, eles(1)%ele%ix_ele, ele%ix_branch, one_turn = .true.)
  do i = 1, 6
    array(1:7) = [mat6(i,:), vec0(i)]
    nl=incr(nl); write (li(nl), ramt) int_str(i), (';', array(j), j = 1, 7)
  enddo

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% merit
!
! Output merit value.
!
! Notes
! -----
! Command syntax:
!   pipe merit
! 
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:

case ('merit')

  nl=incr(nl); write (li(nl), '(es22.14)') tao_merit()

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% orbit_at_s
!
! Output twiss at given s position.
!
! Notes
! -----
! Command syntax:
!   pipe orbit_at_s {ix_uni}@{ele}->{s_offset}|{which}
!
! Where:
!   {ix_uni}   - Universe index. Defaults to s%global%default_universe.
!   {ele}      - Element name or index. 
!                  Default at the Beginning element at start of branch 0.
!   {s_offset} - Offset of the evaluation point from the downstream end of ele. 
!                  Default is 0. If {s_offset} is present, the preceeding "->" sign
!                  must be present. EG: Something like "23|model" will {which} is 
!                  one of: "model", "base" or "design".
!
! Example:
!   pipe orbit_at_s Q10->0.4|model   ! Orbit at 0.4 meters from Q10 element exit end in model lattice.
! 
! Parameters
! ----------
! ix_uni : optional
! ele : optional
! s_offset : optional
! which : default=model
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!    ix_uni: 1
!    ele: 10
!    s_offset: 0.7
!    which: model

case ('orbit_at_s')

  u => point_to_uni(line, .true., err); if (err) return
  tao_lat => point_to_tao_lat(line, u, err); if (err) return
  s_pos = parse_ele_with_s_offset(line, tao_lat, ele, err); if (err) return
  ix_branch = ele%ix_branch

  call twiss_and_track_at_s (tao_lat%lat, s_pos, orb = tao_lat%tao_branch(ix_branch)%orbit, orb_at_s = orb, ix_branch = ix_branch)
  call orbit_out (orb)

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% place_buffer
!
! Output the place command buffer and reset the buffer.
! The contents of the buffer are the place commands that the user has issued.
! See the Tao manual for more details.
!
! Notes
! -----
! Command syntax:
!   pipe place_buffer
! 
! Returns
! -------
! None
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:

case ('place_buffer')

  if (.not. allocated(s%com%plot_place_buffer)) return

  do i = 1, size(s%com%plot_place_buffer)
    nl=incr(nl); li(nl) = trim(s%com%plot_place_buffer(i)%name) // ';' // s%com%plot_place_buffer(i)%plot%name
  enddo

  deallocate(s%com%plot_place_buffer)

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% plot_curve
!
! Output curve information for a plot.
!
! Notes
! -----
! Command syntax:
!   pipe plot_curve {curve_name}
! 
! Parameters
! ----------
! curve_name
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_optics_matching
!  args:
!    curve_name: r13.g.a

case ('plot_curve')

  call tao_find_plots (err, line, 'BOTH', curve = curves, only_visible = .false.)

  if (err .or. size(curves) /= 1) then
    call invalid ('Not a valid curve')
    return
  endif

  c => curves(1)%c
  ix_uni = c%ix_universe

  nl=incr(nl); write (li(nl), amt) 'name;STR;T;',                             trim(c%name)
  nl=incr(nl); write (li(nl), amt) 'data_source;ENUM;T;',                     trim(c%data_source)
  nl=incr(nl); write (li(nl), amt) 'data_type_x;DAT_TYPE_Z;T;',               trim(c%data_type_x)
  nl=incr(nl); write (li(nl), amt) 'data_type;DAT_TYPE;T;',                   trim(c%data_type)
  nl=incr(nl); write (li(nl), amt) 'component;COMPONENT;T;',                  trim(c%component)
  nl=incr(nl); write (li(nl), amt) 'ele_ref_name;STR;T;',                     trim(c%ele_ref_name)
  nl=incr(nl); write (li(nl), amt) 'legend_text;STR;T;',                      trim(c%legend_text)
  nl=incr(nl); write (li(nl), amt) 'message_text;STR;F;',                     trim(c%message_text)
  nl=incr(nl); write (li(nl), amt) 'why_invalid;STR;I;',                      trim(c%why_invalid)
  nl=incr(nl); write (li(nl), rmt) 'y_axis_scale_factor;REAL;T;',             c%y_axis_scale_factor
  nl=incr(nl); write (li(nl), imt) 'ix_universe;INUM;T;',                     c%ix_universe
  nl=incr(nl); write (li(nl), imt) 'symbol_every;INT;T;',                     c%symbol_every
  nl=incr(nl); write (li(nl), jmt) ix_uni, '^ix_branch;INUM;T;',              c%ix_branch
  nl=incr(nl); write (li(nl), jmt) ix_uni, '^ix_bunch;INUM;T;',               c%ix_bunch
  nl=incr(nl); write (li(nl), lmt) 'use_y2;LOGIC;T;',                         c%use_y2
  nl=incr(nl); write (li(nl), lmt) 'draw_line;LOGIC;T;',                      c%draw_line
  nl=incr(nl); write (li(nl), lmt) 'draw_symbols;LOGIC;T;',                   c%draw_symbols
  nl=incr(nl); write (li(nl), lmt) 'draw_symbol_index;LOGIC;T;',              c%draw_symbol_index
  nl=incr(nl); write (li(nl), lmt) 'draw_error_bars;LOGIC;T;',                c%draw_error_bars
  nl=incr(nl); write (li(nl), lmt) 'smooth_line_calc;LOGIC;T;',               c%smooth_line_calc
  nl=incr(nl); write (li(nl), lmt) 'z_color.is_on;LOGIC;I;',                  c%z_color%is_on
  nl=incr(nl); write (li(nl), rmt) 'z_color.min;REAL;T;',                     c%z_color%min
  nl=incr(nl); write (li(nl), rmt) 'z_color.max;REAL;T;',                     c%z_color%max
  nl=incr(nl); write (li(nl), lmt) 'z_color.autoscale;LOGIC;I;',              c%z_color%autoscale
  nl=incr(nl); write (li(nl), amt) 'z_color.data_type;ENUM;T;',               trim(c%z_color%data_type)
  nl=incr(nl); write (li(nl), lmt) 'valid;LOGIC;I;',                          c%valid
  nl=incr(nl); write (li(nl), '(a, i0, 4a)') 'line;STRUCT;T;width;INT;', c%line%width, &
                      ';color;ENUM;', trim(c%line%color), ';line^pattern;ENUM;', c%line%pattern

  nl=incr(nl); write (li(nl), '(9a, i0)')  'symbol;STRUCT;T;symbol^type;ENUM;', trim(c%symbol%type), &
                      ';color;ENUM;', trim(c%symbol%color), ';height;REAL;', to_str(c%symbol%height, 4), &
                      ';fill_pattern;ENUM;', trim(c%symbol%fill_pattern), ';line_width;INT;', c%symbol%line_width

  nl=incr(nl); write (li(nl), imt)  'symbol.line_width;INT;T;',               c%symbol%line_width

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% plot_graph
!
! Output graph info.
!
! Notes
! -----
! Command syntax:
!   pipe plot_graph {graph_name}
!
! {graph_name} is in the form:
!   {p_name}.{g_name}
! where
!   {p_name} is the plot region name if from a region or the plot name if a template plot.
!   This name is obtained from the pipe plot_list command.
!   {g_name} is the graph name obtained from the pipe plot1 command.
! 
! Parameters
! ----------
! graph_name
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_optics_matching
!  args:
!    graph_name: beta.g

case ('plot_graph')

  call tao_find_plots (err, line, 'BOTH', graph = graphs, only_visible = .false.)

  if (err .or. size(graphs) /= 1) then
    call invalid ('Bad graph name')
    return
  endif

  g => graphs(1)%g

  if (g%type == 'floor_plan') then
    call tao_set_floor_plan_axis_label (g, g%x, x_ax, 'X')
    call tao_set_floor_plan_axis_label (g, g%y, y_ax, 'Y')
  else
    x_ax = g%x
    y_ax = g%y
  endif

  n = 0
  if (allocated(g%curve)) n = size(g%curve)

  nl=incr(nl); write (li(nl), imt) 'num_curves;INT;T;',                       n
  do i = 1, n
    nl=incr(nl); write (li(nl), vamt) 'curve[', i, '];STR;T;',                g%curve(i)%name
  enddo

  nl=incr(nl); write (li(nl), amt) 'name;STR;T;',                               trim(g%name)
  nl=incr(nl); write (li(nl), amt) 'graph^type;ENUM;T;',                        trim(g%type)
  nl=incr(nl); write (li(nl), amt) 'title;STR;T;',                              trim(g%title)
  nl=incr(nl); write (li(nl), amt) 'title_suffix;STR;F;',                       trim(g%title_suffix)
  nl=incr(nl); write (li(nl), amt) 'why_invalid;STR;F;',                        trim(g%why_invalid)
  nl=incr(nl); write (li(nl), rmt) 'x_axis_scale_factor;REAL;T;',               g%x_axis_scale_factor
  nl=incr(nl); write (li(nl), rmt) 'symbol_size_scale;REAL;T;',                 g%symbol_size_scale
  nl=incr(nl); write (li(nl), jmt) g%ix_universe, '^ix_branch;INUM;T;',         g%ix_branch
  nl=incr(nl); write (li(nl), imt) 'ix_universe;INUM;T;',                       g%ix_universe
  nl=incr(nl); write (li(nl), lmt) 'clip;LOGIC;T;',                             g%clip
  nl=incr(nl); write (li(nl), lmt) 'is_valid;LOGIC;F;',                         g%is_valid
  nl=incr(nl); write (li(nl), lmt) 'y2_mirrors_y;LOGIC;T;',                     g%y2_mirrors_y
  nl=incr(nl); write (li(nl), lmt) 'limited;LOGIC;F;',                          g%limited
  nl=incr(nl); write (li(nl), lmt) 'draw_axes;LOGIC;T;',                        g%draw_axes
  nl=incr(nl); write (li(nl), lmt) 'draw_curve_legend;LOGIC;T;',                g%draw_curve_legend
  nl=incr(nl); write (li(nl), lmt) 'draw_grid;LOGIC;T;',                        g%draw_grid
  nl=incr(nl); write (li(nl), lmt) 'draw_only_good_user_data_or_vars;LOGIC;T;', g%draw_only_good_user_data_or_vars

  fp => g%floor_plan
  nl=incr(nl); write (li(nl), amt) 'floor_plan.view;ENUM;T;',                   fp%view
  nl=incr(nl); write (li(nl), amt) 'floor_plan.rotation;REAL;T;',               to_str(fp%rotation, 6)
  nl=incr(nl); write (li(nl), amt) 'floor_plan.flip_label_side;LOGIC;T;',       logic_str(fp%flip_label_side)
  nl=incr(nl); write (li(nl), amt) 'floor_plan.size_is_absolute;LOGIC;T;',      logic_str(fp%size_is_absolute)
  nl=incr(nl); write (li(nl), amt) 'floor_plan.draw_building_wall;LOGIC;T;',    logic_str(fp%draw_building_wall)
  nl=incr(nl); write (li(nl), amt) 'floor_plan.draw_only_first_pass;LOGIC;T;',  logic_str(fp%draw_only_first_pass)
  nl=incr(nl); write (li(nl), amt) 'floor_plan.correct_distortion;LOGIC;T;',    logic_str(fp%correct_distortion)
  nl=incr(nl); write (li(nl), amt) 'floor_plan.orbit_scale;REAL;T;',            to_str(fp%orbit_scale, 4)
  nl=incr(nl); write (li(nl), amt) 'floor_plan.orbit_color;ENUM;T;',            trim(fp%orbit_color)
  nl=incr(nl); write (li(nl), amt) 'floor_plan.orbit_lattice;ENUM;T;',          trim(fp%orbit_lattice)
  nl=incr(nl); write (li(nl), amt) 'floor_plan.orbit_width;INT;T;',             int_str(fp%orbit_width)
  nl=incr(nl); write (li(nl), amt) 'floor_plan.orbit_pattern;ENUM;T;',          trim(fp%orbit_pattern)

  nl=incr(nl); write (li(nl), amt) 'x.label;STR;T;',                            trim(x_ax%label)
  nl=incr(nl); write (li(nl), amt) 'x.label_color;ENUM;T;',                     trim(x_ax%label_color)
  nl=incr(nl); write (li(nl), amt) 'x.label_offset;REAL;T;',                    to_str(x_ax%label_offset,6)
  nl=incr(nl); write (li(nl), amt) 'x.max;REAL;T;',                             to_str(x_ax%max,6)
  nl=incr(nl); write (li(nl), amt) 'x.min;REAL;T;',                             to_str(x_ax%min,6)
  nl=incr(nl); write (li(nl), amt) 'x.axis^type;ENUM;T;',                       trim(x_ax%type)
  nl=incr(nl); write (li(nl), amt) 'x.bounds;ENUM;T;',                          trim(x_ax%bounds)
  nl=incr(nl); write (li(nl), amt) 'x.number_offset;REAL;T;',                   to_str(x_ax%number_offset,6)
  nl=incr(nl); write (li(nl), imt) 'x.major_div_nominal;INT;T;',                x_ax%major_div_nominal
  nl=incr(nl); write (li(nl), imt) 'x.minor_div;INT;T;',                        x_ax%minor_div
  nl=incr(nl); write (li(nl), imt) 'x.minor_div_max;INT;T;',                    x_ax%minor_div_max
  nl=incr(nl); write (li(nl), lmt) 'x.draw_label;LOGIC;T;',                     x_ax%draw_label
  nl=incr(nl); write (li(nl), lmt) 'x.draw_numbers;LOGIC;T;',                   x_ax%draw_numbers
  nl=incr(nl); write (li(nl), imt) 'x.tick_side;INUM;T;',                       x_ax%tick_side
  nl=incr(nl); write (li(nl), imt) 'x.number_side;INUM;T;',                     x_ax%number_side
  nl=incr(nl); write (li(nl), amt) 'x.major_tick_len;REAL;T;',                  to_str(x_ax%major_tick_len,6)
  nl=incr(nl); write (li(nl), amt) 'x.minor_tick_len;REAL;T;',                  to_str(x_ax%minor_tick_len,6)

  nl=incr(nl); write (li(nl), amt) 'y.label;STR;T;',                            trim(y_ax%label)
  nl=incr(nl); write (li(nl), amt) 'y.label_color;ENUM;T;',                     trim(y_ax%label_color)
  nl=incr(nl); write (li(nl), amt) 'y.label_offset;REAL;T;',                    to_str(y_ax%label_offset,6)
  nl=incr(nl); write (li(nl), amt) 'y.max;REAL;T;',                             to_str(y_ax%max,6)
  nl=incr(nl); write (li(nl), amt) 'y.min;REAL;T;',                             to_str(y_ax%min,6)
  nl=incr(nl); write (li(nl), amt) 'y.axis^type;ENUM;T;',                       trim(y_ax%type)
  nl=incr(nl); write (li(nl), amt) 'y.bounds;ENUM;T;',                          trim(y_ax%bounds)
  nl=incr(nl); write (li(nl), amt) 'y.number_offset;REAL;T;',                   to_str(y_ax%number_offset,6)
  nl=incr(nl); write (li(nl), imt) 'y.major_div_nominal;INT;T;',                y_ax%major_div_nominal
  nl=incr(nl); write (li(nl), imt) 'y.minor_div;INT;T;',                        y_ax%minor_div
  nl=incr(nl); write (li(nl), imt) 'y.minor_div_max;INT;T;',                    y_ax%minor_div_max
  nl=incr(nl); write (li(nl), lmt) 'y.draw_label;LOGIC;T;',                     y_ax%draw_label
  nl=incr(nl); write (li(nl), lmt) 'y.draw_numbers;LOGIC;T;',                   y_ax%draw_numbers
  nl=incr(nl); write (li(nl), imt) 'y.tick_side;INUM;T;',                       y_ax%tick_side
  nl=incr(nl); write (li(nl), imt) 'y.number_side;INUM;T;',                     y_ax%number_side
  nl=incr(nl); write (li(nl), amt) 'y.major_tick_len;REAL;T;',                  to_str(y_ax%major_tick_len,6)
  nl=incr(nl); write (li(nl), amt) 'y.minor_tick_len;REAL;T;',                  to_str(y_ax%minor_tick_len,6)

  nl=incr(nl); write (li(nl), amt) 'y2.label;STR;T;',                           trim(g%y2%label)
  nl=incr(nl); write (li(nl), amt) 'y2.label_color;ENUM;T;',                    trim(g%y2%label_color)
  nl=incr(nl); write (li(nl), amt) 'y2.label_offset;REAL;T;',                   to_str(g%y2%label_offset,6)
  nl=incr(nl); write (li(nl), amt) 'y2.max;REAL;T;',                            to_str(g%y2%max,6)
  nl=incr(nl); write (li(nl), amt) 'y2.min;REAL;T;',                            to_str(g%y2%min,6)
  nl=incr(nl); write (li(nl), amt) 'y2.axis^type;ENUM;T;',                      trim(g%y2%type)
  nl=incr(nl); write (li(nl), amt) 'y2.bounds;ENUM;T;',                         trim(g%y2%bounds)
  nl=incr(nl); write (li(nl), amt) 'y2.number_offset;REAL;T;',                  to_str(g%y2%number_offset,6)
  nl=incr(nl); write (li(nl), imt) 'y2.major_div_nominal;INT;T;',               g%y2%major_div_nominal
  nl=incr(nl); write (li(nl), imt) 'y2.minor_div;INT;T;',                       g%y2%minor_div
  nl=incr(nl); write (li(nl), imt) 'y2.minor_div_max;INT;T;',                   g%y2%minor_div_max
  nl=incr(nl); write (li(nl), lmt) 'y2.draw_label;LOGIC;T;',                    g%y2%draw_label
  nl=incr(nl); write (li(nl), lmt) 'y2.draw_numbers;LOGIC;T;',                  g%y2%draw_numbers
  nl=incr(nl); write (li(nl), imt) 'y2.tick_side;INUM;T;',                      g%y2%tick_side
  nl=incr(nl); write (li(nl), imt) 'y2.number_side;INUM;T;',                    g%y2%number_side
  nl=incr(nl); write (li(nl), amt) 'y2.major_tick_len;REAL;T;',                 to_str(g%y2%major_tick_len,6)
  nl=incr(nl); write (li(nl), amt) 'y2.minor_tick_len;REAL;T;',                 to_str(g%y2%minor_tick_len,6)

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% plot_histogram
!
! Output plot histogram info.
!
! Notes
! -----
! Command syntax:
!   pipe plot_histogram {curve_name}
! 
! Parameters
! ----------
! curve_name
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_optics_matching
!  args:
!    curve_name: r33.g.x

case ('plot_histogram')

  call tao_find_plots (err, line, 'BOTH', curve = curves, only_visible = .false.)

  if (err .or. size(curves) /= 1) then
    call invalid ('Bad curve name')
    return
  endif

  c => curves(1)%c

  nl=incr(nl); write (li(nl), lmt) 'density_normalized;LOGIC;T;',          c%hist%density_normalized
  nl=incr(nl); write (li(nl), lmt) 'weight_by_charge;LOGIC;T;',            c%hist%weight_by_charge
  nl=incr(nl); write (li(nl), rmt) 'minimum;REAL;T;',                      c%hist%minimum
  nl=incr(nl); write (li(nl), rmt) 'maximum;REAL;T;',                      c%hist%maximum
  nl=incr(nl); write (li(nl), rmt) 'width;REAL;T;',                        c%hist%width
  nl=incr(nl); write (li(nl), rmt) 'center;REAL;T;',                       c%hist%center
  nl=incr(nl); write (li(nl), imt) 'number;INT;T;',                        c%hist%number

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% plot_lat_layout
!
! Output plot Lat_layout info
!
! Notes
! -----
! Command syntax:
!   pipe plot_lat_layout {ix_uni}@{ix_branch}
!
! Note: The returned list of element positions is not ordered in increasing
!       longitudinal position.
! 
! Parameters
! ----------
! ix_uni: 1
! ix_branch: 0
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!    ix_uni: 1
!    ix_branch: 0 

case ('plot_lat_layout')

  u => point_to_uni(line, .true., err); if (err) return
  ix_branch = parse_branch(line, u, .false., err); if (err) return
  lat => u%model%lat
  branch => lat%branch(ix_branch)

  do i = 1, branch%n_ele_track
    ele => branch%ele(i)
    if (ele%slave_status == super_slave$) cycle

    ix_shape_min = 1
    do
      call tao_ele_shape_info (u%ix_uni, ele, s%plot_page%lat_layout%ele_shape, shape, label_name, y1, y2, ix_shape_min)
      if (.not. associated(shape)) exit
      if (.not. shape%draw) cycle
      y1 = y1 * s%plot_page%lat_layout_shape_scale
      y2 = y2 * s%plot_page%lat_layout_shape_scale
      nl=incr(nl); write (li(nl), '(2(i0, a), 2(es22.14, a), (i0, a), 2a, 2(es10.2, a), 4a)') ele%ix_branch, ';', ele%ix_ele, &
                ';', ele%s_start, ';', ele%s, ';', shape%line_width, ';', trim(shape%shape), ';', &
                y1, ';', y2, ';', trim(shape%color), ';', trim(label_name)
    enddo
  enddo

  do i = lat%n_ele_track+1, lat%n_ele_max
    ele => lat%ele(i)
    if (ele%lord_status == multipass_lord$) cycle
    branch2 => pointer_to_branch(ele)
    if (branch2%ix_branch /= branch%ix_branch) cycle

    ix_shape_min = 1
    do
      call tao_ele_shape_info (u%ix_uni, ele, s%plot_page%lat_layout%ele_shape, shape, label_name, y1, y2, ix_shape_min)
      if (.not. associated(shape)) exit
      if (.not. shape%draw) cycle
      y1 = y1 * s%plot_page%lat_layout_shape_scale
      y2 = y2 * s%plot_page%lat_layout_shape_scale
      nl=incr(nl); write (li(nl), '(2(i0, a), 2(es22.14, a), (i0, a), 2a, 2(es10.2, a), 4a)') ele%ix_branch, ';', ele%ix_ele, &
                ';', ele%s_start, ';', ele%s, ';', shape%line_width, ';', trim(shape%shape), ';', &
                y1, ';', y2, ';', trim(shape%color), ';', trim(label_name)
    enddo
  enddo

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% plot_list
!
! Output list of plot templates or plot regions.
!
! Notes
! -----
! Command syntax:
!   pipe plot_list {r_or_g}
!
! where "{r/g}" is:
!   "r"      ! list regions of the form ix;region_name;plot_name;visible;x1;x2;y1;y2
!   "t"      ! list template plots of the form ix;name
! 
! Parameters
! ----------
! r_or_g
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!    r_or_g: r


case ('plot_list')
  if (line == 't') then
    do i = 1, size(s%plot_page%template)
      p => s%plot_page%template(i)
      if (p%phantom) cycle
      if (p%name == '') cycle
      if (p%name == 'scratch') cycle
      nl=incr(nl); write (li(nl), '(i0, 2a)') i, ';', trim(p%name)
    enddo

  elseif (line == 'r') then
    do i = 1, size(s%plot_page%region)
      pr => s%plot_page%region(i)
      if (pr%name == '') cycle
      nl=incr(nl); write (li(nl), '(i0, 5a, l1, 8a)') i, ';', trim(pr%name), ';', trim(pr%plot%name), ';', pr%visible, ';', &
                      re_str(pr%location(1), 4), ';', re_str(pr%location(2), 4), ';', re_str(pr%location(3), 4), ';', re_str(pr%location(4), 4)
    enddo

  else
    call invalid ('Expect "r" or "t"')
  endif

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% plot_template_manage
!
! Template plot creation or destruction.
!
! Notes
! -----
! Command syntax:
!   pipe plot_template_manage {template_location}^^{template_name}^^
!                          {n_graph}^^{graph_names}
!
! Where:
!   {template_location} - Location to place or delete a template plot. 
!                           Use "@Tnnn" syntax for the location.
!   {template_name}     - The name of the template plot. 
!                           If deleting a plot this name is immaterial.
!   {n_graph}           - The number of associated graphs. 
!                           If set to -1 then any existing template plot is deleted.
!   {graph_names}       - Names of the graphs. graph_names should be in the form:
!                             graph1_name^^graph2_name^^...^^graphN_name
!                           where N=n_graph names
!
! Parameters
! ----------
! template_location
! template_name
! n_graph : default=-1
! graph_names : default=
!
! Returns
! -------
! None
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_optics_matching
!  args:
!    template_location: @T1
!    template_name: beta
!    n_graph: 2
!    graph_names: g1^^g2


case ('plot_template_manage')

  call split_this_line (line, name1, -1, err);         if (err) return

  call tao_find_plots (err, name1(1), 'TEMPLATE', plots, only_visible = .false.)
  if (size(plots) == 0) then
    call invalid('No plot template location found for: ' // name1(1))
    return
  endif
  p => plots(1)%p

  if (.not. is_integer(name1(3))) then
    call invalid ('Number of graphs string not an integer: ' // name1(3))
    return
  endif

  read(name1(3), *) n

  !

  if (n == -1) then
    ix = p%ix_plot
    n1 = size(s%plot_page%template)
    s%plot_page%template(ix:n1-1) = s%plot_page%template(ix+1:n1)
    do i = ix, n1
      s%plot_page%template(i)%ix_plot = i
    enddo
    return
  endif

  !

  if (allocated(p%graph)) deallocate(p%graph)
  allocate(p%graph(n))

  p%name = name1(2)

  do i = 1, n
    p%graph(i)%name = name1(i+3)
    p%graph(i)%p => p
  enddo

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% plot_curve_manage
!
! Template plot curve creation/destruction
!
! Notes
! -----
! Command syntax:
!   pipe plot_curve_manage {graph_name}^^{curve_index}^^{curve_name}
!
! If {curve_index} corresponds to an existing curve then this curve is deleted.
! In this case the {curve_name} is ignored and does not have to be present.
! If {curve_index} does not not correspond to an existing curve, {curve_index}
! must be one greater than the number of curves.
! 
! Parameters
! ----------
! graph_name
! curve_index
! curve_name
!
! Returns
! -------
! None
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_optics_matching
!  args:
!    graph_name: beta.g
!    curve_index: 1
!    curve_name: r13.g.a

case ('plot_curve_manage')

  call split_this_line (line, name1, -1, err);         if (err) return
  call tao_find_plots (err, name1(1), 'TEMPLATE', graph = graphs, only_visible = .false.)
  if (size(graphs) /= 1) then
    if (size(graphs) == 0) call invalid('No graph found for: ' // name1(1))
    if (size(graphs) > 1)  call invalid('Multiple graphs found for: ' // name1(1))
    return
  endif

  g => graphs(1)%g

  if (allocated(g%curve)) then
    n1 = size(g%curve)
    call move_alloc(g%curve, curve_temp)
  else
    n1 = 0
  endif

  if (.not. is_integer(name1(2))) then
    call invalid ('Curve index not an integer: ' // name1(2))
    return
  endif
  read(name1(2), *) n
  if (n > n1 + 1) then
    call invalid ('Curve index out of range: ' // name1(2))
    return
  endif

  if (n == n1 + 1) then
    allocate (g%curve(n))
    if (n1 /= 0) g%curve(1:n1) = curve_temp
    g%curve(n)%name = name1(3)
    g%curve(n)%g => g

  else  ! Remove curve
    allocate (g%curve(n1-1))
    g%curve(1:n-1) = curve_temp(1:n-1)
    g%curve(n:n1-1) = curve_temp(n+1:n1)
  endif

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% plot_graph_manage
!
! Template plot graph creation/destruction
!
! Notes
! -----
! Command syntax:
!   pipe plot_graph_manage {plot_name}^^{graph_index}^^{graph_name}
!
! If {graph_index} corresponds to an existing graph then this graph is deleted.
! In this case the {graph_name} is ignored and does not have to be present.
! If {graph_index} does not not correspond to an existing graph, {graph_index}
! must be one greater than the number of graphs.
! 
! Parameters
! ----------
! plot_name
! graph_index
! graph_name
!
! Returns
! -------
! None
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_optics_matching
!  args:
!    plot_name: beta
!    graph_index: 1
!    graph_name: beta.g

case ('plot_graph_manage')

  call split_this_line (line, name1, -1, err);         if (err) return
  call tao_find_plots (err, name1(1), 'TEMPLATE', plots, only_visible = .false.)
  if (size(plots) /= 1) then
    if (size(plots) == 0) call invalid('No plot found for: ' // name1(1))
    if (size(plots) > 1)  call invalid('Multiple plots found for: ' // name1(1))
    return
  endif

  p => plots(1)%p

  if (allocated(p%graph)) then
    n1 = size(p%graph)
    call move_alloc(p%graph, graph_temp)
  else
    n1 = 0
  endif

  if (.not. is_integer(name1(2))) then
    call invalid ('Graph index not an integer: ' // name1(2))
    return
  endif
  read(name1(2), *) n
  if (n > n1 + 1) then
    call invalid ('Graph index out of range: ' // name1(2))
    return
  endif

  if (n == n1 + 1) then
    allocate (p%graph(n))
    if (n1 /= 0) p%graph(1:n1) = graph_temp
    p%graph(n)%name = name1(3)
    p%graph(n)%p => p

  else  ! Remove graph
    allocate (p%graph(n1-1))
    p%graph(1:n-1) = graph_temp(1:n-1)
    p%graph(n:n1-1) = graph_temp(n+1:n1)
  endif

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% plot_line
!
! Output points used to construct the "line" associated with a plot curve.
!
! Notes
! -----
! Command syntax:
!   pipe plot_line {region_name}.{graph_name}.{curve_name} {x_or_y}
!
! Optional {x-or-y} may be set to "x" or "y" to get the smooth line points x or y 
! component put into the tao_c_interface_com%c_real array buffer.
! Note: The plot must come from a region, and not a template, since no template plots 
!       have associated line data.
! Examples:
!   pipe plot_line r13.g.a   ! String array output.
!   pipe plot_line r13.g.a x ! x-component of line points put in array buffer.
!   pipe plot_line r13.g.a y ! y-component of line points put in array buffer.
! 
! Parameters
! ----------
! region_name
! graph_name
! curve_name
! x_or_y : optional
!
! Returns
! -------
! string_list
!   if x_or_y == ''
! real_array
!   if x_or_y != ''
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_plot_line -external_plotting
!  args:
!    region_name: beta
!    graph_name: g
!    curve_name: a
!    x_or_y:
!
! Example: 2
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_plot_line -external_plotting
!  args:
!    region_name: beta
!    graph_name: g
!    curve_name: a
!    x_or_y: y

case ('plot_line')

  call string_trim(line(ix_line+1:), tail_str, ix2)
  line = line(1:ix_line)
  call tao_find_plots (err, line, 'REGION', curve = curves, only_visible = .false.)

  if (size(curves) /= 1) then
    call invalid ('Not a valid curve name')
    return
  endif

  c => curves(1)%c
  if (.not. allocated(c%x_line)) then
    call invalid ('No line associated with curve')
    return
  endif

if (.not. c%valid) then
  call invalid ('Invalid since: ' // c%why_invalid)
  return
endif

  n = size(c%x_line)

  select case (tail_str)
  case ('x', 'y')
    if (.not. allocated(tao_c_interface_com%c_real)) allocate (tao_c_interface_com%c_real(n))
    if (size(tao_c_interface_com%c_real) < n) then
      deallocate (tao_c_interface_com%c_real)
      allocate (tao_c_interface_com%c_real(n))
    endif

    tao_c_interface_com%n_real = n

    if (tail_str == 'x') then
      tao_c_interface_com%c_real(1:n) = c%x_line
    else
      tao_c_interface_com%c_real(1:n) = c%y_line
    endif

  case ('')
    do i = 1, n
      nl=incr(nl); write (li(nl), '(i0, 2(a, es22.14))') i, ';', c%x_line(i), ';', c%y_line(i)
    enddo

  case default
    call invalid ('word after curve name not "x" nor "y"')
  end select


!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% plot_symbol
!
! Output locations to draw symbols for a plot curve.
!
! Notes
! -----
! Command syntax:
!   pipe plot_symbol {region_name}.{graph_name}.{curve_name} {x_or_y}
!
! Optional {x_or_y} may be set to "x" or "y" to get the symbol x or y 
! positions put into the real array buffer.
! Note: The plot must come from a region, and not a template, 
!       since no template plots have associated symbol data.
! Examples:
!   pipe plot_symbol r13.g.a       ! String array output.
!   pipe plot_symbol r13.g.a x     ! x-component of the symbol positions 
!                                      loaded into the real array buffer.
!   pipe plot_symbol r13.g.a y     ! y-component of the symbol positions 
!                                      loaded into the real array buffer.
! 
! Parameters
! ----------
! region_name
! graph_name
! curve_name
! x_or_y
!
! Returns
! -------
! string_list
!   if x_or_y == ''
! real_array
!   if x_or_y != ''
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_plot_line -external_plotting
!  args:
!    region_name: r13
!    graph_name: g
!    curve_name: a
!    x_or_y: 
!
! Example: 2
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_plot_line -external_plotting
!  args:
!    region_name: r13
!    graph_name: g
!    curve_name: a
!    x_or_y: y

case ('plot_symbol')

  call string_trim(line(ix_line+1:), tail_str, ix2)
  line = line(1:ix_line)
  call tao_find_plots (err, line, 'REGION', curve = curves, only_visible = .false.)

  if (size(curves) /= 1) then
    call invalid ('Not a valid curve name.')
    return
  endif

  c => curves(1)%c
  if (.not. allocated(c%x_symb)) then
    call invalid ('No symbol array associated with curve.')
    return
  endif

  n = size(c%x_symb)

  select case (tail_str)
  case ('x', 'y')
    if (.not. allocated(tao_c_interface_com%c_real)) allocate (tao_c_interface_com%c_real(n))
    if (size(tao_c_interface_com%c_real) < n) then
      deallocate (tao_c_interface_com%c_real)
      allocate (tao_c_interface_com%c_real(n))
    endif

    tao_c_interface_com%n_real = n

    if (tail_str == 'x') then
      tao_c_interface_com%c_real(1:n) = c%x_symb
    else
      tao_c_interface_com%c_real(1:n) = c%y_symb
    endif

  case ('')
    do i = 1, size(c%x_symb)
      if (allocated(c%ix_symb)) then
        nl=incr(nl); write (li(nl), '(2(i0, a), 2(es22.14, a))') i, ';', c%ix_symb(i), ';', c%x_symb(i), ';', c%y_symb(i)
      else
        nl=incr(nl); write (li(nl), '(2(i0, a), 2(es22.14, a))') i, ';', 0, ';', c%x_symb(i), ';', c%y_symb(i)
      endif
    enddo

  case default
    call invalid ('word after curve name not "x" nor "y"')
  end select

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% plot_transfer
!
! Output transfer plot parameters from the "from plot" to the "to plot" (or plots).
!
! Notes
! -----
! Command syntax:
!   pipe plot_transfer {from_plot} {to_plot}
!
! To avoid confusion, use "@Tnnn" and "@Rnnn" syntax for {from_plot}.
! If {to_plot} is not present and {from_plot} is a template plot, the "to plots" 
!  are the equivalent region plots with the same name. And vice versa 
!  if {from_plot} is a region plot.
! 
! Parameters
! ----------
! from_plot
! to_plot
!
! Returns
! -------
! None
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_optics_matching
!  args:
!    from_plot: r13
!    to_plot: r23 

case ('plot_transfer')

  call tao_find_plots (err, line(1:ix_line), 'BOTH', plots, only_visible = .false.)
  if (size(plots) /= 1) then
    call invalid ('Number of "from plots" found is not exactly one.')
    return
  endif

  p => plots(1)%p

  call string_trim(line(ix_line+1:), line, ix_line)
  if (line == '') then
    if (associated(p%r)) then
      do i = 1, size(s%plot_page%template)
        if (s%plot_page%template(i)%name /= p%name) cycle
        call tao_plot_struct_transfer (p, s%plot_page%template(i))
      enddo

    else
      do i = 1, size(s%plot_page%region)
        if (s%plot_page%region(i)%plot%name /= p%name) cycle
        call tao_plot_struct_transfer (p, s%plot_page%region(i)%plot)
      enddo
    endif

  else
    call tao_find_plots (err, line(1:ix_line), 'BOTH', plots, only_visible = .false.)
    if (size(plots) == 0) then
      call invalid ('Number of "to plots" is zero.')
      return
    endif

    do i = 1, size(plots)
      call tao_plot_struct_transfer (p, plots(i)%p)
    enddo
  endif

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% plot1
!
! Output info on a given plot.
!
! Notes
! -----
! Command syntax:
!   pipe plot1 {name}
!
! {name} should be the region name if the plot is associated with a region.
! Output syntax is parameter list form. See documentation at the beginning of this file.
! 
! Parameters
! ----------
! name
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_optics_matching
!  args:
!    name: beta

case ('plot1')

  call tao_find_plots (err, line, 'BOTH', plots, print_flag = .false., only_visible = .false.)
  if (err) then
    call invalid ('Expect "r" or "t" at end.')
    return
  endif

  p => plots(1)%p

  n = 0
  if (allocated(p%graph)) n = size(p%graph)

  nl=incr(nl); write (li(nl), imt) 'num_graphs;INT;T;',                       n
  do i = 1, n
    nl=incr(nl); write (li(nl), vamt) 'graph[', i, '];STR;T;',              p%graph(i)%name
  enddo

  nl=incr(nl); write (li(nl), amt) 'name;STR;T;',                             trim(p%name)
  nl=incr(nl); write (li(nl), amt) 'description;STR;T;',                      trim(p%description)
  nl=incr(nl); write (li(nl), amt) 'x_axis_type;ENUM;T;',                     trim(p%x_axis_type)
  nl=incr(nl); write (li(nl), lmt) 'autoscale_x;LOGIC;T;',                    p%autoscale_x
  nl=incr(nl); write (li(nl), lmt) 'autoscale_y;LOGIC;T;',                    p%autoscale_y
  nl=incr(nl); write (li(nl), lmt) 'autoscale_gang_x;LOGIC;T;',               p%autoscale_gang_x
  nl=incr(nl); write (li(nl), lmt) 'autoscale_gang_y;LOGIC;T;',               p%autoscale_gang_y
  nl=incr(nl); write (li(nl), imt) 'n_curve_pts;INT;T;',                      p%n_curve_pts


!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% ptc_com
!
! Output Ptc_com structure components.
!
! Notes
! -----
! Command syntax:
!   pipe ptc_com
! 
! Returns
! -------
! string_list 
!
! Examples
! -------- 
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init 
!  args:

case ('ptc_com')

  nl=incr(nl); write (li(nl), imt) 'max_fringe_order;INT;T;',               ptc_com%max_fringe_order
  nl=incr(nl); write (li(nl), imt) 'old_integrator;INT;T;',                 ptc_com%old_integrator
  nl=incr(nl); write (li(nl), lmt) 'exact_model;LOGIC;T;',                  ptc_com%exact_model
  nl=incr(nl); write (li(nl), lmt) 'exact_misalign;LOGIC;T;',               ptc_com%exact_misalign
  nl=incr(nl); write (li(nl), lmt) 'translate_patch_drift_time;LOGIC;T;',   ptc_com%translate_patch_drift_time


!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% ring_general
!
! Output lattice branch with closed geometry info (emittances, etc.)
!
! Notes
! -----
! Command syntax:
!   pipe ring_general {ix_uni}@{ix_branch}|{which}
!
! where {which} is one of:
!   model
!   base
!   design
! Example:
!   pipe ring_general 1@0|model
! 
! Parameters
! ----------
! ix_uni : optional
! ix_branch : optional
! which : default=model
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!     ix_uni: 1
!     ix_branch: 0
!     which: model
!

case ('ring_general')

  u => point_to_uni(line, .true., err); if (err) return
  tao_lat => point_to_tao_lat(line, u, err); if (err) return
  ix_branch = parse_branch(line, u, .false., err); if (err) return
  tao_branch => tao_lat%tao_branch(ix_branch)
  branch => tao_lat%lat%branch(ix_branch)

  if (branch%param%geometry == open$) then
    call invalid ('Branch geometry must be closed')
    return
  endif

  call chrom_calc (branch%lat, s%global%delta_e_chrom, tao_branch%a%chrom, tao_branch%b%chrom, &
                                  pz = tao_branch%orbit(0)%vec(6), ix_branch = branch%ix_branch, orb0 = tao_branch%orbit(0))
  call calc_z_tune(branch)
  call radiation_integrals (branch%lat, tao_branch%orbit, tao_branch%modes_ri, tao_branch%ix_rad_int_cache, branch%ix_branch)
  call emit_6d(branch%ele(0), .true., tao_branch%modes_6d, mat6, tao_branch%orbit, tao_lat%rad_int_by_ele_6d)
  n = branch%n_ele_track
  time1 = branch%ele(n)%ref_time
  gamma = branch%ele(0)%value(e_tot$) / mass_of(branch%param%particle)

  nl=incr(nl); write (li(nl), rmt) 'param.unstable_factor;REAL;F;',             branch%param%unstable_factor
  nl=incr(nl); write (li(nl), rmt) 'Q_a;REAL;F;',                               branch%ele(n)%a%phi/twopi
  nl=incr(nl); write (li(nl), rmt) 'Q_b;REAL;F;',                               branch%ele(n)%b%phi/twopi
  nl=incr(nl); write (li(nl), rmt) 'Q_z;REAL;F;',                              -branch%z%tune/twopi
  nl=incr(nl); write (li(nl), rmt) 'Q_spin;REAL;F;',                            branch%param%spin_tune/twopi
  nl=incr(nl); write (li(nl), rmt) 'chrom_a;REAL;F;',                           tao_branch%a%chrom
  nl=incr(nl); write (li(nl), rmt) 'chrom_b;REAL;F;',                           tao_branch%b%chrom
  nl=incr(nl); write (li(nl), rmt) 'J_damp_a;REAL;F;',                          tao_branch%modes_6d%a%j_damp
  nl=incr(nl); write (li(nl), rmt) 'J_damp_b;REAL;F;',                          tao_branch%modes_6d%b%j_damp
  nl=incr(nl); write (li(nl), rmt) 'J_damp_z;REAL;F;',                          tao_branch%modes_6d%z%j_damp
  nl=incr(nl); write (li(nl), rmt) 'emit_a;REAL;F;',                            tao_branch%modes_6d%a%emittance
  nl=incr(nl); write (li(nl), rmt) 'emit_b;REAL;F;',                            tao_branch%modes_6d%b%emittance
  nl=incr(nl); write (li(nl), rmt) 'alpha_damp_a;REAL;F;',                      tao_branch%modes_6d%a%alpha_damp
  nl=incr(nl); write (li(nl), rmt) 'alpha_damp_b;REAL;F;',                      tao_branch%modes_6d%b%alpha_damp
  nl=incr(nl); write (li(nl), rmt) 'alpha_damp_z;REAL;F;',                      tao_branch%modes_6d%z%alpha_damp
  nl=incr(nl); write (li(nl), rmt) 'damping_time_a;REAL;F;',                    time1/tao_branch%modes_6d%a%alpha_damp
  nl=incr(nl); write (li(nl), rmt) 'damping_time_b;REAL;F;',                    time1/tao_branch%modes_6d%b%alpha_damp
  nl=incr(nl); write (li(nl), rmt) 'damping_time_z;REAL;F;',                    time1/tao_branch%modes_6d%z%alpha_damp
  nl=incr(nl); write (li(nl), rmt) 'sigE_E;REAL;F;',                            tao_branch%modes_6d%sigE_E
  nl=incr(nl); write (li(nl), rmt) 'sig_z;REAL;F;',                             tao_branch%modes_6d%sig_z
  nl=incr(nl); write (li(nl), rmt) 'energy_loss;REAL;F;',                       tao_branch%modes_6d%e_loss
  nl=incr(nl); write (li(nl), rmt) 'I0;REAL;F;',                                tao_branch%modes_ri%synch_int(0)
  nl=incr(nl); write (li(nl), rmt) 'I1;REAL;F;',                                tao_branch%modes_ri%synch_int(1)
  nl=incr(nl); write (li(nl), rmt) 'I2;REAL;F;',                                tao_branch%modes_ri%synch_int(2)
  nl=incr(nl); write (li(nl), rmt) 'I3;REAL;F;',                                tao_branch%modes_ri%synch_int(3)
  nl=incr(nl); write (li(nl), rmt) 'I4_a;REAL;F;',                              tao_branch%modes_ri%a%synch_int(4)
  nl=incr(nl); write (li(nl), rmt) 'I4_b;REAL;F;',                              tao_branch%modes_ri%b%synch_int(4)
  nl=incr(nl); write (li(nl), rmt) 'I5_a;REAL;F;',                              tao_branch%modes_ri%a%synch_int(5)
  nl=incr(nl); write (li(nl), rmt) 'I5_b;REAL;F;',                              tao_branch%modes_ri%b%synch_int(5)
  nl=incr(nl); write (li(nl), rmt) 'I6_g2_b;REAL;F;',                           tao_branch%modes_ri%b%synch_int(6) / gamma**2
  nl=incr(nl); write (li(nl), lmt) 'twiss_valid;LOGIC;F;',                      tao_branch%twiss_valid

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% shape_list
!
! Output lat_layout or floor_plan shapes list
!
! Notes
! -----
! Command syntax:
!   pipe shape_list {who}
!
! {who} is one of:
!   lat_layout
!   floor_plan
! 
! Parameters
! ----------
! who
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!    who: floor_plan  

case ('shape_list')

  select case (line)
  case ('lat_layout')
    shapes => s%plot_page%lat_layout%ele_shape
  case ('floor_plan')
    shapes => s%plot_page%floor_plan%ele_shape
  case default
    call invalid ('Bad {who}')
    return
  end select

  do i = 1, size(shapes)
    shape => shapes(i)
    if (shape%ele_id == '') cycle
    nl=incr(nl); write (li(nl), '(i0, 7a, es12.4, 3a, 2(l1, a), 2a)') i, ';', &
          trim(shape%ele_id), ';', trim(shape%shape), ';', trim(shape%color), ';', shape%size, ';', &
          trim(shape%label), ';', shape%draw, ';', shape%multi, ';', int_str(shape%line_width)
  enddo


!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% shape_manage
!
! Element shape creation or destruction
!
! Notes
! -----
! Command syntax:
!   pipe shape_manage {who} {index} {add_or_delete}
! 
! {who} is one of:
!   lat_layout
!   floor_plan
! {add_or_delete} is one of:
!   add     -- Add a shape at {index}. 
!              Shapes with higher index get moved up one to make room.
!   delete  -- Delete shape at {index}. 
!              Shapes with higher index get moved down one to fill the gap.
! 
! Example:
!   pipe shape_manage floor_plan 2 add
! Note: After adding a shape use "pipe shape_set" to set shape parameters.
! This is important since an added shape is in a ill-defined state.
! 
! Parameters
! ----------
! who
! index
! add_or_delete
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!    who: floor_plan
!    index: 1
!    add_or_delete: add

case ('shape_manage')

  call tao_next_switch (line, [character(12):: 'lat_layout', 'floor_plan'], .false., switch, err)
  select case (switch)
  case ('lat_layout')
    drawing => s%plot_page%lat_layout
  case ('floor_plan')
    drawing => s%plot_page%floor_plan
  case default
    call invalid ('Expected "floor_plan" or "lat_layout" after "shape_manage"')
    return
  end select

  n = size(drawing%ele_shape)
  ix = parse_int(line, err, 1, n+1); if (err) return

  call tao_next_switch (line, [character(8):: 'add', 'delete'], .false., switch, err)
  select case (switch)
  case ('add')
    call move_alloc (drawing%ele_shape, shapes_temp)
    allocate (drawing%ele_shape(n+1))
    drawing%ele_shape(1:ix-1) = shapes_temp(1:ix-1)
    drawing%ele_shape(ix+1:) = shapes_temp(ix:)
  case ('delete')
    call move_alloc (drawing%ele_shape, shapes_temp)
    allocate (drawing%ele_shape(n-1))
    drawing%ele_shape(1:ix-1) = shapes_temp(1:ix-1)
    drawing%ele_shape(ix:) = shapes_temp(ix+1:)
  case default
    call invalid ('Expected "add" or "delete" after shape index.')
    return
  end select

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% shape_pattern_list
!
! Output list of shape patterns or shape pattern points
!
! Notes
! -----
! Command syntax:
!   pipe shape_pattern_list {ix_pattern}
!
! If optional {ix_pattern} index is omitted then list all the patterns.
! If {ix_pattern} is present, list points of given pattern.
! 
! Parameters
! ----------
! ix_pattern : optional
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_shape
!  args:
!    ix_pattern: 

case ('shape_pattern_list')

  if (line == '') then
    do i = 1, size(s%plot_page%pattern)
      pattern => s%plot_page%pattern(i)
      nl=incr(nl); write (li(nl), '(2a, i0)') trim(pattern%name), ';', pattern%line%width
    enddo

  else
    ix = parse_int (line, err, 1, size(s%plot_page%pattern));  if (err) return
    pattern => s%plot_page%pattern(ix)
    do i = 1, size(pattern%pt)
      nl=incr(nl); write (li(nl), '(3a)') re_str(pattern%pt(i)%s, 6), ';', re_str(pattern%pt(i)%y, 6)
    enddo
  endif

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% shape_pattern_manage
!
! Add or remove shape pattern
!
! Notes
! -----
! Command syntax:
!   pipe shape_pattern_manage {ix_pattern}^^{pat_name}^^{pat_line_width}
!
! Where:
!   {ix_pattern}      -- Pattern index. Patterns with higher indexes will be moved up 
!                                       if adding a pattern and down if deleting.
!   {pat_name}        -- Pattern name.
!   {pat_line_width}  -- Line width. Integer. If set to "delete" then section 
!                                             will be deleted.
! 
! Parameters
! ----------
! ix_pattern
! pat_name
! pat_line_width
!
! Returns
! -------
! None
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_shape
!  args:
!    ix_pattern : 1
!    pat_name : new_pat
!    pat_line_width : 1

case ('shape_pattern_manage')

  n = size(s%plot_page%pattern)
  call split_this_line (line, name1, 3, err); if (err) return

  select case (name1(3))
  case ('delete')
    is = parse_int(name1(1), err, 1, n);  if (err) return
    call move_alloc(s%plot_page%pattern, pat_temp)
    allocate (s%plot_page%pattern(n-1))
    s%plot_page%pattern(1:is-1) = pat_temp(1:is-1)
    s%plot_page%pattern(is:) = pat_temp(is+1:)

  case default
    is = parse_int(name1(1), err, 1, n+1);  if (err) return
    call move_alloc(s%plot_page%pattern, pat_temp)
    allocate (s%plot_page%pattern(n+1))
    s%plot_page%pattern(1:is-1) = pat_temp(1:is-1)
    s%plot_page%pattern(is+1:) = pat_temp(is:)

    pattern => s%plot_page%pattern(is)
    pattern%name       = name1(2)
    pattern%line%width = parse_int(name1(3), err, 1)
  end select

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% shape_pattern_point_manage
!
! Add or remove shape pattern point
!
! Notes
! -----
! Command syntax:
!   pipe shape_pattern_point_manage {ix_pattern}^^{ix_point}^^{s}^^{x}
!
! Where:
!   {ix_pattern}      -- Pattern index.
!   {ix_point}        -- Point index. Points of higher indexes will be moved up
!                                     if adding a point and down if deleting.
!   {s}, {x}          -- Point location. If {s} is "delete" then delete the point.
! 
! Parameters
! ----------
! ix_pattern
! ix_point
! s
! x
!
! Returns
! -------
! None
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_shape
!  args:
!    ix_pattern: 1
!    ix_point: 1
!    s: 0
!    x: 0

case ('shape_pattern_point_manage')

  call split_this_line (line, name1, 4, err); if (err) return

  is = parse_int(name1(1), err, 1, size(s%plot_page%pattern));  if (err) return
  pattern => s%plot_page%pattern(is)
  n = 0
  if (allocated(pattern%pt)) n = size(pattern%pt)

  select case (name1(3))
  case ('delete')
    ip = parse_int(name1(2), err, 1, n)
    call move_alloc(pattern%pt, pat_pt_temp)
    allocate (pattern%pt(n-1))
    pattern%pt(1:ip-1) = pat_pt_temp(1:ip-1)
    pattern%pt(ip:) = pat_pt_temp(ip+1:)

  case default
    ip = parse_int(name1(2), err, 1, n+1)
    call move_alloc(pattern%pt, pat_pt_temp)
    allocate (pattern%pt(n+1))
    if (n > 0) then
      pattern%pt(1:ip-1) = pat_pt_temp(1:ip-1)
      pattern%pt(ip+1:) = pat_pt_temp(ip:)
    endif

    pattern%pt(ip)%s = parse_real(name1(3), err);  if (err) return
    pattern%pt(ip)%y = parse_real(name1(4), err);  if (err) return
  end select

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% shape_set
!
! Set lat_layout or floor_plan shape parameters.
!
! Notes
! -----
! Command syntax:
!   pipe shape_set {who}^^{shape_index}^^{ele_name}^^{shape}^^{color}^^
!                    {shape_size}^^{type_label}^^{shape_draw}^^
!                    {multi_shape}^^{line_width}
!
! {who} is one of:
!   lat_layout
!   floor_plan
! 
! Parameters
! ----------
! who
! shape_index
! ele_name
! shape
! color
! shape_size
! type_label
! shape_draw
! multi_shape
! line_width
!
! Returns
! -------
! None
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!    who: floor_plan
!    shape_index: 1
!    ele_name: Q1
!    shape: circle
!    color:
!    shape_size:
!    type_label:
!    shape_draw:
!    multi_shape:
!    line_width:
 

case ('shape_set')

  allocate (name_arr(10))
  call split_this_line (line, name_arr, -1, err, n_arr);  if (err) return

  select case (name_arr(1))
  case ('lat_layout')
    drawing => s%plot_page%lat_layout
  case ('floor_plan')
    drawing => s%plot_page%floor_plan
  case default
    call invalid ('Expected "floor_plan" or "lat_layout" after "shape_manage"')
    return
  end select

  ix = parse_int(name_arr(2), err, 1, size(drawing%ele_shape)); if (err) return

  shape_input%ele_id     = name_arr(3)
  shape_input%shape      = name_arr(4)
  shape_input%color      = name_arr(5)
  shape_input%size       = real_val(name_arr(6), 0.0_rp, err);   if (err) return
  shape_input%label      = name_arr(7)
  shape_input%draw       = logic_val(name_arr(8), .false., err);  if (err) return
  shape_input%multi      = logic_val(name_arr(9), .false., err);  if (err) return
  shape_input%line_width = int_val(name_arr(10), 1, err);         if (err) return

  drawing%ele_shape(ix) = tao_ele_shape_input_to_struct (shape_input)
  call tao_shape_init(drawing%ele_shape(ix), err, .true.)

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% show
!
! Output the output from a show command.
!
! Notes
! -----
! Command syntax:
!   pipe show {line}
!
! {line} is the string to pass through to the show command.
! Example:
!   pipe show lattice -pipe
! 
! Parameters
! ----------
! line
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!    line: -pipe

case ('show')

  call tao_show_this(trim(line), name, li, nl)

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% space_charge_com
!
! Output space_charge_com structure parameters.
!
! Notes
! -----
! Command syntax:
!   pipe space_charge_com
!
! Output syntax is parameter list form. See documentation at the beginning of this file.
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:

case ('space_charge_com')

  nl=incr(nl); write(li(nl), rmt) 'ds_track_step;REAL;T;',                    space_charge_com%ds_track_step
  nl=incr(nl); write(li(nl), rmt) 'dt_track_step;REAL;T;',                    space_charge_com%dt_track_step
  nl=incr(nl); write(li(nl), rmt) 'cathode_strength_cutoff;REAL;T;',          space_charge_com%cathode_strength_cutoff
  nl=incr(nl); write(li(nl), rmt) 'rel_tol_tracking;REAL;T;',                 space_charge_com%rel_tol_tracking
  nl=incr(nl); write(li(nl), rmt) 'abs_tol_tracking;REAL;T;',                 space_charge_com%abs_tol_tracking
  nl=incr(nl); write(li(nl), rmt) 'beam_chamber_height;REAL;T;',              space_charge_com%beam_chamber_height
  nl=incr(nl); write(li(nl), rmt) 'lsc_sigma_cutoff;REAL;T;',                 space_charge_com%lsc_sigma_cutoff
  nl=incr(nl); write(li(nl), rmt) 'particle_sigma_cutoff;REAL;T;',            space_charge_com%particle_sigma_cutoff

  nl=incr(nl); write(li(nl), '(a, 3(a, i0))') 'space_charge_mesh_size;INT_ARR;T', (';', space_charge_com%space_charge_mesh_size(j), j = 1, 3)
  nl=incr(nl); write(li(nl), '(a, 3(a, i0))') 'csr3d_mesh_size;INT_ARR;T',        (';', space_charge_com%csr3d_mesh_size(j), j = 1, 3)
  nl=incr(nl); write(li(nl), imt) 'n_bin;INT;T;',                             space_charge_com%n_bin
  nl=incr(nl); write(li(nl), imt) 'particle_bin_span;INT;T;',                 space_charge_com%particle_bin_span
  nl=incr(nl); write(li(nl), imt) 'n_shield_images;INT;T;',                   space_charge_com%n_shield_images
  nl=incr(nl); write(li(nl), imt) 'sc_min_in_bin;INT;T;',                     space_charge_com%sc_min_in_bin

  nl=incr(nl); write(li(nl), lmt) 'lsc_kick_transverse_dependence;LOGIC;T;',  space_charge_com%lsc_kick_transverse_dependence

  nl=incr(nl); write(li(nl), amt) 'diagnostic_output_file;STR;T;',            trim(space_charge_com%diagnostic_output_file)

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% species_to_int
!
! Convert species name to corresponding integer
!
! Notes
! -----
! Command syntax:
!   pipe species_to_int {species_str}
!
! Example:
!   pipe species_to_int CO2++
! 
! Parameters
! ----------
! species_str
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!    species_str: electron

case ('species_to_int')

  n = species_id(line)
  if (n == invalid$ .or. line == '') then
    call invalid ('Not a valid species name.')
    return
  endif

  nl=incr(nl); write (li(nl), '(i0)') n

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% species_to_str
!
! Convert species integer id to corresponding
!
! Notes
! -----
! Command syntax:
!   pipe species_to_str {species_int}
!
! Example:
!   pipe species_to_str -1     ! Returns 'Electron'
! 
! Parameters
! ----------
! species_int
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!    species_int: -1

case ('species_to_str')

  n = string_to_int (line, invalid$, err)
  name = species_name(n)

  if (err .or. name == invalid_name) then
    call invalid ('Not a valid species integer id number.')
    return
  endif

  nl=incr(nl); write (li(nl), '(a)') trim(name)

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% spin_invariant
!
! Output closed orbit spin axes n0, l0, or m0 at the ends of all lattice elements in a branch.
! n0, l0, and m0 are solutions of the T-BMT equation.
! n0 is periodic while l0 and m0 are not. At the beginning of the branch, the orientation of the 
! l0 or m0 axes in the plane perpendicular to the n0 axis is chosen a bit arbitrarily.
! See the Bmad manual for more details.
!
! Notes
! -----
! Command syntax:
!   pipe spin_invariant {flags} {who} {ix_uni}@{ix_branch}|{which}
!
! Where:
!   {flags}       - Optional flags (currently there is only one):
!                     -array_out  If present, the output will be available in 
!                                                 the tao_c_interface_com%c_real.
!   {who}         - One of: l0, n0, or m0
!   {ix_uni}      - A universe index. Defaults to s%global%default_universe.
!   {ix_branch}   - A branch index. Defaults to s%global%default_branch.
!   {which}       - Switch which is one of:
!                      model
!                      base
!                      design
!
! Example:
!   pipe spin_invariant 1@0|model
! 
! Note: This command is under development. If you want to use please contact David Sagan.
! 
! Parameters
! ----------
! who
! ix_uni : optional
! ix_branch : optional
! which : default=model
! flags : default=-array_out
!
! Returns
! -------
! string_list
!   if '-array_out' not in flags
! real_array
!   if '-array_out' in flags
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args: 
!    who: l0
!    ix_uni: 1
!    ix_branch: 0
!    which: model

case ('spin_invariant')

  use_real_array_buffer = .false.

  if (index('-array_out', line(1:ix_line)) == 1) then
    call string_trim(line(ix_line+1:), line, ix_line)
    use_real_array_buffer = .true.
  endif

  who = line(:ix_line)
  call string_trim(line(ix_line+1:), line, ix_line) 

  u => point_to_uni(line, .true., err); if (err) return
  tao_lat => point_to_tao_lat(line, u, err); if (err) return
  ix_branch = parse_branch(line, u, .false., err); if (err) return
  tao_branch => tao_lat%tao_branch(ix_branch)
  branch => tao_lat%lat%branch(ix_branch)

  if (.not. bmad_com%spin_tracking_on) call tao_spin_tracking_turn_on()

  orb = tao_lat%tao_branch(ix_branch)%orbit(0)
  n0 = orb%spin

  select case (who)
  case ('n0')
    ! Nothing to do
  case ('l0', 'm0')
    j = maxloc(abs(n0), 1)
    select case (j)
    case (1); l0 = [-n0(3), 0.0_rp, n0(1)]
    case (2); l0 = [n0(2), -n0(1), 0.0_rp]
    case (3); l0 = [0.0_rp, n0(3), -n0(2)]
    end select
    l0 = l0 / norm2(l0)
    m0 = cross_product(l0, n0)
    select case (who)
    case ('l0');  orb%spin = l0
    case ('m0');  orb%spin = m0
    end select
  case default
    call invalid ('BAD {WHO}: ' // who)
    return
  end select

  !

  if (use_real_array_buffer) then
    n = 3*branch%n_ele_track+3
    call re_allocate_c_double(tao_c_interface_com%c_real, n, .false.)
    tao_c_interface_com%n_real = n
  endif

  n = 0
  do ie = 0, branch%n_ele_track
    if (ie /= 0) call track1(orb, branch%ele(ie), branch%param, orb)
    nl=incr(nl); write (li(nl), '(i0, 3(a, es22.14))') ie, (';', orb%spin(j), j = 1, 3)
    if (use_real_array_buffer) then
      tao_c_interface_com%c_real(n+1:n+3) = orb%spin
      n = n + 3
    endif
  enddo

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% spin_polarization
!
! Output spin polarization information
!
! Notes
! -----
! Command syntax:
!   pipe spin_polarization {ix_uni}@{ix_branch}|{which}
!
! Where:
!   {ix_uni} is a universe index. Defaults to s%global%default_universe.
!   {ix_branch} is a branch index. Defaults to s%global%default_branch.
!   {which} is one of:
!     model
!     base
!     design
!
! Example:
!   pipe spin_polarization 1@0|model
! 
! Note: This command is under development. If you want to use please contact David Sagan.
! 
! Parameters
! ----------
! ix_uni : optional
! ix_branch : optional
! which : default=model
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args: 
!    ix_uni: 1
!    ix_branch: 0
!    which: model

case ('spin_polarization')

  u => point_to_uni(line, .true., err); if (err) return
  tao_lat => point_to_tao_lat(line, u, err); if (err) return
  ix_branch = parse_branch(line, u, .false., err); if (err) return
  tao_branch => tao_lat%tao_branch(ix_branch)
  branch => tao_lat%lat%branch(ix_branch)

  call tao_spin_polarization_calc (branch, tao_branch)

  z = anomalous_moment_of(branch%param%particle) * branch%ele(0)%value(e_tot$) / mass_of(branch%param%particle)
  nl=incr(nl); write (li(nl), rmt) 'anom_moment_times_gamma;REAL;F;',           z
  nl=incr(nl); write (li(nl), rmt) 'spin_tune;REAL;F;',                         tao_branch%spin%tune/twopi
  nl=incr(nl); write (li(nl), rmt) 'polarization_limit_st;REAL;F;',             tao_branch%spin%pol_limit_st
  nl=incr(nl); write (li(nl), rmt) 'polarization_limit_dk;REAL;F;',             tao_branch%spin%pol_limit_dk
  nl=incr(nl); write (li(nl), rmt) 'polarization_limit_dk_partial_a;REAL;F;',   tao_branch%spin%pol_limit_dk_partial(1)
  nl=incr(nl); write (li(nl), rmt) 'polarization_limit_dk_partial_b;REAL;F;',   tao_branch%spin%pol_limit_dk_partial(2)
  nl=incr(nl); write (li(nl), rmt) 'polarization_limit_dk_partial_c;REAL;F;',   tao_branch%spin%pol_limit_dk_partial(3)
  nl=incr(nl); write (li(nl), rmt) 'polarization_limit_dk_partial2_a;REAL;F;',  tao_branch%spin%pol_limit_dk_partial2(1)
  nl=incr(nl); write (li(nl), rmt) 'polarization_limit_dk_partial2_b;REAL;F;',  tao_branch%spin%pol_limit_dk_partial2(2)
  nl=incr(nl); write (li(nl), rmt) 'polarization_limit_dk_partial2_c;REAL;F;',  tao_branch%spin%pol_limit_dk_partial2(3)
  nl=incr(nl); write (li(nl), rmt) 'polarization_rate_bks;REAL;F;',             tao_branch%spin%pol_rate_bks
  nl=incr(nl); write (li(nl), rmt) 'depolarization_rate;REAL;F;',               tao_branch%spin%depol_rate
  nl=incr(nl); write (li(nl), rmt) 'depolarization_rate_partial_a;REAL;F;',     tao_branch%spin%depol_rate_partial(1)
  nl=incr(nl); write (li(nl), rmt) 'depolarization_rate_partial_b;REAL;F;',     tao_branch%spin%depol_rate_partial(2)
  nl=incr(nl); write (li(nl), rmt) 'depolarization_rate_partial_c;REAL;F;',     tao_branch%spin%depol_rate_partial(3)
  nl=incr(nl); write (li(nl), rmt) 'depolarization_rate_partial2_a;REAL;F;',    tao_branch%spin%depol_rate_partial2(1)
  nl=incr(nl); write (li(nl), rmt) 'depolarization_rate_partial2_b;REAL;F;',    tao_branch%spin%depol_rate_partial2(2)
  nl=incr(nl); write (li(nl), rmt) 'depolarization_rate_partial2_c;REAL;F;',    tao_branch%spin%depol_rate_partial2(3)

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% spin_resonance
!
! Output spin resonance information
!
! Notes 
! -----
! Command syntax:
!   pipe spin_resonance {ix_uni}@{ix_branch}|{which} {ref_ele}
!
! Where:
!   {ix_uni} is a universe index. Defaults to s%global%default_universe.
!   {ix_branch} is a lattice branch index. Defaults to s%global%default_branch.
!   {which} is one of: "model", "base" or "design"
!   {ref_ele} is an element name or index.
! This will return a string_list with the following fields:
!   spin_tune                   -- Spin tune
!   dq_X_sum, dq_X_diff         -- Tune sum Q_spin+Q_mode and tune difference 
!                                    Q_spin-Q_mode for modes X = a, b, and c.
!   xi_res_X_sum, xi_res_X_diff -- The linear spin/orbit sum and difference resonance 
!                                    strengths for X = a, b, and c modes.  
!
! Parameters
! ----------
! ix_uni : optional
! ix_branch : optional
! which : default=model
! ref_ele : default=0
!   Reference element to calculate at.
!
! Returns
! -------
! string_list  
!
! Examples
! --------
!
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args: 
!    ix_uni: 1
!    ix_branch: 0
!    which: model

case ('spin_resonance')

  u => point_to_uni(line, .true., err); if (err) return
  tao_lat => point_to_tao_lat(line, u, err, which, tail_str); if (err) return
  ix_branch = parse_branch(line, u, .false., err); if (err) return
  tao_branch => tao_lat%tao_branch(ix_branch)
  branch => tao_lat%lat%branch(ix_branch)

  if (tail_str == '') then
    ele => branch%ele(0)
  else
    call lat_ele_locator (tail_str, u%model%lat, eles, n_loc, err, ix_dflt_branch = branch%ix_branch)
    if (err .or. n_loc /= 1) then
      call invalid('UNIQUE LATTICE START ELEMENT NOT FOUND FOR: ' // quote(tail_str))
      return
    endif
    ele => eles(1)%ele
  endif

  datum%ix_branch = branch%ix_branch
  sm => datum%spin_map
  call tao_spin_matrix_calc (datum, u, ele, ele)
  call spin_mat_to_eigen (sm%map1%orb_mat, sm%map1%spin_q, eval, evec, n0, n_eigen, err)

  qs = branch%param%spin_tune/twopi
  nl=incr(nl); write (li(nl), rmt) 'spin_tune;REAL;F;',   qs

  do i = 1, 3
    j = 2 * i - 1
    q = atan2(aimag(eval(j)), real(eval(j),rp)) / twopi
    call spin_quat_resonance_strengths(evec(j,:), sm%map1%spin_q, xi_sum, xi_diff)
    nl=incr(nl); write (li(nl), amt) 'dq_', mode(i), '_sum;REAL;F;', re_str(modulo2(qs+q, 0.5_rp), 6)
    nl=incr(nl); write (li(nl), amt) 'dq_', mode(i), '_diff;REAL;F;', re_str(modulo2(qs-q, 0.5_rp), 6)
    nl=incr(nl); write (li(nl), amt) 'xi_res_', mode(i), '_sum;REAL;F;', re_str(xi_sum, 6)
    nl=incr(nl); write (li(nl), amt) 'xi_res_', mode(i), '_diff;REAL;F;', re_str(xi_diff, 6)
    nl=incr(nl); write (li(nl), rmt) 'n0;REAL_ARR;F;', n0(1), ';', n0(2), ';', n0(3)
  enddo

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% super_universe
!
! Output super_Universe parameters.
!
! Notes
! -----
! Command syntax:
!   pipe super_universe
! 
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args: 

case ('super_universe')

  nl=incr(nl); write (li(nl), imt) 'n_universe;INT;F;',                ubound(s%u, 1)
  nl=incr(nl); write (li(nl), imt) 'n_v1_var_used;INT;F;',             s%n_v1_var_used
  nl=incr(nl); write (li(nl), imt) 'n_var_used;INT;F;',                s%n_var_used

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% taylor_map
!
! Output Taylor map between two points.
!
! Notes
! -----
! Command syntax:
!   pipe taylor_map {ele1_id} {ele2_id} {order}
!
! Where:
!   {ele1_id}   - The start element.
!   {ele2_id}   - The end element.
!   {order}     - The map order. Default is order set in the lattice file. 
!                   {order} cannot be larger than what is set by the lattice file. 
!
! If {ele2_id} = {ele1_id}, the 1-turn transfer map is computed.
! Note: {ele2_id} should just be an element name or index without universe, 
!       branch, or model/base/design specification.
! Example:
!   pipe taylor_map 2@1>>q01w|design q02w  2
! 
! Parameters
! ----------
! ele1_id 
! ele2_id 
! order : default=1
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!    ele1_id: 1@0>>q01w|design
!    ele2_id: q02w
!    order: 1

case ('taylor_map')

  u => point_to_uni(line, .true., err); if (err) return
  tao_lat => point_to_tao_lat(line, u, err, which, tail_str); if (err) return
  ele1 => point_to_ele(line, tao_lat%lat, err); if (err) return
  ix_branch = ele1%ix_branch

  call string_trim(tail_str, tail_str, ix)
  call lat_ele_locator (tail_str(:ix), tao_lat%lat, eles, n_loc, err, ix_dflt_branch = ix_branch)
  if (err .or. n_loc == 0) then
    call invalid ('Bad ele2_id: ' // line)
    return
  endif
  if (n_loc > 1) then
    call invalid ('More than one element matches name: ' // line)
    return
  endif
  ele2 => eles(1)%ele

  n_order = string_to_int(tail_str(ix+1:), -1, err)
  if (err) then
    call invalid ('Invalid integer order: ' // quote(tail_str(ix+1:)))
    return
  endif
  if (n_order > ptc_private%taylor_order_ptc) then
    call invalid ('Taylor order cannot be above order set in lattice: ' // int_str(ptc_private%taylor_order_ptc))
    return
  endif

  call transfer_map_calc (tao_lat%lat, taylor, err, ele1%ix_ele, ele2%ix_ele, &
          tao_lat%tao_branch(ix_branch)%orbit(ele1%ix_ele), one_turn = .true., concat_if_possible = s%global%concatenate_maps)
  if (n_order > 0) call truncate_taylor_to_order (taylor, n_order, taylor)


  do i = 1, 6
    do j = 1, size(taylor(i)%term)
      tt => taylor(i)%term(j)
      nl=incr(nl); write (li(nl), '(i0, a, es22.14, 6(a, i0))') i, ';term;', tt%coef, (';', tt%expn(k), k = 1, 6)
    enddo
  enddo

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% twiss_at_s
!
! Output twiss parameters at given s position.
!
! Notes
! -----
! Command syntax:
!   pipe twiss_at_s {ix_uni}@{ele}->{s_offset}|{which}
!
! Where:
!   {ix_uni}    - A universe index. Defaults to s%global%default_universe.
!   {ele}       - An element name or index. Default is the Beginning element of branch 0.
!   {s_offset}  - Evaluation point offset from the downstream end of ele. Default is 0.
!                   If {s_offset} is present, "->" must also be present. 
!   {which}     - One of: "model", "base" or "design".
! 
! Parameters
! ----------
! ix_uni : optional
! ele : optional
! s_offset : optional
! which : default=model
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args: 
!    ix_uni: 1
!    ele: 10
!    s_offset: 0.7
!    which: model 

case ('twiss_at_s')

  u => point_to_uni(line, .true., err); if (err) return
  tao_lat => point_to_tao_lat(line, u, err); if (err) return
  s_pos = parse_ele_with_s_offset(line, tao_lat, ele, err); if (err) return
  ix_branch = ele%ix_branch

  call twiss_and_track_at_s (tao_lat%lat, s_pos, this_ele, tao_lat%tao_branch(ix_branch)%orbit, ix_branch = ix_branch)
  call twiss_out (this_ele%a, '', 'a')
  call twiss_out (this_ele%b, '', 'b')
  nl=incr(nl); write (li(nl), rmt) 'c_mat11;REAL;F;',            ele%c_mat(1,1)
  nl=incr(nl); write (li(nl), rmt) 'c_mat12;REAL;F;',            ele%c_mat(1,2)
  nl=incr(nl); write (li(nl), rmt) 'c_mat21;REAL;F;',            ele%c_mat(2,1)
  nl=incr(nl); write (li(nl), rmt) 'c_mat22;REAL;F;',            ele%c_mat(2,2)
  nl=incr(nl); write (li(nl), rmt) 'gamma_c;REAL;F;',            ele%gamma_c

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% universe
!
! Output universe info.
!
! Notes
! -----
! Command syntax:
!   pipe universe {ix_uni}
!
! Use "pipe global" to get the number of universes.
! 
! Parameters
! ----------
! ix_uni
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args: 
!    ix_uni: 1

case ('universe')

  u => point_to_uni(line, .false., err); if (err) return

  nl=incr(nl); write (li(nl), imt) 'ix_universe;INUM;F;',                     u%ix_uni
  nl=incr(nl); write (li(nl), imt) 'n_d2_data_used;INT;F;',                   u%n_d2_data_used
  nl=incr(nl); write (li(nl), imt) 'n_data_used;INT;F;',                      u%n_data_used
  nl=incr(nl); write (li(nl), lmt) 'is_on;LOGIC;T;',                          u%is_on

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% var
!
! Output parameters of a given variable.
!
! Notes
! -----
! Command syntax:
!   pipe var {var} {slaves}
!
! Note: use "pipe var_general" to get a list of variables.
! 
! Parameters
! ----------
! var
! slaves : optional
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_optics_matching
!  args: 
!    var: quad[1]
!    slaves:
!
! Example: 2
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_optics_matching
!  args: 
!    var: quad[1]
!    slaves: slaves

case ('var')

  call tao_find_var (err, line(:ix_line), v_array = v_array)

  if (.not. allocated(v_array) .or. size(v_array) /= 1) then
    call invalid ('Not a valid variable name')
    return
  endif

  call string_trim(line(ix_line+1:), line, ix_line)

  v_ptr => v_array(1)%v

  select case (line)
  case ('')
    nl=incr(nl); write (li(nl), rmt)  'model_value;REAL;T;',          v_ptr%model_value
    nl=incr(nl); write (li(nl), rmt)  'base_value;REAL;T;',           v_ptr%base_value

    nl=incr(nl); write (li(nl), amt) 'ele_name;STR;F;',                         trim(v_ptr%ele_name)
    nl=incr(nl); write (li(nl), amt) 'attrib_name;STR;F;',                      trim(v_ptr%attrib_name)
    nl=incr(nl); write (li(nl), imt) 'ix_v1;INT;F;',                            v_ptr%ix_v1
    nl=incr(nl); write (li(nl), imt) 'ix_var;INT;F;',                           v_ptr%ix_var
    nl=incr(nl); write (li(nl), imt) 'ix_dvar;INT;F;',                          v_ptr%ix_dvar
    nl=incr(nl); write (li(nl), imt) 'ix_attrib;INT;F;',                        v_ptr%ix_attrib
    nl=incr(nl); write (li(nl), imt) 'ix_key_table;INT;T;',                     v_ptr%ix_key_table
    nl=incr(nl); write (li(nl), rmt) 'design_value;REAL;F;',                    v_ptr%design_value
    nl=incr(nl); write (li(nl), rmt) 'scratch_value;REAL;F;',                   v_ptr%scratch_value
    nl=incr(nl); write (li(nl), rmt) 'old_value;REAL;F;',                       v_ptr%old_value
    nl=incr(nl); write (li(nl), rmt) 'meas_value;REAL;T;',                      v_ptr%meas_value
    nl=incr(nl); write (li(nl), rmt) 'ref_value;REAL;T;',                       v_ptr%ref_value
    nl=incr(nl); write (li(nl), rmt) 'correction_value;REAL;F;',                v_ptr%correction_value
    nl=incr(nl); write (li(nl), rmt) 'high_lim;REAL;T;',                        v_ptr%high_lim
    nl=incr(nl); write (li(nl), rmt) 'low_lim;REAL;T;',                         v_ptr%low_lim
    nl=incr(nl); write (li(nl), rmt) 'step;REAL;T;',                            v_ptr%step
    nl=incr(nl); write (li(nl), rmt) 'weight;REAL;T;',                          v_ptr%weight
    nl=incr(nl); write (li(nl), rmt) 'delta_merit;REAL;F;',                     v_ptr%delta_merit
    nl=incr(nl); write (li(nl), rmt) 'merit;REAL;F;',                           v_ptr%merit
    nl=incr(nl); write (li(nl), rmt) 'dmerit_dvar;REAL;F;',                     v_ptr%dMerit_dVar
    nl=incr(nl); write (li(nl), rmt) 'key_val0;REAL;F;',                        v_ptr%key_val0
    nl=incr(nl); write (li(nl), rmt) 'key_delta;REAL;T;',                       v_ptr%key_delta
    nl=incr(nl); write (li(nl), rmt) 's;REAL;F;',                               v_ptr%s
    nl=incr(nl); write (li(nl), amt) 'var^merit_type;ENUM;T;',                  trim(v_ptr%merit_type)
    nl=incr(nl); write (li(nl), lmt) 'exists;LOGIC;F;',                         v_ptr%exists
    nl=incr(nl); write (li(nl), lmt) 'good_var;LOGIC;F;',                       v_ptr%good_var
    nl=incr(nl); write (li(nl), lmt) 'good_user;LOGIC;T;',                      v_ptr%good_user
    nl=incr(nl); write (li(nl), lmt) 'good_opt;LOGIC;T;',                       v_ptr%good_opt
    nl=incr(nl); write (li(nl), lmt) 'good_plot;LOGIC;T;',                      v_ptr%good_plot
    nl=incr(nl); write (li(nl), lmt) 'useit_opt;LOGIC;F;',                      v_ptr%useit_opt
    nl=incr(nl); write (li(nl), lmt) 'useit_plot;LOGIC;F;',                     v_ptr%useit_plot
    nl=incr(nl); write (li(nl), lmt) 'key_bound;LOGIC;T;',                      v_ptr%key_bound

  case ('slaves')
    do i = 1, size(v_ptr%slave)
      nl=incr(nl); write (li(nl), '(3(i0, a))') v_ptr%slave(i)%ix_uni, ';', &
                                                       v_ptr%slave(i)%ix_branch, ';', v_ptr%slave(i)%ix_ele
    enddo

  case default
    call invalid ('BAD SWITCH: ' // line)
    return
  end select

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% var_create
!
! Create a single variable
!
! Notes
! -----
! Command syntax:
!   pipe var_create {var_name}^^{ele_name}^^{attribute}^^{universes}^^
!                     {weight}^^{step}^^{low_lim}^^{high_lim}^^{merit_type}^^
!                     {good_user}^^{key_bound}^^{key_delta}
!
! {var_name} is something like "kick[5]".
! Before using var_create, setup the appropriate v1_var array using 
! the "pipe var_v1_create" command.
! 
! Parameters
! ----------
! var_name
! ele_name
! attribute
! universes
! weight
! step
! low_lim
! high_lim
! merit_type
! good_user
! key_bound
! key_delta
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/tao.init_optics_matching
!  args:
!    var_name: quad[1]
!    ele_name: Q1
!    attribute: L
!    universes: 1
!    weight: 0.001
!    step: 0.001
!    low_lim: -10
!    high_lim: 10
!    merit_type: 
!    good_user: T
!    key_bound: T
!    key_delta: 0.01 

case ('var_create')

  allocate(name_arr(12))
  call split_this_line (line, name_arr, 12, err);         if (err) return

  call tao_find_var (err, name_arr(1), v_array = v_array);
  if (err .or. size(v_array) /= 1) then
    call invalid('BAD VARIABLE NAME')
    return
  endif
  v_ptr => v_array(1)%v

  call tao_pick_universe ('[' // trim(name_arr(4)) // ']@', name, picked, err, dflt_uni = -1);
  if (err .or. name /= '') then
    call invalid('INVALID UNIVERSE SPECIFICATION')
    return
  endif

  ele_name = upcase(name_arr(2))
  attrib_name = upcase(name_arr(3))

  v_ptr%exists      = .true.
  v_ptr%merit_type  = str_val(name_arr(9), 'limit')
  v_ptr%good_var    = .true.
  v_ptr%good_opt    = .true.
  v_ptr%good_user   = logic_val(name_arr(10), .true., err);    if (err) return
  v_ptr%ele_name    = ele_name
  v_ptr%attrib_name = attrib_name
  v_ptr%low_lim     = real_val(name_arr(7), -1d30, err);       if (err) return
  v_ptr%high_lim    = real_val(name_arr(8), 1d30, err);        if (err) return
  v_ptr%weight      = real_val(name_arr(5), 0.0_rp, err);      if (err) return
  v_ptr%step        = real_val(name_arr(6), 0.0_rp, err);      if (err) return
  v_ptr%key_bound   = logic_val(name_arr(11), .false., err);   if (err) return
  v_ptr%key_delta   = real_val(name_arr(12), 0.0_rp, err);     if (err) return

  call tao_var_stuffit2(picked, v_ptr, 'pipe var_create command')

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% var_general
!
! Output list of all variable v1 arrays
!
! Notes
! -----
! Command syntax:
!   pipe var_general
!
! Output syntax:
!   {v1_var name};{v1_var%v lower bound};{v1_var%v upper bound}
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:

case ('var_general')

  do i = 1, s%n_v1_var_used
    v1_ptr => s%v1_var(i)
    if (v1_ptr%name == '') cycle
    call location_encode (line, v1_ptr%v%useit_opt, v1_ptr%v%exists, lbound(v1_ptr%v, 1))
    nl=incr(nl); write (li(nl), '(4a, 2(i0, a))') trim(v1_ptr%name), ';', trim(line), ';', lbound(v1_ptr%v, 1), ';', ubound(v1_ptr%v, 1)
  enddo

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% var_v_array
!
! Output list of variables for a given data_v1.
!
! Notes
! -----
! Command syntax:
!   pipe var_v_array {v1_var}
!
! Example:
!   pipe var_v_array quad_k1
! 
! Parameters
! ----------
! v1_var
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!    v1_var: quad_k1

case ('var_v_array')

  call tao_find_var (err, line, v_array = v_array)

  if (.not. allocated(v_array)) then
    call invalid ('Not a valid v1_var name')
    return
  endif

  do i = 1, size(v_array)
    v_ptr => v_array(i)%v
    if (.not. v_ptr%exists) cycle
    nl=incr(nl); write(li(nl), '(i0, 3a, 3(es22.14, a), 2(l1, a), es22.14)') &
                  v_ptr%ix_v1, ';', trim(tao_var_attrib_name(v_ptr)), ';', v_ptr%meas_value, ';', &
                  v_ptr%model_value, ';', v_ptr%design_value, ';', v_ptr%useit_opt, ';', v_ptr%good_user, ';', v_ptr%weight
  enddo


!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% var_v1_array
!
! Output list of variables in a given variable v1 array
!
! Notes
! -----
! Command syntax:
!   pipe var_v1_array {v1_var}
! 
! Parameters
! ----------
! v1_var
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!    v1_var: quad_k1 

case ('var_v1_array')

  call tao_find_var (err, line, v1_array = v1_array)

  if (err .or. .not. allocated(v1_array)) then
    call invalid ('Not a valid v1 name')
    return
  endif

  v1_ptr => v1_array(1)%v1

  do i = lbound(v1_ptr%v, 1), ubound(v1_ptr%v, 1)
    v_ptr => v1_ptr%v(i)
    if (.not. v_ptr%exists) cycle
    nl=incr(nl); write (li(nl), '(2a, i0, 5a, 3(es22.14, a), 2 (l1, a))') trim(v1_ptr%name), '[', &
                     v_ptr%ix_v1, '];', trim(v_ptr%ele_name), ';', trim(v_ptr%attrib_name), ';', &
                     v_ptr%meas_value, ';', v_ptr%model_value, ';', &
                     v_ptr%design_value, ';', v_ptr%good_user, ';', v_ptr%useit_opt
  enddo

  nl=incr(nl); write (li(nl), imt) 'ix_v1_var;INT;F;',                       v1_ptr%ix_v1_var

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% var_v1_create
!
! Create a v1 variable structure along with associated var array.
!
! Notes
! -----
! Command syntax:
!   pipe var_v1_create {v1_name} {n_var_min} {n_var_max}
!
! {n_var_min} and {n_var_max} are the lower and upper bounds of the var
! Example:
!   pipe var_v1_create quad_k1 0 45
! This example creates a v1 var structure called "quad_k1" with an associated
! variable array that has the range [0, 45].
! 
! Use the "pipe var_create" and "set variable" commands to set variable parameters.
! Note: When setting multiple variable parameters, first set
!   set global lattice_calc_on = F")
! to prevent Tao trying to evaluate the 
! partially created variable and generating unwanted error messages.
! 
! Parameters
! ----------
! v1_name
! n_var_min
! n_var_max
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!    v1_name: quad_k1 
!    n_var_min: 0 
!    n_var_max: 45 

case ('var_v1_create')

  call tao_cmd_split (line, 3, name1, .true., err)

  call tao_find_var(err, name1(1), v1_array = v1_array, print_err = .false.)
  if (size(v1_array) == 1) call destroy_this_var_v1(name1(1))

  if (.not. is_integer(name1(2)) .or. .not. is_integer(name1(3))) then
    call invalid ('Is Malformed')
    return
  endif

  read (name1(2), *) ix_min(1)
  read (name1(3), *) ix_max(1)

  n1 = size(s%v1_var)
  if (s%n_v1_var_used + 1 > n1) then
    call move_alloc (s%v1_var, v1_temp)
    allocate (s%v1_var(s%n_v1_var_used + 1))
    s%v1_var(1:n1) = v1_temp
  endif

  n = size(s%var)
  n_delta = ix_max(1) + 1 - ix_min(1)
  if (s%n_var_used + n_delta  > n) then
    call move_alloc (s%var, v_temp)
    allocate (s%var(s%n_var_used+n_delta))
    s%var(1:n) = v_temp
    do k = s%n_var_used+1, size(s%var)
      s%var(k)%ix_var = k
    enddo
  endif

  i2 = 0   ! In case there are no v1 structures yet defined.

  do i = 1, s%n_v1_var_used
    v1_ptr => s%v1_var(i)
    i1 = lbound(v1_ptr%v, 1)
    i1 = v1_ptr%v(i1)%ix_var
    i2 = ubound(v1_ptr%v, 1)
    i2 = v1_ptr%v(i2)%ix_var
    call tao_point_v1_to_var (v1_ptr, s%var(i1:i2), s%var(i1)%ix_v1)
  enddo

  nn = s%n_v1_var_used + 1
  s%n_v1_var_used = nn
  s%v1_var(nn)%ix_v1_var = nn
  s%v1_var(nn)%name = name1(1)
  s%n_var_used = s%n_var_used + n_delta
  i1 = i2 + 1
  i2 = i2 + n_delta
  call tao_point_v1_to_var (s%v1_var(nn), s%var(i1:i2), ix_min(1))

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% var_v1_destroy
!
! Destroy a v1 var structure along with associated var sub-array.
!
! Notes
! -----
! Command syntax:
!   pipe var_v1_destroy {v1_datum}
! 
! Parameters
! ----------
! v1_datum
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!    v1_datum: quad_k1

case ('var_v1_destroy')

  call destroy_this_var_v1(line)

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% wall3d_radius
!
! Output vaccum chamber wall radius for given s-position and angle in (x,y) plane.
! The radius is with respect to the local wall origin which may not be the (x,y) = (0,0) origin.
!
! Notes
! -----
! Command syntax:
!   pipe wall3d_radius {ix_uni}@{ix_branch} {s_position} {angle}
!
! Where:
!   {ix_uni} is a universe index. Defaults to s%global%default_universe.
!   {ix_branch} is a lattice branch index. 
!   {s_position} is the s-position to evaluate at.
!   {angle} is the angle to evaluate at.
!
! Parameters
! ----------
! ix_uni : ""
! ix_branch : ""
! s_position
! angle
!
! Returns
! -------
! string_list
!
! Examples
! --------

case ('wall3d_radius')

  u => s%u(s%global%default_universe)
  if (index(line, '@') /= 0) then
    u => point_to_uni(line, .true., err); if (err) return
  endif

  ix_branch = parse_branch(line, u, .false., err); if (err) return
  branch => u%model%lat%branch(ix_branch)
  call split_this_line (line, name1, 2, err, space_sep = .true.); if (err) return
  s_here = parse_real(name1(1), err); if (err) return
  angle  = parse_real(name1(2), err); if (err) return

  !

  ele => pointer_to_element_at_s (branch, s_here, .true., err)
  if (err) then
    call invalid ('Bad s-position')
    return
  endif

  vec0 = [cos(angle), 0.0_rp, sin(angle), 0.0_rp, s_here - ele%s_start, 1.0_rp]

  value =  wall3d_d_radius(vec0, ele, 1, perp, ix, not, origin, r_wall, err)
  if (not) then
    call invalid ('No Wall3D here')
    return
  elseif (err) then
    call invalid ('Error in radius calc. Please report this.')
    return
  endif

  nl=incr(nl); write (li(nl), rmt)  'wall_radius;REAL;F;',            r_wall
  nl=incr(nl); write (li(nl), ramt) 'origin;REAL_ARR;F',              (';', origin(i), i = 1, 3)
  nl=incr(nl); write (li(nl), ramt) 'perpendicular;REAL_ARR;F',       (';', perp(i), i = 1, 3)

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!%% wave
!
! Output Wave analysis info.
!
! Notes
! -----
! Command syntax:
!   pipe wave {who}
!
! Where {who} is one of:
!   params
!   loc_header
!   locations
!   plot1, plot2, plot3
! 
! Parameters
! ----------
! who
!
! Returns
! -------
! string_list
!
! Examples
! --------
! Example: 1
!  init: -init $ACC_ROOT_DIR/regression_tests/pipe_test/cesr/tao.init
!  args:
!    who: params

case ('wave')

  select case (line)
  case ('params')

    nl=incr(nl); write (li(nl), amt) 'wave_data_type;ENUM;T;',          trim(s%wave%data_type)
    nl=incr(nl); write (li(nl), imt) 'ix_a1;INT;T;',                    s%wave%ix_a1
    nl=incr(nl); write (li(nl), imt) 'ix_a2;INT;T;',                    s%wave%ix_a2
    nl=incr(nl); write (li(nl), imt) 'ix_b1;INT;T;',                    s%wave%ix_b1
    nl=incr(nl); write (li(nl), imt) 'ix_b2;INT;T;',                    s%wave%ix_b2

    select case (s%wave%data_type)
    case ('orbit.x', 'orbit.y', 'eta.x', 'eta.y', 'beta.a', 'beta.b', 'ping_a.amp_x', 'ping_b.amp_y')
      nl=incr(nl); write(li(nl), '(a, f8.3)') 'A Region Sigma_Fit/Amp_Fit;REAL;F;', s%wave%rms_rel_a
      nl=incr(nl); write(li(nl), '(a, f8.3)') 'B Region Sigma_Fit/Amp_Fit;REAL;F;', s%wave%rms_rel_b
      nl=incr(nl); write(li(nl), '(a, f8.3)') 'Sigma_Kick/Kick;REAL;F;', s%wave%rms_rel_k
      nl=incr(nl); write(li(nl), '(a, f8.3)') 'Sigma_phi;REAL;F;', s%wave%rms_phi

    case ('phase.a', 'phase.b', 'ping_a.phase_x', 'ping_b.phase_y')
      nl=incr(nl); write(li(nl), '(a, f8.3)') 'A Region Sigma_Fit/Amp_Fit;REAL;F;', s%wave%rms_rel_a
      nl=incr(nl); write(li(nl), '(a, f8.3)') 'B Region Sigma_Fit/Amp_Fit;REAL;F;', s%wave%rms_rel_b
      nl=incr(nl); write(li(nl), '(a, f8.3)') 'Sigma_Kick/Kick;REAL;F;', s%wave%rms_rel_k
      nl=incr(nl); write(li(nl), '(a, f8.3)') 'Sigma_phi;REAL;F;', s%wave%rms_phi
      nl=incr(nl); write(li(nl), '(a, f8.3, a)') 'Chi_C [Figure of Merit];REAL;F;', s%wave%chi_c

    case ('ping_a.amp_sin_rel_y', 'ping_a.amp_cos_rel_y', 'ping_b.amp_sin_rel_x', 'ping_b.amp_cos_rel_x', &
          'ping_a.amp_sin_y', 'ping_a.amp_cos_y', 'ping_b.amp_sin_x', 'ping_b.amp_cos_x', 'cbar.11', 'cbar.12', 'cbar.22')
      nl=incr(nl); write(li(nl), '(a, f8.3)') 'A Region Sigma_+/Amp_+;REAL;F;', s%wave%rms_rel_as
      nl=incr(nl); write(li(nl), '(a, f8.3)') 'A Region Sigma_-/Amp_-;REAL;F;', s%wave%rms_rel_ar
      nl=incr(nl); write(li(nl), '(a, f8.3)') 'B Region Sigma_+/Amp_+;REAL;F;', s%wave%rms_rel_bs
      nl=incr(nl); write(li(nl), '(a, f8.3)') 'B Region Sigma_-/Amp_-;REAL;F;', s%wave%rms_rel_br
      nl=incr(nl); write(li(nl), '(a, f8.3)') 'Kick |K+|', 2*s%wave%amp_ba_s
      nl=incr(nl); write(li(nl), '(a, f8.3)') 'Sigma_K+/K+', 2*s%wave%rms_rel_ks
      nl=incr(nl); write(li(nl), '(a, f8.3)') 'Kick |K-|', 2*s%wave%amp_ba_r
      nl=incr(nl); write(li(nl), '(a, f8.3)') 'Sigma_K-/K-', 2*s%wave%rms_rel_kr
      nl=incr(nl); write(li(nl), '(a, f8.3)') 'Sigma_phi+;REAL;F;', s%wave%rms_phi_s
      nl=incr(nl); write(li(nl), '(a, f8.3)') 'Sigma_phi-;REAL;F;', s%wave%rms_phi_r
      nl=incr(nl); write(li(nl), '(a, f8.3, a)') 'Chi_a [Figure of Merit];REAL;F;', s%wave%chi_a
    end select

  !

  case ('loc_header')
    select case (s%wave%data_type)
    case ('beta.a', 'beta.b')
      nl=incr(nl); li(nl) = 'header1;STR;F;Normalized Kick = kick * beta  [urad * meter]'
      nl=incr(nl); li(nl) = 'columns;After Dat#;Norm_Kick;s;ix_ele;ele@kick;phi'

    case ('orbit.x', 'orbit.y', 'eta.x', 'eta.y', 'ping_a.amp_x', 'ping_b.amp_y')
      nl=incr(nl); li(nl) = 'header1;STR;F;Normalized Kick = kick * sqrt(beta)  [urad * sqrt(meter)]'
      nl=incr(nl); li(nl) = 'columns;After Dat#;Norm_Kick;s;ix_ele;ele@kick;phi'

    case ('phase.a', 'phase.b', 'ping_a.phase_x', 'ping_b.phase_y')
      nl=incr(nl); li(nl) = 'header1;STR;F;Normalized Kick = k * l * beta [dimensionless]'
      nl=incr(nl); li(nl) = 'header2;STR;F;where k = quadrupole gradient [rad/m^2].'
      nl=incr(nl); li(nl) = 'columns;After Dat#;Norm_Kick;s;ix_ele;ele@kick;phi'

    case ('ping_a.amp_sin_rel_y', 'ping_a.amp_cos_rel_y', 'ping_b.amp_sin_rel_x', 'ping_b.amp_cos_rel_x', &
          'ping_a.amp_sin_y', 'ping_a.amp_cos_y', 'ping_b.amp_sin_x', 'ping_b.amp_cos_x', 'cbar.11', 'cbar.12', 'cbar.22')
      nl=incr(nl); li(nl) = 'columns;After Dat#;Norm_K;phi+;phi-;phi_a;phi_b'
    end select

  !

  case ('locations')
    select case (s%wave%data_type)
    case ('orbit.x', 'orbit.y', 'eta.x', 'eta.y', 'beta.a', 'beta.b', 'ping_a.amp_x', 'ping_b.amp_y')
      do i = 1, s%wave%n_kick
        wk => s%wave%kick(i)
        nl=incr(nl); write(li(nl), '(i0, a, f0.2, a, f0.2, a, i0, 3a, f0.3)') wk%ix_dat_before_kick, ';', 1e6*wk%amp, ';', wk%s, ';', wk%ele%ix_ele, ';', trim(wk%ele%name), ';', wk%phi
      enddo

    case ('phase.a', 'phase.b', 'ping_a.phase_x', 'ping_b.phase_y')
      do i = 1, s%wave%n_kick
        wk => s%wave%kick(i)
        nl=incr(nl); write(li(nl), '(i0, a, f0.4, a, f0.2, a, i0, 3a, f0.3)') wk%ix_dat_before_kick, ';', wk%amp, ';', wk%s, ';', wk%ele%ix_ele, ';', trim(wk%ele%name), ';', wk%phi
      enddo

    case ('ping_a.amp_sin_rel_y', 'ping_a.amp_cos_rel_y', 'ping_b.amp_sin_rel_x', 'ping_b.amp_cos_rel_x', &
          'ping_a.amp_sin_y', 'ping_a.amp_cos_y', 'ping_b.amp_sin_x', 'ping_b.amp_cos_x', 'cbar.11', 'cbar.12', 'cbar.22')
      do i = 1, s%wave%n_kick
        wk => s%wave%kick(i)
        nl=incr(nl); write(li(nl), '(i0, a, f0.4, a, f0.2, a, i0, 3a, 4(f0.3, a))') wk%ix_dat_before_kick, ';', &
                  wk%amp, ';', wk%s, ';', wk%ele%ix_ele, ';', trim(wk%ele%name), ';', wk%phi_s, ';', wk%phi_r, ';', (wk%phi_s+wk%phi_r)/2, ';', (wk%phi_s-wk%phi_r)/2
      enddo
    end select

  case ('plot1', 'plot2', 'plot3')
    if (.not. associated(s%wave%region)) then
      call invalid ('Wave plot regions not yet setup.')
      return
    endif

    select case (line)
    case ('plot1')
      g => s%wave%region%plot%graph(1)
    case ('plot2')
      g => s%wave%region%plot%graph(2)
    case ('plot3')
      g => s%wave%region%plot%graph(3)
    end select

    c => g%curve(1)
    do i = 1, size(c%x_symb)
      nl=incr(nl); write (li(nl), '(i0, 2(a, es14.6))') i, ';', c%x_symb(i), ';', c%y_symb(i)
    enddo

  case default
    call invalid ('Bad {who}: ' // line)
  end select

!----------------------------------------------------------------------

case default

  call out_io (s_error$, r_name, "pipe command internal error, shouldn't be here!")

end select

call end_stuff(li, nl)

!----------------------------------------------------------------------
! return through scratch

contains

subroutine end_stuff(li, nl)

type (out_io_output_direct_struct) out_dir_state
character(n_char_show), allocatable :: li(:)
integer nl, i

!

call output_direct (get = out_dir_state)
call output_direct(print_and_capture = .true.)

if (doprint) then
  call out_io (s_blank$, r_name, li(1:nl))
endif

if (opened) then
  do i = 1, nl
    write (iu_write, '(a)') trim(li(i))
  enddo
  close (iu_write)
endif

call output_direct (get = out_dir_state)

end subroutine

!----------------------------------------------------------------------
! contains

function point_to_uni (line, compound_word, err) result (u)

type (tao_universe_struct), pointer :: u
integer ix, ix_universe, ios
logical compound_word, err
character(*) line

! A compound_word is something like "2@q10w" or "q10w". A non-compound word is something like "2" which
! just represents a universe index.

nullify(u)
err = .false.

if (compound_word) then
  ix = tao_uni_ampersand_index(line)
  if (ix == 0) then
    ix_universe = s%global%default_universe
  elseif (line(1:ix-1) == '') then
    ix_universe = s%global%default_universe
    line = line(ix+1:)
  else
    read (line(1:ix-1), *,  iostat = ios)  ix_universe
    if (ios /= 0) ix_universe = -999
    line = line(ix+1:)
  endif
else
  if (len_trim(line) == 0) then
    ix_universe = s%global%default_universe
  else
    ! In this case line is just a universe number
    read (line, *,  iostat = ios)  ix_universe
    if (ios /= 0) ix_universe = -999
  endif
endif

u => tao_pointer_to_universe(ix_universe, .true.)

if (.not. associated(u)) then
  call invalid ('bad universe index')
  err = .true.
endif

end function point_to_uni

!----------------------------------------------------------------------
! contains

function incr(n) result (n1)

integer n, n1

n1 = n + 1
if (n1 > size(li)) call re_allocate_lines (li, int(1.5 * n1))

end function

!----------------------------------------------------------------------
! contains

subroutine re_allocate_lines (li, n_lines)

character(n_char_show), allocatable :: li(:)
integer n_lines

!

if (.not. allocated(li)) allocate (li(n_lines))
if (size(li) < n_lines) call re_allocate (li, n_lines)

end subroutine re_allocate_lines

!----------------------------------------------------------------------
! contains

function point_to_tao_lat (line, u, err, which, tail_str) result (tao_lat)

type (tao_lattice_struct), pointer :: tao_lat
type (tao_universe_struct) u
integer ix
logical err
character(*) line
character(*), optional :: which, tail_str


err = .true.
nullify(tao_lat)

call string_trim(line, line, ix)
if (present(tail_str)) call string_trim(line(ix+1:), tail_str, i)
line = line(1:ix)

ix = index(line, '|')
if (ix == 0) then
  tao_lat => u%model
  err = .false.
  return
endif

select case (line(ix+1:))
case ('model')
  tao_lat => u%model
case ('base')
  tao_lat => u%base
case ('design')
  tao_lat => u%design
case ('')
  tao_lat => u%model
case default
  call invalid ('Expecting "|{which}" where {which} must be one of "model", "base", or "design"')
  return
end select

if (present(which)) then
  which = line(ix+1:)
  if (which == '') which = 'model'
endif

line = line(1:ix-1)
err = .false.

end function point_to_tao_lat

!----------------------------------------------------------------------
! contains

function point_to_ele (line, lat, err) result (ele)

type (lat_struct) lat
type (ele_struct), pointer :: ele
type (ele_pointer_struct), allocatable :: eles(:)
integer n_loc
character(*) line
logical err

!

err = .true.
nullify(ele)
call lat_ele_locator (line, lat, eles, n_loc)

select case (n_loc)
case (0) 
  call invalid ('Cannot locate element.')
  return
case (1)
  ! Good
case default
  call invalid ('Multiple matches to element.')
  return
end select

ele => eles(1)%ele
err = .false.

end function point_to_ele

!----------------------------------------------------------------------
! contains

function parse_branch (line, u, has_separator, err) result (ix_branch)

type (tao_universe_struct) u
integer ix, ios, ix_branch
logical has_separator, err
character(*) line
character(40) str

!

err = .false.
ix_branch = s%global%default_branch
if (line(1:1) == ' ') return

if (has_separator) then
  ix = index(line, '>>')

  if (ix == 0) then
    call invalid ('Missing ">>"')
    err = .true.
    return
  endif

  if (ix /= 1) then
    read (line(1:ix-1), *, iostat = ios) ix_branch
    if (ios /= 0) ix_branch = -999
  endif
  line = line(ix+2:)

elseif (len_trim(line) /= 0) then
  read (line, *, iostat = ios) ix_branch
  if (ios /= 0) ix_branch = -999
  call string_trim(line, line, ix)
  call string_trim(line(ix+1:), line, ix)
endif

if (ix_branch < 0 .or. ix_branch > ubound(u%design%lat%branch, 1)) then
  call invalid ('Out of range branch index')
  err = .true.
  return
endif

end function parse_branch

!----------------------------------------------------------------------
! contains

function parse_real (line, err) result (a_real)

real(rp) a_real
logical err
character(*) line

a_real = string_to_real (line, real_garbage$, err)
if (err .or. a_real == real_garbage$) then
  call invalid ('Bad real number')
  return
endif

end function parse_real

!----------------------------------------------------------------------
! contains

function parse_int (line, err_flag, min_bound, max_bound, dflt_val) result (a_int)

integer a_int
integer, optional :: min_bound, max_bound, dflt_val
logical err, err_flag
character(*) line

!

err_flag = .true.


a_int = string_to_int (line, integer_option(int_garbage$, dflt_val), err)

if (err .or. a_int == int_garbage$) then
  call invalid ('Bad int number')
  return
endif

if (present(min_bound)) then
  if (a_int < min_bound) then
    call invalid ('Integer below lower bound')
    return
  endif
endif

if (present(max_bound)) then
  if (a_int > max_bound) then
    call invalid ('Integer above upper bound')
    return
  endif
endif

err_flag = .false.

end function parse_int

!----------------------------------------------------------------------
! contains

subroutine orbit_out (orbit)

type (coord_struct) orbit

nl=incr(nl); write (li(nl), rmt) 'x;REAL;F;',                                orbit%vec(1)
nl=incr(nl); write (li(nl), rmt) 'px;REAL;F;',                               orbit%vec(2)
nl=incr(nl); write (li(nl), rmt) 'y;REAL;F;',                                orbit%vec(3)
nl=incr(nl); write (li(nl), rmt) 'py;REAL;F;',                               orbit%vec(4)
nl=incr(nl); write (li(nl), rmt) 'z;REAL;F;',                                orbit%vec(5)
nl=incr(nl); write (li(nl), rmt) 'pz;REAL;F;',                               orbit%vec(6)

nl=incr(nl); write (li(nl), ramt) 'spin;REAL_ARR;F',                         (';', orbit%spin(i), i = 1, 3)
nl=incr(nl); write (li(nl), ramt) 'field;REAL_ARR;F',                        (';', orbit%field(i), i = 1, 2)
nl=incr(nl); write (li(nl), ramt) 'phase;REAL_ARR;F',                        (';', orbit%phase(i), i = 1, 2)

nl=incr(nl); write (li(nl), rmt) 's;REAL;F;',                                orbit%s
nl=incr(nl); write (li(nl), rmt) 't;REAL;F;',                                orbit%t
nl=incr(nl); write (li(nl), rmt) 'charge;REAL;F;',                           orbit%charge
nl=incr(nl); write (li(nl), rmt) 'dt_ref;REAL;F;',                           orbit%dt_ref
nl=incr(nl); write (li(nl), rmt) 'p0c;REAL;F;',                              orbit%p0c
nl=incr(nl); write (li(nl), rmt) 'beta;REAL;F;',                             orbit%beta
nl=incr(nl); write (li(nl), imt) 'ix_ele;INT;F;',                            orbit%ix_ele
nl=incr(nl); write (li(nl), amt) 'state;STR;F;',                             trim(coord_state_name(orbit%state))
nl=incr(nl); write (li(nl), imt) 'direction;INT;F;',                         orbit%direction
nl=incr(nl); write (li(nl), amt) 'species;SPECIES;F;',                       trim(species_name(orbit%species))
nl=incr(nl); write (li(nl), amt) 'location;STR;F;',                          trim(location_name(orbit%location))

end subroutine orbit_out

!----------------------------------------------------------------------
! contains

subroutine coord_out(bunch, coordinate)
type (bunch_struct) :: bunch
character(40) coordinate
integer :: i_vec, n

! Allocate scratch
n = size(bunch%particle)
call reallocate_c_real_scratch(n)

! Add data
select case (coordinate)
case ('x')
  tao_c_interface_com%c_real(1:n) = bunch%particle(:)%vec(1)
case ('px')
  tao_c_interface_com%c_real(1:n) = bunch%particle(:)%vec(2)
case ('y')
  tao_c_interface_com%c_real(1:n) = bunch%particle(:)%vec(3)
case ('py')
  tao_c_interface_com%c_real(1:n) = bunch%particle(:)%vec(4)
case ('z')
  tao_c_interface_com%c_real(1:n) = bunch%particle(:)%vec(5)
case ('pz')
  tao_c_interface_com%c_real(1:n) = bunch%particle(:)%vec(6)
case ('s')
  tao_c_interface_com%c_real(1:n) = bunch%particle(:)%s
case ('t')
  tao_c_interface_com%c_real(1:n) = bunch%particle(:)%t
case ('charge')
  tao_c_interface_com%c_real(1:n) = bunch%particle(:)%charge
case ('p0c')
  tao_c_interface_com%c_real(1:n) = bunch%particle(:)%p0c
case ('state')
  call reallocate_c_integer_scratch(n)
  tao_c_interface_com%c_integer(1:n) = bunch%particle(:)%state
case ('ix_ele')
  call reallocate_c_integer_scratch(n)
  tao_c_interface_com%c_integer(1:n) = bunch%particle(:)%ix_ele
case default
  call invalid ('coordinate not "x", "px", etc. ')
  return
end select

end subroutine coord_out

!----------------------------------------------------------------------
! contains

subroutine reallocate_c_real_scratch(n)
integer :: n
if (.not. allocated(tao_c_interface_com%c_real)) allocate (tao_c_interface_com%c_real(n))
if (size(tao_c_interface_com%c_real) < n) then
  deallocate (tao_c_interface_com%c_real)
  allocate (tao_c_interface_com%c_real(n))
endif
tao_c_interface_com%n_real = n
end subroutine

subroutine reallocate_c_integer_scratch(n)
integer :: n
if (.not. allocated(tao_c_interface_com%c_integer)) allocate (tao_c_interface_com%c_integer(n))
if (size(tao_c_interface_com%c_integer) < n) then
  deallocate (tao_c_interface_com%c_integer)
  allocate (tao_c_interface_com%c_integer(n))
endif
tao_c_interface_com%n_int = n
end subroutine

!----------------------------------------------------------------------
! contains

subroutine twiss_out (twiss, prefix, suffix, emit_out, can_vary)

type (twiss_struct) twiss
character(*) prefix, suffix
character(40) fmt
character(8) v_str
logical, optional :: emit_out, can_vary

if (logic_option(.false., can_vary)) then
  v_str = ';REAL;T;'
else
  v_str = ';REAL;F;'
endif

fmt = '(4a, es22.14)'

nl=incr(nl); write (li(nl), fmt) prefix, 'beta_', suffix, v_str,                twiss%beta
nl=incr(nl); write (li(nl), fmt) prefix,  'alpha_', suffix, v_str,              twiss%alpha
nl=incr(nl); write (li(nl), fmt) prefix,  'gamma_', suffix, ';REAL;F;',         twiss%gamma
nl=incr(nl); write (li(nl), fmt) prefix,  'phi_', suffix, v_str,                twiss%phi
nl=incr(nl); write (li(nl), fmt) prefix,  'eta_', suffix, v_str,                twiss%eta
nl=incr(nl); write (li(nl), fmt) prefix,  'etap_', suffix, v_str,               twiss%etap

if (logic_option(.false., emit_out)) then
  nl=incr(nl); write (li(nl), fmt) prefix, 'sigma_', suffix, ';REAL;F;',       twiss%sigma
  nl=incr(nl); write (li(nl), fmt) prefix, 'sigma_p_', suffix, ';REAL;F;',     twiss%sigma_p
  nl=incr(nl); write (li(nl), fmt) prefix, 'emit_', suffix, ';REAL;F;',        twiss%emit
  nl=incr(nl); write (li(nl), fmt) prefix, 'norm_emit_', suffix, ';REAL;F;',   twiss%norm_emit
endif

end subroutine twiss_out

!----------------------------------------------------------------------
! contains

subroutine xy_disp_out (xy_disp, suffix, can_vary)
! Similar to twiss_out
type (xy_disp_struct) xy_disp
character(*) suffix
character(40) fmt
character(8) v_str
logical, optional ::  can_vary

if (logic_option(.false., can_vary)) then
  v_str = ';REAL;T;'
else
  v_str = ';REAL;F;'
endif

fmt = '(3a, es22.14)'

nl=incr(nl); write (li(nl), fmt) 'eta_', suffix, v_str,                           xy_disp%eta
nl=incr(nl); write (li(nl), fmt) 'etap_', suffix, v_str,                          xy_disp%etap

end subroutine xy_disp_out

!----------------------------------------------------------------------
! contains

subroutine destroy_this_data_d2 (d2_name)

type (tao_d2_data_struct), pointer :: d2_ptr
type (tao_d1_data_struct), pointer :: d1_ptr
type (tao_d2_data_array_struct), allocatable :: d2_array(:)
type (tao_d1_data_array_struct), allocatable :: d1_array(:)
type (tao_universe_struct), pointer :: u

integer i, j, ix_d2, i1, i2, n1, n_delta
logical err

character(*) d2_name

!

call tao_find_data (err, d2_name, d2_array = d2_array)
if (err .or. .not. allocated(d2_array)) then
  call invalid ('Not a valid d2 data name')
  return
endif

d2_ptr => d2_array(1)%d2
u => s%u(d2_ptr%ix_universe)
ix_d2 = d2_ptr%ix_d2_data

d1_ptr => d2_ptr%d1(1)
i1 = lbound(d1_ptr%d, 1)
i1 = d1_ptr%d(i1)%ix_data

n1 = size(d2_ptr%d1)
d1_ptr => d2_ptr%d1(n1)
i2 = ubound(d1_ptr%d, 1)
i2 = d1_ptr%d(i2)%ix_data

n_delta = i2 + 1 - i1

! Squeeze u%d2_data and u%data arrays

do i = ix_d2, u%n_d2_data_used - 1
  u%d2_data(i) = u%d2_data(i+1)
  u%d2_data(i)%ix_d2_data = i
  do j = 1, size(u%d2_data(i)%d1)
    d1_ptr => u%d2_data(i)%d1(j)
    d1_ptr%d2 => u%d2_data(i)
    i1 = d1_ptr%d(lbound(d1_ptr%d,1))%ix_data - n_delta
    i2 = d1_ptr%d(ubound(d1_ptr%d,1))%ix_data - n_delta
    u%data(i1:i2) = u%data(i1+n_delta:i2+n_delta)
    call tao_point_d1_to_data(d1_ptr, u%data(i1:i2), u%data(i1)%ix_d1)
    do k = i1, i2
      u%data(k)%ix_data = k
    enddo
  enddo
enddo

u%n_d2_data_used = u%n_d2_data_used - 1
u%n_data_used = u%n_data_used - n_delta

end subroutine destroy_this_data_d2

!----------------------------------------------------------------------
! contains

subroutine destroy_this_var_v1(v1_name)

type (tao_v1_var_array_struct), allocatable, target :: v1_array(:)
type (tao_v1_var_struct), pointer :: v1_ptr

integer j, k, i1, i2, n_delta, n

character(*) v1_name

logical err

!

call tao_find_var (err, v1_name, v1_array = v1_array)
if (err .or. .not. allocated(v1_array)) then
  call invalid ('Not a valid v1 var name')
  return
endif

v1_ptr => v1_array(1)%v1
i1 = lbound(v1_ptr%v, 1)
i1 = v1_ptr%v(i1)%ix_var
i2 = ubound(v1_ptr%v, 1)
i2 = v1_ptr%v(i2)%ix_var

n_delta = i2 + 1 - i1

n = s%n_var_used

do j = v1_ptr%ix_v1_var, s%n_v1_var_used - 1
  s%v1_var(j) = s%v1_var(j+1)
  v1_ptr => s%v1_var(j)
  v1_ptr%ix_v1_var = j
  i1 = v1_ptr%v(lbound(v1_ptr%v,1))%ix_var
  i2 = v1_ptr%v(ubound(v1_ptr%v,1))%ix_var
  s%var(i1-n_delta:i2-n_delta) = s%var(i1:i2)
  call tao_point_v1_to_var(s%v1_var(j), s%var(i1-n_delta:i2-n_delta), s%var(i1-n_delta)%ix_v1)
  do k = i1-n_delta, i2-n_delta
    s%var(k)%ix_var = k
  enddo
enddo

s%n_v1_var_used = s%n_v1_var_used - 1
s%n_var_used = s%n_var_used - n_delta

end subroutine destroy_this_var_v1

!----------------------------------------------------------------------
! contains

function match_ele_name (match_str, ele, err) result (is_a_match)

type (ele_struct) ele

integer ix, key

character(*) match_str
character(60) string

logical err, is_a_match

!

is_a_match = .false.
err = .false.

if (match_str == '*') then
  is_a_match = .true.
  return
endif

! key::name construct

string = match_str
ix = index(string, '::')
if (ix /= 0) then
  if (string(:ix-1) /= "*") then
    key = key_name_to_key_index (string(:ix-1), .true.)
    if (key < 1) then
      call invalid ('BAD ELEMENT KEY: ' // string(:ix-1))
      err = .true.
      return
    endif

    if (ele%key /= key) then
      is_a_match = .false.
      return
    endif
  endif

  string = string(ix+2:)
endif

!

if (index(string, "*") /= 0 .or. index(string, "%") /= 0) then
  is_a_match = match_wild(ele%name, string)
else
  is_a_match = (ele%name == string)
endif

end function match_ele_name

!----------------------------------------------------------------------
! contains

subroutine invalid (why_invalid, err)

character(*) why_invalid
logical, optional :: err

nl=incr(nl); li(nl) = 'INVALID'
call out_io (s_error$, r_name, '"pipe ' // trim(input_str) // '": ' // why_invalid)
call end_stuff(li, nl)
if (present(err)) err = .true.

end subroutine invalid

!----------------------------------------------------------------------
! contains

function re_str(r, n_signif) result (str)

real(rp) r
integer n_signif
character(:), allocatable :: str
character(40) string

string = real_to_string(r, 20, n_signif = n_signif)
allocate (character(len_trim(adjustl(string))):: str)
str = trim(adjustl(string))

end function re_str

!----------------------------------------------------------------------
! contains

function real_part_str(z) result (str)

complex(rp) z
character(23) str

write (str, '(a, es22.14)') ';', real(z)

end function real_part_str

!----------------------------------------------------------------------
! contains

function cmplx_str(z) result (str)

complex(rp) z
character(46) str

write (str, '(2(a, es22.14))') ';', real(z), ';', aimag(z)

end function cmplx_str

!----------------------------------------------------------------------
! contains

subroutine split_this_line (line, array, target_num, err, actual_num, space_sep)

character(*) line
character(len(line)) str
character(*) :: array(:)

integer target_num
integer, optional :: actual_num
integer i, ix
logical err, space
logical, optional :: space_sep

! For input, "^^" is used as the separator instead of ";" since the Tao code that
! calls pipe_cmd will interpret ";" as a command separator and will thus mangle
! the input_str argument.

str = line
err = .true.
array = ''
space = logic_option(.false., space_sep)

do i = 1, 1000
  if (i > size(array)) then
    call invalid('LINE SPLITTING ARRAY OVERFLOW.')
    return
  endif

  if (space) then
    call string_trim(str, str, ix)
    if (ix == 0) exit
    array(i) = str(1:ix)
    str = str(ix+1:)
  else
    ix = index(str, '^^')
    if (ix == 0) then
      array(i) = str
      exit
    endif
    array(i) = str(1:ix-1)
    str = str(ix+2:)
  endif
enddo

if (space) i = i - 1
if (present(actual_num)) actual_num = i

err = (target_num > 0 .and. i /= target_num)
if (err) then
  call invalid('NUMBER OF COMPONENTS ON LINE NOT CORRECT.')
endif

end subroutine split_this_line

!----------------------------------------------------------------------
! contains

function str_val (str_in, str_dflt) result (str_out)

character(*) str_in, str_dflt
character(max(len(str_in), len(str_dflt))) str_out

if (str_in == '') then
  str_out = str_dflt
else
  str_out = str_in
endif

end function str_val

!----------------------------------------------------------------------
! contains

function logic_val (str_in, logic_dflt, err) result (logic_out)

character(*) str_in
logical logic_dflt, err, logic_out
integer ios

err = .false.

if (str_in == '') then
  logic_out = logic_dflt
  return
endif

read (str_in, *, iostat = ios) logic_out
err = (ios /= 0)
if (err) then
  call invalid ('Not a logical: ' // str_in)
endif

end function logic_val

!----------------------------------------------------------------------
! contains

function real_val (str_in, real_dflt, err) result (real_out)

character(*) str_in
real(rp) real_dflt, real_out
integer ios
logical err

err = .false.

if (str_in == '') then
  real_out = real_dflt
  return
endif

read (str_in, *, iostat = ios) real_out
err = (ios /= 0)
if (err) then
  call invalid ('Not a real: ' // str_in)
endif

end function real_val

!----------------------------------------------------------------------
! contains

function int_val (str_in, int_dflt, err) result (int_out)

character(*) str_in
integer int_dflt, int_out
integer ios
logical err

err = .false.

if (str_in == '') then
  int_out = int_dflt
  return
endif

read (str_in, *, iostat = ios) int_out
err = (ios /= 0)
if (err) then
  call invalid ('Not an integer: ' // str_in)
endif

end function int_val

!----------------------------------------------------------------------
! contains

function ele_param_value(name, ele, orbit, data_type, err) result (value)

type (ele_struct) ele
type (coord_struct) orbit
real(rp) value
integer data_type, ix
logical err
character(*) name
character(40) attrib_name

!

err = .true.
data_type = is_real$

select case (name)
case ('orbit.floor.x', 'orbit.floor.y', 'orbit.floor.z')
  floor%r = [orbit%vec(1), orbit%vec(3), ele%value(l$)]
  floor1 = coords_local_curvilinear_to_floor (floor, ele, .false.)
  select case (name)
  case ('orbit.floor.x')
    value = floor1%r(1)
  case ('orbit.floor.y')
    value = floor1%r(2)
  case ('orbit.floor.z')
    value = floor1%r(3)
  end select
case ('orbit.spin.1')
  value = orbit%spin(1)
case ('orbit.spin.2')
  value = orbit%spin(2)
case ('orbit.spin.3')
  value = orbit%spin(3)
case ('orbit.vec.1')
  value = orbit%vec(1)
case ('orbit.vec.2')
  value = orbit%vec(2)
case ('orbit.vec.3')
  value = orbit%vec(3)
case ('orbit.vec.4')
  value = orbit%vec(4)
case ('orbit.vec.5')
  value = orbit%vec(5)
case ('orbit.vec.6')
  value = orbit%vec(6)
case ('orbit.t')
  value = orbit%t
case ('orbit.beta')
  value = orbit%beta
case ('orbit.state')
  value = orbit%state
  data_type = is_integer$
case ('orbit.pc')  
  value = (1 + orbit%vec(6)) * orbit%p0c
case ('orbit.energy', 'orbit.e_tot')  ! orbit.e_tot is old style
  call convert_pc_to ((1 + orbit%vec(6)) * orbit%p0c, orbit%species, E_tot = value)
case ('ele.ix_ele')
  value = ele%ix_ele
  data_type = is_integer$        
case ('ele.ix_branch')
  value = ele%ix_branch
case ('ele.a.beta')
  value = ele%a%beta
case ('ele.a.alpha')
  value = ele%a%alpha
case ('ele.a.eta')
  value = ele%a%eta
case ('ele.a.etap')
  value = ele%a%etap
case ('ele.a.gamma')
  value = ele%a%gamma
case ('ele.a.phi')
  value = ele%a%phi
case ('ele.b.beta')
  value = ele%b%beta
case ('ele.b.alpha')
  value = ele%b%alpha
case ('ele.b.eta')
  value = ele%b%eta
case ('ele.b.etap')
  value = ele%b%etap
case ('ele.b.gamma')
  value = ele%b%gamma
case ('ele.b.phi')
  value = ele%b%phi
case ('ele.e_tot')
  value = ele%value(e_tot$)
case ('ele.p0c')
  value = ele%value(p0c$)
case ('ele.ref_time')
  value = ele%ref_time
case ('ele.ref_time_start')
  value = ele%value(ref_time_start$)
case ('ele.x.eta')
  value = ele%x%eta
case ('ele.x.etap')
  value = ele%x%etap
case ('ele.y.eta')
  value = ele%y%eta
case ('ele.y.etap')
  value = ele%y%etap
case ('ele.s')
  value = ele%s
case ('ele.l')
  value = ele%value(l$)
case ('ele.gamma_c')
  value = ele%gamma_c
case default
  call str_upcase (attrib_name, name)
  ix = index(attrib_name, '.')

  if (attrib_name(1:ix-1) == 'ELE') then
    attrib_name = attrib_name(ix+1:)
    call pointer_to_attribute (ele, attrib_name, .true., a_ptr, err, .false.)
  else
    err = .true.
  endif

  if (err) then
    call invalid ('Bad {who}: ' // name); return
  endif

  if (associated(a_ptr%r)) then
    value = a_ptr%r
  elseif (associated(a_ptr%i)) then
    data_type = is_integer$
    value = a_ptr%i
  else
    call invalid ('{who} does not evaluate to an integer or real: ' // name); return
    err = .true.
    return
  endif
end select

err = .false.

end function ele_param_value

!----------------------------------------------------------------------
! contains

!+
! Function parse_ele_with_s_offset(line, tao_lat, ele, err) result (s_pos)
!
! Parse something like:  "{ix_branch}>>{ele}@{s_offset}".
!
! Input:
!   line      -- character(*): String to parse
!   tao_lat   -- tao_lattice_struct: Tao structure containing lattice
!
! Output:
!   line      -- character(*): String with parsed stuff removed.
!   ele       -- ele_struct, pointer: Pointer to element.
!   err       -- logical: Set True if there is an error.
!   s_pos     -- real(rp): S-position -> ele%s + offset.
!-

function parse_ele_with_s_offset(line, tao_lat, ele, err) result (s_pos)

type (tao_lattice_struct) tao_lat
type (ele_struct), pointer :: ele
type (ele_pointer_struct), allocatable :: eles(:)
type (tao_expression_info_struct), allocatable :: info(:)

real(rp) s_pos, s_offset
real(rp), allocatable :: values(:)

integer ix, ixa, ix_branch, n_loc

character(*) line
character(40) ele_name
logical err

!

err = .false.
s_pos = 0

call string_trim (line, line, ix)
if (ix == 0) then
  ele => tao_lat%lat%ele(0)
  return
endif

ixa = index(line, '->')
if (ixa == 0) then
  ele_name = line(:ix)

else
  ele_name = line(:ixa-1)
  if (ix > ixa+1) then
    call tao_evaluate_expression (line(ixa+2:ix), 1, .false., values, err, .true., info, dflt_uni = tao_lat%u%ix_uni)
    if (err) return
    s_pos = values(1)
  endif
endif

call string_trim(line(ix+1:), line, ix)


if (ele_name == '') then
  ele => tao_lat%lat%ele(0)
else
  call lat_ele_locator (ele_name, tao_lat%lat, eles, n_loc, err)

  if (err) return
  if (n_loc == 0) then
    call invalid('No element found matching: ' // ele_name, err)
    return
  elseif (n_loc > 1) then
    call invalid('Multiple elements found matching: ' // ele_name, err)
    return
  endif

  ele => eles(1)%ele
endif

s_pos = s_pos + ele%s

end function parse_ele_with_s_offset

!----------------------------------------------------------------------
! contains

subroutine bunch_params_out (bunch_params)

type (bunch_params_struct) bunch_params

!

call twiss_out(bunch_params%x, 'twiss_', 'x', emit_out = .true.)
call twiss_out(bunch_params%y, 'twiss_', 'y', emit_out = .true.)
call twiss_out(bunch_params%z, 'twiss_', 'z', emit_out = .true.)
call twiss_out(bunch_params%a, 'twiss_', 'a', emit_out = .true.)
call twiss_out(bunch_params%b, 'twiss_', 'b', emit_out = .true.)
call twiss_out(bunch_params%c, 'twiss_', 'c', emit_out = .true.)

! Sigma matrix
do i = 1, 6
  do j = 1,6
    nl=incr(nl); write (li(nl), '(a, i0, i0, a, es22.14)') 'sigma_', i, j, ';REAL;F;', bunch_params%sigma(i,j)
  enddo
enddo

! Relative min, max, centroid
do i = 1, 6
  nl=incr(nl); write (li(nl), '(a, i0, a, es22.14)') 'rel_min_', i, ';REAL;F;',      bunch_params%rel_min(i)
  nl=incr(nl); write (li(nl), '(a, i0, a, es22.14)') 'rel_max_', i, ';REAL;F;',      bunch_params%rel_max(i)
  nl=incr(nl); write (li(nl), '(a, i0, a, es22.14)') 'centroid_vec_', i, ';REAL;F;', bunch_params%centroid%vec(i)
enddo

nl=incr(nl); write (li(nl), rmt) 'centroid_t;REAL;F;',                       bunch_params%centroid%t
nl=incr(nl); write (li(nl), rmt) 'centroid_p0c;REAL;F;',                     bunch_params%centroid%p0c
nl=incr(nl); write (li(nl), rmt) 'centroid_beta;REAL;F;',                    bunch_params%centroid%beta
nl=incr(nl); write (li(nl), imt) 'ix_ele;INT;F;',                            bunch_params%centroid%ix_ele
nl=incr(nl); write (li(nl), imt) 'direction;INT;F;',                         bunch_params%centroid%direction
nl=incr(nl); write (li(nl), amt) 'species;SPECIES;F;',                       trim(species_name(bunch_params%centroid%species))
nl=incr(nl); write (li(nl), amt) 'location;ENUM;F;',                         trim(location_name(bunch_params%centroid%location))
nl=incr(nl); write (li(nl), rmt) 's;REAL;F;',                                bunch_params%s
nl=incr(nl); write (li(nl), rmt) 't;REAL;F;',                                bunch_params%t
nl=incr(nl); write (li(nl), rmt) 'sigma_t;REAL;F;',                          bunch_params%sigma_t
nl=incr(nl); write (li(nl), rmt) 'charge_live;REAL;F;',                      bunch_params%charge_live
nl=incr(nl); write (li(nl), imt) 'n_particle_tot;INT;F;',                    bunch_params%n_particle_tot
nl=incr(nl); write (li(nl), imt) 'n_particle_live;INT;F;',                   bunch_params%n_particle_live
nl=incr(nl); write (li(nl), imt) 'n_particle_lost_in_ele;INT;F;',            bunch_params%n_particle_lost_in_ele

end subroutine bunch_params_out

!----------------------------------------------------------------------
! contains

subroutine real_array_out(val_arr, use_buffer, ix0, ix1)

real(rp) val_arr(:)
integer, optional :: ix0, ix1
integer i, j, n_arr
logical use_buffer

!

n_arr = integer_option(size(val_arr), ix1) - integer_option(1, ix0) + 1

if (use_buffer) then
  call re_allocate_c_double(tao_c_interface_com%c_real, n_arr, .false.)
  tao_c_interface_com%n_real = n_arr
  tao_c_interface_com%c_real(1:n_arr) = val_arr(1:n_arr)

else  ! string_list
  do i = 1, n_arr
    j = i + integer_option(1, ix0) - 1
    nl=incr(nl); write (li(nl), '(i0, a, es22.14)') j, ';', val_arr(i)
  enddo
endif

end subroutine real_array_out

!----------------------------------------------------------------------
! contains

subroutine write_this_floor(floor, name, can_vary)

type(floor_position_struct) :: floor
character(*) name
logical can_vary
integer i, j

!

nl=incr(nl); write (li(nl), rmt2) trim(name) // ';REAL_ARR;', can_vary, (';', floor%r(i), i = 1, 3), ';', floor%theta, ';', floor%phi, ';', floor%psi
nl=incr(nl); write (li(nl), rmt2) trim(name) // '-W;REAL_ARR;', .false., ((';', floor%w(i,j), i = 1, 3), j = 1, 3)

end subroutine write_this_floor

!----------------------------------------------------------------------
! contains

recursive subroutine write_this_ele_floor(ele, loc, can_vary, suffix)

type (ele_struct), target :: ele
type (ele_struct), pointer :: ele0

integer n, ie, loc
logical can_vary
character(*) suffix

!


if (ele%lord_status /= super_lord$ .and. ele%lord_status /= girder_lord$ .and. ele%n_slave > 0) then
  do ie = 1, ele%n_slave
    ele0 => pointer_to_slave(ele, ie)
    call write_this_ele_floor(ele0, loc, can_vary, '-Slave' // int_str(ie))
  enddo
  return
endif

!

if (ele%ix_ele == 0) then
  ele0 => ele
else
  ele0 => pointer_to_next_ele(ele, -1)
endif

select case (loc)
case (1)
  call write_this_floor(ele0%floor, 'Reference' // suffix, can_vary)
  call write_this_floor(ele_geometry_with_misalignments (ele, 0.0_rp), 'Actual' // suffix, .false.)
case (2)
  call ele_geometry(ele0%floor, ele, floor, 0.5_rp)
  call write_this_floor(floor, 'Reference' // suffix, can_vary)
  call write_this_floor(ele_geometry_with_misalignments (ele, 0.5_rp), 'Actual' // suffix, .false.)
case (3)
  call write_this_floor(ele%floor, 'Reference' // suffix, can_vary)
  call write_this_floor(ele_geometry_with_misalignments (ele), 'Actual' // suffix, .false.)
end select

end subroutine write_this_ele_floor

!----------------------------------------------------------------------
! contains

subroutine this_floor_plan (iu, graph)

type (tao_graph_struct) :: graph
type (lat_struct), pointer :: lat
type (tao_ele_shape_struct), pointer :: ele_shape, ele_shape2
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele, slave

real(rp) y1, y2
integer iu, n, i, j, ix_shape_min, ix_pass, n_links
character(40) label_name

!

lat => s%u(iu)%model%lat

do n = 0, ubound(lat%branch, 1)
  branch => lat%branch(n)
  branch%ele%logic = .false.  ! Used to mark as drawn.
  do i = 0, branch%n_ele_max
    ele => branch%ele(i)
    if (ele%slave_status == super_slave$) cycle

    ix_shape_min = 1
    do
      call tao_ele_shape_info(iu, ele, s%plot_page%floor_plan%ele_shape, ele_shape, label_name, y1, y2, ix_shape_min)
      if (.not. associated(ele_shape) .and. (ele%key == overlay$ .or. &
                                             ele%key == group$ .or. ele%key == girder$)) exit   ! Nothing to draw

      if (graph%floor_plan%draw_only_first_pass .and. ele%slave_status == multipass_slave$) then
        call multipass_chain (ele, ix_pass, n_links)
        if (ix_pass > 1) exit
      endif

      if (ele%lord_status == multipass_lord$) then
        do j = 1, ele%n_slave
          if (graph%floor_plan%draw_only_first_pass .and. j > 1) exit
          slave => pointer_to_slave(ele, j)
          ele_shape2 => tao_pointer_to_ele_shape (iu, slave, s%plot_page%floor_plan%ele_shape)
          if (associated(ele_shape2)) cycle ! Already drawn. Do not draw twice
          call this_floor_plan2 (graph, slave, ele_shape, label_name, y1, y2)
        enddo
      else
        call this_floor_plan2 (graph, ele, ele_shape, label_name, y1, y2)
      endif
      if (.not. associated(ele_shape)) exit
      if (.not. ele_shape%multi) exit
    enddo

  enddo
enddo

end subroutine this_floor_plan

!----------------------------------------------------------------------
! contains

subroutine this_floor_plan2 (graph, ele, ashape, label_name, y1, y2)

type (tao_graph_struct) :: graph
type (ele_struct) ele
type (tao_ele_shape_struct), pointer :: ashape
type (ele_struct), pointer :: ele1, ele2
type (floor_position_struct) floor, floor1, floor2

real(rp) y1, y2, off1, off2
integer line_width
character(40) color, label_name, shape_shape

!

call find_element_ends (ele, ele1, ele2)
if (.not. associated(ele1)) return

if (.not. associated(ashape)) then
  color = ''
  label_name = ''
  shape_shape = ''
  line_width = 0
else
  color = ashape%color
  shape_shape = ashape%shape
  line_width = ashape%line_width
endif

off1 = y1 * s%plot_page%floor_plan_shape_scale
off2 = y2 * s%plot_page%floor_plan_shape_scale

floor%r = [0.0_rp, 0.0_rp, 0.0_rp]
floor1 = coords_local_curvilinear_to_floor (floor, ele, .true.)

floor%r = [0.0_rp, 0.0_rp, ele%value(l$)]
floor2 = coords_local_curvilinear_to_floor (floor, ele, .true.)

call tao_floor_to_screen_coords (graph, floor1, end1)
call tao_floor_to_screen_coords (graph, floor2, end2)

if (ele%key == sbend$) then
  nl=incr(nl); write (li(nl), '(2(i0, a), 2a, 6(es14.7, a), (i0, a), 2a, 2(es10.2, a), 4a, 4(es14.7, a))') &
              ele%ix_branch, ';', ele%ix_ele, ';', &
              trim(key_name(ele%key)), ';', end1%r(1), ';', end1%r(2), ';', end1%theta, ';', &
              end2%r(1), ';', end2%r(2), ';', end2%theta, ';', &
              line_width, ';', trim(shape_shape), ';', off1, ';', off2, ';', trim(color), ';', trim(label_name), ';', &
              ele%value(l$), ';', ele%value(angle$), ';', ele%value(e1$), ';', ele%value(e2$)
else
  nl=incr(nl); write (li(nl), '(2(i0, a), 2a, 6(es14.7, a), (i0, a), 2a, 2(es10.2, a), 4a)') &
              ele%ix_branch, ';', ele%ix_ele, ';', &
              trim(key_name(ele%key)), ';', end1%r(1), ';', end1%r(2), ';', end1%theta, ';', &
              end2%r(1), ';', end2%r(2), ';', end2%theta, ';', &
              line_width, ';', trim(shape_shape), ';', off1, ';', off2, ';', trim(color), ';', trim(label_name)
endif

end subroutine this_floor_plan2

end subroutine tao_pipe_cmd
