!+
! Module tao_interface
!
! Module to define the interfaces for the tao routines.
!-

module tao_interface

use tao_struct

!+
! Function tao_pointer_to_universe (...) result (u)
!
! Routine to set a pointer to a universe.
!
! This is an overloaded routine for the:
!  tao_pointer_to_universe_int (ix_uni, neg2_to_default) result (u)
!  tao_pointer_to_universe_str (string, neg2_to_default) result (u)
!
! Note: With a string argument, this routine can only handle single universe picks. 
! That is, it cannot handlle something like "[1,3,4]@...". To handle multiple universe picks, use:
!   tao_pointer_to_universes
!
! Input:
!   ix_uni          -- Integer: Index to the s%u(:) array
!                        If ix_uni is -1 -> u(s%global%default_universe) will be used.
!   string          -- character(*): String in the form "<ix_uni>@..." or, if 
!                        no "@" is present, u will point to the default universe.
!   neg2_to_default -- logical, optional: i_uni = -2 (all universes) maps to the default uni?
!                             Default if False.
!
! Output:
!   string      -- character(*): String with universe prefix stripped off.
!   u           -- Tao_universe_struct, pointer: Universe pointer.
!                     u will be nullified if there is an error and an error message will be printed.
!-

interface tao_pointer_to_universe
  module procedure tao_pointer_to_universe_int
  module procedure tao_pointer_to_universe_str
end interface

private tao_pointer_to_universe_int, tao_pointer_to_universe_str

interface

subroutine tao_abort_command_file(force_abort)
  implicit none
  logical, optional :: force_abort
end subroutine

subroutine tao_alias_cmd (alias, string)
  implicit none
  character(*) :: alias
  character(*) :: string
end subroutine

function tao_beam_emit_calc (plane, emit_type, ele, bunch_params) result (emit)
  import
  implicit none
  integer plane, emit_type
  type (ele_struct) ele
  type (bunch_params_struct) bunch_params
  real(rp) emit
end function

function tao_beam_track_endpoint (ele_id, lat, branch_str, where, u) result (ele)
  import
  implicit none
  type (ele_struct), pointer :: ele
  type (lat_struct), target :: lat
  character(*) ele_id, where, branch_str
  type (tao_universe_struct), target :: u
end function

function tao_branch_index (ix_branch) result (ix_this)
  import
  implicit none
  integer ix_branch, ix_this
end function

subroutine tao_call_cmd (file_name, cmd_arg)
  implicit none
  character(*) :: file_name
  character(*), optional :: cmd_arg(:)
end subroutine

subroutine tao_var_check (eles, attribute, silent)
  import
  implicit none
  type (ele_pointer_struct), allocatable :: eles(:)
  logical silent
  character(*) attribute
end subroutine

function tao_chrom_calc_needed (data_type, data_source) result (do_chrom)
  import
  implicit none
  character(*) data_type, data_source
  logical do_chrom
end function

subroutine tao_clear_cmd (cmd_line)
  implicit none
  character(*) cmd_line
end subroutine

subroutine tao_clip_cmd (gang, where, value1, value2)
  import
  implicit none
  logical gang
  character(*) :: where
  real(rp) value1, value2
end subroutine

subroutine tao_close_command_file()
end subroutine

subroutine tao_command (command_line, err, err_is_fatal)
  implicit none
  character(*) :: command_line
  logical err, err_is_fatal
end subroutine

subroutine tao_control_tree_list (ele, tree)
  import
  implicit none
  type (ele_struct) ele
  type (ele_pointer_struct), allocatable :: tree(:)
end subroutine

function tao_constraint_type_name(datum) result (datum_name)
  import
  implicit none
  type (tao_data_struct) datum
  character(400) datum_name
end function

subroutine tao_count_strings (string, pattern, num)
  import
  implicit none
  character(*) string, pattern
  integer num
end subroutine

function tao_curve_ele_ref (curve, point_to_ele_ref) result (ele_track)
  import
  implicit none
  type (tao_curve_struct) curve
  type (ele_struct), pointer :: ele_track
  logical point_to_ele_ref
end function

function tao_curve_ix_uni (curve) result (ix_uni)
  import
  implicit none
  type (tao_curve_struct) curve
  integer ix_uni
end function

function tao_curve_name(curve, use_region) result (curve_name)
  import
  implicit none
  type (tao_curve_struct) curve
  character(60) curve_name
  logical, optional :: use_region
end function

subroutine tao_curve_rms_calc (curve, who, rms, mean)
  import
  implicit none
  type (tao_curve_struct) curve
  real(rp) rms, mean, ys, dx
  integer i, n
  character(*) who
end subroutine

function tao_d2_d1_name(d1, show_universe) result (d2_d1_name)
  import
  implicit none
  type (tao_d1_data_struct) d1
  character(60) d2_d1_name
  logical, optional :: show_universe
end function

subroutine tao_data_check (err)
  import
  implicit none
  logical err
end subroutine

function tao_data_sanity_check (datum, print_err, default_data_type, uni) result (is_valid)
  import
  type (tao_data_struct) datum
  type (tao_universe_struct), optional, target :: uni
  logical print_err, is_valid
  character(*) default_data_type
end function

subroutine tao_data_show_use (d2_data, lines, nl)
  import
  implicit none
  type (tao_d2_data_struct), intent(in), target :: d2_data
  character(*), optional, allocatable :: lines(:)
  integer, optional :: nl
end subroutine

function tao_datum_has_associated_ele (data_type, branch_geometry) result (has_associated_ele)
  implicit none
  character(*) data_type
  integer has_associated_ele
  integer, optional :: branch_geometry
end function

function tao_datum_name(datum, show_universe) result (datum_name)
  import
  implicit none
  type (tao_data_struct) datum
  character(60) datum_name
  logical, optional :: show_universe
end function

subroutine tao_de_optimizer (abort)
  implicit none
  logical abort
end subroutine

subroutine tao_ele_shape_info (ix_uni, ele, ele_shapes, e_shape, label_name, y1, y2, ix_shape_min)
  import
  implicit none
  type (ele_struct) ele
  type (tao_ele_shape_struct) ele_shapes(:)
  type (tao_ele_shape_struct), pointer :: e_shape
  real(rp) y1, y2
  integer ix_uni
  integer, optional :: ix_shape_min
  character(*) label_name
end subroutine

recursive subroutine tao_evaluate_a_datum (datum, u, tao_lat, datum_value, valid_value, &
                                                            why_invalid, called_from_lat_calc)
  import
  implicit none
  type (tao_data_struct) datum
  type (tao_universe_struct), target :: u
  type (tao_lattice_struct), target :: tao_lat
  real(rp) datum_value
  logical valid_value
  logical, optional :: called_from_lat_calc
  character(*), optional :: why_invalid
end subroutine

recursive &
subroutine tao_evaluate_expression (expression, n_size, use_good_user, value, err_flag, print_err, &
                      info, stack, dflt_component, dflt_source, dflt_ele_ref, dflt_ele_start, dflt_ele, &
                      dflt_dat_or_var_index, dflt_uni, dflt_eval_point, dflt_s_offset, dflt_orbit, datum)
  import
  implicit none
  character(*) :: expression
  character(*), optional :: dflt_component, dflt_source
  character(*), optional :: dflt_dat_or_var_index
  type (tao_eval_stack1_struct), allocatable, optional :: stack(:)
  type (ele_struct), optional, pointer :: dflt_ele_ref, dflt_ele_start, dflt_ele
  type (coord_struct), optional :: dflt_orbit
  type (tao_expression_info_struct), allocatable, optional :: info(:)
  type (tao_data_struct), optional :: datum
  real(rp), allocatable :: value(:)
  real(rp), optional :: dflt_s_offset
  integer n_size
  integer, optional :: dflt_uni, dflt_eval_point
  logical use_good_user, err_flag
  logical, optional :: print_err
end subroutine

function tao_evaluate_tune (q_str, q0, delta_input) result (q_val)
  import
  implicit none
  real(rp) q0, q_val
  real(rp), allocatable :: set_val(:)
  character(*) q_str
  logical delta_input
end function

subroutine tao_evaluate_element_parameters (err, param_name, values, print_err, dflt_ele, &
                                              dflt_source, dflt_component, dflt_uni, eval_point, info)
  import
  implicit none
  type (ele_struct), pointer, optional :: dflt_ele
type (tao_expression_info_struct), allocatable, optional :: info(:)
  character(*) param_name
  character(*) dflt_source
  character(*), optional :: dflt_component
  real(rp), allocatable :: values(:)
  integer, optional :: dflt_uni, eval_point
  logical err
  logical :: print_err
end subroutine


subroutine tao_find_data (err, data_name, d2_array, d1_array, d_array, re_array, &
                           log_array, str_array, int_array, ix_uni, dflt_index, print_err, component)
  import
  implicit none
  type (tao_d2_data_array_struct), allocatable, optional :: d2_array(:)
  type (tao_d1_data_array_struct), allocatable, optional :: d1_array(:)
  type (tao_data_array_struct), allocatable, optional    :: d_array(:)
  type (tao_real_pointer_struct), allocatable, optional    :: re_array(:)
  type (tao_integer_array_struct), allocatable, optional :: int_array(:)
  type (tao_logical_array_struct), allocatable, optional :: log_array(:)
  type (tao_string_array_struct), allocatable, optional  :: str_array(:)
  character(*) :: data_name
  character(*), optional :: component
  character(*), optional :: dflt_index
  integer, optional :: ix_uni
  logical err
  logical, optional :: print_err
end subroutine


subroutine tao_find_var (err, var_name, v1_array, v_array, re_array, log_array, &
                                               str_array, print_err, component, dflt_var_index)
  import
  implicit none
  type (tao_v1_var_array_struct), allocatable, optional  :: v1_array(:)
  type (tao_var_array_struct), allocatable, optional     :: v_array(:)
  type (tao_real_pointer_struct), allocatable, optional  :: re_array(:)
  type (tao_logical_array_struct), allocatable, optional :: log_array(:)
  type (tao_string_array_struct), allocatable, optional  :: str_array(:)
  character(*) :: var_name
  character(*), optional :: component, dflt_var_index
  logical, optional :: print_err
  logical err, print_error
end subroutine

subroutine tao_find_plot_region (err, where, region, print_flag)
  import
  implicit none
  type (tao_plot_region_struct), pointer :: region
  character(*) where
  logical, optional :: print_flag
  logical err
end subroutine

subroutine tao_find_plots (err, name, where, plot, graph, curve, print_flag, blank_means_all, only_visible)
  import
  implicit none
  type (tao_plot_array_struct), allocatable, optional :: plot(:)
  type (tao_graph_array_struct), allocatable, optional :: graph(:)
  type (tao_curve_array_struct), allocatable, optional :: curve(:)
  character(*) name, where
  logical, optional :: print_flag, blank_means_all, only_visible
  logical err
end subroutine
 
subroutine tao_floor_to_screen (graph, r_floor, x_screen, y_screen)
  import
  implicit none
  type (tao_graph_struct) graph
  real(rp) r_floor(3), x_screen, y_screen 
end subroutine

subroutine tao_floor_to_screen_coords (graph, floor, screen)
  import
  implicit none
  type (tao_graph_struct) graph
  type (floor_position_struct) floor, screen
end subroutine

subroutine tao_get_opt_vars (var_value, var_step, var_delta, var_weight, var_ix, &
                                    ignore_if_weight_is_zero, ignore_if_not_limited)
  import
  implicit none
  real(rp), allocatable, optional :: var_value(:), var_delta(:)
  real(rp), allocatable, optional :: var_step(:), var_weight(:)
  integer, allocatable, optional :: var_ix(:)
  logical, optional :: ignore_if_weight_is_zero, ignore_if_not_limited
  logical ignore_weight_is_zero, ignore_not_limited
end subroutine

function tao_graph_name(graph, use_region) result (graph_name)
  import
  implicit none
  type (tao_graph_struct) graph
  character(60) graph_name
  logical, optional :: use_region
end function

subroutine tao_has_been_created ()
end subroutine
 
subroutine tao_help (what1, what2, lines, n_lines)
  implicit none
  character(*) what1, what2
  character(*), optional, allocatable :: lines(:)
  integer, optional :: n_lines
end subroutine

subroutine tao_init (err_flag)
  implicit none
  logical :: err_flag
end subroutine

subroutine tao_init_find_elements (u, search_string, eles, attribute, found_one)
  import
  implicit none
  type (tao_universe_struct), target :: u
  type (ele_pointer_struct), allocatable :: eles(:)
  character(*) search_string
  character(*), optional :: attribute
  logical, optional :: found_one
end subroutine

subroutine tao_init_lattice (lat_file, err_flag)
  implicit none
  character(*) lat_file
  logical err_flag
end subroutine

subroutine tao_init_plotting (plot_file)
  implicit none
  character(*) plot_file
end subroutine

subroutine tao_init_single_mode (single_mode_file)
  implicit none
  character(*) single_mode_file
end subroutine

function tao_is_valid_name (name, why_invalid) result (is_valid)
  implicit none
  character(*) name, why_invalid
  logical is_valid
end function

subroutine tao_json_cmd (input_str)
  import
  implicit none
  character(*) input_str
end subroutine

subroutine tao_key_info_to_str (ix_key, ix_min_key, ix_max_key, key_str, header_str)
  import
  implicit none
  integer ix_key, ix_min_key, ix_max_key
  character(*) key_str
  character(*) header_str
end subroutine

subroutine tao_lat_bookkeeper (u, tao_lat, err_flag)
  import
  implicit none
  type (tao_universe_struct), target :: u
  type (tao_lattice_struct) :: tao_lat
  logical err_flag
end subroutine

function tao_lat_emit_calc (plane, emit_type, ele, modes) result (emit)
  import
  implicit none
  integer plane, emit_type
  type (ele_struct) ele
  type (normal_modes_struct) modes
  real(rp) emit
end function

function tao_lat_sigma_calc_needed (data_type, data_source) result (do_lat_sigma)
  import
  implicit none
  character(*) data_type, data_source
  logical do_lat_sigma
end function
 
subroutine tao_lattice_calc (calc_ok, print_err)
  implicit none
  logical calc_ok
  logical, optional :: print_err
end subroutine

subroutine tao_limit_calc (limited)
  implicit none
  logical limited
end subroutine

subroutine tao_lmdif_optimizer (abort)
  implicit none
  logical abort
end subroutine

subroutine tao_locate_all_elements (ele_list, eles, err, ignore_blank)
  import
  implicit none
  type (ele_pointer_struct), allocatable :: eles(:)
  character(*) ele_list
  logical err
  logical, optional :: ignore_blank
end subroutine

subroutine tao_locate_elements (ele_list, ix_universe, eles, err, lat_type, ignore_blank, &
                                       err_stat_level, above_ubound_is_err, ix_branch, multiple_eles_is_err)
  import
  implicit none
  character(*) ele_list
  integer ix_universe
  type (ele_pointer_struct), allocatable :: eles(:)
  logical err
  integer, optional :: lat_type, err_stat_level, ix_branch
  logical, optional :: ignore_blank, above_ubound_is_err, multiple_eles_is_err
end subroutine


subroutine tao_mark_lattice_ele (lat)
  import
  implicit none
  type (lat_struct), target :: lat
end subroutine

function tao_merit (calc_ok) result (this_merit)
  import
  implicit none
  real(rp) this_merit
  logical, optional :: calc_ok
end function

function tao_one_turn_map_calc_needed (data_type, data_source) result (do_one_turn_map)
  import
  implicit none
  character(*) data_type, data_source
  logical do_one_turn_map
end function

subroutine tao_open_file (file, iunit, file_name, error_severity, binary)
  implicit none
  character(*) file
  character(*) file_name
  integer iunit, error_severity
  logical, optional :: binary
end subroutine

Function tao_open_scratch_file (err) result (iu)
  implicit none
  integer iu
  logical err
end function

function tao_optimization_status (datum) result (why_str)
  import
  implicit none
  type (tao_data_struct) :: datum
  character(60) why_str
end function

function tao_oreint_building_wall_pt(pt_in) result (pt_out)
  import
  implicit none
  type (tao_building_wall_point_struct) pt_in, pt_out
end function
 
function tao_pointer_to_datum (d1, ele_name) result (datum_ptr)
  import
  implicit none
  type (tao_d1_data_struct), target :: d1
  type (tao_data_struct), pointer :: datum_ptr
  character(*) ele_name
end function

subroutine tao_pointer_to_branches (branch_str, branches, unis, err)
  import
  implicit none
  character(*) branch_str
  type (branch_pointer_struct), allocatable :: branches(:)
  type (tao_universe_pointer_struct), allocatable, target :: unis(:)
  logical err
end subroutine

subroutine tao_pointer_to_universes (name_in, unis, err, name_out, explicit_uni, dflt_uni)
  import
  implicit none
  type (tao_universe_pointer_struct), allocatable :: unis(:)
  character(*) name_in
  character(*), optional :: name_out
  integer, optional :: dflt_uni
  logical err
  logical, optional :: explicit_uni
end subroutine

function tao_param_value_at_s (dat_name, ele_to_s, ele_here, orbit, err_flag, why_invalid, print_err, bad_datum) result (value)
  import
  implicit none
  type (ele_struct) ele_to_s, ele_here
  type (coord_struct) orbit
  real(rp) value
  character(*) dat_name
  character(*), optional :: why_invalid
  logical err_flag
  logical, optional :: print_err, bad_datum
end function

subroutine tao_parse_command_args (error, cmd_line)
  import
  implicit none
  character(*), optional :: cmd_line
  logical error
end subroutine

subroutine tao_parse_element_param_str (err, in_str, uni, element, parameter, where, component)
  import
  implicit none
  character(*) in_str, uni, element, parameter, component
  integer where
  logical err
end subroutine

subroutine tao_pause_cmd (time)
  import
  implicit none
  real(rp) time
end subroutine

subroutine tao_pick_universe (name_in, name_out, picked, err, ix_uni, explicit_uni, dflt_uni, pure_uni)
  import
  implicit none
  character(*) name_in, name_out
  integer, optional :: ix_uni, dflt_uni
  logical, allocatable :: picked(:)
  logical err
  logical, optional :: explicit_uni, pure_uni
end subroutine
 
subroutine tao_pipe_cmd (input_str)
  import
  implicit none
  character(*) input_str
end subroutine

subroutine tao_place_cmd (where, who, no_buffer)
  implicit none
  character(*) who
  character(*) where
  logical, optional :: no_buffer
end subroutine
 
subroutine tao_plot_cmd (where, component)
  implicit none
  character(*) :: where
  character(*) :: component
end subroutine

subroutine tao_plot_setup ()
  implicit none
end subroutine

subroutine tao_plot_struct_transfer (plot_in, plot_out)
  import
  implicit none
  type (tao_plot_struct), target :: plot_in, plot_out
end subroutine

function tao_pointer_to_building_wall_shape (wall_name) result (e_shape)
  import
  implicit none
  type (tao_ele_shape_struct), pointer :: e_shape
  character(*) wall_name
end function

function tao_pointer_to_ele_shape (ix_uni, ele, ele_shape, dat_var_name, dat_var_value, ix_shape_min) result (e_shape)
  import
  implicit none
  integer ix_uni
  type (ele_struct), target :: ele
  type (tao_ele_shape_struct), target :: ele_shape(:)
  character(*), optional :: dat_var_name
  real(rp), optional :: dat_var_value
  integer, optional :: ix_shape_min
  type (tao_ele_shape_struct), pointer :: e_shape
end function

function tao_pointer_to_tao_lat (u, lat_type) result (tao_lat)
  import
  implicit none
  type (tao_universe_struct), target :: u
  type (tao_lattice_struct), pointer :: tao_lat
  integer, optional :: lat_type
end function

subroutine tao_print_command_line_info
  import
  implicit none
end subroutine

subroutine tao_ptc_normal_form (do_calc, tao_lat, ix_branch, rf_on)
  import
  type (tao_lattice_struct), target :: tao_lat
  integer ix_branch
  logical do_calc
  integer, optional :: rf_on
end subroutine

subroutine tao_python_cmd (input_str)
  import
  implicit none
  character(*) input_str
end subroutine

function tao_rad_int_calc_needed (data_type, data_source) result (do_rad_int)
  import
  implicit none
  character(*) data_type, data_source
  logical do_rad_int
end function

subroutine tao_re_allocate_expression_info (info, n, exact)
  import
  implicit none
  type (tao_expression_info_struct), allocatable :: info(:)
  integer, intent(in) :: n
  logical, optional :: exact
end subroutine

subroutine tao_regression_test ()
  implicit none
end subroutine

subroutine tao_remove_blank_characters (str)
  implicit none
  character(*) str
end subroutine

subroutine tao_read_cmd (which, unis, file, silent)
  implicit none
  character(*) which, unis, file
  logical silent
end subroutine

function tao_read_phase_space_index (name, ixc, print_err) result (ix_ps)
  import
  implicit none
  character(*) name
  integer ix_ps, ixc
  logical, optional :: print_err
end function
 
subroutine tao_run_cmd (which, abort)
  implicit none
  character(*) which
  logical abort
end subroutine

subroutine tao_scale_ping_data (u)
  import
  implicit none
  type (tao_universe_struct) u
end subroutine

subroutine tao_set_data_useit_opt (data)
  import
  implicit none
  type (tao_data_struct), optional :: data(:)
end subroutine

subroutine tao_set_flags_for_changed_attribute (u, ele_name, ele_ptr, val_ptr, who)
  import
  implicit none
  type (tao_universe_struct), target :: u
  type (ele_struct), pointer, optional :: ele_ptr
  type (all_pointer_struct), optional :: val_ptr
  character(*) ele_name
  character(*), optional :: who
end subroutine

subroutine tao_set_invalid (datum, message, why_invalid, exterminate, err_level)
  import
  implicit none
  type (tao_data_struct) datum
  logical, optional :: exterminate
  integer, optional :: err_level
  character(*) message
  character(*), optional :: why_invalid
end subroutine

subroutine tao_set_var_model_value (var, value, print_limit_warning)
  import
  implicit none
  type (tao_var_struct), target :: var
  real(rp) value
  logical, optional :: print_limit_warning
end subroutine

subroutine tao_set_var_useit_opt ()
end subroutine

subroutine tao_set_opt_vars (var_vec, print_limit_warning)
  import
  implicit none
  real(rp) var_vec(:)
  logical, optional :: print_limit_warning
end subroutine

subroutine tao_setup_key_table ()
  import
  implicit none
end subroutine

function tao_srdt_calc_needed (data_type, data_source) result (do_srdt)
  import
  implicit none
  character(*) data_type, data_source
  integer do_srdt
end function

subroutine tao_symbol_import_from_lat(lat)
  import
  implicit none
  type (lat_struct) lat
end subroutine

subroutine tao_quiet_set (set)
  import
  implicit none
  character(*) set
end subroutine

subroutine tao_single_mode (char)
  implicit none
  character(1) :: char
end subroutine

subroutine tao_split_component (comp_str, comp, err)
  import
  implicit none
  character(*) comp_str
  type (tao_data_var_component_struct), allocatable :: comp(:)
  logical err
end subroutine

subroutine tao_spin_matrix_calc (datum, u, ele_ref, ele, excite_zero)
  import
  implicit none
  type (tao_data_struct), target :: datum
  type (tao_universe_struct), target :: u
  type (ele_struct), pointer :: ele_ref, ele
  character(*), optional :: excite_zero(3)
end subroutine

subroutine tao_spin_polarization_calc (branch, tao_branch, excite_zero, ignore_kinetic, err_flag)
  import
  implicit none
  type (branch_struct), target :: branch
  type (tao_lattice_branch_struct), target :: tao_branch
  character(*), optional :: excite_zero(3), ignore_kinetic
  logical, optional :: err_flag
end subroutine

function tao_spin_matrices_calc_needed (data_type, data_source) result (do_calc)
  import
  implicit none
  character(*) data_type, data_source
  logical do_calc
end function

subroutine tao_spin_tracking_turn_on()
end subroutine

subroutine tao_shape_init (shape, err, print_err)
  import
  implicit none
  type (tao_ele_shape_struct) shape
  logical err
  logical, optional :: print_err
end subroutine

subroutine tao_show_cmd (what)
  implicit none
  character(*) what
end subroutine

subroutine tao_show_this (what, result_id, lines, nl)
  implicit none
  character(*) :: what
  character(*) result_id
  character(*), allocatable :: lines(:)
  integer nl
end subroutine

function tao_subin_uni_number (name_in, ix_uni, name_out) result (ok)
  import
  implicit none
  character(*) name_in, name_out
  integer ix_uni
  logical ok
end function

subroutine tao_taper_cmd(except, uni_names)
  import
  implicit none
  character(*) except, uni_names
end subroutine

subroutine tao_top_level (command, errcode)
  implicit none
  character(*), optional :: command
  integer, optional :: errcode
end subroutine
 
subroutine tao_turn_on_special_calcs_if_needed_for_plotting ()
  import
  implicit none
end subroutine

function tao_universe_index (i_uni, neg2_to_default) result (i_this_uni)
  import
  implicit none
  integer i_uni, i_this_uni
  logical, optional :: neg2_to_default
end function

subroutine tao_use_data (action, data_name)
  implicit none
  character(*) :: action
  character(*) :: data_name
end subroutine

subroutine tao_use_var (action, var_name)
  implicit none
  character(*) :: action
  character(*) :: var_name
end subroutine

function tao_user_is_terminating_optimization () result (is_terminating)
  implicit none
  logical is_terminating
end function

function tao_var1_name(var) result (var1_name)
  import
  implicit none
  type (tao_var_struct) var
  character(60) var1_name
end function

function tao_var_attrib_name(var) result (var_attrib_name)
  import
  implicit none
  type (tao_var_struct) var
  character(60) var_attrib_name
end function
 
subroutine tao_var_repoint ()
end subroutine

subroutine tao_var_show_use (v1_var, lines, nl)
  import
  implicit none
  type (tao_v1_var_struct), intent(in) :: v1_var
  character(*), optional, allocatable :: lines(:)
  integer, optional :: nl
end subroutine

subroutine tao_var_target_calc ()
  import
  implicit none
end subroutine

subroutine tao_var_useit_plot_calc (graph, var)
  import
  implicit none
  type (tao_graph_struct) graph
  type (tao_var_struct) var(:)
end subroutine

subroutine tao_write_cmd (what)
  implicit none
  character(*) :: what
end subroutine
 
subroutine tao_x_axis_cmd (where, what)
  implicit none
  character(*) where
  character(*) what
end subroutine

end interface

!-----------------------------------------------------------------------

abstract interface

subroutine tao_hook_branch_calc_def (u, tao_lat, branch)
  import
  implicit none
  type (tao_universe_struct), target :: u
  type (tao_lattice_struct), target :: tao_lat
  type (branch_struct), target :: branch
end subroutine
 
subroutine tao_hook_command_def (command_line, found)
  implicit none
  character(*) command_line
  logical found
end subroutine
 
function tao_hook_curve_s_pt_def (s_default, ix_now, x1, x2, n_pts, tao_lat, curve) result (s_pt)
  import
  implicit none
  type (tao_curve_struct) curve
  type (tao_lattice_struct) tao_lat
  real(rp) s_default, x1, x2, s_pt
  integer ix_now, n_pts
end function

function tao_hook_data_sanity_check_def (found, datum, print_err, default_data_type, uni) result (is_valid)
  import
  implicit none
  type (tao_data_struct) datum
  type (tao_universe_struct), optional, target :: uni
  logical found, print_err, is_valid
  character(*) default_data_type
end function

subroutine tao_hook_draw_floor_plan_def (plot, graph)
  import
  implicit none
  type (tao_plot_struct) plot
  type (tao_graph_struct) graph
end subroutine
 
subroutine tao_hook_draw_graph_def (plot, graph, found)
  import
  implicit none
  type (tao_plot_struct) plot
  type (tao_graph_struct) graph
  logical found
end subroutine

subroutine tao_hook_evaluate_a_datum_def (found, datum, u, tao_lat, datum_value, valid_value, why_invalid)
  import
  implicit none
  type (tao_data_struct) datum
  type (tao_universe_struct), target :: u
  type (tao_lattice_struct), target :: tao_lat
  real(rp) datum_value
  logical found, valid_value
  character(*), optional :: why_invalid
end subroutine

subroutine tao_hook_graph_postsetup_def (plot, graph)
  import
  implicit none
  type (tao_plot_struct) plot
  type (tao_graph_struct) graph
end subroutine
 
subroutine tao_hook_graph_setup_def (plot, graph, found)
  import
  implicit none
  type (tao_plot_struct) plot
  type (tao_graph_struct) graph
  logical found
end subroutine
 
subroutine tao_hook_init_beam_def ()
  implicit none
end subroutine

subroutine tao_hook_init_data_def ()
  implicit none
end subroutine

subroutine tao_hook_init_global_def (init_file, global)
  import
  implicit none
  type (tao_global_struct) global
  character(*) init_file
end subroutine
 
subroutine tao_hook_init_lattice_post_parse_def (u)
  import
  implicit none
  type (tao_universe_struct) u
end subroutine

subroutine tao_hook_init_plotting_def ()
  import
  implicit none
end subroutine

subroutine tao_hook_init_read_lattice_info_def (lat_file)
  implicit none
  character(*) lat_file
end subroutine

subroutine tao_hook_init1_def (init_file_name)
  implicit none
  character(*) init_file_name
end subroutine

subroutine tao_hook_init2_def ()
  implicit none
end subroutine

subroutine tao_hook_init_var_def()
  implicit none
end subroutine

subroutine tao_hook_lattice_calc_def (calc_ok)
  implicit none
  logical calc_ok
end subroutine

subroutine tao_hook_merit_data_def (i_uni, j_data, data, valid_value_set)
  import
  implicit none
  type (tao_data_struct) data
  integer, intent(in) :: i_uni, j_data
  logical valid_value_set
end subroutine

subroutine tao_hook_merit_var_def (i_uni, j_var, var)
  import
  implicit none
  type (tao_var_struct) var
  integer, intent(in) :: i_uni, j_var
end subroutine

subroutine tao_hook_optimizer_def (abort)
  implicit none
  logical abort
end subroutine
 
subroutine tao_hook_parse_command_args_def()
  implicit none
end subroutine

subroutine tao_hook_plot_setup_def()
  import
  implicit none
end subroutine

subroutine tao_hook_post_process_data_def ()
  implicit none
end subroutine
 
subroutine tao_hook_show_cmd_def (what, result_id, lines, nl)
  implicit none
  character(*) what, result_id
  character(*), allocatable :: lines(:)
  integer nl
end subroutine

end interface  ! abstract

! Function pointers

procedure(tao_hook_branch_calc_def), pointer :: tao_hook_branch_calc_ptr => null()
procedure(tao_hook_command_def), pointer :: tao_hook_command_ptr => null()
procedure(tao_hook_curve_s_pt_def), pointer :: tao_hook_curve_s_pt_ptr => null()
procedure(tao_hook_data_sanity_check_def), pointer :: tao_hook_data_sanity_check_ptr => null()
procedure(tao_hook_draw_floor_plan_def), pointer :: tao_hook_draw_floor_plan_ptr => null()
procedure(tao_hook_draw_graph_def), pointer :: tao_hook_draw_graph_ptr => null()
procedure(tao_hook_evaluate_a_datum_def), pointer :: tao_hook_evaluate_a_datum_ptr => null()
procedure(tao_hook_graph_postsetup_def), pointer :: tao_hook_graph_postsetup_ptr => null()
procedure(tao_hook_graph_setup_def), pointer :: tao_hook_graph_setup_ptr => null()
procedure(tao_hook_init_beam_def), pointer :: tao_hook_init_beam_ptr => null()
procedure(tao_hook_init_data_def), pointer :: tao_hook_init_data_ptr => null()
procedure(tao_hook_init_global_def), pointer :: tao_hook_init_global_ptr => null()
procedure(tao_hook_init_lattice_post_parse_def), pointer :: tao_hook_init_lattice_post_parse_ptr => null()
procedure(tao_hook_init_plotting_def), pointer :: tao_hook_init_plotting_ptr => null()
procedure(tao_hook_init_read_lattice_info_def), pointer :: tao_hook_init_read_lattice_info_ptr => null()
procedure(tao_hook_init1_def), pointer :: tao_hook_init1_ptr => null()
procedure(tao_hook_init2_def), pointer :: tao_hook_init2_ptr => null()
procedure(tao_hook_init_var_def), pointer :: tao_hook_init_var_ptr => null()
procedure(tao_hook_lattice_calc_def), pointer :: tao_hook_lattice_calc_ptr => null()
procedure(tao_hook_merit_data_def), pointer :: tao_hook_merit_data_ptr => null()
procedure(tao_hook_merit_var_def), pointer :: tao_hook_merit_var_ptr => null()
procedure(tao_hook_optimizer_def), pointer :: tao_hook_optimizer_ptr => null()
procedure(tao_hook_parse_command_args_def), pointer :: tao_hook_parse_command_args_ptr => null()
procedure(tao_hook_plot_setup_def), pointer :: tao_hook_plot_setup_ptr => null()
procedure(tao_hook_post_process_data_def), pointer :: tao_hook_post_process_data_ptr => null()
procedure(tao_hook_show_cmd_def), pointer :: tao_hook_show_cmd_ptr => null()

contains

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Function tao_pointer_to_universe_int (ix_uni, neg2_to_default) result (u)
!
! Overloaded by tao_pointer_to_universe. See this routine for more details.
!-

function tao_pointer_to_universe_int (ix_uni, neg2_to_default) result(u)

implicit none

type (tao_universe_struct), pointer :: u
integer ix_uni, ix_u
logical, optional :: neg2_to_default
character(*), parameter :: r_name = 'tao_pointer_to_universe_int'

!

ix_u = tao_universe_index(ix_uni, neg2_to_default)

if (ix_u < lbound(s%u, 1) .or. ix_u > ubound(s%u, 1)) then
  call out_io (s_fatal$, r_name, 'UNIVERSE INDEX OUT OF RANGE: \I0\ ', ix_u)
  nullify (u)
  return
endif

u => s%u(ix_u)

end function tao_pointer_to_universe_int

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Function tao_pointer_to_universe_str (string, neg2_to_default) result (u)
!
! Overloaded by tao_pointer_to_universe. See this routine for more details.
!-

function tao_pointer_to_universe_str (string, neg2_to_default) result(u)

implicit none

type (tao_universe_struct), pointer :: u
integer ix, ix_u
logical, optional :: neg2_to_default
character(*) string
character(*), parameter :: r_name = 'tao_pointer_to_universe_str'

!

nullify(u)

ix = tao_uni_ampersand_index(string)
if (ix == 0) then
  u => s%u(tao_universe_index(-1))
  return
elseif (string(1:ix-1) == '') then
  u => s%u(tao_universe_index(-1))
  string = string(ix+1:)
  return
endif

!

if (.not. is_integer(string(1:ix-1))) then
  call out_io (s_fatal$, r_name, 'MALFORMED UNIVERSE STRING')
  return
endif
read (string(1:ix-1), *) ix_u
string = string(ix+1:)

ix_u = tao_universe_index(ix_u, neg2_to_default)

if (ix_u < lbound(s%u, 1) .or. ix_u > ubound(s%u, 1)) then
  call out_io (s_fatal$, r_name, 'UNIVERSE INDEX OUT OF RANGE: \I0\ ', ix_u)
  return
endif

u => s%u(ix_u)

end function tao_pointer_to_universe_str

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Function tao_uni_ampersand_index(string) result (ix_amp)
!
! Routine to return the index of an ampersand ("@") sign in a string if the ampersand is
! being used as a separator between a universe spec and the rest of the string.
!
! For example:
!   string = "[1:3]@orbit.x[5] => ix_amp = 6
!   string = "orbit.x[5@0.2]   => ix_amp = 0 (no universe "@" present)
!
! Input:
!   string      -- character(*): String to parse
!
! Output:
!   ix_amp      -- integer: Index of universe "@". Set to zero if no universe "@" found.
!-

function tao_uni_ampersand_index(string) result (ix_amp)

implicit none

integer ix_amp, i
character(*) string

! Any characters before a uni "@" must be in the set '-0123456789:,[]*'.

ix_amp = 0

do i = 1, len(string)
  if (index('-0123456789:,[]*', string(i:i)) == 0 .and. string(i:i) /= '@') return
  if (string(i:i) /= '@') cycle
  if (index(string(1:i), '::') /= 0) return
  ix_amp = i
  return
enddo

end function tao_uni_ampersand_index

end module
