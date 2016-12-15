!+
! Module tao_interface
!
! Module to define the interfaces for the tao routines.
!-

! Note: To overcome a problem with the intel compiler the following routine interfaces
! Have been deleted:
!   tao_command
!   tao_command_end_calc
!   tao_plot_setup

module tao_interface

use tao_struct

interface

subroutine tao_alias_cmd (alias, string)
  implicit none
  character(*) :: alias
  character(*) :: string
end subroutine
 
subroutine tao_call_cmd (file_name, cmd_arg)
  implicit none
  character(*) :: file_name
  character(*), optional :: cmd_arg(:)
end subroutine
 
subroutine tao_top_level (command, errcode)
  implicit none
  character(*), optional :: command
  integer, optional :: errcode
end subroutine
 
subroutine tao_clip_cmd (gang, where, value1, value2)
  import
  implicit none
  logical gang
  character(*) :: where
  real(rp) value1, value2
end subroutine
 
subroutine tao_data_show_use (d2_data, lines, nl)
  import
  implicit none
  type (tao_d2_data_struct) :: d2_data
  character(*), optional, allocatable :: lines(:)
  integer, optional :: nl
end subroutine

subroutine tao_de_optimizer (abort)
  implicit none
  logical abort
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

subroutine tao_find_plots (err, name, where, plot, graph, curve, print_flag, always_allocate)
  import
  implicit none
  type (tao_plot_array_struct), allocatable, optional :: plot(:)
  type (tao_graph_array_struct), allocatable, optional :: graph(:)
  type (tao_curve_array_struct), allocatable, optional :: curve(:)
  character(*) name, where
  logical, optional :: print_flag, always_allocate
  logical err
end subroutine
 
subroutine tao_get_user_input (cmd_out, prompt_str, wait_flag, cmd_in, will_need_input)
  implicit none
  character(*) :: cmd_out
  character(*), optional :: prompt_str, cmd_in
  logical, optional :: wait_flag, will_need_input
end subroutine
 
subroutine tao_has_been_created ()
end subroutine
 
subroutine tao_help (what1, what2, lines, n_lines)
  implicit none
  character(*) what1, what2
  character(*), optional, allocatable :: lines(:)
  integer, optional :: n_lines
end subroutine

subroutine tao_hook_branch_calc (u, tao_lat, branch)
  import
  implicit none
  type (tao_universe_struct), target :: u
  type (tao_lattice_struct), target :: tao_lat
  type (branch_struct), target :: branch
end subroutine
 
subroutine tao_hook_command (command_line, found)
  implicit none
  character(*) command_line
  logical found
end subroutine
 
subroutine tao_hook_init_var() 
  implicit none
end subroutine

subroutine tao_hook_parse_command_args()
  implicit none
end subroutine

subroutine tao_hook_draw_floor_plan (plot, graph)
  import
  implicit none
  type (tao_plot_struct) plot
  type (tao_graph_struct) graph
end subroutine
 
subroutine tao_hook_draw_graph (plot, graph, found)
  import
  implicit none
  type (tao_plot_struct) plot
  type (tao_graph_struct) graph
  logical found
end subroutine
 
subroutine tao_hook_graph_setup (plot, graph, found)
  import
  implicit none
  type (tao_plot_struct) plot
  type (tao_graph_struct) graph
  logical found
end subroutine
 
subroutine tao_hook_graph_postsetup (plot, graph)
  import
  implicit none
  type (tao_plot_struct) plot
  type (tao_graph_struct) graph
end subroutine
 
subroutine tao_hook_evaluate_a_datum (found, datum, u, tao_lat, datum_value, valid_value, why_invalid)
  import
  implicit none
  type (tao_data_struct) datum
  type (tao_universe_struct) u
  type (tao_lattice_struct) tao_lat
  real(rp) datum_value
  logical found, valid_value
  character(*), optional :: why_invalid
end subroutine

subroutine tao_hook_merit_data (i_uni, j_data, data, valid_value_set)
  import
  implicit none
  type (tao_data_struct) data
  integer i_uni, j_data
  logical valid_value_set
end subroutine

subroutine tao_hook_merit_var (i_uni, j_var, var)
  import
  implicit none
  type (tao_var_struct) var
  integer i_uni, j_var
end subroutine

subroutine tao_hook_optimizer (abort)
  implicit none
  logical abort
end subroutine
 
subroutine tao_hook_post_process_data ()
  implicit none
end subroutine
 
subroutine tao_hook_show_cmd (what, stuff, result_id, lines, nl)
  implicit none
  character(*) what, stuff, result_id
  character(*), allocatable :: lines(:)
  integer nl
end subroutine

subroutine tao_hook_lattice_calc (calc_ok)
  implicit none
  logical calc_ok
end subroutine

subroutine tao_hook_init_data () 
  implicit none
end subroutine

subroutine tao_hook_init_beam ()
  implicit none
end subroutine

subroutine tao_hook_init_global (init_file, global)
  import
  implicit none
  type (tao_global_struct) global
  character(*) init_file
end subroutine
 
subroutine tao_hook_init_lattice_post_process (u)
  import
  implicit none
  type (tao_universe_struct) u
end subroutine

subroutine tao_hook_init_read_lattice_info (lat_file)
  implicit none
  character(*) lat_file
end subroutine

subroutine tao_hook_init1 (init_file_name)
  implicit none
  character(*) init_file_name
end subroutine

subroutine tao_hook_init2 ()
  implicit none
end subroutine

subroutine tao_init (err_flag)
  implicit none
  logical :: err_flag
end subroutine
 
subroutine tao_init_lattice (lat_file)
  implicit none
  character(*) lat_file
end subroutine

subroutine tao_init_plotting (plot_file)
  implicit none
  character(*) plot_file
end subroutine

subroutine tao_init_single_mode (single_mode_file)
  implicit none
  character(*) single_mode_file
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

function tao_merit (calc_ok) result (this_merit)
  import
  implicit none
  real(rp) this_merit
  logical, optional :: calc_ok
end function
 
subroutine tao_open_file (file, iunit, file_name, error_severity)
  implicit none
  character(*) file
  character(*) file_name
  integer iunit, error_severity
end subroutine

function tao_pointer_to_datum (d1, ele_name) result (datum_ptr)
  import
  implicit none
  type (tao_d1_data_struct), target :: d1
  type (tao_data_struct), pointer :: datum_ptr
  character(*) ele_name
end function

 
subroutine tao_write_cmd (what)
  implicit none
  character(*) :: what
end subroutine
 
subroutine tao_pause_cmd (time)
  import
  implicit none
  real(rp) time
end subroutine
 
subroutine tao_place_cmd (where, who)
  implicit none
  character(*) who
  character(*) where
end subroutine
 
subroutine tao_plot_cmd (where, component)
  implicit none
  character(*) :: where
  character(*) :: component
end subroutine
 
subroutine tao_plot_struct_transfer (plot_in, plot_out)
  import
  implicit none
  type (tao_plot_struct) plot_in
  type (tao_plot_struct) plot_out
end subroutine
 
subroutine tao_read_cmd (which, file)
  implicit none
  character(*) which, file
end subroutine
 
subroutine tao_run_cmd (which, abort)
  implicit none
  character(*) which
  logical abort
end subroutine
 
subroutine tao_set_data_useit_opt (data)
  import
  implicit none
  type (tao_data_struct), optional :: data(:)
end subroutine

subroutine tao_set_var_useit_opt ()
end subroutine

subroutine tao_single_mode (char)
  implicit none
  character(1) :: char
end subroutine

function tao_universe_number (i_uni) result (i_this_uni)
  import
  implicit none
  integer i_uni, i_this_uni
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
 
subroutine tao_var_show_use (v1_var, lines, nl)
  import
  implicit none
  type (tao_v1_var_struct) :: v1_var
  character(*), optional, allocatable :: lines(:)
  integer, optional :: nl
end subroutine

subroutine tao_x_axis_cmd (where, what)
  implicit none
  character(*) where
  character(*) what
end subroutine

end interface

integer, private, save :: dummy = 0 ! So ranlib will not complain about no symbols

end module


