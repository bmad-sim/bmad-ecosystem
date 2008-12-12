!+
! Module tao_interface
!
! Module to define the interfaces for the tao routines.
!-

module tao_interface

interface
  function merit_wrapper (var_vec, status, iter_count) result (merit)
    use precision_def
    real(rp) var_vec(:)           ! Input: trial solution.
    integer status
    integer iter_count
    real(rp) merit                ! Output: Merit value corresponting to vec.
  end function
end interface

interface
  subroutine tao_alias_cmd (alias, string)
    implicit none
    character(*) :: alias
    character(*) :: string
  end subroutine
end interface
 
interface
  subroutine tao_call_cmd (file_name, cmd_arg)
    implicit none
    character(*) :: file_name
    character(*), optional :: cmd_arg(:)
  end subroutine
end interface
 
interface
  subroutine tao_cl (prompt_string)
    implicit none
    character(*), optional :: prompt_string
  end subroutine
end interface
 
interface
  subroutine tao_clip_cmd (gang, where, value1, value2)
    use precision_def
    implicit none
    logical gang
    character(*) :: where
    real(rp) value1, value2
  end subroutine
end interface
 
interface
  subroutine tao_cmd_end_calc ()
    implicit none
  end subroutine
end interface
 
interface
  subroutine tao_command (cmd_line, err)
    implicit none
    character(*) :: cmd_line
    logical err
  end subroutine
end interface
 
interface
  subroutine tao_data_show_use (d2_data)
    use tao_struct, only: tao_d2_data_struct
    implicit none
    type (tao_d2_data_struct) :: d2_data
  end subroutine
end interface

interface
  subroutine tao_de_optimizer (abort)
    implicit none
    logical abort
  end subroutine
end interface
 
interface
  subroutine tao_get_user_input (cmd_line, prompt_str)
    implicit none
    character(*) :: cmd_line
    character(*), optional :: prompt_str 
  end subroutine
end interface
 
interface
  subroutine tao_has_been_created ()
  end subroutine
end interface
 
interface
  subroutine tao_help (what1, what2)
    implicit none
    character(*) what1, what2
  end subroutine
end interface

interface
  subroutine tao_hook_command (command_line, found)
    implicit none
    character(*) command_line
    logical found
  end subroutine
end interface
 
interface
  subroutine tao_hook_init_var(is_set) 
    implicit none
    logical is_set
  end subroutine
end interface

interface
  subroutine tao_hook_parse_command_args(is_set)
    implicit none
    logical is_set
  end subroutine
end interface

interface
  subroutine tao_hook_draw_graph (plot, graph, found)
    use tao_struct, only: tao_plot_struct, tao_graph_struct
    implicit none
    type (tao_plot_struct) plot
    type (tao_graph_struct) graph
    logical found
  end subroutine
end interface
 
interface
  subroutine tao_hook_graph_setup (plot, graph, found)
    use tao_struct, only: tao_plot_struct, tao_graph_struct
    implicit none
    type (tao_plot_struct) plot
    type (tao_graph_struct) graph
    logical found
  end subroutine
end interface
 
interface
  subroutine tao_hook_graph_postsetup (plot, graph)
    use tao_struct, only: tao_plot_struct, tao_graph_struct
    implicit none
    type (tao_plot_struct) plot
    type (tao_graph_struct) graph
  end subroutine
end interface
 
interface
  subroutine tao_hook_evaluate_a_datum (found, datum, u, tao_lat, datum_value, valid_value)
    use tao_struct, only: tao_data_struct, tao_universe_struct, tao_lattice_struct
    use precision_def, only: rp
    implicit none
    type (tao_data_struct) datum
    type (tao_universe_struct) u
    type (tao_lattice_struct) tao_lat
    real(rp) datum_value
    logical found, valid_value
  end subroutine
end interface

interface
  subroutine tao_hook_merit_data (i_uni, j_data, data)
    use tao_struct, only: tao_data_struct
    implicit none
    type (tao_data_struct) data
    integer i_uni, j_data
  end subroutine
end interface

interface
  subroutine tao_hook_merit_var (i_uni, j_var, var)
    use tao_struct, only: tao_var_struct
    implicit none
    type (tao_var_struct) var
    integer i_uni, j_var
  end subroutine
end interface

interface
  subroutine tao_hook_optimizer (abort)
    implicit none
    logical abort
  end subroutine
end interface
 
interface
  subroutine tao_hook_post_process_data ()
    implicit none
  end subroutine
end interface
 
interface
  subroutine tao_hook_show_cmd (what, stuff, result_id, lines, nl)
    implicit none
    character(*) what, stuff, result_id
    character(*), allocatable :: lines(:)
    integer nl
  end subroutine
end interface

interface
  subroutine tao_hook_lattice_calc (calc_ok)
    implicit none
    logical calc_ok
  end subroutine
end interface

interface
  subroutine tao_hook_init_data (is_set) 
    implicit none
    logical is_set
  end subroutine
end interface

interface
  subroutine tao_hook_init_beam (is_set)
    implicit none
    logical is_set
  end subroutine
end interface

interface
  subroutine tao_hook_init_global (init_file, global)
    use tao_struct, only: tao_global_struct
    implicit none
    type (tao_global_struct) global
    character(*) init_file
  end subroutine
end interface
 
interface
  subroutine tao_hook_init_connected_uni (is_set)
    implicit none
    logical is_set
  end subroutine
end interface
 
interface
  subroutine tao_hook_init_lattice_post_process (u)
    use tao_struct, only: tao_universe_struct
    implicit none
    type (tao_universe_struct) u
  end subroutine
end interface

interface
  subroutine tao_hook_init_read_lattice_info (lat_file, is_set)
    implicit none
    character(*) lat_file
    logical is_set
  end subroutine
end interface

interface
  subroutine tao_hook_init1 (init_file_name)
    implicit none
    character(*) init_file_name
  end subroutine
end interface

interface
  subroutine tao_hook_init2 ()
    implicit none
  end subroutine
end interface

interface
  subroutine tao_init (err_flag)
    implicit none
    logical, optional :: err_flag
  end subroutine
end interface
 
interface
  subroutine tao_init_lattice (lat_file)
    implicit none
    character(*) lat_file
  end subroutine
end interface

interface
  subroutine tao_init_plotting (plot_file)
    implicit none
    character(*) plot_file
  end subroutine
end interface

interface
  subroutine tao_init_single_mode (single_mode_file)
    implicit none
    character(*) single_mode_file
  end subroutine
end interface
 
interface
  subroutine tao_limit_calc (limited)
    implicit none
    logical limited
  end subroutine
end interface

interface
  subroutine tao_lmdif_optimizer (abort)
    implicit none
    logical abort
  end subroutine
end interface
 
interface
  function tao_merit (calc_ok) result (this_merit)
    use precision_def
    implicit none
    real(rp) this_merit
    logical, optional :: calc_ok
  end function
end interface
 
interface
  subroutine tao_open_file (logical_dir, file, iunit, file_name)
    implicit none
    character(*) logical_dir
    character(*) file
    character(*) file_name
    integer iunit
  end subroutine
end interface
 
interface
  subroutine tao_output_cmd (what)
    implicit none
    character(*) :: what
  end subroutine
end interface
 
interface
  subroutine tao_pause_cmd (time)
    use precision_def
    implicit none
    real(rp) time
  end subroutine
end interface
 
interface
  subroutine tao_place_cmd (where, who)
    implicit none
    character(*) who
    character(*) where
  end subroutine
end interface
 
interface
  subroutine tao_plot_cmd (where, component)
    implicit none
    character(*) :: where
    character(*) :: component
  end subroutine
end interface
 
interface
  subroutine tao_plot_setup ()
    implicit none
  end subroutine
end interface
 
interface
  subroutine tao_plot_struct_transfer (plot_in, plot_out)
    use tao_struct, only: tao_plot_struct
    implicit none
    type (tao_plot_struct) plot_in
    type (tao_plot_struct) plot_out
  end subroutine
end interface
 
interface
  subroutine tao_read_cmd (which, file)
    implicit none
    character(*) which, file
  end subroutine
end interface
 
interface
  subroutine tao_run_cmd (which, abort)
    implicit none
    character(*) which
    logical abort
  end subroutine
end interface
 
interface
  subroutine tao_set_data_useit_opt (data)
    use tao_struct, only : tao_data_struct
    implicit none
    type (tao_data_struct), optional :: data(:)
  end subroutine
end interface

interface
  subroutine tao_set_var_useit_opt ()
  end subroutine
end interface

interface
  subroutine tao_single_mode (char)
    implicit none
    character(1) :: char
  end subroutine
end interface

interface
  subroutine tao_use_data (action, data_name)
    implicit none
    character(*) :: action
    character(*) :: data_name
  end subroutine
end interface

interface
  subroutine tao_use_var (action, var_name)
    implicit none
    character(*) :: action
    character(*) :: var_name
  end subroutine
end interface
 
interface
  subroutine tao_var_show_use (v1_var)
    use tao_struct, only: tao_v1_var_struct
    implicit none
    type (tao_v1_var_struct) :: v1_var
  end subroutine
end interface

interface
  subroutine tao_view_cmd (i_universe)
    implicit none
    integer i_universe
  end subroutine
end interface
 
interface
  subroutine tao_x_axis_cmd (where, what)
    implicit none
    character(*) where
    character(*) what
  end subroutine
end interface

end module


