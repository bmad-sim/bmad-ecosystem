!+
! Module tao_interface
!
! Module to define the interfaces for the tao routines.
!-

module tao_interface

interface
  function merit_wrapper (var_vec, end_flag) result (merit)
    use precision_def
    real(rp) var_vec(:)           ! Input: trial solution.
    logical, optional :: end_flag ! Output: Set True to terminate opti_de.
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
  subroutine tao_change_cmd (who, name, where, num_str)
    implicit none
    character(*) :: who, name, where, num_str
  end subroutine
end interface
 
interface
  subroutine tao_clip_cmd (where, value1, value2)
    use precision_def
    implicit none
    character(*) :: where
    real(rp) value1, value2
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
  subroutine tao_de_optimizer ()
    implicit none
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
  subroutine tao_help (help_what)
    implicit none
    character(*) help_what
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
  subroutine tao_hook_plot_graph (plot, graph, found)
    use tao_struct, only: tao_plot_struct, tao_graph_struct
    implicit none
    type (tao_plot_struct) plot
    type (tao_graph_struct) graph
    logical found
  end subroutine
end interface
 
interface
  subroutine tao_hook_evaluate_a_datum (found, datum, u, tao_lat, datum_value)
    use tao_struct, only: tao_data_struct, tao_universe_struct, tao_lattice_struct
    use precision_def, only: rp
    implicit none
    type (tao_data_struct) datum
    type (tao_universe_struct) u
    type (tao_lattice_struct) tao_lat
    real(rp) datum_value
    logical found
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
  subroutine tao_hook_optimizer ()
    implicit none
  end subroutine
end interface
 
interface
  subroutine tao_hook_post_process_data ()
    implicit none
  end subroutine
end interface
 
interface
  subroutine tao_hook_lattice_calc (universe, tao_lat, used)
    use tao_struct, only: tao_universe_struct, tao_lattice_struct
    implicit none
    type (tao_universe_struct) :: universe
    type (tao_lattice_struct) :: tao_lat
    logical used
  end subroutine
end interface
 
interface
  subroutine tao_hook_init_design_lattice (design_lattice_file, custom_init)
    use tao_input_struct, only: tao_design_lat_input
    implicit none
    type (tao_design_lat_input)  design_lattice_file(:)
    logical custom_init
  end subroutine
end interface

interface
  subroutine tao_hook_init (init_file_name)
    implicit none
    character(*) init_file_name
  end subroutine
end interface

interface
  subroutine tao_init (init_file)
    implicit none
    character(*) init_file
  end subroutine
end interface
 
interface
  subroutine tao_init_design_lattice (lat_file)
    implicit none
    character(*) lat_file
  end subroutine
end interface

interface
  subroutine tao_init_global_and_universes (init_file, data_file, var_file)
    implicit none
    character(*) init_file, data_file, var_file
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
  subroutine tao_lm_optimizer ()
    implicit none
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
  subroutine tao_mrq_func (x, a, y_fit, dy_da)
    use precision_def
    implicit none
    real(rp), intent(in) :: x(:), a(:)
    real(rp), intent(out) :: y_fit(:)
    real(rp), intent(out) :: dy_da(:, :)
  end subroutine
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
  subroutine tao_output_cmd (what, who)
    implicit none
    character(*) :: what, who
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
  subroutine tao_plot_cmd (where, who)
    implicit none
    character(*) :: where
    character(*) :: who(:)
  end subroutine
end interface
 
interface
  subroutine tao_plot_data_setup ()
    implicit none
  end subroutine
end interface
 
interface
  subroutine tao_cmd_end_calc ()
    implicit none
  end subroutine
end interface
 
interface
  subroutine tao_plot_struct_transfer (plot_in, plot_out)
    use tao_struct, only: tao_plot_struct
    type (tao_plot_struct) plot_in
    type (tao_plot_struct) plot_out
  end subroutine
end interface
 
interface
  subroutine tao_run_cmd (which)
    implicit none
    character(*) which
  end subroutine
end interface
 
interface
  subroutine tao_set_data_useit_opt ()
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


