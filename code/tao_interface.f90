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
  subroutine tao_cl (prompt_string)
    implicit none
    character(*), optional :: prompt_string
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
  subroutine tao_hook_count_d2_data (n_d2_data, is_set) 
    use tao_struct, only: s
    implicit none
    integer n_d2_data(lbound(s%u, 1):)
    logical is_set
  end subroutine
end interface

interface
  subroutine tao_hook_count_var(n_var, is_set) 
    implicit none
    integer n_var
    logical is_set
  end subroutine
end interface

interface
  subroutine tao_hook_read_d1_data (d1_data, data, ix_d1_data, ix_min_data, ix_max_data, &
                   default_weight, default_data_type, default_data_source, &
                   use_same_lat_eles_as, search_for_lat_eles, is_set)
    use tao_input_struct, only: tao_d1_data_input, tao_data_input, n_data_minn, &
                                n_data_maxx, rp
    implicit none
    type (tao_d1_data_input) d1_data
    type (tao_data_input) data(n_data_minn:n_data_maxx)
    integer ix_d1_data, ix_min_data, ix_max_data
    real(rp) default_weight        ! default merit function weight
    character(*) default_data_type, default_data_source
    character(*) use_same_lat_eles_as, search_for_lat_eles
    logical is_set
  end subroutine
end interface

interface
  subroutine tao_hook_read_d2_data (d2_data, n_d1_data, &
                                 default_merit_type, universe, is_set, ended)
    use tao_input_struct, only: tao_d2_data_input
    implicit none
    type (tao_d2_data_input) d2_data
    integer n_d1_data
    character(*) default_merit_type, universe
    logical is_set, ended
  end subroutine
end interface
 
interface
  subroutine tao_hook_read_var (v1_var, var, default_weight, default_step, default_key_delta, &
                    ix_min_var, ix_max_var, default_universe, default_attribute, &
                    default_low_lim, default_high_lim, default_merit_type, &
                    use_same_lat_eles_as, search_for_lat_eles, default_key_bound, &
                    is_set, ended)
    use tao_input_struct, only: rp, n_var_minn, n_var_maxx, tao_v1_var_input, tao_var_input
    implicit none
    type (tao_v1_var_input) v1_var
    type (tao_var_input) var(n_var_minn:n_var_maxx)
    real(rp) default_weight, default_step, default_key_delta, &
                    default_low_lim, default_high_lim
    integer ix_min_var, ix_max_var
    character(*) default_universe, default_attribute, default_merit_type
    character(*) use_same_lat_eles_as, search_for_lat_eles
    character(*) default_key_bound
    logical is_set, ended
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
  subroutine tao_hook_lattice_calc (calc_ok)
    implicit none
    logical calc_ok
  end subroutine
end interface

interface
  subroutine tao_hook_init_beam (ix_universe, beam0_file, ix_track_start, ix_track_end, &
                beam_all_file, beam_init, save_beam_at, is_set, ended)
    use tao_input_struct, only: beam_init_struct
    implicit none
    type (beam_init_struct) beam_init
    integer ix_universe, ix_track_start, ix_track_end
    character(*) beam0_file, beam_all_file
    character(*) save_beam_at(:)
    logical is_set, ended
  end subroutine
end interface

interface
  subroutine tao_hook_init_params (global, is_set)
    use tao_input_struct, only: tao_global_struct
    implicit none
    type (tao_global_struct) global
    logical is_set
  end subroutine
end interface
 
interface
  subroutine tao_hook_init_connected_uni (ix_universe, connect, is_set, ended)
    use tao_input_struct, only: tao_connected_uni_input
    implicit none
    integer ix_universe
    type (tao_connected_uni_input) connect
    logical is_set, ended
  end subroutine
end interface
 
interface
  subroutine tao_hook_init_lattice (lattice_file, custom_init)
    use tao_input_struct, only: tao_design_lat_input
    implicit none
    type (tao_design_lat_input)  lattice_file(:)
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
  subroutine tao_plot_struct_transfer (plot_in, plot_out)
    use tao_struct, only: tao_plot_struct
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


