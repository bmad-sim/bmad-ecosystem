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
  function tao_merit (s) result (this_merit)
    use tao_struct, only: tao_super_universe_struct, rp
    implicit none
    type (tao_super_universe_struct), target :: s
    real(rp) this_merit
  end function
end interface
 
interface
  subroutine tao_plot_cmd (s, where, who)
    use tao_struct, only: tao_super_universe_struct
    implicit none
    type (tao_super_universe_struct) :: s
    character(*) :: where
    character(*) :: who(:)
  end subroutine
end interface
 
interface
  subroutine tao_view_cmd (s, i_universe)
    use tao_struct, only: tao_super_universe_struct
    implicit none
    type (tao_super_universe_struct) :: s
    integer i_universe
  end subroutine
end interface
 
interface
  subroutine tao_call_cmd (s, file_name, cmd_arg)
    use tao_struct, only: tao_super_universe_struct
    implicit none
    type (tao_super_universe_struct) :: s
    character(*) :: file_name
    character(*) :: cmd_arg(:)
  end subroutine
end interface
 
interface
  subroutine tao_alias_cmd (s, alias, string)
    use tao_struct, only: tao_super_universe_struct
    implicit none
    type (tao_super_universe_struct) :: s
    character(*) :: alias
    character(*) :: string
  end subroutine
end interface
 
interface
  subroutine tao_output_cmd (s, what)
    use tao_struct, only: tao_super_universe_struct
    implicit none
    type (tao_super_universe_struct) :: s
    character(*) :: what
  end subroutine
end interface
 
interface
  subroutine tao_clip_cmd (s, where, value1, value2)
    use tao_struct, only: tao_super_universe_struct, rp
    implicit none
    type (tao_super_universe_struct) :: s
    character(*) :: where
    real(rp) value1, value2
  end subroutine
end interface
 
interface
  subroutine tao_scale_cmd (s, where, value1, value2)
    use tao_struct, only: tao_super_universe_struct, rp
    implicit none
    type (tao_super_universe_struct) :: s
    character(*) :: where
    real(rp) value1, value2
  end subroutine
end interface
 
interface
  subroutine tao_change_cmd (s, who, name, where, num_str)
    use tao_struct, only: tao_super_universe_struct
    implicit none
    type (tao_super_universe_struct) :: s
    character(*) :: who, name, where, num_str
  end subroutine
end interface
 
interface
  subroutine tao_set_data_cmd (s, who, component, set_value, list)
    use tao_struct, only: tao_super_universe_struct
    implicit none
    type (tao_super_universe_struct) :: s
    character(*) who, component, set_value, list
  end subroutine
end interface
 
interface
  subroutine tao_set_var_cmd (s, who, component, set_value, list)
    use tao_struct, only: tao_super_universe_struct
    implicit none
    type (tao_super_universe_struct) :: s
    character(*) who, component, set_value, list
  end subroutine
end interface
 
interface
  subroutine tao_set_global_cmd (s, who, set_value)
    use tao_struct, only: tao_super_universe_struct
    implicit none
    type (tao_super_universe_struct) :: s
    character(*) who, set_value
  end subroutine
end interface
 
interface
  subroutine tao_command (s, cmd_line, err)
    use tao_struct, only: tao_super_universe_struct
    implicit none
    type (tao_super_universe_struct) :: s
    character(*) :: cmd_line
    logical err
  end subroutine
end interface
 
interface
  subroutine tao_plot_data_setup (s)
    use tao_struct, only: tao_super_universe_struct
    implicit none
    type (tao_super_universe_struct), target :: s
  end subroutine
end interface
 
interface
  subroutine tao_de_optimizer (s)
    use tao_struct, only: tao_super_universe_struct
    implicit none
    type (tao_super_universe_struct), target :: s
  end subroutine
end interface
 
interface
  subroutine tao_get_user_input (s, cmd_line)
    use tao_struct, only: tao_super_universe_struct
    implicit none
    type (tao_super_universe_struct) :: s
    character(*) :: cmd_line
  end subroutine
end interface
 
interface
  subroutine tao_get_vars (s, var_value, var_del, var_weight, var_data_value)
    use tao_struct, only: tao_super_universe_struct, rp
    implicit none
    type (tao_super_universe_struct), target :: s
    real(rp), allocatable :: var_value(:)
    real(rp), allocatable, optional :: var_del(:)
    real(rp), allocatable, optional :: var_weight(:)
    real(rp), allocatable, optional :: var_data_value(:)
  end subroutine
end interface
 
interface
  subroutine tao_init (s, init_file)
    use tao_struct, only: tao_super_universe_struct
    implicit none
    type (tao_super_universe_struct) s
    character(*) init_file
  end subroutine
end interface
 
interface
  subroutine tao_init_global_and_universes (s, data_and_var_file)
    use tao_struct, only: tao_super_universe_struct
    implicit none
    type (tao_super_universe_struct) s
    character(*) data_and_var_file
  end subroutine
end interface

interface
  subroutine tao_init_plotting (s, plot_file)
    use tao_struct, only: tao_super_universe_struct
    implicit none
    type (tao_super_universe_struct), target :: s
    character(*) plot_file
  end subroutine
end interface
 
interface
  subroutine tao_lm_optimizer (s)
    use tao_struct, only: tao_super_universe_struct
    implicit none
    type (tao_super_universe_struct), target :: s
  end subroutine
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
  subroutine tao_place_cmd (s, where, who)
    use tao_struct, only: tao_super_universe_struct
    implicit none
    type (tao_super_universe_struct), target :: s
    character(*) who
    character(*) where
  end subroutine
end interface
 
interface
  subroutine tao_plot_out (s)
    use tao_struct, only: tao_super_universe_struct
    implicit none
    type (tao_super_universe_struct), target :: s
  end subroutine
end interface
 
interface
  subroutine tao_plot_struct_transfer (plot_in, plot_out, preserve_region)
    use tao_struct, only: tao_plot_struct
    type (tao_plot_struct) plot_in
    type (tao_plot_struct) plot_out
    logical preserve_region
  end subroutine
end interface
 
interface
  subroutine tao_run_cmd (s, which)
    use tao_struct, only: tao_super_universe_struct
    implicit none
    type (tao_super_universe_struct), target :: s
    character(*) which
  end subroutine
end interface
 
interface
  subroutine tao_set_vars (s, var_vec)
    use tao_struct, only: tao_super_universe_struct, rp
    implicit none
    type (tao_super_universe_struct), target :: s
    real(rp) var_vec(:)
  end subroutine
end interface
 
interface
  subroutine tao_show_cmd (s, show_word1, show_word2, show_word3, show_word4)
    use tao_struct, only: tao_super_universe_struct
    implicit none
    type (tao_super_universe_struct) :: s
    character(*) :: show_word1
    character(*) :: show_word2
    character(*) :: show_word3
    character(*) :: show_word4
  end subroutine
end interface
 
interface
  subroutine tao_hook_optimizer (s)
    use tao_struct, only: tao_super_universe_struct
    implicit none
    type (tao_super_universe_struct), target :: s
  end subroutine
end interface
 
interface
  subroutine tao_hook_command (s, cmd_line, found)
    use tao_struct, only: tao_super_universe_struct
    implicit none
    type (tao_super_universe_struct), target :: s
    character(*) cmd_line
    logical found
  end subroutine
end interface
 
interface
  subroutine tao_use_data (s, action, data_name, range)
    use tao_struct, only: tao_super_universe_struct
    implicit none
    type (tao_super_universe_struct) :: s
    character(*)                :: action
    character(*)                :: data_name
    character(*)                :: range
  end subroutine
end interface

interface
  subroutine tao_set_data_useit_opt (s)
    use tao_struct, only: tao_super_universe_struct
    implicit none
    type (tao_super_universe_struct) s
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
  subroutine tao_use_var (s, action, var_name, range)
    use tao_struct, only: tao_super_universe_struct
    implicit none
    type (tao_super_universe_struct) :: s
    character(*)                :: action
    character(*)                :: var_name
    character(*)                :: range
  end subroutine
end interface
 
interface
  subroutine tao_set_var_useit_opt (s)
    use tao_struct, only: tao_super_universe_struct
    implicit none
    type (tao_super_universe_struct) s
  end subroutine
end interface

interface
  subroutine tao_hook_merit_var (s, i_uni, j_var, var)
    use tao_struct, only: tao_super_universe_struct, tao_var_struct
    implicit none
    type (tao_super_universe_struct) s
    type (tao_var_struct) var
    integer i_uni, j_var
  end subroutine
end interface

interface
  subroutine tao_hook_merit_data (s, i_uni, j_data, data)
    use tao_struct, only: tao_super_universe_struct, tao_data_struct
    implicit none
    type (tao_super_universe_struct) s
    type (tao_data_struct) data
    integer i_uni, j_data
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
  subroutine tao_limit_calc (s, limited)
    use tao_struct, only: tao_super_universe_struct
    implicit none
    type (tao_super_universe_struct) s
    logical limited
  end subroutine
end interface

interface
  subroutine tao_hook_init_design_lattice (s, lat_file)
    use tao_struct, only: tao_super_universe_struct
    implicit none
    type (tao_super_universe_struct) s
    character(*) lat_file
  end subroutine
end interface

interface
  subroutine tao_lattice_calc (s)
    use tao_struct, only: tao_super_universe_struct
    implicit none
    type (tao_super_universe_struct) s
  end subroutine
end interface

interface
  subroutine tao_load_data_array (s)
    use tao_struct, only: tao_super_universe_struct
    implicit none
    type (tao_super_universe_struct) s
  end subroutine
end interface

interface
  subroutine tao_hook_load_data_array (s, data, found)
    use tao_struct, only: tao_super_universe_struct, tao_data_struct
    implicit none
    type (tao_super_universe_struct) s
    type (tao_data_struct) data
    logical found
  end subroutine
end interface

interface
  subroutine tao_help (s, help_what)
    use tao_struct, only: tao_super_universe_struct
    implicit none
    type (tao_super_universe_struct) s
    character(*) help_what
  end subroutine
end interface

end module


