!+
! Module tao_hook_mod
!
! Module to define the custom variables used for a particular machine.
! See tao_struct.f90 for where the hooks go.
!
! This particular file is a dummy module so the compiler will
! not complain if machine specific variables are not defined.
!-

module tao_hook_mod

  type tao_data_hook
    integer dummy
  end type

  type tao_d1_data_hook
    integer dummy
  end type

  type tao_d2_data_hook
    integer dummy
  end type

  type tao_var_hook
    integer dummy
  end type

  type tao_v1_var_hook
    integer dummy
  end type

  type tao_global_hook
    integer dummy
  end type

  type tao_curve_hook
    integer dummy
  end type

  type tao_graph_hook
    integer dummy
  end type

  type tao_plot_hook
    integer dummy
  end type

  type tao_plot_page_hook
    integer dummy
  end type

end module
