!+
! Module tao_hook_mod
!
!  Module to define  custom variables used in a custom version of Tao.
! The following structures are referenced in tao_struct.f90. See this file to find where the
! "hooks" go.
!
!   Keep in mind that the following type structures are used in the tao_struct
!  module. So, even if they are not used a "dummy" structure must remain here so
!  that the compiler doesn't complain. 
!
!   More structures can be added to this module for use in custom code.
!
!  Note: because tao_struct uses the structures defined in here, the tao library
!  must be compiled *after* this module is created.
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
