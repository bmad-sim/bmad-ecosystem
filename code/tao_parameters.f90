!+
! module tao_parameters
!
! File of parameters like the size of the u%data array.
!-

module tao_parameters

  integer, parameter :: n_legend_maxx = 20   ! number of lines for a legend.
  integer, parameter :: n_descrip_maxx = 10  ! d2_data descriptive info lines.
  integer, parameter :: n_template_maxx = 10 ! number of plot templates.
  integer, parameter :: n_data_maxx = 1000   ! max index of datum per d1_data
  integer, parameter :: n_data_minn = -100   ! min index of datum per d1_data
  integer, parameter :: n_var_maxx = 500     ! max index of datum per v1_var 
  integer, parameter :: n_var_minn = -100    ! min index of datum per v1_var 
  integer, parameter :: n_region_maxx = 10   ! number of plotting regions.
  integer, parameter :: n_curve_maxx = 20    ! number of curves per graph

! the name struct
  
  integer, parameter :: data$ = 1, variable$ = 2
  integer, parameter :: use$ = 1, veto$ = 2, restore$ = 3
  
  type  name_struct
    character(16) :: xy_plane(2) = (/ "X", "Y" /)
    character(16) :: data_or_var(2) = (/ "data    ", "variable"/)
    character(16) :: use_veto_restore(3) = (/ "use    ", "veto   ", "restore"/)
  end type

  type (name_struct), save :: name$

end module

