!+
! Subroutine init_lat (lat, n, init_beginning_ele)
!
! Subroutine to initialize a Bmad lat_struct.
! 
! Input:
!   n                  -- integer, optional: Upper bound lat%ele(0:) array is initialized to.
!                           Default is 10.
!   init_beginning_ele -- logical, optional: Init lat%ele(0)? Default is False.
!
! Output:
!   lat -- lat_struct: Initialized lat.
!-

subroutine init_lat (lat, n, init_beginning_ele)

use equal_mod, dummy => init_lat

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele

integer, optional :: n
logical, optional :: init_beginning_ele

!

call deallocate_lat_pointers (lat)
call allocate_branch_array (lat, 0)

if (present(n)) call allocate_lat_ele_array(lat, n)
call init_ele (lat%ele_init)

call reallocate_control (lat, 100)

lat%title = ' '
lat%use_name = ' '
lat%lattice = ' '
lat%input_file_name = ' '

lat%param = lat_param_struct()
call set_status_flags (lat%param%bookkeeping_state, stale$)

call init_mode_info (lat%a)
call init_mode_info (lat%b)
call init_mode_info (lat%z)

lat%n_ele_track = 0
lat%n_ele_max = 0
lat%n_control_max = 0
lat%n_ic_max = 0
lat%input_taylor_order = 0
lat%version = -1

call nametable_init(lat%nametable, 0)

if (logic_option(.false., init_beginning_ele)) then
  ele => lat%ele(0)
  call init_ele(ele, beginning_ele$, 0, 0, lat%branch(0))
  ele%name = 'BEGINNING'
  call set_ele_defaults (ele)   ! Defaults for beginning_ele element
  call nametable_add (lat%nametable, ele%name, 0)
endif

!----------------------------------------
contains

subroutine init_mode_info (t)
type (mode_info_struct) t
t%tune = 0
t%emit = 0
t%chrom = 0
end subroutine init_mode_info

end subroutine init_lat

