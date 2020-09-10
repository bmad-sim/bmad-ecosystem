!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine init_lat (lat, n)
!
! Subroutine to initialize a Bmad lat.
! 
! Input:
!   n    -- Integer, optional: Upper bound lat%ele(0:) array is initialized to.
!
! Output:
!   lat -- lat_struct: Initialized lat.
!-

subroutine init_lat (lat, n)

use equal_mod, dummy => init_lat

implicit none

type (lat_struct)  lat

integer, optional :: n

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
lat%absolute_time_tracking   = bmad_com%absolute_time_tracking_default
lat%nametable%n_max = -1

!----------------------------------------
contains

subroutine init_mode_info (t)
type (mode_info_struct) t
t%tune = 0
t%emit = 0
t%chrom = 0
end subroutine init_mode_info

end subroutine init_lat

