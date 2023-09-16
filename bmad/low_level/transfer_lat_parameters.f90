!+
! Subroutine transfer_lat_parameters (lat_in, lat_out)
!
! Subroutine to transfer the lat parameters (such as lat%name, lat%param, etc.)
! from one lat to another. The only stuff that is not transfered are things
! that are (or have) pointers.
!
! Input:
!   lat_in -- lat_struct: Input lat.
!
! Output:
!   lat_out -- lat_struct: Output lat with parameters set.
!-

subroutine transfer_lat_parameters (lat_in, lat_out)

use bmad_struct

implicit none

type (lat_struct), intent(in) :: lat_in
type (lat_struct) :: lat_out

!

lat_out%use_name                  = lat_in%use_name
lat_out%lattice                   = lat_in%lattice
lat_out%machine                   = lat_in%machine
lat_out%input_file_name           = lat_in%input_file_name
lat_out%title                     = lat_in%title
lat_out%a                         = lat_in%a
lat_out%b                         = lat_in%b
lat_out%z                         = lat_in%z
lat_out%param                     = lat_in%param
lat_out%lord_state                = lat_in%lord_state
lat_out%particle_start            = lat_in%particle_start
lat_out%beam_init                 = lat_in%beam_init
lat_out%pre_tracker               = lat_in%pre_tracker
lat_out%nametable                 = lat_in%nametable
lat_out%version                   = lat_in%version
lat_out%n_ele_track               = lat_in%n_ele_track
lat_out%n_ele_max                 = lat_in%n_ele_max
lat_out%n_control_max             = lat_in%n_control_max
lat_out%n_ic_max                  = lat_in%n_ic_max
lat_out%input_taylor_order        = lat_in%input_taylor_order
lat_out%photon_type               = lat_in%photon_type
lat_out%creation_hash             = lat_in%creation_hash

! Allocatable components

if (allocated(lat_in%print_str)) then
  lat_out%print_str = lat_in%print_str
else
  if (allocated(lat_out%print_str)) deallocate(lat_out%print_str)
endif

if (allocated(lat_in%constant)) then
  lat_out%constant = lat_in%constant
else
  if (allocated(lat_out%constant)) deallocate(lat_out%constant)
endif

if (allocated(lat_in%control)) then
  lat_out%control = lat_in%control
else
  if (allocated(lat_out%control)) deallocate(lat_out%control)
endif

if (allocated(lat_in%custom)) then
  lat_out%custom = lat_in%custom
else
  if (allocated(lat_out%custom)) deallocate(lat_out%custom)
endif

if (allocated(lat_in%ic)) then
  lat_out%ic = lat_in%ic
else
  if (allocated(lat_out%ic)) deallocate(lat_out%ic)
endif

end subroutine transfer_lat_parameters
