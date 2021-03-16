!+
! Subroutine new_control (lat, ix_ele, ele_name)
!
! Routine to create a new control element.
!
! Note: This routine may reallocate the lat%branch arrays so any pointers to elements will need to be repointed.
! Note: If ele_name is missing, set_ele_name should later be called to properly set lat%nametable.
!
! Input:
!     lat       -- lat_struct: Lat used
!     ele_name  -- character(*), optional: Name of the new element.
!
! Output
!     ix_ele -- Integer: Index of the new control element
!-

subroutine new_control (lat, ix_ele, ele_name)

use bmad_interface, except_dummy => new_control

implicit none

type (lat_struct)  lat
integer ix_ele
character(*), optional :: ele_name
character(40) name

! If the actual argument for ele_name is something like lat%ele(j)%name, a call to
! allocate_lat_ele_array will make the value of ele_name undefined. To get around this
! save ele_name.

lat%n_ele_max = lat%n_ele_max + 1
ix_ele = lat%n_ele_max
if (present(ele_name)) name = ele_name
if (ix_ele > ubound(lat%ele, 1)) call allocate_lat_ele_array (lat)
if (present(ele_name)) lat%ele(ix_ele)%name = name
call nametable_add(lat%nametable, lat%ele(ix_ele)%name, ix_ele)
lat%lord_state%attributes = stale$
lat%lord_state%control = stale$
lat%lord_state%s_position = stale$
lat%lord_state%floor_position = stale$

end subroutine
