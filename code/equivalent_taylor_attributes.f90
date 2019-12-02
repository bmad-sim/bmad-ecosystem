!+ 
! Function equivalent_taylor_attributes (ele_taylor, ele2) result (equiv)
!
! Subroutine to see if two elements are equivalent in terms of attributes so
! that their Taylor Maps would be the same. 
!
! This routine is used to see if a taylor map from one element may be 
! used for another and thus save some computation time. elements of type taylor
! are considered *never* to be equivalent since their maps are never computed.
!
! Input: 
!   ele_taylor -- Ele_struct: Element with a Taylor map
!   ele2       -- Ele_struct: Element that might receive the Taylor map from ele_taylor.
!
! Output:
!   equiv -- logical: True if elements are equivalent.
!-

function equivalent_taylor_attributes (ele_taylor, ele2) result (equiv)

use bmad_interface, dummy => equivalent_taylor_attributes

implicit none

type (ele_struct) :: ele_taylor, ele2

integer it

logical equiv
logical vmask(num_ele_attrib$), vnot(num_ele_attrib$)

!

equiv = .false.

if (ele_taylor%key == taylor$) return  
if (ele_taylor%key /= ele2%key) return
if (ele_taylor%sub_key /= ele2%sub_key) return
if (ele_taylor%taylor_map_includes_offsets .neqv. ele2%taylor_map_includes_offsets) return
if (ele_taylor%value(integrator_order$) /= ele2%value(integrator_order$)) return

vmask = .true.
vmask(delta_ref_time$) = .false.
vmask(ref_time_start$) = .false.
if ((ele_taylor%key == wiggler$ .or. ele_taylor%key == undulator$) .and. ele_taylor%field_calc == fieldmap$) then
  vmask( [k1$, g_max$, b_max$] ) = .false.  ! These are dependent attributes.
endif
if (.not. ele_taylor%taylor_map_includes_offsets) then
  vmask( [x_offset$, y_offset$, z_offset$, tilt$, x_pitch$, &
            y_pitch$, x_offset_tot$, y_offset_tot$, z_offset_tot$, &
            tilt_tot$, x_pitch_tot$, y_pitch_tot$] ) = .false.
endif

vnot = (ele_taylor%value /= ele2%value)
vnot = vnot .and. vmask
if (any(vnot)) return

if (associated(ele_taylor%cartesian_map) .neqv. associated(ele2%cartesian_map)) return
if (associated(ele_taylor%cartesian_map)) then
  if (.not. all(ele_taylor%cartesian_map == ele2%cartesian_map)) return
endif

equiv = .true.

end function equivalent_taylor_attributes 
