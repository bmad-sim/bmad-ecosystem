!+
! Subroutine combine_consecutive_elements (lat)
!
! Routine to combine consecutive elements in the lattice that have the same name.
! This allows simplification, for example, of lattices where elements have been 
! split to compute the beta function at the center.
!
! Note: Lattice_bookkeeper is not called by this routine.
!
! Input:
!   lat -- Lat_struct: Lattice.
!
! Output:
!   lat -- Lat_struct: Lattice with elements combined.
!-

subroutine combine_consecutive_elements (lat)

use bookkeeper_mod, except => combine_consecutive_elements
use ptc_interface_mod, except2 => combine_consecutive_elements

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele1, ele2

integer i, jv

character(*), parameter :: r_name = 'combine_consecutive_elements'

! loop over all elements...

ele_loop: do i = 1, lat%n_ele_track - 1

  ele1 => lat%ele(i)
  ele2 => lat%ele(i+1)

  if (ele1%name /= ele2%name .or. ele1%key /= ele2%key) cycle ! ignore if not a matching pair

  ! Error checks

  select case (ele1%key)
  case (patch$, match$)
    call out_io (s_error$, r_name, 'ELEMENTS OF TYPE ' // trim(key_name(ele1%key)) // ' CANNOT BE COMBINED: ' // ele1%name)
  end select

  if (lat%ele(i-1)%name == ele1%name .and. ele2%key /= marker$) then
    call out_io (s_error$, r_name, 'TRIPLE CONSECUTIVE ELEMENTS HAVE SAME NAME! ' // ele1%name)
    cycle
  endif

  if (ele1%key == sbend$) then
    if (ele1%value(e2$) /= 0 .or. ele2%value(e1$) /= 0) then
      call out_io (s_error$, r_name, 'CONSECUTIVE BENDS HAVE "INTERNAL" FACE ANGLES: ' // ele1%name)
      cycle
    endif

    if (nint(ele1%value(fringe_at$)) == exit_end$ .or. nint(ele1%value(fringe_at$)) == both_ends$ .or. &
        nint(ele2%value(fringe_at$)) == entrance_end$ .or. nint(ele2%value(fringe_at$)) == both_ends$) then
      call out_io (s_error$, r_name, 'CONSECUTIVE BENDS HAVE "INTERNAL" FRINGE FIELDS: ' // ele1%name)
      cycle
    endif
  endif

  do jv = 1, size(ele1%value)
    if (attribute_name(ele1, jv) == 'REF_TIME_START' .or. attribute_name(ele1, jv) == null_name$) cycle
    if (ele1%key == sbend$) then
      select case (jv)
      case (e1$, e2$, hgap$, hgapx$, fint$, fintx$, h1$, h2$, fringe_at$); cycle
      end select
    end if

    if (abs(ele1%value(jv) - ele2%value(jv)) > 1d-14 * (abs(ele1%value(jv)) + abs(ele2%value(jv)))) then
      call out_io (s_error$, r_name, 'ELEMENT PARAMETERS DO NOT MATCH FOR: ' // ele1%name)
      cycle ele_loop
    endif
  enddo

  if (ele1%value(x_pitch_tot$) /= 0 .or. ele1%value(y_pitch_tot$) /= 0) then
    call out_io (s_error$, r_name, 'ELEMENT HAS NON-ZERO PITCH: ' // ele1%name)
    cycle
  endif

  ! Now combine

  ele2%ix_ele = -1   ! mark for deletion
  call set_flags_for_changed_attribute (ele1)

  ele1%value(l$) = 2 * ele1%value(l$)

  if (has_attribute(ele1, 'HKICK')) then
    ele1%value(hkick$) = 2 * ele1%value(hkick$)
    ele1%value(vkick$) = 2 * ele1%value(vkick$)
    ele1%value(BL_hkick$) = 2 * ele1%value(BL_hkick$)
    ele1%value(BL_vkick$) = 2 * ele1%value(BL_vkick$)
  endif

  if (has_attribute(ele1, 'KICK')) then
    ele1%value(kick$) = 2 * ele1%value(kick$)
    ele1%value(BL_kick$) = 2 * ele1%value(BL_kick$)
  endif

  if (associated(ele1%a_pole)) then
    ele1%a_pole = 2 * ele1%a_pole
    if (ele1%key /= multipole$) ele1%b_pole = 2 * ele1%b_pole
  endif

  if (associated(ele1%a_pole_elec)) then
    ele1%a_pole_elec = 2 * ele1%a_pole_elec
    ele1%b_pole_elec = 2 * ele1%b_pole_elec
  endif

  select case (ele1%key)
  case (sbend$)
    ele1%value(e2$) = ele2%value(e2$)
    ele1%value(h2$) = ele2%value(h2$)
    ele1%value(hgapx$) = ele2%value(hgapx$)
    ele1%value(fintx$) = ele2%value(fintx$)
    ele1%value(angle$) = 2 * ele1%value(angle$)
    if (nint(ele2%value(fringe_at$)) == exit_end$) then
      if (nint(ele1%value(fringe_at$)) == no_end$) then
        ele1%value(fringe_at$) = exit_end$
      else
        ele1%value(fringe_at$) = both_ends$
      endif
    endif

  case (rfcavity$, lcavity$)
    ele1%value(voltage$) = 2 * ele1%value(voltage$)

  case (taylor$)
    call concat_ele_taylor (ele1%taylor, ele2, ele1%taylor)

  end select

enddo ele_loop

call remove_eles_from_lat (lat)     ! Remove all null_ele elements

end subroutine
