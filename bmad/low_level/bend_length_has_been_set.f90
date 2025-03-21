!+
! Function bend_length_has_been_set(ele) result (is_set)
!
! Routine to determine if a bend element (rf_bend, rbend, or sbend) has its (arc) length set.
!
! Example: If an rbend element has b_field set along with l_chord, the arc length is not set until the reference energy has been computed.
! But superpositions happen before the reference energy is set (since superimposed elements may affect the reference energy.
! This chicken-or-egg problem needs to be flagged and resolved by the User.
!
! Input:
!   ele           -- ele_struct: Element to be checked.
!
! Ouput:
!   is_set        -- logical: Note: will be set True for non-bend elements.

function bend_length_has_been_set(ele) result (is_set)

use bmad_struct

implicit none

type (ele_struct) ele
logical is_set

!

is_set = .true.

select case(ele%key)
case (rf_bend$, sbend$, rbend$)
case default
  return
end select

if (is_false(ele%value(init_needed$))) return
if (ele%value(b_field$) == 0 .and. ele%value(db_field$) == 0) return
if (ele%value(l_chord$) == 0 .and. ele%value(l_rectangle$) == 0 .and. ele%value(angle$) == 0) return

is_set = .false.

end function
