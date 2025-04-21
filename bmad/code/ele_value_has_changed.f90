!+
! Function ele_value_has_changed (ele, list, abs_tol, set_old) result (has_changed)
!
! Routine to see if a parameter value in a lattice element has changed significantly.
! A relative tolerance of small_rel_change$ = 1d-14 is also added to abs_tol.
!
! Input:
!   ele         -- ele_struct: Element under consideration.
!   list(:)     -- integer: List of indexes of ele%value(:) array to check.
!   abs_tol(:)  -- real(rp): List of values such that if the change in parameter value is
!                    less than this it is not considered to have changed significantly.
!   set_old     -- logical: If True then set ele%old_value(j) = ele%value(j) for j in list
!
! Output:
!   ele         -- ele_struct: ele%old_value may be set depending upon setting of set_old
!   has_changed -- logical: Set True if a value has changed significantly.
!-

function ele_value_has_changed (ele, list, abs_tol, set_old) result (has_changed)

use bmad_routine_interface, dummy => ele_value_has_changed

implicit none

type (ele_struct) ele
integer list(:)
integer i, j
real(rp) abs_tol(:)
logical set_old, has_changed

!

has_changed = .false.
do i = 1, size(list)
  j = list(i)
  if (.not. significant_difference(ele%value(j), ele%old_value(j), abs_tol(i), small_rel_change$)) cycle
  has_changed = .true.
  exit
enddo

if (set_old) ele%old_value(list) = ele%value(list)

end function ele_value_has_changed
