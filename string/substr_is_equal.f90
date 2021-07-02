!+
! Function substr_is_equal(var_str, n1, n2, str) result (is_equal)
!
! Routine to compare var_str(n1:n2) to str to see if the two are equal.
!
! Input:
!   var_str     -- character(:), allocatable: Variable length string.
!   n1          -- integer: Start index of var_str(n1:n2) substring.
!   n2          -- integer: End index of var_str(n1:n2) substring.
!   str         -- character(*) string to compare to var_str(n1:n2).
!
! Output:
!   is_equal    -- logical: True if var_str(n1:n2) and str are the same.
!-

function substr_is_equal(var_str, n1, n2, str) result (is_equal)

implicit none

character(*) str
character(:), allocatable :: var_str
integer n1, n2
logical is_equal

!

is_equal = .false.

if (.not. allocated(var_str)) return
if (n2 > len(var_str)) return

is_equal = (var_str(n1:n2) == str)

end function
