!+
! Function substr(var_str, n1, n2) result (sub_str)
!
! Routine to extract a substring from a variable length string
!
! Input:
!   var_str     -- character(:), allocatable: Variable length string.
!   n1          -- integer: Start index of var_str(n1:n2) substring.
!   n2          -- integer: End index of var_str(n1:n2) substring.
!
! Output:
!   sub_str     -- character(n2-n1+1): var_str(n1:n2) substring. 
!                   Blanks are added if len(var_str) < n2
!-

function substr(var_str, n1, n2) result (sub_str)

implicit none

character(:), allocatable :: var_str
integer n1, n2, n
character(n2-n1+1) sub_str

!

if (.not. allocated(var_str)) then
  sub_str = ''

elseif (n2 > len(var_str)) then
  n = len(var_str)
  sub_str = ''
  sub_str(:n1-n) = var_str(n1:n)

else
  sub_str = var_str(n1:n2)
endif

end function
