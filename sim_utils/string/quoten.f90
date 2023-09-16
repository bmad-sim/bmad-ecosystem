!+
! Function quoten(str, delim) result (q_str)
!
! Function to put double quote marks (") around each string in an array and
! return the concatenated string with a delimitor between the strings.
!
! Input:
!   str(:)  -- character(*): Input string array
!   delim   -- character(*), optional: Delimitor between strings. Default is two blank spaces.
! Output:
!   q_str   -- character(:), allocatable: String with quote marks.
!-

function quoten(str, delim) result (q_str)

character(*) str(:)
character(*), optional :: delim
character(:), allocatable :: q_str
integer i, n, ns

!

ns = size(str)
n = 0
do i = 1, ns
  n = n + len_trim(str(i)) + 2
enddo

if (present(delim)) then
  n = n + (ns-1) * len(delim)
else
  n = n + (ns-1) * 2
endif

allocate(character(n) :: q_str)

n = 0
do i = 1, ns
  q_str = q_str(:n) // '"' // trim(str(i)) // '"'
  n = n + len_trim(str(i)) + 2
  if (i == ns) exit
  if (present(delim)) then
    q_str = q_str(:n) // delim
    n = n + len(delim)
  else
    q_str = q_str(:n) // '  '
    n = n + 2
  endif
enddo

end function quoten
