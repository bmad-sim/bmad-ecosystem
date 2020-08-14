!+
! Function djb_hash (str) result (hash)
!
! Routine to return an integer hash of a string.
! Routine by Daniel J. Bernstein.
!
! Input:
!   str -- character(*) Input string. Trailing blanks are ignored in hashing
!
! Output:
!   hash  -- integer: Hash of string.
!-

function djb_hash (str) result (hash)

implicit none

character(*) :: str
integer :: hash
integer :: i

!

hash = 5381

do i = 1, len_trim(str)
  hash = (ishft(hash,5) + hash) + ichar(str(i:i))
end do

end function djb_hash
