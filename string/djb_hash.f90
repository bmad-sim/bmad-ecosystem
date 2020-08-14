!+
! Function djb_hash (str, old_hash) result (hash)
!
! Routine to return an integer hash of a string.
! Routine by Daniel J. Bernstein.
!
! The old_hash optional argument is used when the string to be hashed is fed into
! this routine in pieces. 
!
! Input:
!   str       -- character(*) Input string. Trailing blanks are ignored in hashing.
!   old_hash  -- integer, optional: Old value of the hash. 
!
! Output:
!   hash      -- integer: Hash of string.
!-

function djb_hash (str, old_hash) result (hash)

implicit none

character(*) :: str
integer :: hash
integer, optional :: old_hash
integer :: i

!

if (present(old_hash)) then
  hash = old_hash
else
  hash = 5381
endif

do i = 1, len_trim(str)
  hash = (ishft(hash,5) + hash) + ichar(str(i:i))
end do

end function djb_hash
