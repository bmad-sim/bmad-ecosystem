!+
! Function djb_str_hash (in_str) result (hash_str)
!
! Routine to return a string hash of a string.
! The characters in the hash string are drawn from the set [0-9, a-z, A-Z].
!
! This routine is modeled after the algorithm by Daniel J. Bernstein.
! Also see djb_hash
!
! Input:
!   in_str    -- character(*) Input string. Trailing blanks are ignored in hashing.
!
! Output:
!   hash      -- character(6): Hash of in_str
!-

function djb_str_hash (in_str) result (hash_str)

implicit none

character(*) :: in_str
character(6) :: hash_str
character(62), parameter :: tran = '0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'

integer i, ix, hash

!

hash = 5381

do i = 1, len_trim(in_str)
  hash = (ishft(hash,5) + hash) + ichar(in_str(i:i))
end do

do i = 1, 6
  if (i == 6 .and. hash < 0) hash = 62 + hash
  ix = modulo(hash, 62) + 1
  hash_str(i:i) = tran(ix:ix)
  hash = hash / 62
enddo

end function 
