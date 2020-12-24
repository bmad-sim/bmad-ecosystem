module bit_mod

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine bit_set (word, pos, set_to_1)
!
! Routine to set a bit in a word.
!
! Input:
!   word     -- Integer: Input word
!   pos      -- Integer: position to set.
!   set_to_1 -- Logical: If True then bit is set to 1. If False bit is set to 0.
!
! Output:
!   word   -- Integer: Word with bit set.
!-

subroutine bit_set (word, pos, set_to_1)

implicit none

integer word, pos, result
logical set_to_1

!

if (set_to_1) then
  result = ibset(word, pos)
else
  result = ibclr(word, pos)
endif

end subroutine

end module
