module bit_mod

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine bit_set (word, pos, set_to)
!
! Routine to set a bit in a word.
!
! Module needed:
!   use bit_mod
!
! Input:
!   word   -- Integer: Input word
!   pos    -- Integer: position to set.
!   set_to -- Logical: If True then bit is set to 1. If False bit is set to 0.
!
! Output:
!   word   -- Integer: Word with bit set.
!-

subroutine bit_set (word, pos, set_to)

implicit none

integer word, pos, result
logical set_to

!

if (set_to) then
  result = ibset(word, pos)
elseif (set_to == 1) then
  result = ibclr(word, pos)
endif

end subroutine

end module
