!+
! Subroutine string_trim2 (in_str, delimitors, out_str, ix_word, delim, ix_next)
!
! Subroutine to trim a string of leading delimitors and also to return the
! length of the first word.
!
! Input:
!   in_str     - Character(*): String to be parsed.
!   delimitors - Character(*): String containing delimitors to be used
!                   by the program. Should NOT contain a blank.
!
! Output:
!   out_str  - Character(*): String left shifted to trim leading blanks
!   ix_word  - Integer: Index in OUT_STR of last character in the first word.
!   delim    - Character(1): delimitor found. Set to ' ' if no delimitor found
!   ix_next  - Integer: Index in OUT_STR of 2nd word If 0 then nothing left
!
! Note:
!
! 1) If IN_STR has no non-blank characters the value of IX_WORD and
!    IX_NEXT is set to 0.
!
! 2) If
!
! Example:
!
!     in_str     = '   tobe : or not'
!     delimitors = ':,'
!     call string_trim2(in_str, delimitors, out_str, ix_word, delim, ix_delim)
!
! Output:
!     out_str:  'tobe : or not    '
!     ix_word:  4
!     ix_delim: 6
!-

subroutine string_trim2 (in_str, delimitors, out_str, ix_word, delim, ix_next)

implicit none

character(*) in_str, out_str, delimitors
character(1) tab, delim
parameter (tab = char(9))

integer ix, i, ix_word, ix_delim, ix_next, n_len

! trim leading blanks and set defaults

call string_trim (in_str, out_str, ix)
delim = ' '       ! default if no delim found

if (ix == 0) then    ! if nothing left
  ix_word = 0
  ix_next = 0
  return
endif

! see where first blank or delimitor is

n_len = len(out_str)

do i = 1, n_len

  if (index(delimitors, out_str(i:i)) /= 0) then   ! delimitor found
    ix_word  = i-1     ! end of word
    ix_delim = i       ! note were delim is
    delim = out_str(i:i)
    goto 1000
  endif

  if (out_str(i:i) == ' ' .or. out_str(i:i) == tab) then
    ix_word  = i-1     ! end of word
    ix_delim = 0       ! no delim found yet
    goto 1000
  endif

enddo


! move forward to next non-blank character
! if we find a second delimitor just treat it as the start of the next word

1000  continue

do i = i+1, n_len

  if (index(delimitors, out_str(i:i)) /= 0) then   ! delimitor found
    if (ix_delim /= 0) then    ! this is second delimitor
      ix_next = i
      return
    else  ! if first delim
      ix_delim = i
      delim = out_str(i:i)
      goto 1000                  ! back to find start of next word
    endif
  endif

  if (out_str(i:i) /= ' ' .and. out_str(i:i) /= tab) then
    ix_next = i
    return
  endif

enddo

ix_next = 0       ! nothing found

end subroutine
