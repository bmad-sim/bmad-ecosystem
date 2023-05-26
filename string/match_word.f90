!+
! Subroutine match_word (string, names, ix, exact_case, can_abbreviate, matched_name)
!
! Subroutine to match the first word in a string against a list of names.
!
! Input:
!   string          -- Character(*): String whose first word is to be matched
!   names(:)        -- Character(*): Array of names. 
!   exact_case      -- Logical, optional: If present and True then the match must
!                       be exact. Default: Match is case insensitive.
!   can_abbreviate  -- Logical, optional: If present and False then abbreviations
!                       are not permitted. Default is True.
! Output:
!   ix              -- Integer: Index in names that matched.
!                       = 0 if no match.
!                       < 0 if multiple matches
!   matched_name    -- Character(*), optional: Name that is matched to.
!                       Set to '' if no or multiple matches.
!
! Note: If an exact match is found then it superceeds all abreviated matches
!-

subroutine match_word (string, names, ix, exact_case, can_abbreviate, matched_name)

use sim_utils, except_dummy => match_word

implicit none

character(*) string, names(:)
character(*), optional :: matched_name
character(60) str, name
integer ix, i, ixs, ixm
integer :: match, exact$ = 1, no_match$ = 2, abrev$ = 3, mult_abrev$ = 4
logical, optional :: exact_case, can_abbreviate
logical exact, can_abbrev

! Init

exact = logic_option (.false., exact_case)
can_abbrev = logic_option (.true., can_abbreviate)

ix = 0                         ! if no match found
call string_trim (string, str, ixs)
if (.not. exact) call str_upcase(str, str)

if (ixs == 0) then
  if (present(matched_name)) matched_name = ''
  return           ! blank match string
endif

! loop over all match names

match = no_match$

do i = 1, size(names)

  if (names(i) == ' ') cycle

  call string_trim (names(i), name, ixm)

  if (.not. can_abbrev .and. len_trim(name) /= ixs) cycle

  if (.not. exact) call str_upcase(name, name)

  if (str(:ixs) == name(:ixs)) then
    if (ixm == ixs) then
      if (match == exact$) then
        ix = -ix
        return
      endif
      match = exact$
      ix = i
    else
      if (match == abrev$) then
        match = mult_abrev$
        ix = -ix
      elseif (match == no_match$) then
        match = abrev$
        ix = i
      endif
    endif
  endif

enddo

if (present(matched_name)) then
  if (ix > 0) then
    matched_name = names(ix)
  else
    matched_name = ''
  endif
endif

end subroutine
