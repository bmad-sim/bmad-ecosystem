!+
! Function str_first_in_set (line, set, ignore_clauses) result (ix_match)
!
! Function to locate the first character in "line" that matches
! a character in "set". If no match is found, function returns 0.
!
! Input:
!   line            -- character(*): Line to parse
!   set             -- character(*): Set of characters to match to
!   ignore_clauses  -- logical, optional: If present and True, ignore characters in any clause.
!                       A clause is any substring delimited with parentheses "(...)", "[...]" "{...}" plus
!                       any substring delimited with quotes (single or double). Note parenteses in quotes are ignored.
!
! Output:
!   ix_match        -- integer: Index of the first matched character in "line"
!-

function str_first_in_set(line, set, ignore_clauses) result (ix_match)

use sim_utils, only: logic_option

implicit none
character(*) line, set
character(1) ch
integer ix_match, ix1, ix2, ix3
logical, optional :: ignore_clauses
logical ignore, in_single, in_double, match

! Find first character in "line" that is in "set"

ix1 = 0; ix2 = 0; ix3 = 0
in_single = .false.; in_double = .false.
ignore = logic_option(.false., ignore_clauses)

do ix_match = 1, len(line)
  ch = line(ix_match:ix_match)
  match = (index(set, ch) /= 0)

  if (.not. ignore) then
    if (match) return
    cycle
  endif

  !

  if (in_double) then
    if (ch == '"') in_double = .false.
  elseif (in_single) then
    if (ch == "'") in_single = .false.
  else
    select case (ch)
    case ('"');  in_double = .true.
    case ("'");  in_single = .true.
    case ('(');  if (match .and. ix1 == 0 .and. ix2 == 0 .and. ix3 == 0) return;  ix1 = ix1 + 1
    case ('[');  if (match .and. ix1 == 0 .and. ix2 == 0 .and. ix3 == 0) return;  ix2 = ix2 + 1
    case ('{');  if (match .and. ix1 == 0 .and. ix2 == 0 .and. ix3 == 0) return;  ix3 = ix3 + 1
    case (')');  if (match .and. ix1 == 0 .and. ix2 == 0 .and. ix3 == 0) return;  ix1 = ix1 - 1
    case (']');  if (match .and. ix1 == 0 .and. ix2 == 0 .and. ix3 == 0) return;  ix2 = ix2 - 1
    case ('}');  if (match .and. ix1 == 0 .and. ix2 == 0 .and. ix3 == 0) return;  ix3 = ix3 - 1
    end select
  endif

  if (match .and. .not. in_single .and. .not. in_double .and. &
                ix1 == 0 .and. ix2 == 0 .and. ix3 == 0 .and. index(')]}', ch) == 0) exit
enddo

ix_match = 0 ! no match found

end function
