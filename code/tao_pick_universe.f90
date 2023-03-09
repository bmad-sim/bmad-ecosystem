!+
! Subroutine tao_pick_universe (name_in, name_out, picked, err, ix_uni, explicit_uni, dflt_uni, pure_uni)
!
! Subroutine to pick what universe the data name is comming from.
! Examples:
!   "*@..."           -- Choose all universes.
!   "3@..."           -- Choose universe 3. 
!   "[1:30,34]@..."   -- Choose universes 1 through 30 and 34.
!   "-1@"             -- Default universe.
!   No "@" in name and universe not specified
!                     -- Choose dflt_uni if it is present. Otherwise pick s%global%default_universe.
!
! Normally a universe prefix must be followed by a "@" sign. For example, "3@". 
! Something like "3" without a "@" sign is not considered a universe prefix. 
! Exception: If pure_uni is set to False, a universe prefix does not need a "@" sign to follow.
!
! Also see:
!   tao_pointer_to_universe
!
! Input:
!   name_in      -- character(*): data name with possible universe spec.
!   dflt_uni     -- integer, optional: Default universe to use. Set to -1 if explicit universe is required.
!   pure_uni     -- logical, optional: Default is False. See above
!
! Output:
!   name_out     -- character(*): name_in without any "n@" beginning.
!   picked(:)    -- logical, allocatable: Array showing picked universes.
!                     The array will be resized if necessary.
!   err          -- logical: Set True if an error is detected.
!   ix_uni       -- integer, optional: Set to the picked universe with the highest index.
!   explicit_uni -- logical, optional: Set True if name_in has explicit universe "n@" specification.
!-

subroutine tao_pick_universe (name_in, name_out, picked, err, ix_uni, explicit_uni, dflt_uni, pure_uni)

use tao_interface, except_dummy => tao_pick_universe

implicit none

character(*) name_in, name_out
character(*), parameter :: r_name = 'tao_pick_universe'
character(40) uni
character(1) ch

integer, optional :: ix_uni, dflt_uni
integer i, j, ix, n, ios, iu, num, ic, iu_dflt

logical, allocatable :: picked(:)
logical, allocatable :: p(:)
logical err
logical, optional :: explicit_uni, pure_uni

! Init

call re_allocate2 (picked, lbound(s%u, 1), ubound(s%u, 1))
call re_allocate2 (p, -1, ubound(s%u, 1))

err = .false.
name_out = trim(name_in)
picked = .false.
p = .false.
if (present(ix_uni)) ix_uni = -1

! Look for universe substring.

if (logic_option(.false., pure_uni)) then
  do i = 1, len(name_out)
    ch = name_out(i:i)
    if (ch == ' ' .and. name_out(1:1) /= '[') exit  ! Blanks allowed in bracket construct.
    if (index('[],-*:0123456789', ch) == 0) exit
    j = max(1, i-1)
    if (name_out(j:j) == ']') exit
  enddo

  uni = name_out(1:i-1)

  if (ch == '@') then
    name_out = name_out(i+1:)
  elseif (i == len(name_out) + 1) then
    name_out = ''
  else
    name_out = name_out(i:)
  endif

  if (ch /= '@' .and. uni /= ' ' .and. name_out /= ' ') then
    call out_io (s_error$, r_name, 'MALFORMED UNIVERSE SPEC')
    err = .true.
    return
  endif

else
  ix = tao_uni_ampersand_index(name_out)
  if (ix == 0) then
    uni = ''
  else
    uni = name_out(1:ix-1)
    name_out = name_out(ix+1:)
  endif
endif

!

if (present(explicit_uni)) explicit_uni = (uni /= '')

if (uni == '') then
  iu_dflt = integer_option(s%global%default_universe, dflt_uni)
  if (iu_dflt < 0) then
    call out_io (s_error$, r_name, 'NO UNIVERSE NUMBER GIVEN')
    err = .true.
    return
  endif

  picked (iu_dflt) = .true.
  if (present(ix_uni)) ix_uni = iu_dflt
  return
endif


! Strip off '[' and ']'

ix = len_trim(uni)
if (uni(1:1) == '[' .and. uni(ix:ix) == ']') uni = uni(2:ix-1)

if (uni == '*') then
  picked = .true.
  if (present(ix_uni)) ix_uni = lbound(s%u, 1)
  return
endif

! "show var" uses a blank universe to correspond to the default universe.

if (uni == '') then
  iu_dflt = integer_option(s%global%default_universe, dflt_uni)
  picked (iu_dflt) = .true.
  if (present(ix_uni)) ix_uni = iu_dflt
  return
endif

call location_decode (uni, p, lbound(p, 1), num)
if (num == -1) then
  call out_io (s_error$, r_name, 'BAD UNIVERSE NUMBER: ' // uni)
  err = .true.
  return
endif

do i = lbound(p, 1), ubound(p, 1)
  if (.not. p(i)) cycle

  iu = tao_universe_index(i)
  if (iu < lbound(s%u, 1) .or. iu > ubound(s%u, 1)) then
    call out_io (s_error$, r_name, 'NUMBER DOES NOT CORRESPOND TO A UNIVERSE: ' // uni)
    err = .true.
    return
  endif
  picked(iu) = .true.
  if (present(ix_uni)) ix_uni = iu
enddo

end subroutine tao_pick_universe
