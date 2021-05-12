!+
! Subroutine tao_pick_universe (name_in, name_out, picked, err, ix_uni, explicit_uni, dflt_uni)
!
! Subroutine to pick what universe the data name is comming from.
! Examples:
!   "*@..."           -- Choose all universes.
!   "3@..."           -- Choose universe 3. 
!   "[1:30,34]@..."   -- Choose universes 1 through 30 and 34
!   No "@" in name    -- Choose dflt_uni if it is present. Otherwise pick s%global%default_universe.
!
! Also see:
!   tao_pointer_to_universe
!
! Input:
!   name_in      -- character(*): data name with possible universe spec.
!   dflt_uni     -- integer, optional: Default universe to use. Set to -1 if explicit universe is required.
!
! Output:
!   name_out     -- character(*): name_in without any "n@" beginning.
!   picked(:)    -- logical, allocatable: Array showing picked universes.
!                     The array will be resized if necessary.
!   err          -- logical: Set True if an error is detected.
!   ix_uni       -- integer, optional: Set to the picked universe with the highest index.
!   explicit_uni -- logical, optional: Set True if name_in has explicit universe "n@" specification.
!-

subroutine tao_pick_universe (name_in, name_out, picked, err, ix_uni, explicit_uni, dflt_uni)

use tao_interface, except_dummy => tao_pick_universe

implicit none

character(*) name_in, name_out
character(*), parameter :: r_name = 'tao_pick_universe'
character(40) uni

integer, optional :: ix_uni, dflt_uni
integer i, ix, n, ios, iu, num, ic, iu_dflt

logical, allocatable :: picked(:)
logical, allocatable :: p(:)
logical err
logical, optional :: explicit_uni

! Init

call re_allocate2 (picked, lbound(s%u, 1), ubound(s%u, 1))
call re_allocate2 (p, -1, ubound(s%u, 1))

err = .false.
picked = .false.
p = .false.
if (present(ix_uni)) ix_uni = -1

! No "@" then simply choose s%global%default_universe.

ix = index (name_in, '@')
ic = index (name_in, '::')

if (present(explicit_uni)) explicit_uni = (ix /= 0)

if (ix == 0 .or. (ic /= 0 .and. ix > ic)) then
  iu_dflt = integer_option(s%global%default_universe, dflt_uni)
  if (iu_dflt < 0) then
    call out_io (s_error$, r_name, 'NO UNIVERSE NUMBER GIVEN')
    err = .true.
    return
  endif

  picked (iu_dflt) = .true.
  name_out = name_in
  if (present(ix_uni)) ix_uni = iu_dflt
  return
endif

! Here when "@" is found...

uni = name_in(:ix-1)
name_out = name_in(ix+1:)

! Strip off '[' and ']'

if (ix > 2) then
  if (uni(1:1) == '[' .and. uni(ix-1:ix-1) == ']') uni = uni(2:ix-2)
endif

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
  iu = tao_universe_number (i)
  if (iu < lbound(s%u, 1) .or. iu > ubound(s%u, 1)) then
    call out_io (s_error$, r_name, 'NUMBER DOES NOT CORRESPOND TO A UNIVERSE: ' // uni)
    err = .true.
    return
  endif
  picked(iu) = .true.
  if (present(ix_uni)) ix_uni = iu
enddo

end subroutine tao_pick_universe
