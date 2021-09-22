!+
! Subroutine tao_universe_locator (name_in, name_out, unis, err, explicit_uni, dflt_uni)
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
!   unis(:)      -- tao_universe_pointer_struct, allocatable: Array of pointers to picked universes.
!                     The array will be resized if necessary.
!   err          -- logical: Set True if an error is detected.
!   explicit_uni -- logical, optional: Set True if name_in has explicit universe "n@" specification.
!-

subroutine tao_universe_locator (name_in, name_out, unis, err, explicit_uni, dflt_uni)

use tao_interface, except_dummy => tao_universe_locator

implicit none

type (tao_universe_pointer_struct), allocatable :: unis(:)

character(*) name_in, name_out
character(*), parameter :: r_name = 'tao_pick_universe'
character(40) uni

integer, optional :: dflt_uni
integer i, ix, n, ios, iu, num, ic, iu_dflt
integer, allocatable :: arr(:)

logical err
logical, optional :: explicit_uni

! Init

err = .false.
name_out = name_in

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

  call tao_allocate_uni_pointers(unis, 1)
  unis(1)%u => s%u(iu_dflt)
  return
endif

! Here when "@" is found...

uni = name_in(:ix-1)
name_out = name_out(ix+1:)

! Strip off '[' and ']'

if (ix > 2) then
  if (uni(1:1) == '[' .and. uni(ix-1:ix-1) == ']') uni = uni(2:ix-2)
endif

if (uni == '*') then
  call tao_allocate_uni_pointers (unis, size(s%u))
  do i = 1, ubound(s%u, 1)
    unis(i)%u => s%u(i)
  enddo
  return
endif

! "show var" uses a blank universe to correspond to the default universe.

if (uni == '') then
  iu_dflt = integer_option(s%global%default_universe, dflt_uni)
  call tao_allocate_uni_pointers(unis, 1)
  unis(1)%u => s%u(iu_dflt)
  return
endif

call pointer_to_locations (uni, arr, num, 1, ubound(s%u, 1))
if (num == -1) then
  call out_io (s_error$, r_name, 'BAD UNIVERSE NUMBER: ' // uni)
  err = .true.
  return
endif

call tao_allocate_uni_pointers (unis, num)

do i = 1, num
	unis(i)%u => s%u(arr(i))
enddo

!-------------------------------------
contains

subroutine tao_allocate_uni_pointers (unis, n)

type (tao_universe_pointer_struct), allocatable :: unis(:)
integer n

!

if (.not. allocated(unis)) then
  allocate(unis(n))
  return
endif

if (size(unis) /= n) then
  deallocate(unis)
  allocate(unis(n))
endif

end subroutine tao_allocate_uni_pointers

end subroutine tao_universe_locator
