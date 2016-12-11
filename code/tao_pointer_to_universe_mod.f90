module tao_pointer_to_universe_mod

use tao_struct
use tao_interface

!+
! Function tao_pointer_to_universe (...) result (u)
!
! Routine to set a pointer to a universe.
!
! This is an overloaded routine for the:
!  tao_pointer_to_universe_int (ix_uni) result (u)
!  tao_pointer_to_universe_str (string) result (u)
!
! Note: With a string argument, this routine can only handle single universe picks. 
! That is, it cannot handlle something like "[1,3,4]@...". To handle multiple universe picks, use:
!   tao_pick_universe
!
! Input:
!   ix_uni -- Integer: Index to the s%u(:) array
!               If ix_uni is -1 then u(s%com%default_universe) will be used.
!   string -- character(*): String in the form "<ix_uni>@..." or just "<ix_uni>".
!
! Output:
!   string -- character(*): String with universe prefix stripped off.
!   u      -- Tao_universe_struct, pointer: Universe pointer.
!               u will be nullified if there is an error and an error message will be printed.
!-

interface tao_pointer_to_universe
  module procedure tao_pointer_to_universe_int
  module procedure tao_pointer_to_universe_str
end interface

private tao_pointer_to_universe_int, tao_pointer_to_universe_str

contains

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Function tao_pointer_to_universe_int (ix_uni) result (u)
!
! Overloaded by tao_pointer_to_universe. See this routine for more details.
!-

function tao_pointer_to_universe_int (ix_uni) result(u)

implicit none

type (tao_universe_struct), pointer :: u
integer ix_uni, ix_u
character(*), parameter :: r_name = 'tao_pointer_to_universe_int'

!

ix_u = tao_universe_number(ix_uni)

if (ix_u < lbound(s%u, 1) .or. ix_u > ubound(s%u, 1)) then
  call out_io (s_fatal$, r_name, 'UNIVERSE INDEX OUT OF RANGE: \I0\ ', ix_u)
  nullify (u)
  return
endif

u => s%u(ix_u)

end function tao_pointer_to_universe_int

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Function tao_pointer_to_universe_str (string) result (u)
!
! Overloaded by tao_pointer_to_universe. See this routine for more details.
!-

function tao_pointer_to_universe_str (string) result(u)

implicit none

type (tao_universe_struct), pointer :: u
integer ix, ix_u
character(*) string
character(*), parameter :: r_name = 'tao_pointer_to_universe_str'

!

nullify(u)

ix = index(string, '@')
if (ix == 0) then
  if (.not. is_integer(string)) then
    call out_io (s_fatal$, r_name, 'MALFORMED UNIVERSE STRING')
    return
  endif
  read (string, *) ix_u
  string = ''
else
  if (.not. is_integer(string(1:ix-1))) then
    call out_io (s_fatal$, r_name, 'MALFORMED UNIVERSE STRING')
    return
  endif
  read (string(1:ix-1), *) ix_u
  string = string(ix+1:)
endif

ix_u = tao_universe_number(ix_u)

if (ix_u < lbound(s%u, 1) .or. ix_u > ubound(s%u, 1)) then
  call out_io (s_fatal$, r_name, 'UNIVERSE INDEX OUT OF RANGE: \I0\ ', ix_u)
  return
endif

u => s%u(ix_u)

end function tao_pointer_to_universe_str

end module
