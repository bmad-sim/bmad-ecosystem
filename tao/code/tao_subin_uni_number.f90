!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function tao_subin_uni_number (name_in, ix_uni, name_out) result (ok)
!
! Routine to mangle a file name based upon the universe number ix_uni.
! If there is a "#" in the name and the number of universes is more than one,
! the universe number will be substituted.
! If the number of universes is 1, the "#" will be stripped from the name.
!
! Input:
!   name_in -- Character(*): Input name with "#" character
!   ix_uni  -- Integer: Universe index.
!
! Output:
!   name_out -- Character(*): Output name.
!   ok       -- Logical: False if multiple universes and no "#" in name_in.
!                 True otherwise.
!-

function tao_subin_uni_number (name_in, ix_uni, name_out) result (ok)

use tao_struct

implicit none

character(*) name_in, name_out
character(28) :: r_name = 'tao_subin_uni_number'
integer ix, ix_uni
logical ok

!

ok = .true.
name_out = name_in

ix = index(name_out, '#')
if (size(s%u) > 1 .and. ix == 0) then
  call out_io (s_info$, r_name, 'FILE NAME DOES NOT HAVE A "#" CHARACTER!', &
    ' YOU NEED THIS TO GENERATE A UNIQUE FILE NAME FOR EACH UNIVERSE!')
  ok = .false.
endif

if (ix /= 0) write (name_out, '(a, i0, a)') name_out(1:ix-1), ix_uni, trim(name_out(ix+1:))

end function tao_subin_uni_number
