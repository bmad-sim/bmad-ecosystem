!+
! Subroutine ele_to_fibre_hook (ele, ptc_fibre, param, use_offsets, err_flag)
!
! Routine that can be customized for creating a PTC fibre from a Bmad element.
! When setup, this routine is called by ele_to_fibre.
! To use, see the Bmad manual.
!
! Input:
!   ele             -- Ele_struct: Bmad element.
!   param           -- lat_param_struct: 
!   use_offsets     -- Logical: Does ptc_fibre include element offsets, pitches and tilt?
!
! Output:
!   err_flag        -- logical: Set True if setup OK. False otherwise.
!   ptc_fibre       -- Fibre: PTC fibre element.
!-

subroutine ele_to_fibre_hook (ele, ptc_fibre, param, use_offsets, err_flag)

use bmad
use s_family, only: work, suntao, assignment(=)  ! PTC

implicit none

type (ele_struct) ele
type (fibre) ptc_fibre
type (lat_param_struct) param

logical use_offsets, err_flag

character(*), parameter :: r_name = 'ele_to_fibre_hook'

! Stuff specific to Suntao custom tracking

type (work) wrk
real(rp) sfactor

! This is for Suntao's custom tracking through a Cesr wiggler.

if (ele%key /= wiggler$ .or. .not. associated(ele%r)) return

if (ubound(ele%r, 2) /= 24) then
  call out_io (s_fatal$, r_name, 'BAD R_CUSTOM ARRAY SETUP!')
  if (global_com%exit_on_error) call err_exit
  return
endif

ptc_fibre%mag%wi%w%ex(1:24) = ele%r(1, 1:24, 0) * ele%value(polarity$)
ptc_fibre%mag%wi%w%ey(1:24) = ele%r(2, 1:24, 0) * ele%value(polarity$)

wrk = ptc_fibre
sfactor = wrk%gamma0i*suntao/ptc_fibre%mag%l
ptc_fibre%mag%wi%w%ex(1:22:3) = ptc_fibre%mag%wi%w%ex(1:22:3) * sfactor
ptc_fibre%mag%wi%w%ey(1:22:3) = ptc_fibre%mag%wi%w%ey(1:22:3) * sfactor

ptc_fibre%magp%wi%w%ex = ptc_fibre%mag%wi%w%ex
ptc_fibre%magp%wi%w%ey = ptc_fibre%mag%wi%w%ey

ptc_fibre%mag%p%permfringe = 2
ptc_fibre%magp%p%permfringe = 2

end subroutine ele_to_fibre_hook
