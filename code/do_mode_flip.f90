!+
! Subroutine do_mode_flip (ele)
!
! Subroutine to mode flip the twiss_parameters of an element.
! That is, the normal mode associated with ele%a is transfered to ele%b
! and the normal mode associated with ele%b is trnsfered to ele%a.
!
! Modules Needed:
!   use bmad
!
! Input:
!   ele          -- Ele_struct: Starting Element
!   bmad_status  -- Status_struct:
!      %type_out   -- If true then type out when an error occurs.
!
! Output:
!     ele         -- Ele_struct: Flipped element
!     bmad_status -- Status_struct: Status (in common block)
!            %ok      -- Set true if flip was done.
!-

#include "CESR_platform.inc"

subroutine do_mode_flip (ele)

use bmad_struct
use bmad_interface, except_dummy => do_mode_flip

implicit none

type (ele_struct)  ele
type (twiss_struct) a

real(rp) cg_mat(2,2), cg_conj(2,2), gamma_flip

! Check that a flip can be done.

if (ele%gamma_c >= 1.0) then
bmad_status%ok = .false.
if (bmad_status%type_out)  &
            print *, 'ERROR IN DO_MODE_FLIP: CANNOT MODE FLIP ELEMENT'
if (bmad_status%exit_on_error) call err_exit
return
endif

! Do the flip

gamma_flip = sqrt(1 - ele%gamma_c**2)

cg_mat = ele%c_mat / gamma_flip
call mat_symp_conj (cg_mat, cg_conj)

a = ele%a
call twiss1_propagate (ele%b, cg_mat, 0.0_rp, ele%a)
call twiss1_propagate (a, -cg_conj, 0.0_rp, ele%b)

ele%mode_flip = .not. ele%mode_flip
ele%c_mat = -cg_mat * ele%gamma_c
ele%gamma_c = gamma_flip

end subroutine
