!+
! Subroutine do_mode_flip (ele, ele_flip)
!
! Subroutine to mode flip the twiss_parameters of an element
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
!     ele_flip    -- Ele_struct: Flipped element
!     bmad_status -- Status_struct: Status (in common block)
!            %ok      -- Set true if flip was done.
!-

!$Id$
!$Log$
!Revision 1.4  2002/02/23 20:32:14  dcs
!Double/Single Real toggle added
!
!Revision 1.3  2002/01/08 21:44:38  dcs
!Aligned with VMS version  -- DCS
!
!Revision 1.2  2001/09/27 18:31:51  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"



subroutine do_mode_flip (ele, ele_flip)

  use bmad

  implicit none

  type (ele_struct)  ele, ele2, ele_flip

  real(rdef) c_conj(2,2), c_mat(2,2), gamma_flip, v_mat(4,4), v_inv_mat(4,4)

  logical init_needed / .true. /

!

  if (init_needed) then
    ele2%mat6 = 0
    init_needed = .false.
  endif

  if (ele%gamma_c >= 1.0) then
    bmad_status%ok = .false.
    if (bmad_status%type_out)  &
            type *, 'ERROR IN DO_MODE_FLIP: CANNOT MODE FLIP ELEMENT'
    if (bmad_status%exit_on_error) call err_exit
    return
  endif

  gamma_flip = sqrt(1 - ele%gamma_c**2)

  call mat_symp_conj (ele%c_mat, c_conj, 2, 2)
  ele2%mat6(1:2,1:2) = -c_conj / gamma_flip
  ele2%mat6(3:4,3:4) = ele%c_mat/gamma_flip

  call twiss_decoupled_propagate (ele, ele2)

  ele_flip = ele
  ele_flip%mode_flip = .not. ele%mode_flip
  ele_flip%c_mat = -ele%c_mat * ele%gamma_c / gamma_flip
  ele_flip%gamma_c = gamma_flip
  ele_flip%x = ele2%x
  ele_flip%y = ele2%y

  call make_v_mats (ele_flip, v_mat, v_inv_mat)
  call mobius_twiss_calc (ele_flip, v_mat)

end subroutine
