!+
! Subroutine dispersion_to_orbit (ele, disp_orb)
!
! Subroutine to make an orbit vector proportional to the dispersion.
!
! Note that, despite the names, ele%x and ele%y are the "a" and "b" mode
! Twiss parameters. Thus, for example, with coupling ele%x%eta is not simply
! proportional to disp_orb%vec(1).
!
! Note: to calculate ele%x%eta, ele%c_mat, etc. you typically need to have 
! called twiss_propagate_all.
!
! Modules needed:
!   use bamd
!
! Input:
!   ele      -- Ele_struct: Element containing the dispersion info.
!     %x%eta  -- "a" mode eta.
!     %x%etap -- "a" mode etap.
!     %y%eta  -- "b" mode eta.
!     %y%etap -- "b" mode etap.
!     %c_mat   -- coupling matrix needed by the routine. 
!     %gamma_c -- coupling gamma factor.
!
! Output:
!   disp_orb -- Coord_struct: dispersive orbit normalized so that
!                 disp_orb%vec(6) = 1.
!-

#include "CESR_platform.inc"

Subroutine dispersion_to_orbit (ele, disp_orb)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ele_struct), intent(in) :: ele
  type (coord_struct), intent(out) :: disp_orb

  real v_mat(4,4), v_inv_mat(4,4)
  
!

  disp_orb%vec = (/ ele%x%eta, ele%x%etap, ele%y%eta, ele%y%etap, &
                                                       0.0_rdef, 1.0_rdef /)

  if (all(ele%c_mat == 0)) return

  call ele_to_v_mats (ele, v_mat, v_inv_mat)
  disp_orb%vec(1:4) = matmul (v_mat, disp_orb%vec(1:4))

end subroutine
