!+
! Subroutine orbit_to_dispersion (orb_diff, ele)
!
! Subroutine to take an orbit vector difference and calculate the dispersion.
!
! Note that, despite the names, ele%x and ele%y are the "a" and "b" mode
! Twiss parameters. Thus, for example, with coupling ele%x%eta is not simply
! proportional to orb_diff%x%pos.
!
! Modules needed:
!   use bamd
!
! Input:
!   orb_diff -- Coord_struct: orbit difference between 2 orbits with differing
!               energies. That is, orb_diff%vec(6) must be nonzero.
!   ele      -- Ele_struct: Element containing the coupling info.
!     %c_mat   -- coupling matrix needed by the routine. You need to have
!                 called twiss_propagate_all to have computed this.
!     %gamma_c -- coupling gamma factor (gets computed with %c_mat)
!
! Output:
!   ele   -- Ele_struct: Element containing the dispersion
!     %x%eta  -- "a" mode eta.
!     %x%etap -- "a" mode etap.
!     %y%eta  -- "b" mode eta.
!     %y%etap -- "b" mode etap.
!-

#include "CESR_platform.inc"

Subroutine orbit_to_dispersion (orb_diff, ele)

  use bmad_struct
  use bmad_interface

  implicit none

  type (coord_struct), intent(in) :: orb_diff
  type (ele_struct) :: ele

  real(rdef) v_mat(4,4), v_inv_mat(4,4), disp_vec(4)
  
!

  if (orb_diff%vec(6) == 0) then
    print *, 'ERROR IN ORBIT_TO_DISPERSION: ORB_DIFF%VEC(6) = 0 !'
    call err_exit
  endif

  call ele_to_v_mats (ele, v_mat, v_inv_mat)
  disp_vec = matmul (v_inv_mat, orb_diff%vec(1:4)) / orb_diff%vec(6)
  ele%x%eta  = disp_vec(1)
  ele%x%etap = disp_vec(2)
  ele%y%eta  = disp_vec(3)
  ele%y%etap = disp_vec(4)

end subroutine
