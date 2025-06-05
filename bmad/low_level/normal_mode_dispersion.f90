!+
! Subroutine normal_mode_dispersion(ele, reverse) 
!
! Routine to calculate the normal mode dispersion from the x,y dispersions.
! Or vice versa if reverse = True.
!
! Input:
!   ele       -- ele_struct: Element whose dispersions are to be adjusted.
!   reverse   -- logical, optional: Default is False. If True, calculate the x,y dispersions 
!                   from the normal mode ones.
!
! Output:
!   ele       -- ele_struct: Element with adjusted dispersions.
!-

subroutine normal_mode_dispersion(ele, reverse) 

use bmad_interface, dummy => normal_mode_dispersion

implicit none

type (ele_struct) ele

real(rp) v_mat(4,4), v_inv_mat(4,4), eta_vec(4), orb_vec(4), rel_p
logical, optional :: reverse

! Normal mode to x,y
! Note: rel_p = 0 can happen at or near an e_gun element.

if (logic_option(.false., reverse)) then

  call make_v_mats (ele, v_mat = v_mat)
  eta_vec = [ele%a%eta, ele%a%etap, ele%b%eta, ele%b%etap]
  eta_vec = matmul (v_mat, eta_vec)
  orb_vec = ele%map_ref_orb_out%vec(1:4)
  rel_p = 1 + ele%map_ref_orb_out%vec(6)

  ele%x%eta     = eta_vec(1)
  ele%x%etap    = eta_vec(2)

  ele%y%eta     = eta_vec(3)
  ele%y%etap    = eta_vec(4)

  if (rel_p < 1e-3) then
    ele%x%deta_ds = 0
    ele%y%deta_ds = 0
  else
    ele%x%deta_ds = eta_vec(2) / rel_p - orb_vec(2) / rel_p**2
    ele%y%deta_ds = eta_vec(4) / rel_p - orb_vec(4) / rel_p**2
  endif

! x,y to normal mode

else
  call make_v_mats (ele, v_inv_mat = v_inv_mat)
  eta_vec = [ele%x%eta, ele%x%etap, ele%y%eta, ele%y%etap]
  eta_vec = matmul (v_inv_mat, eta_vec)
  orb_vec = matmul(v_inv_mat, ele%map_ref_orb_out%vec(1:4))
  rel_p = 1 + ele%map_ref_orb_out%vec(6)

  ele%a%eta     = eta_vec(1)
  ele%a%etap    = eta_vec(2)
  ele%b%eta     = eta_vec(3)
  ele%b%etap    = eta_vec(4)

  if (rel_p < 1e-3) then
    ele%a%deta_ds = 0
    ele%b%deta_ds = 0
  else
    ele%a%deta_ds = eta_vec(2) / rel_p - orb_vec(2) / rel_p**2
    ele%b%deta_ds = eta_vec(4) / rel_p - orb_vec(4) / rel_p**2
  endif

  eta_vec = [ele%x%deta_dpz, ele%x%detap_dpz, ele%y%deta_dpz, ele%y%detap_dpz]
  eta_vec = matmul (v_inv_mat, eta_vec)

  ele%a%deta_dpz     = eta_vec(1)
  ele%a%detap_dpz    = eta_vec(2)
  ele%b%deta_dpz     = eta_vec(3)
  ele%b%detap_dpz    = eta_vec(4)
endif

end subroutine
