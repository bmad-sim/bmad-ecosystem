!+
! Subroutine transfer_twiss (ele_in, ele_out, reverse)
!
! Routine to transfer the twiss parameters from one element to another.
!
!
! Input:
!   ele_in   -- ele_struct: Element with existing Twiss parameters.
!   reverse  -- logical, optional: Reverse alpha and coupling as if particle is going in 
!                 the reversed direction? Default is False.
!
! Output:
!   ele_out  -- ele_struct: Element receiving the Twiss parameters.
!-

subroutine transfer_twiss (ele_in, ele_out, reverse)

use bmad_struct

implicit none

type (ele_struct) ele_in, ele_out
logical, optional :: reverse

!

ele_out%x         = ele_in%x
ele_out%y         = ele_in%y
ele_out%a         = ele_in%a
ele_out%b         = ele_in%b
ele_out%z         = ele_in%z
ele_out%c_mat     = ele_in%c_mat
ele_out%gamma_c   = ele_in%gamma_c
ele_out%mode_flip = ele_in%mode_flip

if (logic_option(.false., reverse)) then
  ele_out%x%etap    = -ele_in%x%etap
  ele_out%x%deta_ds = -ele_in%x%deta_ds
  ele_out%y%etap    = -ele_in%y%etap
  ele_out%y%deta_ds = -ele_in%y%deta_ds
  ele_out%a%alpha   = -ele_in%a%alpha
  ele_out%a%etap    = -ele_in%a%etap
  ele_out%a%deta_ds = -ele_in%a%deta_ds
  ele_out%b%alpha   = -ele_in%b%alpha
  ele_out%b%etap    = -ele_in%b%etap
  ele_out%b%deta_ds = -ele_in%b%deta_ds
  ele_out%z%alpha   = -ele_in%z%alpha
  ele_out%z%etap    = -ele_in%z%etap
  ele_out%z%deta_ds = -ele_in%z%deta_ds
  ele_out%c_mat     = -mat_symp_conj(ele_in%c_mat)
endif

end subroutine transfer_twiss

