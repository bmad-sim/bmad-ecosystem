!+
! Subroutine do_mode_flip (ele, err_flag)
!
! Subroutine to mode flip the twiss_parameters of an element.
! That is, the normal mode associated with ele%a is transfered to ele%b
! and the normal mode associated with ele%b is trnsfered to ele%a.
!
! Input:
!   ele          -- Ele_struct: Starting Element
!
! Output:
!     ele         -- Ele_struct: Flipped element
!     err_flag    -- Logical, optional: Set True if there is an error. False otherwise.
!-

subroutine do_mode_flip (ele, err_flag)

use bmad_interface, except_dummy => do_mode_flip

implicit none

type (ele_struct)  ele
type (twiss_struct) a

real(rp) cg_mat(2,2), cg_conj(2,2), gamma_flip
logical, optional :: err_flag
logical err

character(12), parameter :: r_name = 'do_mode_flip'

! Check that a flip can be done.

if (present(err_flag)) err_flag = .true.

if (ele%gamma_c >= 1.0) then
  call out_io (s_fatal$, r_name, 'CANNOT MODE FLIP ELEMENT')
  if (global_com%exit_on_error) call err_exit
  return
endif

! Do the flip

gamma_flip = sqrt(1 - ele%gamma_c**2)

cg_mat = ele%c_mat / gamma_flip
cg_conj = mat_symp_conj (cg_mat)

a = ele%a
call twiss1_propagate (ele%b, cg_mat, drift$, 0.0_rp, ele%a, err)
if (err) return
call twiss1_propagate (a, -cg_conj, drift$, 0.0_rp, ele%b, err)
if (err) return

ele%mode_flip = .not. ele%mode_flip
ele%c_mat = -cg_mat * ele%gamma_c
ele%gamma_c = gamma_flip

if (present(err_flag)) err_flag = .false.

end subroutine
