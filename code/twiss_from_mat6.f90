!+
! Subroutine twiss_from_mat6 (mat6, map0, ele, stable, growth_rate)
!
! Subroutine to calculate the Twiss parameters from a 1-turn matrix.
! Note: The 1-turn matrix needs to be formed with the RF turned off.
!
! Modules needed:
!   use bmad
!
! Input:
!   mat6(6,6)   -- Real(rp): 6x6 matrix (linear part) of the 1-turn map.
!   map0(6)     -- Real(rp): 0th order part of the 1-turn map.
!
! Output:
!   ele         -- Ele_struct: Structure holding the Twiss parameters.
!     %a           -- X Twiss parameters at the start of the lat.
!     %a%phi       -- Fractional part of the tune in radians.
!     %b           -- Y Twiss parameters at the start of the lat.
!     %b%phi       -- Fractional part of the tune in radians.
!     %c_mat       -- Coupling matrix.
!   stable      -- Set true or false.
!   growth_rate -- unstable growth rate (= 0 if stable)
! 
!   bmad_status -- BMAD Common block status structure
!     %ok          -- Logical: .True. if everything is OK,
!     %status      -- Integer: Calculation status.
!                        See MAT_SYMP_DECOUPLE for for more info
!-

#include "CESR_platform.inc"

subroutine twiss_from_mat6 (mat6, map0, ele, stable, growth_rate)

  use bmad_struct
  use bmad_interface, except => twiss_from_mat6

  implicit none

  type (ele_struct) :: ele
  real(rp), intent(in) :: mat6(:,:), map0(:)
  real(rp), intent(out) :: growth_rate
  logical, intent(out) :: stable

  real(rp) mat4(4,4), eta_vec(4), vec(4), orb(4)
  real(rp) u(4,4), v(4,4), ubar(4,4), vbar(4,4), g(4,4), v_inv(4,4)
  real(rp) det, rate1, rate2
  real(rp) :: tol = 1.0e-3

  integer i

  character(20) :: r_name = 'twiss_from_mat6'

! init

  mat4 = mat6(1:4, 1:4)

!

  call mat_symp_decouple (mat4, tol, bmad_status%status, u, v, ubar, &
                                             vbar, g, ele%a, ele%b, .false.)

  if (bmad_status%status /= ok$) then
    if (bmad_status%type_out) then
       call out_io (s_error$, r_name, 'ERROR IN TWISS_FROM_MAT6: BAD 1-TURN MATRIX: ' // &
            status_name(bmad_status%status),'TWISS PARAMETERS NOT COMPUTED')
    endif
    if (bmad_status%status == non_symplectic$) then
      rate1 = 10.0
      rate2 = 10.0
      rate1 = max(rate1, maxval(abs(mat4)))
    else
      rate1 = sqrt(max(abs(u(1,1) + u(2,2)) - 2, 0.0_rp))
      rate2 = sqrt(max(abs(u(3,3) + u(4,4)) - 2, 0.0_rp))
    endif
    growth_rate = max(rate1, rate2)
    stable = .false.
    return
  else
    growth_rate = 0          ! no growth
    stable = .true.          ! stable lat
  endif

! here if everything normal so load twiss parameters

  if (ele%a%beta /= 0 .and. ele%b%beta /= 0) then
    ele%mode_flip = .false.
    ele%c_mat = v(1:2,3:4)
    call mat_det (ele%c_mat, det)
    ele%gamma_c = sqrt(1-det)
  endif

! Compute normal mode and lab dispersion.

  mat4 = -mat4
  forall (i = 1:4) mat4(i,i) = mat4(i,i) + 1
  call mat_inverse (mat4, mat4)

  orb = matmul(mat4, map0(1:4))
  ele%ref_orb(1:4) = orb

  vec(1) = mat6(1,6) + mat6(1,2) * orb(2) + mat6(1,4) * orb(4)
  vec(2) = mat6(2,6) - mat6(2,1) * orb(1) - mat6(2,3) * orb(3) - map0(2)
  vec(3) = mat6(3,6) + mat6(3,2) * orb(2) + mat6(3,4) * orb(4)
  vec(4) = mat6(4,6) - mat6(4,1) * orb(1) - mat6(4,3) * orb(3) - map0(4)
  eta_vec = matmul(mat4, vec)

  ele%x%eta  = eta_vec(1)
  ele%x%etap = eta_vec(2)
  ele%y%eta  = eta_vec(3)
  ele%y%etap = eta_vec(4)

  call mat_symp_conj (v, v_inv)
  eta_vec = matmul(v_inv, eta_vec)
  ele%a%eta  = eta_vec(1)
  ele%a%etap = eta_vec(2)
  ele%b%eta  = eta_vec(3)
  ele%b%etap = eta_vec(4)

  bmad_status%ok = .true.

end subroutine
