!+
! Subroutine twiss_from_mat6 (mat6, orb0, ele, stable, growth_rate, status, type_out)
!
! Subroutine to calculate the Twiss parameters from a 1-turn matrix.
! Note: The 1-turn matrix needs to be formed with the RF turned off.
!
! Modules needed:
!   use bmad
!
! Input:
!   mat6(6,6)   -- Real(rp): 6x6 matrix (linear part) of the 1-turn map.
!   orb0(6)     -- Real(rp): Initial orbit around which the map is made.
!   type_out    -- Logical: If True then print a message when calculation fails.
!
! Output:
!   ele         -- Ele_struct: Structure holding the Twiss parameters.
!     %a           -- X Twiss parameters at the start of the lat.
!     %a%phi       -- Fractional part of the tune in radians.
!     %b           -- Y Twiss parameters at the start of the lat.
!     %b%phi       -- Fractional part of the tune in radians.
!     %c_mat       -- Coupling matrix.
!   stable      -- Logical: Set true or false.
!   growth_rate -- Real(rp): Unstable growth rate (= 0 if stable)
!   status      -- Integer: Calculation status: ok$, in_stop_band$, or unstable$.
!-

subroutine twiss_from_mat6 (mat6, orb0, ele, stable, growth_rate, status, type_out)

use bmad_interface, except_dummy => twiss_from_mat6

implicit none

type (ele_struct) :: ele

real(rp) :: mat6(:,:), orb0(:)
real(rp) :: growth_rate
real(rp) mat4(4,4), eta_vec(4), vec(4), rel_p
real(rp) u(4,4), v(4,4), ubar(4,4), vbar(4,4), g(4,4), v_inv(4,4)
real(rp) rate1, rate2, symp_err
real(rp) :: symp_tol = 1.0d-3

integer :: status
integer i

logical :: stable, type_out

character(20) :: r_name = 'twiss_from_mat6'

! init

mat4 = mat6(1:4, 1:4)

symp_err = mat_symp_error(mat4)

if (mat_symp_error(mat4) > 1) then
  call out_io (s_error$, r_name, 'BAD 1-TURN MATRIX. NON_SYMPLECTIC WITH SYMPLECTIC ERROR OF: \f8.1\ ', &
                                 'TWISS PARAMETERS NOT COMPUTED', r_array = [symp_err])
  status = non_symplectic$
  return
endif

if (mat_symp_error(mat4) > symp_tol) then
  call out_io (s_warn$, r_name, '1-TURN MATRIX MARGINALLY SYMPLECTIC WITH SYMPLECTIC ERROR OF: \f8.1\ ', &
                                r_array = [symp_err])
endif

!

call mat_symp_decouple (mat4, status, u, v, ubar, vbar, g, ele%a, ele%b, ele%gamma_c, .false.)

if (status /= ok$) then
  if (type_out) then
     call out_io (s_error$, r_name, 'BAD 1-TURN MATRIX: ' // &
                      matrix_status_name(status), 'TWISS PARAMETERS NOT COMPUTED')
  endif
  if (status == non_symplectic$) then
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
endif

! Compute normal mode and lab dispersion.

mat4 = -mat4
forall (i = 1:4) mat4(i,i) = mat4(i,i) + 1
call mat_inverse (mat4, mat4)

rel_p = 1 + orb0(6)

vec(1) = mat6(1,6) + (mat6(1,2) * orb0(2) + mat6(1,4) * orb0(4)) / rel_p
vec(2) = mat6(2,6) + (mat6(2,2) * orb0(2) + mat6(2,4) * orb0(4) - orb0(2)) / rel_p
vec(3) = mat6(3,6) + (mat6(3,2) * orb0(2) + mat6(3,4) * orb0(4)) / rel_p
vec(4) = mat6(4,6) + (mat6(4,2) * orb0(2) + mat6(4,4) * orb0(4) - orb0(4)) / rel_p
eta_vec = matmul(mat4, vec)

eta_vec(2) = eta_vec(2) / rel_p
eta_vec(4) = eta_vec(4) / rel_p

ele%x%eta  = eta_vec(1)
ele%x%etap = eta_vec(2)
ele%y%eta  = eta_vec(3)
ele%y%etap = eta_vec(4)
ele%z%eta  = 0
ele%z%etap = 1

call mat_symp_conj (v, v_inv)
eta_vec = matmul(v_inv, eta_vec)
ele%a%eta  = eta_vec(1)
ele%a%etap = eta_vec(2)
ele%b%eta  = eta_vec(3)
ele%b%etap = eta_vec(4)

end subroutine
