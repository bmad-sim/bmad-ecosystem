!+
! Subroutine twiss_from_mat6 (mat6, orb0, ele, stable, growth_rate, status, type_out)
!
! Subroutine to calculate the Twiss parameters from a 1-turn matrix.
! Note: The 1-turn matrix needs to be formed with the RF turned off.
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
real(rp) m6(6,6), mat4(4,4), eta_vec(4), vec(4), rel_p
real(rp) u(4,4), v(4,4), ubar(4,4), vbar(4,4), g(4,4)
real(rp) rate1, rate2, symp_err
real(rp) :: symp_tol = 3.0d-3

integer :: status
integer i

logical :: stable, type_out

character(20) :: r_name = 'twiss_from_mat6'

! init

stable = .false.
m6 = mat6
rel_p = 1 + orb0(6)

mat4 = m6(1:4, 1:4)

if (maxval(abs(mat4)) > 1d10) then
  if (type_out) call out_io (s_error$, r_name, 'BAD 1-TURN MATRIX: UNSTABLE.', &
                                               'TWISS PARAMETERS NOT COMPUTED')
  status = unstable$
  growth_rate = maxval(abs(mat4))
  return
endif

symp_err = mat_symp_error(mat4)

if (symp_err > 1) then
  if (type_out) call out_io (s_error$, r_name, 'BAD 1-TURN MATRIX. NON_SYMPLECTIC WITH SYMPLECTIC ERROR OF: \f8.1\ ', &
                                               'TWISS PARAMETERS NOT COMPUTED', r_array = [symp_err])
  status = non_symplectic$
  growth_rate = min(1d5 * symp_err, maxval(abs(mat4)))

  return
endif

if (symp_err > symp_tol) then
  if (type_out) call out_io (s_warn$, r_name, '1-TURN MATRIX MARGINALLY SYMPLECTIC WITH SYMPLECTIC ERROR OF: \f10.4\ ', &
                                r_array = [symp_err])
endif

!

call mat_symp_decouple (mat4, status, u, v, ubar, vbar, g, ele%a, ele%b, ele%gamma_c, .false.)

if (status /= ok$) then
  if (type_out) call out_io (s_error$, r_name, 'BAD 1-TURN MATRIX: ' // &
                                                matrix_status_name(status), 'TWISS PARAMETERS NOT COMPUTED')
  rate1 = sqrt(max(abs(u(1,1) + u(2,2)) - 2, 0.0_rp))
  rate2 = sqrt(max(abs(u(3,3) + u(4,4)) - 2, 0.0_rp))
  growth_rate = max(rate1, rate2)
  if (growth_rate > 1d4) growth_rate = 1d4 + 1d2 * log10(growth_rate - 1d4 + 1.0_rp)
  return
endif

! Here if everything normal so load twiss parameters

growth_rate = 0          ! no growth
stable = .true.          ! stable lat

if (ele%a%beta /= 0 .and. ele%b%beta /= 0) then
  ele%mode_flip = .false.
  ele%c_mat = v(1:2,3:4)
endif

! Compute normal mode and lab dispersion.

mat4 = -mat4
forall (i = 1:4) mat4(i,i) = mat4(i,i) + 1
call mat_inverse (mat4, mat4)

eta_vec = matmul(mat4, mat6(1:4,6))

ele%x%eta     = eta_vec(1)
ele%x%etap    = eta_vec(2)
ele%x%deta_ds = ele%x%etap / rel_p - orb0(2) / rel_p**2 

ele%y%eta     = eta_vec(3)
ele%y%etap    = eta_vec(4)
ele%y%deta_ds = eta_vec(4) / rel_p - orb0(4) / rel_p*2

ele%z%eta     = 0
ele%z%etap    = 1
ele%z%deta_ds = 1

!

eta_vec = matmul(mat_symp_conj(v), eta_vec)
vec = matmul(mat_symp_conj(v), orb0(1:4))

ele%a%eta     = eta_vec(1)
ele%a%etap    = eta_vec(2)
ele%a%deta_ds = ele%a%etap / rel_p - vec(2) / rel_p**2 

ele%b%eta     = eta_vec(3)
ele%b%etap    = eta_vec(4)
ele%b%deta_ds = eta_vec(4) / rel_p - vec(4) / rel_p*2

end subroutine
