!+
! Subroutine twiss_propagate1 (ele1, ele2, err_flag)
!
! Subroutine to propagate the twiss, coupling, and dispersion parameters from 
! the exit end of ele1 to the exit end of ele2.
!
! Input:
!   ele1        -- Ele_struct: Structure holding the starting parameters.
!   ele2        -- Ele_struct: Structure holding the transfer matrix.
!     %key                -- Needed since, for example, Match element are handled 
!                              differently from other elements.
!     %map_ref_orb_in     -- Important for the dispersion calc.
!     %map_ref_orb_out    -- Important for the dispersion calc.
!
! Output:
!   ele2     -- Ele_struct: Structure for the ending parameters.
!   err_flag -- Logical, optional: Set True if there is an error. False otherwise.
!-

subroutine twiss_propagate1 (ele1, ele2, err_flag)

use bmad_interface, except_dummy => twiss_propagate1

implicit none
type (ele_struct), target :: ele1, ele2
type (twiss_struct) twiss_a
type (lat_param_struct) param
type (coord_struct), pointer :: orb, orb_out
integer key2

real(rp) :: mat6(6,6)
real(rp) v_mat(4,4), v_inv_mat(4,4), det, mat2_a(2,2), mat2_b(2,2)
real(rp) big_M(2,2), small_m(2,2), big_N(2,2), small_n(2,2)
real(rp) c_conj_mat(2,2), E_inv_mat(2,2), F_inv_mat(2,2)
real(rp) mat2(2,2), eta1_vec(6), eta_vec(6), vec(6), dpz2_dpz1, rel_p1, rel_p2
real(rp) det_factor, deriv_rel, gamma2_c, df

logical error
logical, optional :: err_flag

character(*), parameter :: r_name = 'twiss_propagate1'

! Beginning element bookkeeping

if (ele1%key == beginning_ele$) then
  ele1%map_ref_orb_out = ele2%map_ref_orb_in
  rel_p1 = 1 + ele1%map_ref_orb_out%vec(6)

  if (ele1%x%etap == real_garbage$ .and. ele1%x%deta_ds == real_garbage$) then  ! Not set so use default
    ele1%x%etap = 0
    ele1%x%deta_ds = 0
  elseif (ele1%x%etap == real_garbage$) then
    ele1%x%etap = ele1%x%deta_ds * rel_p1 + ele1%map_ref_orb_out%vec(2) / rel_p1
  elseif (ele1%x%deta_ds == real_garbage$) then
    ele1%x%deta_ds = ele1%x%etap / rel_p1 - ele1%map_ref_orb_out%vec(2) / rel_p1**2
  endif

  if (ele1%y%etap == real_garbage$ .and. ele1%y%deta_ds == real_garbage$) then  ! Not set so use default
    ele1%y%etap = 0
    ele1%y%deta_ds = 0
  elseif (ele1%y%etap == real_garbage$) then
    ele1%y%etap = ele1%y%deta_ds * rel_p1 + ele1%map_ref_orb_out%vec(4) / rel_p1
  elseif (ele1%y%deta_ds == real_garbage$) then
    ele1%y%deta_ds = ele1%y%etap / rel_p1 - ele1%map_ref_orb_out%vec(4) / rel_p1**2
  endif
endif

!

if (present(err_flag)) err_flag = .true.

if (ele2%key == match$ .and. is_true(ele2%value(recalc$)) .and. &
    (nint(ele2%value(matrix$)) == match_twiss$ .or. nint(ele2%value(matrix$)) == phase_trombone$)) then
  call match_ele_to_mat6(ele2, ele1%map_ref_orb_out, ele2%mat6, ele2%vec0, error, set_trombone = .true.)
  if (error) return
endif

!---------------------------------------------------------------------
! init

if (ele1%a%beta <= 0 .or. ele1%b%beta <= 0) then

  ! For x-ray lines assume that since beta = 0 there is no interest in calculating the Twiss parameters.
  ! So this is not treated as an error.
  if (associated(ele1%branch)) then
    if (ele1%branch%param%particle == photon$) then
      if (present(err_flag)) err_flag = .false.
      return
    endif
  endif

  call out_io (s_error$, r_name, 'NON-POSITIVE BETA DETECTED AT ELEMENT: ' // &
                                  trim(ele1%name) // '  ' // ele_loc_name(ele1, .true., '()'))
  return
endif

ele2%mode_flip = ele1%mode_flip          ! assume no flip
key2 = ele2%key

!---------------------------------------------------------------------
! markers are easy

if (key2 == marker$ .or. key2 == photon_fork$ .or. key2 == fork$) then
  call transfer_twiss (ele1, ele2)
  if (present(err_flag)) err_flag = .false.
  return
endif

!

orb  => ele2%map_ref_orb_in
orb_out => ele2%map_ref_orb_out
rel_p1 = 1 + orb%vec(6)               ! reference energy 
rel_p2 = 1 + orb_out%vec(6)

mat6 = ele2%mat6
if (ele2%key /= e_gun$) then   ! Energy change normalization is not applied to an e-gun
  if (bmad_private%normalize_twiss) then
    mat6(:, 2:6:2) = mat6(:, 2:6:2) * rel_p1
    mat6(2:6:2, :) = mat6(2:6:2, :) / rel_p2
  endif
  rel_p2 = rel_p2 / rel_p1
  rel_p1 = 1
endif

!---------------------------------------------------------------------
! det_factor is a renormalization factor to handle non-symplectic errors.

det_factor = sqrt(determinant(mat6(1:4,1:4)))
if (det_factor == 0) return  ! Can happen if matrix was never computed.

!---------------------------------------------------------------------
! if transfer matrix is not coupled...
! propagate c_mat coupling matrix and setup temporary element for propagation

if (all(mat6(1:2,3:4) == 0)) then
  mat2_a = mat6(1:2,1:2)
  mat2_b = mat6(3:4,3:4) 

  ele2%c_mat = matmul(matmul(mat2_a, ele1%c_mat), mat_symp_conj(mat2_b)) / det_factor ! conj == inverse
  ele2%gamma_c = ele1%gamma_c

!---------------------------------------------------------------------
! here if we are dealing with a coupled transfer matrix

else
  df = (1.0_rp/det_factor)**0.25
  big_M   = df * mat6(1:2,1:2)
  small_m = df * mat6(1:2,3:4)
  big_N   = df * mat6(3:4,3:4)
  small_n = df * mat6(3:4,1:2)

  c_conj_mat = mat_symp_conj (ele1%c_mat)
  mat2 = ele1%gamma_c * big_M - matmul(small_m, c_conj_mat)
  det = determinant(mat2)

  ! we demand that gamma_c > 0.3 (ie det > 0.1)
  ! try to make it so that there is no net mode flip here

  if (det > 0.9 .or. (det > 0.1 .and. .not. ele1%mode_flip)) then

    gamma2_c = sqrt(det)
    mat2_a = mat2 / gamma2_c
    mat2_b = (ele1%gamma_c * big_N + matmul(small_n, ele1%c_mat)) / gamma2_c
    F_inv_mat = mat_symp_conj (mat2_b)
    ele2%c_mat = matmul(matmul(big_M, ele1%c_mat) + ele1%gamma_c * small_m, F_inv_mat)
    ele2%gamma_c = sqrt(1.0_rp - determinant(ele2%c_mat))

  ! else we flip the modes

  else
    mat2 = matmul(big_M, ele1%c_mat) + ele1%gamma_c * small_m
    det = determinant(mat2)
    if (det < 0) then
      call out_io (s_error$, r_name,  '||mat2|| < 0! (Due to roundoff?) \f12.4\ ', &
                                      'When propagating through: [\i0\]  ' // trim(ele2%name), &
                                      r_array = [det], i_array = [ele2%ix_ele])
    endif

    gamma2_c = sqrt(abs(det))
    mat2_a = (ele1%gamma_c * small_n - matmul(big_N, c_conj_mat)) / gamma2_c
    mat2_b = mat2 / gamma2_c

    E_inv_mat = mat_symp_conj (mat2_a)
    ele2%c_mat = matmul(ele1%gamma_c * big_M - matmul(small_m, c_conj_mat), E_inv_mat)
    ele2%mode_flip = .not. ele1%mode_flip
    ele2%gamma_c = sqrt(1.0_rp - determinant(ele2%c_mat))
  endif

endif

!---------------------------------------------------------------------
! Propagate Twiss.
! If there is a mode flip (ele%mode_flip = T), ele%a is the "b"-mode. That is, ele%a is associated with
! the lower right block of the U matrix (See Sagan & Rubin: Linear Analysis of Coupled Lattices).
! Another way of saying this: ele%a always represents the same physical mode.

if (ele1%mode_flip) then
  call twiss1_propagate (ele1%a, mat2_b, ele2%key, ele2%value(l$), ele2%a, error); if (error) return
  call twiss1_propagate (ele1%b, mat2_a, ele2%key, ele2%value(l$), ele2%b, error); if (error) return
else
  call twiss1_propagate (ele1%a, mat2_a, ele2%key, ele2%value(l$), ele2%a, error); if (error) return
  call twiss1_propagate (ele1%b, mat2_b, ele2%key, ele2%value(l$), ele2%b, error); if (error) return
endif

! Comming out of a flipped state, the phase is often off by a factor of twopi from what one would expect.
! Factors of twopi are not physically meaningful so this does not affect any calculations.
! Unfortunately there is no definitive way to calcuate what the "right" way to handle this is.
! Subtracting off a factor of twopi is a bit of a kludge but it is better than nothing. 

if (ele1%mode_flip .and. .not. ele2%mode_flip) then
  ele2%a%phi = ele2%a%phi - twopi
  ele2%b%phi = ele2%b%phi - twopi
endif

!----------------------------------------------------
! Dispersion calc.
! p_z2 is p_z at end of ele2 assuming p_z = 1 at end of ele1.
! This is just 1.0 (except for RF cavities).

eta1_vec = [ele1%x%eta, ele1%x%deta_ds * rel_p1, ele1%y%eta, ele1%y%deta_ds * rel_p1, ele1%z%eta, 1.0_rp]

! For a circular ring, defining the dependence of z on pz is problematical.
! With the RF off, z is not periodic so dz/dpz is dependent upon turn number.
! With the RF on, dz/dpz, along with the other dispersion components, is not well defined.
! This being the case, eta1_vec(5) is just treated as zero for an rfcavity.

if (key2 == rfcavity$) eta1_vec(5) = 0

! Must avoid 0/0 divide at zero reference momentum. 
! If rel_p1 = 0 then total momentum is zero and orb%vec(2) and orb%vec(4) must be zero.

! Also for elements with a static electric field, pz should include the potential energy and so the 
! mat6(6,:) terms should be all zero. However Bmad is not using proper canonical coords so the
! mat6(6,:) terms are forced to zero.

if (rel_p1 == 0) then
  dpz2_dpz1 = dot_product(mat6(6,:), eta1_vec) 
  eta_vec(1:5) = matmul (mat6(1:5,:), eta1_vec) / dpz2_dpz1
else
  if (key2 == rfcavity$ .or. key2 == lcavity$) then
    dpz2_dpz1 = dot_product(mat6(6,:), eta1_vec) 
    dpz2_dpz1 = dpz2_dpz1 + (mat6(6,2) * orb%vec(2) + mat6(6,4) * orb%vec(4)) / rel_p1
  else
    dpz2_dpz1 = 1 / rel_p2
  endif
endif

!

eta1_vec = [ele1%x%eta, ele1%x%etap, ele1%y%eta, ele1%y%etap, ele1%z%eta, 1.0_rp]
eta_vec(1:5) = matmul (mat6(1:5,:), eta1_vec) / dpz2_dpz1

ele2%x%eta     = eta_vec(1)
ele2%x%etap    = eta_vec(2)
ele2%x%deta_ds = eta_vec(2) / rel_p2 - orb_out%vec(2) / rel_p2**2

ele2%y%eta     = eta_vec(3)
ele2%y%etap    = eta_vec(4)
ele2%y%deta_ds = eta_vec(4) / rel_p2 - orb_out%vec(4) / rel_p2**2

ele2%z%eta     = eta_vec(5)
ele2%z%etap    = ele1%z%etap * dpz2_dpz1
ele2%z%deta_ds = ele1%z%deta_ds * dpz2_dpz1

call make_v_mats (ele2, v_mat, v_inv_mat)
eta_vec(1:4) = matmul (v_inv_mat, eta_vec(1:4))
vec(1:4) = matmul(v_inv_mat, orb_out%vec(1:4))

ele2%a%eta     = eta_vec(1)
ele2%a%etap    = eta_vec(2)
ele2%a%deta_ds = eta_vec(2) / rel_p2 - vec(2) / rel_p2**2

ele2%b%eta     = eta_vec(3)
ele2%b%etap    = eta_vec(4)
ele2%b%deta_ds = eta_vec(4) / rel_p2 - vec(4) / rel_p2**2

!

if (present(err_flag)) err_flag = .false.

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine twiss1_propagate (twiss1, mat2, ele_key, length, twiss2, err)
!
! Subroutine to propagate the twiss parameters of a single mode.
!
! The betatron phase advance del_phi is only determined up to a factor of 2pi. 
! The length and ele_key argument are used to determine what range
! of angles del_phi should be in:
!   length > 0:                           0 <= del_phi <  twopi
!   length = 0 or ele_key = patch$:     -pi <  del_phi <= pi
!   length < 0:                      -twopi <  del_phi <= 0
! The patch element is exceptional in that its length is defined in a somewhat 
! arbitrary manner and thus is not a good reference as to what the phase advance should be.
!
! Note: The soft edge multipole fringe may give a slightly negative phase shift.
! So to avoid unwanted 2pi phase shifts if del_phi is small, the above ranges are 
! only enforced if If |del_phi| > 0.1.
!
! Input:
!   twiss1    -- Twiss_struct: Input Twiss parameters.
!   mat2(2,2) -- Real(rp): The transfer matrix.
!   ele_key   -- Integer: quadrupole$, etc.
!   length    -- Real(rp): Determines whether the phase is 
!                            increasing or decreasing.
!
! Output:
!   twiss2    -- Twiss_struct: Output Twiss parameters.
!   err       -- Logical: Set True if there is an error, false otherwise.
!-

subroutine twiss1_propagate (twiss1, mat2, ele_key, length, twiss2, err)

use bmad_interface, except_dummy => twiss1_propagate

implicit none

type (twiss_struct)  twiss1, twiss2, temp

real(rp) m11, m12, m21, m22, del_phi, length
real(rp) a1, b1, g1, a2, b2, g2, mat2(2,2), det
integer ele_key
logical err

!----------------------------------------------------
! Basic equation is given by Bovet 2.5.b page 16
! Linac rf matrices need to be renormalized.

err = .true.

det = determinant (mat2)

if (det == 0 .or. twiss1%beta == 0 .or. twiss1%beta > 1d100) return  ! Limit max beta to prevent numerical overflow.

m11 = mat2(1,1)
m12 = mat2(1,2)
m21 = mat2(2,1)
m22 = mat2(2,2)

a1 = twiss1%alpha
b1 = twiss1%beta
g1 = (1 + a1**2) / b1

b2 =       (m11**2  * b1 - 2*m11*m12 * a1 + m12**2  * g1) / det
a2 = a1 + (-m21*m11 * b1 + 2*m12*m21 * a1 - m12*m22 * g1) / det
g2 =  (1 + a2**2) /b2

del_phi = atan2(m12, m11*b1 - m12*a1)

! The soft edge multipole fringe may give a slightly negative phase shift.
! So avoid unwanted 2pi phase shifts if del_phi is small.

if (ele_key /= patch$ .and. abs(del_phi) > 0.1) then
  if (del_phi < 0 .and. length > 0) del_phi = del_phi + twopi
  if (del_phi > 0 .and. length < 0) del_phi = del_phi - twopi
endif

twiss2%beta = b2
twiss2%alpha = a2
twiss2%gamma = g2
twiss2%phi = twiss1%phi + del_phi

err = .false.

end subroutine
