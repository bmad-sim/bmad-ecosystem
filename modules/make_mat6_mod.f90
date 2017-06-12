module make_mat6_mod

use bmad_utils_mod

contains

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine mat6_add_pitch (x_pitch_tot, y_pitch_tot, orientation, mat6)
!
! Subroutine to modify a first order transfer matrix to include the affect
! of an element pitch. Note that this routine does not correct the 0th order
! part of the map. It is assumed that on input the transfer map
! does not include the affect of any pitches.
!
! Modules needed:
!   use bmad
!
! Input:
!   x_pitch_tot -- Real(rp): Horizontal pitch
!   y_pitch_tot -- Real(rp): Vertical pitch
!   orientation -- integer: Element longitudinal orientation. +1 or -1.
!   mat6(6,6)   -- Real(rp): 1st order part of the transfer map (Jacobian).
!
! Output:
!   mat6(6,6) -- Real(rp): 1st order xfer map with pitches.
!-

subroutine mat6_add_pitch (x_pitch_tot, y_pitch_tot, orientation, mat6)

implicit none

real(rp) mat6(6,6), x_pitch_tot, y_pitch_tot
integer orientation

!

if (x_pitch_tot == 0 .and. y_pitch_tot == 0) return

! The equations below are performing matrix multiplication. The original matrix
! is being multiplied from left and right by matrices that correspond to the pitches. 
! The pitch matrices are obtained by differentiating the corresponding equations in   
! the offset_particle subroutine. The (i,j) numbers mentioned as comments refer to  
! the non-zero elements present in the pitch matrices. 

mat6(:,6) = mat6(:,6) - mat6(:,2) * orientation * x_pitch_tot ! (2,6)
mat6(:,1) = mat6(:,1) + mat6(:,5) * orientation * x_pitch_tot ! (5,1)

mat6(:,6) = mat6(:,6) - mat6(:,4) * orientation * y_pitch_tot ! (4,6)
mat6(:,3) = mat6(:,3) + mat6(:,5) * orientation * y_pitch_tot ! (5,3)

mat6(2,:) = mat6(2,:) + orientation * x_pitch_tot * mat6(6,:) ! (2,6)
mat6(5,:) = mat6(5,:) - orientation * x_pitch_tot * mat6(1,:) ! (5,1)

mat6(4,:) = mat6(4,:) + orientation * y_pitch_tot * mat6(6,:) ! (4,6)
mat6(5,:) = mat6(5,:) - orientation * y_pitch_tot * mat6(3,:) ! (5,3)

end subroutine mat6_add_pitch

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine quad_mat2_calc (k1, length, rel_p, mat2, z_coef, dz_dpz_coef)
!
! Subroutine to calculate the 2x2 transfer matrix for a quad for one plane. 
!
! Modules needed:
!   use bmad
!
! Input:
!   k1       -- Real(rp): Quad strength: k1 > 0 ==> defocus.
!   length   -- Real(rp): Quad length
!   rel_p    -- Real(rp), Relative momentum P/P0.      
!
! Output:
!   mat2(2,2)      -- Real(rp): Transfer matrix.
!   z_coef(3)      -- Real(rp), optional: Coefficients for calculating the
!                       the change in z position:
!                          z = Integral [-(px/(1+pz))^2/2 ds]
!                            = c(1) * x_0^2 + c(2) * x_0 * px_0 + c(3) * px_0^2 
!   dz_dpz_coef(3) -- Real(rp), optional: Coefficients for calculating the
!                       the mat6(5,6) Jacobian matrix element:
!                         dz_dpz = c(1) * x_0^2 + c(2) * x_0 * px_0 + c(3) * px_0^2 
!-

subroutine quad_mat2_calc (k1, length, rel_p, mat2, z_coef, dz_dpz_coef)

implicit none

real(rp) length, mat2(:,:), cx, sx
real(rp) k1, sqrt_k, sk_l, k_l2, zc(3), dsx, dcx, rel_p
real(rp), optional :: z_coef(3), dz_dpz_coef(3)

!

sqrt_k = sqrt(abs(k1))
sk_l = sqrt_k * length

if (abs(sk_l) < 1d-10) then
  k_l2 = k1 * length**2
  cx = 1 + k_l2 / 2
  sx = (1 + k_l2 / 6) * length
elseif (k1 < 0) then       ! focus
  cx = cos(sk_l)
  sx = sin(sk_l) / sqrt_k
else                       ! defocus
  cx = cosh(sk_l)
  sx = sinh(sk_l) / sqrt_k
endif

mat2(1,1) = cx
mat2(1,2) = sx / rel_p
mat2(2,1) = k1 * sx * rel_p
mat2(2,2) = cx

!

if (present(z_coef) .or. present(dz_dpz_coef)) then
  zc(1) = k1 * (-cx * sx + length) / 4
  zc(2) = -k1 * sx**2 / (2 * rel_p)
  zc(3) = -(cx * sx + length) / (4 * rel_p**2)
  if (present(z_coef)) z_coef = zc
endif

! dz_dpz_coef

if (present(dz_dpz_coef)) then

  if (abs(sk_l) < 1d-10) then
    dcx = -k_l2 / (2 * rel_p)
    dsx = -k_l2 * length / (6 * rel_p)
  else
    dcx = -k1 * sx * length / (2 * rel_p)
    dsx = (sx - length * cx) / (2 * rel_p)
  endif

  dz_dpz_coef(1) =   -zc(1)/rel_p - k1 * (cx * dsx + dcx * sx) / 4
  dz_dpz_coef(2) = -2*zc(2)/rel_p - k1 * sx * dsx/rel_p
  dz_dpz_coef(3) = -2*zc(3)/rel_p - (cx * dsx + dcx * sx) / (4 * rel_p**2)

endif

end subroutine quad_mat2_calc

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine multipole_kick_mat (knl, tilt, ref_species, ele, orbit, factor, mat6)
!
! Subroutine to return the multipole kick components needed to
! construct the transfer matrix.
! This routine is not meant for general use.
!
! Input:
!   knl(0:)     -- Real(rp): Strength of multipoles
!   tilt(0:)    -- Real(rp): Tilt of multipoles
!   ref_species -- integer: Reference species.
!   ele         -- ele_struct: Lattice element containing multipoles.
!   orbit       -- Coord_struct: coordinates of particle around which the
!                    multipole kick matrix is computed.
!   factor      -- real(rp): Factor to scale knl by.
!
! Output:
!   mat6(6,6)   -- Real(rp): matrix with kick values at mat6(2:4:2, 1:3:2).
!                   The rest of the matrix is untouched.
!-

subroutine multipole_kick_mat (knl, tilt, ref_species, ele, orbit, factor, mat6)

implicit none

type (coord_struct) orbit
type (ele_struct) ele

real(rp) mat6(6,6), kmat1(4,4), factor, charge_dir
real(rp) knl(0:), tilt(0:)

integer ref_species, n

!                        

mat6(2:4:2, 1:3:2) = 0
if (orbit%vec(1) == 0 .and. orbit%vec(3) == 0 .and. knl(1) == 0) return
charge_dir = orbit%direction * ele%orientation * charge_to_mass_of(orbit%species) / charge_to_mass_of(ref_species)

do n = 1, ubound(knl, 1)
  if (knl(n) /= 0) then
    call mat4_multipole (knl(n)*charge_dir, tilt(n), n, orbit, kmat1)
    mat6(2:4:2, 1:3:2) = mat6(2:4:2, 1:3:2) + factor * kmat1(2:4:2, 1:3:2)
  endif
enddo

end subroutine multipole_kick_mat

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine mat4_multipole (knl, tilt, n, orbit, kick_mat)
!
! Subroutine to find the kick matrix (Jacobian) due to a multipole.
! This routine is not meant for general use.
!
! Input:
!   orbit   -- Coord_struct: coordinates of particle
!   knl     -- Real(rp): Strength of multipole
!   tilt    -- Real(rp): Tilt of multipole
!
! Output:
!   kick_mat(4,4) -- Real(rp): Kick matrix (Jacobian) at orbit.
!-


subroutine mat4_multipole (knl, tilt, n, orbit, kick_mat)

implicit none

type (coord_struct) orbit

real(rp) x_pos, y_pos, x, y, knl, tilt
real(rp) sin_ang, cos_ang, mat(2,2), rot(2,2)
real(rp) kick_mat(4,4)

integer m, n

! init

kick_mat = 0
forall (m = 1:4) kick_mat(m,m) = 1

x_pos = orbit%vec(1)
y_pos = orbit%vec(3)
         
! simple case

if (knl == 0 .or. (x_pos == 0 .and. y_pos == 0 .and. n /= 1)) then
  kick_mat(2:4:2, 1:3:2) = 0
  return
endif

! get position of particle in frame of multipole

if (tilt == 0) then
  x = x_pos
  y = y_pos
else
  sin_ang = sin(tilt)
  cos_ang = cos(tilt)
  x =  x_pos * cos_ang + y_pos * sin_ang
  y = -x_pos * sin_ang + y_pos * cos_ang
endif

! compute kick matrix

mat = 0

do m = 0, n, 2
  mat(1,1) = mat(1,1) +  knl * (n-m) * c_multi(n, m) * mexp(x, n-m-1) * mexp(y, m)
  mat(1,2) = mat(1,2) +  knl * m * c_multi(n, m) * mexp(x, n-m) * mexp (y, m-1)
enddo

do m = 1, n, 2
  mat(2,1) = mat(2,1) +  knl * (n-m) * c_multi(n, m) * mexp(x, n-m-1) * mexp(y, m)
  mat(2,2) = mat(2,2) +  knl * m * c_multi(n, m) * mexp(x, n-m) * mexp(y, m-1)
enddo

! transform back to lab frame

if (tilt /= 0) then
  rot(1,1) =  cos_ang
  rot(1,2) = -sin_ang
  rot(2,1) =  sin_ang
  rot(2,2) =  cos_ang
  mat = matmul(rot, mat)
  rot(1,2) =  sin_ang
  rot(2,1) = -sin_ang
  mat = matmul (mat, rot)
endif

kick_mat(2,1) = mat(1,1)
kick_mat(2,3) = mat(1,2)
kick_mat(4,1) = mat(2,1)
kick_mat(4,3) = mat(2,2)

end subroutine mat4_multipole

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine bbi_slice_calc (ele, n_slice, z_slice)
!
! Routine to compute the longitudinal positions of the slices of
! a beambeam element.
!
! Input:
!   ele        -- ele_struct: beambeam element
!   n_slice    -- integer: Number of slices
!
! Output:
!   z_slice(:) -- real(rp): Array of slice positions 
!-

subroutine bbi_slice_calc (ele, n_slice, z_slice)

implicit none

type (ele_struct) ele

integer :: i, n_slice
real(rp) z_slice(:), y
real(rp) :: z_norm

!

if (n_slice == 1) then
  z_slice(1) = ele%value(z_offset$)
elseif (n_slice > 1) then
  do i = 1, n_slice
    y = (i - 0.5) / n_slice - 0.5
    z_norm = inverse(probability_funct, y, -5.0_rp, 5.0_rp, 1.0e-5_rp)
    z_slice(i) = ele%value(sig_z$) * z_norm + ele%value(z_offset$)
  enddo
else
  print *, 'ERROR IN BBI_SLICE_CALC: N_SLICE IS NEGATIVE:', n_slice
  if (global_com%exit_on_error) call err_exit
endif

z_slice(n_slice+1) = 0

end subroutine bbi_slice_calc

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+      
! Subroutine tilt_mat6 (mat6, tilt)
!
! Subroutine to transform a 6x6 transfer matrix to a new reference frame
! that is tilted in (x, Px, y, Py) with respect to the old reference frame.
!     mat6 -> tilt_mat * mat6 * tilt_mat_inverse
!
! Modules needed:
!   use bmad
!
! Input:
!   mat6(6,6) -- Real(rp): Untilted matrix.
!   tilt      -- Real(rp): Tilt angle.
!
! Output:
!   mat6(6,6) -- Real(rp): Tilted matrix.
!-

subroutine tilt_mat6 (mat6, tilt)

implicit none

real(rp) tilt, mat6(6,6), mm(6,6)
real(rp) c, s

!

if (tilt == 0) return

c = cos(tilt)
s = sin(tilt)

mm(1,:) = c * mat6(1,:) - s * mat6(3,:)
mm(2,:) = c * mat6(2,:) - s * mat6(4,:)
mm(3,:) = c * mat6(3,:) + s * mat6(1,:)
mm(4,:) = c * mat6(4,:) + s * mat6(2,:)
mm(5,:) =     mat6(5,:)
mm(6,:) =     mat6(6,:)

mat6(:,1) = mm(:,1) * c - mm(:,3) * s
mat6(:,2) = mm(:,2) * c - mm(:,4) * s
mat6(:,3) = mm(:,3) * c + mm(:,1) * s
mat6(:,4) = mm(:,4) * c + mm(:,2) * s
mat6(:,5) = mm(:,5)
mat6(:,6) = mm(:,6)

end subroutine tilt_mat6

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine mat6_coord_transformation (mat6, orbit, param, dr, w_mat)
!
! Subroutine to calculate the matrix associated with a coordinate transformation.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele        -- real(rp): Element to drift through.
!   orbit      -- coord_struct: Starting coords
!   param      -- lat_param_struct:
!   dr(3)      -- real(rp): Transformation offset.
!   w_mat(3,3) -- real(rp): Transformation rotation matrix. 
!
! Output:
!   mat6(6,6) -- Real(rp): Transfer matrix
!-

subroutine mat6_coord_transformation (mat6, ele, param, orbit, dr, w_mat)

implicit none

type (ele_struct) ele
type (lat_param_struct) param
type (coord_struct) orbit


real(rp) dr(3), w_mat(3,3), mat6(6,6), m6_drift(6,6), rel_len
real(rp) rel_p, w_inv(3,3), beta_ref, vec(3), pvec(3), px, py, pz

!

w_inv = transpose(w_mat)

rel_p = 1 + orbit%vec(6)
pz = sqrt(rel_p**2 - orbit%vec(2)**2 - orbit%vec(4)**2)
beta_ref = ele%value(p0c$) / ele%value(e_tot$)

mat6(1,:) = [w_inv(1,1), 0.0_rp, w_inv(1,2), 0.0_rp, 0.0_rp, 0.0_rp]
mat6(3,:) = [w_inv(2,1), 0.0_rp, w_inv(2,2), 0.0_rp, 0.0_rp, 0.0_rp]

mat6(2,:) = [0.0_rp, (w_inv(1,1) - w_inv(1,3) * orbit%vec(2) / pz), &
             0.0_rp, (w_inv(1,2) - w_inv(1,3) * orbit%vec(4) / pz), 0.0_rp, w_inv(1,3) * rel_p / pz]
mat6(4,:) = [0.0_rp, (w_inv(2,1) - w_inv(2,3) * orbit%vec(2) / pz), &
             0.0_rp, (w_inv(2,2) - w_inv(2,3) * orbit%vec(4) / pz), 0.0_rp, w_inv(2,3) * rel_p / pz]

mat6(5,:) = [0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 0.0_rp]
mat6(6,:) = [0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp]

! Drift

call mat_make_unit (m6_drift)

vec  = matmul (w_inv, [orbit%vec(1), orbit%vec(3), 0.0_rp] - dr)
pvec = matmul (w_inv, [orbit%vec(2), orbit%vec(4), pz])
px = pvec(1); py = pvec(2); pz = pvec(3)
rel_len = -vec(3) / pz

m6_drift(1,2) = rel_len * (px**2 / pz**2 + rel_p**2)
m6_drift(3,4) = rel_len * (py**2 / pz**2 + rel_p**2)
m6_drift(1,4) = rel_len * px*py / pz**2
m6_drift(3,2) = rel_len * px*py / pz**2
m6_drift(1,6) = - rel_len * rel_p * px / pz**2
m6_drift(3,6) = - rel_len * rel_p * py / pz**2
m6_drift(5,2) = - rel_len * rel_p * px / pz**2
m6_drift(5,4) = - rel_len * rel_p * py / pz**2

! Reference particle drift length is zero.
m6_drift(5,6) = -vec(3) * (px**2 + py**2) / pz**3 

! These matrix terms are due to the variation of vec(3) drift length
m6_drift(1,1) = m6_drift(1,1) + (w_inv(1,3) / w_inv(3,3)) * px / pz
m6_drift(1,3) = m6_drift(1,3) + (w_inv(2,3) / w_inv(3,3)) * px / pz

m6_drift(3,1) = m6_drift(3,1) + (w_inv(1,3) / w_inv(3,3)) * py / pz
m6_drift(3,3) = m6_drift(3,3) + (w_inv(2,3) / w_inv(3,3)) * py / pz

m6_drift(5,1) = m6_drift(5,1) - (w_inv(1,3) / w_inv(3,3)) * rel_p / pz
m6_drift(5,3) = m6_drift(5,3) - (w_inv(2,3) / w_inv(3,3)) * rel_p / pz

mat6 = matmul (m6_drift, mat6)

end subroutine mat6_coord_transformation

end module
