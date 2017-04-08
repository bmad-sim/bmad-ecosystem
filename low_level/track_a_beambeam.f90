!+
! Subroutine track_a_beambeam (orbit, ele, param, mat6, make_matrix)
!
! Particle tracking through a beambeam element. 
!
! Modules Needed:
!   use bmad
!
! Input:
!   orbit       -- Coord_struct: Starting position.
!   ele         -- Ele_struct: Beambeam element.
!   param       -- lat_param_struct: Lattice parameters.
!   make_matrix -- logical, optional: Propagate the transfer matrix? Default is false.
!
! Output:
!   orbit      -- Coord_struct: End position.
!   mat6(6,6)  -- Real(rp), optional: Transfer matrix through the element.
!-

subroutine track_a_beambeam (orbit, ele, param, mat6, make_matrix)

use track1_mod, except_dummy => track_a_beambeam

implicit none

type (coord_struct) :: orbit, start_orb
type (ele_struct), target :: ele
type (lat_param_struct) :: param

real(rp), optional :: mat6(6,6)
real(rp) sig_x0, sig_y0, beta_a0, beta_b0, alpha_a0, alpha_b0, sig_x, sig_y
real(rp) z_slice(100), beta, s_pos, s_pos_old, k0_x, k0_y, k_xx, k_xy, k_yx, k_yy, coef, del
real(rp) mat21, mat23, mat41, mat43, del_s, x_pos, y_pos, ratio, bbi_const

integer i, n_slice

logical, optional :: make_matrix

!

if (ele%value(charge$) == 0 .or. param%n_part == 0) return

if (make_matrix) call mat_make_unit(mat6)

del = 0.001
sig_x0 = ele%value(sig_x$)
sig_y0 = ele%value(sig_y$)
if (sig_x0 == 0 .or. sig_y0 == 0) return

if (ele%value(beta_a$) == 0) then
  beta_a0 = ele%a%beta
  alpha_a0 = ele%a%alpha
else
  beta_a0 = ele%value(beta_a$)
  alpha_a0 = ele%value(alpha_a$)
endif

if (ele%value(beta_b$) == 0) then
  beta_b0 = ele%b%beta
  alpha_b0 = ele%b%alpha
else
  beta_b0 = ele%value(beta_b$)
  alpha_b0 = ele%value(alpha_b$)
endif

call offset_particle (ele, param, set$, orbit)
call canonical_to_angle_coords (orbit)

n_slice = max(1, nint(ele%value(n_slice$)))
call bbi_slice_calc (ele, n_slice, z_slice)
s_pos = 0    ! end at the ip

do i = 1, n_slice
  s_pos_old = s_pos
  s_pos = (orbit%vec(5) + z_slice(i)) / 2
  del_s = s_pos - s_pos_old

  orbit%vec(1) = orbit%vec(1) + orbit%vec(2) * del_s
  orbit%vec(3) = orbit%vec(3) + orbit%vec(4) * del_s

  if (make_matrix) then
    mat6(1,:) = mat6(1,:) + del_s * mat6(2,:)
    mat6(3,:) = mat6(3,:) + del_s * mat6(4,:)
  endif

  if (beta_a0 == 0) then
    sig_x = sig_x0
    sig_y = sig_y0
  else
    beta = beta_a0 - 2 * alpha_a0 * s_pos + (1 + alpha_a0**2) * s_pos**2 / beta_a0
    sig_x = sig_x0 * sqrt(beta / beta_a0)
    beta = beta_b0 - 2 * alpha_b0 * s_pos + (1 + alpha_b0**2) * s_pos**2 / beta_b0
    sig_y = sig_y0 * sqrt(beta / beta_b0)
  endif

  bbi_const = -param%n_part * ele%value(charge$) * classical_radius_factor /  &
                                        (2 * pi * ele%value(p0c$) * (sig_x + sig_y))

  x_pos = orbit%vec(1) / sig_x  ! this has offset in it
  y_pos = orbit%vec(3) / sig_y
  ratio = sig_y / sig_x

  call bbi_kick (x_pos, y_pos, ratio, k0_x, k0_y)

  coef = bbi_const / (n_slice * (1 + orbit%vec(6)))
  orbit%vec(2) = orbit%vec(2) + k0_x * coef
  orbit%vec(4) = orbit%vec(4) + k0_y * coef

  if (make_matrix) then
    call bbi_kick (x_pos+del, y_pos, ratio, k_xx, k_yx)
    call bbi_kick (x_pos, y_pos+del, ratio, k_xy, k_yy)

    coef = bbi_const / (ele%value(n_slice$) * del * (1 + orbit%vec(6)))
    mat21 = coef * (k_xx - k0_x) / sig_x
    mat23 = coef * (k_xy - k0_x) / sig_x
    mat41 = coef * (k_yx - k0_y) / sig_y
    mat43 = coef * (k_yy - k0_y) / sig_y

    mat6(2,:) = mat6(2,:) + mat21 * mat6(1,:) + mat23 * mat6(3,:)
    mat6(4,:) = mat6(4,:) + mat41 * mat6(1,:) + mat43 * mat6(3,:)
  endif
enddo

orbit%vec(1) = orbit%vec(1) - orbit%vec(2) * s_pos
orbit%vec(3) = orbit%vec(3) - orbit%vec(4) * s_pos

call angle_to_canonical_coords (orbit)
call offset_particle (ele, param, unset$, orbit)  

end subroutine

