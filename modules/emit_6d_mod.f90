! Note: A negative emittance is possible and just means that the beam is
! unstable. That is, the corresponding damping partition number is negative.

!+
! Module emit_6d_mod
!
! Module for calculating the equilibrium sigma matrix for a closed geometry lattice.
!-

module emit_6d_mod

use mode3_mod

implicit none

type rad1_map_struct
  real(rp) :: damp_mat(6,6) = 0    ! Transfer matrix = damp_mat + nodamp_mat.
  real(rp) :: nodamp_mat(6,6) = 0  !
  real(rp) :: stoc_mat(6,6) = 0    ! Stochastic matrix.
end type

contains

!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!+
! Subroutine emit_6d (ele_ref, sigma_mat, emit)
!
! Routine to calculate the three normal mode emittances. Since the emattances are
! only an invariant in the limit of zero damping, the calculated emittance will
! vary depending upon the reference element.
!
! Input:
!   ele_ref         -- ele_struct: Origin of the 1-turn maps used to evaluate the emittances.
!
! Output:
!   sigma_mat(6,6)  -- real(rp): Sigma matrix
!   emit(3)         -- real(rp): The three normal mode emittances.
!-

subroutine emit_6d (ele_ref, sigma_mat, emit)

use f95_lapack, only: dgesv_f95

type (ele_struct) ele_ref

real(rp) sigma_mat(6,6), emit(3)
real(rp) damp_mat(6,6), stoc_mat(6,6), nodamp_mat(6,6)
real(rp) mt(21,21), v_sig(21,1)

integer i, j, k, ipev(21), info

integer, parameter :: w1(21) = [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 5, 5, 6]
integer, parameter :: w2(21) = [1, 2, 3, 4, 5, 6, 2, 3, 4, 5, 6, 3, 4, 5, 6, 4, 5, 6, 5, 6, 6]
integer, parameter :: v(6,6) = reshape( &
            [1,  2,  3,  4,  5,  6,   2,  7,  8,  9, 10, 11,   3,  8, 12, 13, 14, 15, &
             4,  9, 13, 16, 17, 18,   5, 10, 14, 17, 19, 20,   6, 11, 15, 18, 20, 21], [6,6])

logical err

! Analysis is documented in the Bmad manual

call damping_and_stochastic_rad_mats (ele_ref, ele_ref, damp_mat, stoc_mat, nodamp_mat)

damp_mat = damp_mat + nodamp_mat
call mat_make_unit(mt)

do i = 1, 21
  v_sig(i,1) = stoc_mat(w1(i), w2(i))

  do j = 1, 6
  do k = 1, 6
    mt(i,v(j,k)) = mt(i,v(j,k)) - damp_mat(w1(i),j) * damp_mat(w2(i),k)
  enddo
  enddo
enddo

call dgesv_f95(mt, v_sig, ipev, info)

if (info /= 0) then
  sigma_mat = -1
  emit = -1
  return
endif

do j = 1, 6
do k = 1, 6
  sigma_mat(j,k) = v_sig(v(j,k), 1)
enddo
enddo

call get_emit_from_sigma_mat(sigma_mat, emit, err_flag = err)

end subroutine emit_6d

!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!+
! Subroutine damping_and_stochastic_rad_mats (ele1, ele2, damp_mat, stoc_mat, nodamp_mat)
!
! Routine to calculate the damping and stochastic variance matrices from exit end of ele1
! to the exit end of ele2. Use ele1 = ele2 to get 1-turn matrices.
!
! If ele2 is before ele1 the integration range if from ele1 to the branch end plus 
! from the beginning to ele2.
!
! The damping matrix is defined to be:
!   M_transport = nodamp_mat + damp_mat
! where M_transport is the full transport matrix including damping.
!
! Note: The ele%mat6 matrices will be remade. By convention, these matrices
! do not include damping.
!
! Input:
!   ele1          -- ele_struct: Start element of integration range.
!   ele2          -- ele_struct: End element of integration range.
!
! Output:
!   damp_mat(6,6)   -- real(rp): Damping part of the transfer matrix.
!   stoc_mat(6,6)   -- real(rp): Stochastic variance matrix.
!   nodamp_mat(6,6) -- real(rp): Nondamping (symnplectic) part of the transfer matrix.
!-

subroutine damping_and_stochastic_rad_mats (ele1, ele2, damp_mat, stoc_mat, nodamp_mat)

use super_recipes_mod

type (ele_struct), target :: ele1, ele2
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele1_track, ele2_track, ele3
type (bmad_common_struct) bmad_com_save
type (rad1_map_struct), allocatable :: ds(:)
type (coord_struct), allocatable :: closed_orb(:)

real(rp) sig_mat(6,6)
real(rp) damp_mat(6,6), stoc_mat(6,6), nodamp_mat(6,6), mt(6,6), damp_coef, stoc_coef, radi
real(rp) :: kd_coef, kf_coef, g2_ave, g3_ave, gamma

integer ie

logical err_flag

!

call find_element_ends(ele1, ele3, ele1_track)
call find_element_ends(ele2, ele3, ele2_track)

branch => ele1_track%branch
allocate (ds(branch%n_ele_track))

bmad_com_save = bmad_com
bmad_com%radiation_fluctuations_on = .false.
bmad_com%radiation_damping_on = .true.

call closed_orbit_calc(branch%lat, closed_orb, 6, +1, branch%ix_branch, err_flag)
bmad_com = bmad_com_save
if (err_flag) return

call lat_make_mat6 (ele1_track%branch%lat, -1, closed_orb, branch%ix_branch, err_flag)
if (err_flag) return

! Calc typical radiation values to get an error tolerance.

g2_ave = 0; g3_ave = 0
do ie = 1, branch%n_ele_track
  ele3 => branch%ele(ie)
  if (ele3%key /= sbend$) cycle
  g2_ave = g2_ave + ele3%value(l$) * ele3%value(g$)**2
  g3_ave = g3_ave + ele3%value(l$) * ele3%value(g$)**3
enddo

call convert_pc_to (ele1_track%value(p0c$), ele1_track%ref_species, gamma = gamma)
radi = classical_radius(ele1_track%ref_species)
kd_coef = 2.0_rp * radi * gamma**3 / 3.0_rp
kf_coef = 55.0_rp * radi * h_bar_planck * gamma**5 * c_light / (24.0_rp * sqrt_3 * mass_of(ele1_track%ref_species))
damp_coef = kd_coef * g2_ave / branch%param%total_length
stoc_coef = kf_coef * g3_ave / branch%param%total_length

! Calculate element-by-element damping and stochastic mats.

do ie = 1, branch%n_ele_track
  call qromb_rad_mat_int(branch%ele(ie), closed_orb(ie-1), closed_orb(ie), damp_coef, stoc_coef, ds(ie))
enddo

!

call mat_make_unit(nodamp_mat)
damp_mat = 0
stoc_mat = 0

ie = ele1_track%ix_ele
do 
  ie = ie + 1
  if (ie > branch%n_ele_track) ie = 0
  if (ie /= 0) then
    ele3 => branch%ele(ie)
    mt = ds(ie)%nodamp_mat + ds(ie)%damp_mat
    nodamp_mat = matmul(ds(ie)%nodamp_mat, nodamp_mat)
    damp_mat = matmul(mt, damp_mat) + ds(ie)%damp_mat
    stoc_mat = matmul(matmul(mt, stoc_mat), transpose(mt)) + ds(ie)%stoc_mat
  endif
  if (ie == ele2%ix_ele) exit
enddo

!---------------------------------------------------------------------------------
contains

! This is adapted from qromb and trapzd from Numerical Recipes.

subroutine qromb_rad_mat_int (ele, orb_in, orb_out, damp_coef, stoc_coef, ds)

type qromb_pt1_struct
  real(rp) xmat(6,6)  ! Transfer map without damping from beginning of element
end type

type qromb_int_struct
  real(rp) h
  real(rp) damp_mat(6,6)  
  real(rp) stoc_mat(6,6) 
end type

type (ele_struct) ele
type (coord_struct) orb_in, orb_out, orb0, orb1
type (rad1_map_struct) :: ds
type (fringe_field_info_struct) fringe_info
type (qromb_pt1_struct) pt1(0:16)
type (qromb_int_struct) qi(0:4)
type (bmad_common_struct) bmad_com_save

real(rp) damp_mat1(6,6), stoc_mat1(6,6), damp_coef, stoc_coef
real(rp) mat0(6,6), mat1(6,6), ddamp(6,6), dstoc(6,6), damp_mat_sum(6,6), stoc_mat_sum(6,6), mat0_inv(6,6)
real(rp) del_z, l_ref, rel_tol, eps_damp, eps_stoc, z_pos, d_max

integer j, j1, i1, i2, n, j_min_test, ll, n_pts
integer :: j_max = 10

! No radiation cases

ds%nodamp_mat = ele%mat6
ds%damp_mat = 0
ds%stoc_mat = 0

if (ele%value(l$) == 0 .or. (orb_out%vec(2) == orb_in%vec(2) .and. &
                         orb_out%vec(4) == orb_in%vec(4) .and. ele%key /= sbend$)) return

!

bmad_com_save = bmad_com
bmad_com%radiation_fluctuations_on = .false.
bmad_com%radiation_damping_on = .false.

! Offsets and fringes at upstream end.

orb0 = orb_in

call mat_make_unit(mat0)
call offset_particle (ele, set$, orb0, set_hvkicks = .false., mat6 = mat0, make_matrix = .true.)
call init_fringe_info (fringe_info, ele)
if (fringe_info%has_fringe) then
  fringe_info%particle_at = first_track_edge$
  call apply_element_edge_kick(orb0, fringe_info, ele, ele%branch%param, .false., mat0, .true.)
endif

! Integrate through body

qi(0)%h = 4
qi(0)%damp_mat = 0
qi(0)%stoc_mat = 0

rel_tol = 1d-4
eps_damp = rel_tol * ele%value(l$) * damp_coef 
eps_stoc = rel_tol * ele%value(l$) * stoc_coef 

j_min_test = 3

if (ele%key == wiggler$ .or. ele%key == undulator$) then
  j_min_test = 5
endif

do j = 1, j_max
  if (j == 1) then
    n_pts = 2
    del_z = ele%value(l$)
    l_ref = 0
  else
    n_pts = 2**(j-2)
    del_z = ele%value(l$) / n_pts
    l_ref = del_z / 2
  endif

  damp_mat_sum = 0
  stoc_mat_sum = 0

  do n = 1, n_pts
    z_pos = l_ref + (n-1) * del_z
    call calc_rad_at_pt(ele, orb_in, z_pos, damp_mat1, stoc_mat1)
    damp_mat_sum = damp_mat_sum + damp_mat1 
    stoc_mat_sum = stoc_mat_sum + stoc_mat1 
  enddo

  j1 = min(j, 4)
  if (j > 4) qi(0:3) = qi(1:4)
  qi(j1)%h = 0.25_rp * qi(j1-1)%h
  qi(j1)%damp_mat = 0.5_rp * (qi(j1-1)%damp_mat + del_z * damp_mat_sum)
  qi(j1)%stoc_mat = 0.5_rp * (qi(j1-1)%stoc_mat + del_z * stoc_mat_sum)

  if (j < j_min_test) cycle

  do i1 = 1, 6
    do i2 = 1, 6
      call super_polint(qi(1:j1)%h, qi(1:j1)%damp_mat(i1,i2), 0.0_rp, damp_mat(i1,i2), ddamp(i1,i2))
      call super_polint(qi(1:j1)%h, qi(1:j1)%stoc_mat(i1,i2), 0.0_rp, stoc_mat(i1,i2), dstoc(i1,i2))
    enddo
  enddo
  d_max = max(eps_damp*maxval(abs(ddamp)), eps_stoc*maxval(abs(dstoc)))
  if (d_max < 1) exit

enddo

! Reference position for matrices are the element exit end

mat0_inv = mat_symp_conj(mat0)   ! mat0 is transport matrix through the upstream edge
ds%damp_mat = matmul(ele%mat6, matmul(mat0_inv, damp_mat)) 

ds%stoc_mat = matmul(matmul(mat0_inv, stoc_mat), transpose(mat0_inv))
ds%stoc_mat = matmul(matmul(ele%mat6, stoc_mat), transpose(ele%mat6))

bmad_com = bmad_com_save

end subroutine qromb_rad_mat_int

!---------------------------------------------------------------------------------
! contains

subroutine calc_rad_at_pt (ele, orb0, z_pos, damp_mat1, stoc_mat1)

type (ele_struct) ele, runt
type (coord_struct) orb0, orbz  ! Orbit at start, orbit at z.

real(rp) z_pos, damp_mat1(6,6), stoc_mat1(6,6), g(3), dg(3,3), g2, dg2_dx, dg2_dy
real(rp) mb(6,6), kd, kf, rel_p, v(6)

integer i, j

logical err_flag

! Note: g from g_bending_strength_from_em_field is g of the actual particle and not the ref particle
! so g has a factor of 1/(1 + pz) in it.

call create_element_slice (runt, ele, z_pos, 0.0_rp, ele%branch%param, .false., .false., err_flag, pointer_to_next_ele(ele, -1))
call track1(orb0, runt, ele%branch%param, orbz)
call make_mat6(runt, ele%branch%param, orb0, orbz)
mb = mat_symp_conj(runt%mat6)   ! matrix from z_pos back to 0

call g_bending_strength_from_em_field (ele, ele%branch%param, z_pos, orbz, .true., g, dg)
v = orbz%vec
rel_p = 1 + v(6)
g = g * rel_p
dg = dg * rel_p

g2 = sum(g*g)
kd = kd_coef
dg2_dx = 2 * dot_product(g, dg(:,1))
dg2_dy = 2 * dot_product(g, dg(:,2))

damp_mat1(:,1) = -kd_coef * (mb(:,2)*dg2_dx*v(2)*rel_p + mb(:,4)*dg2_dx*v(4)*rel_p + mb(:,6)*dg2_dx*rel_p**2)
damp_mat1(:,2) = -kd_coef * (mb(:,2)*g2*rel_p) 
damp_mat1(:,3) = -kd_coef * (mb(:,2)*dg2_dy*v(2)*rel_p + mb(:,4)*dg2_dy*v(4)*rel_p + mb(:,6)*dg2_dy*rel_p**2)
damp_mat1(:,4) = -kd_coef * (mb(:,4)*g2*rel_p)
damp_mat1(:,5) = -kd_coef * (0)
damp_mat1(:,6) = -kd_coef * (mb(:,2)*g2*v(2) + mb(:,4)*g2*v(4) + mb(:,6)*2*g2*rel_p)

kf = kf_coef * rel_p**2 * sqrt(g2)**3
forall (i = 1:6)
  stoc_mat1(i,:) = kf * (rel_p**2 * mb(i,6)*mb(:,6) + v(2)**2 * mb(i,2)*mb(:,2) + v(4)**2 * mb(i,4)*mb(:,4) + &
                                rel_p * v(2) * (mb(i,2) * mb(:,6) + mb(i,6) * mb(:,2)) + &
                                rel_p * v(4) * (mb(i,4) * mb(:,6) + mb(i,6) * mb(:,4)) + &
                                v(2) * v(4) * (mb(i,2) * mb(:,4) + mb(i,4) * mb(:,2)))
end forall

end subroutine calc_rad_at_pt

end subroutine damping_and_stochastic_rad_mats

end module
