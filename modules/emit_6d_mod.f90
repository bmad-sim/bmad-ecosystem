! Note: A negative emittance is possible and just means that the beam is
! unstable. That is, adamping partition number is negative.
!



!+
! Module emit_6d_mod
!
! Module for calculating the equilibrium sigma matrix for a closed geometry lattice.
!-

module emit_6d_mod

use bmad

implicit none

type rad1_map_struct
  real(rp) :: damp_mat(6,6) = 0    ! Transfer matrix with damping.
  real(rp) :: stoc_mat(6,6) = 0    ! Stochastic matrix.
end type

contains

!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!+
! Subroutine damping_and_stochastic_rad_mats (ele1, ele2, damp_mat, stoc_mat)
!
! Routine to calculate the damping and stochastic variance matrices from exit end of ele1
! to the exit end of ele2. Use ele1 = ele2 to get 1-turn matrices.
! If ele2 is before ele1 the integration range if from ele1 to the branch end plus 
! from the beginning to ele2.
!
! Note: The branch%mat6 matrices will be remade with damping on.
!
! Input:
!   ele1          -- ele_struct: Start element of integration range.
!   ele2          -- ele_struct: End element of integration range.
!
! Output:
!   damp_mat(6,6) -- real(rp): One-turn transfer matrix with damping.
!   stoc_mat(6,6) -- real(rp): One-turn stochastic variance matrix.
!-

subroutine damping_and_stochastic_rad_mats (ele1, ele2, damp_mat, stoc_mat)

use super_recipes_mod

type (ele_struct), target :: ele1, ele2
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele1_track, ele2_track, ele3
type (bmad_common_struct) bmad_com_save
type (rad1_map_struct), allocatable :: ds(:)
type (coord_struct), allocatable :: closed_orb(:)

real(rp) sig_mat(6,6)
real(rp) damp_mat(6,6), stoc_mat(6,6), damp_typ, stoc_typ, r_e
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
r_e = classical_radius(ele1_track%ref_species)
kd_coef = 2.0_rp * r_e * gamma**3 / 3.0_rp
kf_coef = 55.0_rp * r_e * h_bar_planck * gamma**5 * c_light / (24.0_rp * sqrt_3 * mass_of(ele1_track%ref_species))
damp_typ = kd_coef * g2_ave / branch%param%total_length
stoc_typ = kf_coef * g3_ave / branch%param%total_length

! Calculate element-by-element damping and stochastic mats.

do ie = 1, branch%n_ele_track
  call qromb_rad_mat_int(branch%ele(ie), closed_orb(ie-1), closed_orb(ie), damp_typ, stoc_typ, ds(ie)%damp_mat, ds(ie)%stoc_mat)
enddo

!

call mat_make_unit(damp_mat)
stoc_mat = 0

ie = ele1_track%ix_ele
do 
  ie = ie + 1
  if (ie > branch%n_ele_track) ie = 0
  if (ie /= 0) then
    ele3 => branch%ele(ie)
    damp_mat = matmul(ds(ie)%damp_mat, damp_mat)
    stoc_mat = matmul(matmul(ele3%mat6, stoc_mat), mat_symp_conj(ele3%mat6)) + ds(ie)%stoc_mat
  endif
  if (ie == ele2%ix_ele) exit
enddo

!---------------------------------------------------------------------------------
contains

! This is adapted from qromb and trapzd from Numerical Recipes.

subroutine qromb_rad_mat_int (ele, orb_in, orb_out, damp_typ, stoc_typ, damp_mat, stoc_mat)

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
type (fringe_field_info_struct) fringe_info
type (qromb_pt1_struct) pt1(0:16)
type (qromb_int_struct) qi(0:4)
type (bmad_common_struct) bmad_com_save

real(rp) damp_mat(6,6), damp_mat1(6,6), stoc_mat1(6,6), stoc_mat(6,6), damp_typ, stoc_typ
real(rp) mat0(6,6), mat1(6,6), ddamp(6,6), dstoc(6,6), damp_mat_sum(6,6), stoc_mat_sum(6,6), mat0_inv(6,6)
real(rp) del_z, l_ref, rel_tol, eps_damp, eps_stoc, z_pos, d_max

integer j, j1, i1, i2, n, j_min_test, ll, n_pts
integer :: j_max = 10

! No radiation cases

if (ele%value(l$) == 0 .or. (orb_out%vec(2) == orb_in%vec(2) .and. &
                         orb_out%vec(4) == orb_in%vec(4) .and. ele%key /= sbend$)) then
  damp_mat = ele%mat6
  stoc_mat = 0
  return
endif

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
eps_damp = rel_tol * ele%value(l$) * damp_typ 
eps_stoc = rel_tol * ele%value(l$) * stoc_typ 

j_min_test = 3

if (ele%key == wiggler$ .or. ele%key == undulator$) then
  j_min_test = 5
endif

damp_mat = 0
stoc_mat = 0

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

!

mat0_inv = mat_symp_conj(mat0)
stoc_mat = matmul(matmul(mat0_inv, stoc_mat), mat0)
stoc_mat = matmul(matmul(ele%mat6, stoc_mat), mat_symp_conj(ele%mat6))
damp_mat = ele%mat6 + matmul(ele%mat6, matmul(mat0_inv, damp_mat)) 

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

!

call create_element_slice (runt, ele, z_pos, 0.0_rp, ele%branch%param, .false., .false., err_flag, pointer_to_next_ele(ele, -1))
call track1(orb0, runt, ele%branch%param, orbz)
call make_mat6(runt, ele%branch%param, orb0, orbz)
mb = mat_symp_conj(runt%mat6)   ! matrix from z_pos back to 0

call g_bending_strength_from_em_field (ele, ele%branch%param, z_pos, orbz, .true., g, dg)
v = orbz%vec
rel_p = 1 + v(6)
g2 = sum(g*g)
kd = kd_coef
dg2_dx = 2 * dot_product(g, dg(:,1))
dg2_dy = 2 * dot_product(g, dg(:,2))

damp_mat1(:,2) = -kd_coef * (dg2_dx*v(2)*rel_p*mb(:,1) + g2*rel_p*mb(:,2) + dg2_dy*v(2)*rel_p*mb(:,3) + g2*v(2)*mb(:,6))
damp_mat1(:,4) = -kd_coef * (dg2_dx*v(4)*rel_p*mb(:,1) + dg2_dy*v(4)*rel_p*mb(:,3) + g2*rel_p*mb(:,4) + g2*v(4)*mb(:,6))
damp_mat1(:,6) = -kd_coef * (dg2_dx*rel_p**2*mb(:,1) + dg2_dy*v(4)*rel_p**2*mb(:,3) + 2*g2*rel_p*mb(:,6))

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
