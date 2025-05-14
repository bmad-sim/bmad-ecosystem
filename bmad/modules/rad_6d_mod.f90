! Note: A negative emittance is possible and just means that the beam is unstable.
! That is, the corresponding damping partition number is negative.

!+
! Module rad_6d_mod
!
! Module for 6D radiation calculations. EG: the equilibrium sigma matrix for a closed geometry lattice.
!-

module rad_6d_mod

use bmad_routine_interface

implicit none

type private_stash_struct
  real(rp) :: gamma0            ! Relativistic gamma factor.
  real(rp) :: eta_x_coef(0:3)   ! Dispersion interpolation coefs.
  real(rp) :: eta_y_coef(0:3)   ! Dispersion interpolation coefs.
end type

private :: private_stash_struct

contains

!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!+
! Subroutine emit_6d (ele_ref, include_opening_angle, mode, sigma_mat, closed_orbit, rad_int_by_ele) 
!
! Routine to calculate the three normal mode emittances, damping partition numbers, radiation integrals, etc. 
! Since the emattances, etc. are only an invariant in the limit of zero damping, the calculated
! values will vary depending upon the reference element.
!
! If the lattice geometry is open, only the radiation integrals is computed.
!
! Input:
!   ele_ref               -- ele_struct: Origin of the 1-turn maps used to evaluate the emittances.
!   include_opening_angle -- logical: If True include the effect of the vertical opening angle of emitted radiation.
!                             Generally use True unless comparing against other codes.
!   closed_orbit(0:)      -- coord_struct, optional: Closed orbit. If not present this routine will calculate it.
!
! Output:
!   mode            -- normal_modes_struct: Emittance and other info.
!   sigma_mat(6,6)  -- real(rp): Sigma matrix.
!   rad_int_by_ele  -- rad_int_all_ele_struct, optional: Radiation integrals element-by-element.
!-

subroutine emit_6d (ele_ref, include_opening_angle, mode, sigma_mat, closed_orbit, rad_int_by_ele)

use f95_lapack, only: dgesv_f95

type (ele_struct), target :: ele_ref
type (ele_struct), pointer :: ele
type (coord_struct) orbit
type (lat_struct), pointer :: lat
type (branch_struct), pointer :: branch
type (normal_modes_struct) mode
type (rad_map_struct) rmap
type (coord_struct), optional, target :: closed_orbit(0:)
type (rad_int_all_ele_struct), target, optional :: rad_int_by_ele
type (rad_int_branch_struct), pointer :: ri_branch

real(rp) sigma_mat(6,6), rf65, sig_s(6,6), mat6(6,6), xfer_nodamp_mat(6,6), mt(21,21), v_sig(21,1)

complex(rp) eval(6), evec(6,6), cmat(6,6)

integer i, j, k, ib, ipev(21), info

integer, parameter :: w1(21) = [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 5, 5, 6]
integer, parameter :: w2(21) = [1, 2, 3, 4, 5, 6, 2, 3, 4, 5, 6, 3, 4, 5, 6, 4, 5, 6, 5, 6, 6]
integer, parameter :: v(6,6) = reshape( &
            [1,  2,  3,  4,  5,  6,   2,  7,  8,  9, 10, 11,   3,  8, 12, 13, 14, 15, &
             4,  9, 13, 16, 17, 18,   5, 10, 14, 17, 19, 20,   6, 11, 15, 18, 20, 21], [6,6])

logical include_opening_angle, err, rf_off

! Allocate rad_int_by_ele

branch => pointer_to_branch(ele_ref)
lat => branch%lat

if (present(rad_int_by_ele)) then
  if (allocated(rad_int_by_ele%branch)) then
    if (ubound(rad_int_by_ele%branch, 1) /= ubound(lat%branch, 1)) deallocate (rad_int_by_ele%branch)
    if (allocated(rad_int_by_ele%branch)) then
      do ib = 0, ubound(lat%branch, 1)
        ri_branch => rad_int_by_ele%branch(ib)
        if (allocated(ri_branch%ele)) then
          if (ubound(ri_branch%ele, 1) == lat%branch(ib)%n_ele_max) cycle
          deallocate (ri_branch%ele)
        endif
        if (.not. allocated(ri_branch%ele)) allocate (ri_branch%ele(0:lat%branch(ib)%n_ele_max))
      enddo
    endif
  endif

  if (.not. allocated(rad_int_by_ele%branch)) then
    allocate(rad_int_by_ele%branch(0:ubound(lat%branch, 1)))
    do ib = 0, ubound(lat%branch, 1)
      allocate (rad_int_by_ele%branch(ib)%ele(0:lat%branch(ib)%n_ele_max))
    enddo
  endif
endif

! Analysis is documented in the Bmad manual.

mode = normal_modes_struct()
sigma_mat = 0

if (branch%param%total_length == 0) return   ! Test lattices can have zero length.

if (present(rad_int_by_ele)) then
  ri_branch => rad_int_by_ele%branch(branch%ix_branch)
  call rad_damp_and_stoc_mats (ele_ref, ele_ref, include_opening_angle, rmap, mode, xfer_nodamp_mat, &
                                                            err, closed_orbit, ri_branch)
  mode%synch_int(0) = sum(ri_branch%ele%i0)
  mode%synch_int(1) = sum(ri_branch%ele%i1)
  mode%synch_int(2) = sum(ri_branch%ele%i2)
  mode%synch_int(3) = sum(ri_branch%ele%i3)
else
  call rad_damp_and_stoc_mats (ele_ref, ele_ref, include_opening_angle, rmap, mode, xfer_nodamp_mat, &
                                                            err, closed_orbit)
endif

if (err) return
if (branch%param%geometry == open$) return

! If there is no RF then add a small amount to enable the calculation to proceed.
! The RF is modeled as a unit matrix with M(6,5) = 1d-4.

rf_off = (rmap%xfer_damp_mat(6,5) == 0)
if (rf_off) then
  rf65 = 1e-4
  rmap%xfer_damp_mat(6,:) = rmap%xfer_damp_mat(6,:) + rf65 * rmap%xfer_damp_mat(5,:)
  rmap%stoc_mat(6,:) = rmap%stoc_mat(6,:) + rf65 * rmap%stoc_mat(5,:)
  rmap%stoc_mat(:,6) = rmap%stoc_mat(:,6) + rmap%stoc_mat(:,5) * rf65
endif

! The 6x6 sigma matrix equation is recast as a linear 21x21 matrix equation and solved.

call mat_make_unit(mt)

do i = 1, 21
  v_sig(i,1) = rmap%stoc_mat(w1(i), w2(i))

  do j = 1, 6
  do k = 1, 6
    mt(i,v(j,k)) = mt(i,v(j,k)) - rmap%xfer_damp_mat(w1(i),j) * rmap%xfer_damp_mat(w2(i),k)
  enddo
  enddo
enddo

call dgesv_f95(mt, v_sig, ipev, info)

if (info /= 0) then
  sigma_mat = -1
  if (include_opening_angle) then
    mode%a%emittance = -1; mode%b%emittance = -1; mode%z%emittance = -1
  else
    mode%a%emittance_no_vert = -1; mode%b%emittance_no_vert = -1; mode%z%emittance_no_vert = -1
  endif
  return
endif

do j = 1, 6
do k = 1, 6
  sigma_mat(j,k) = v_sig(v(j,k), 1)
enddo
enddo

sig_s(:,1) = -sigma_mat(:,2)
sig_s(:,2) =  sigma_mat(:,1)
sig_s(:,3) = -sigma_mat(:,4)
sig_s(:,4) =  sigma_mat(:,3)
sig_s(:,5) = -sigma_mat(:,6)
sig_s(:,6) =  sigma_mat(:,5)

call mat_eigen(sig_s, eval, evec, err)

if (include_opening_angle) then
  mode%a%emittance = aimag(eval(1))
  mode%b%emittance = aimag(eval(3))
  mode%z%emittance = aimag(eval(5))
else
  mode%a%emittance_no_vert = aimag(eval(1))
  mode%b%emittance_no_vert = aimag(eval(3))
  mode%z%emittance_no_vert = aimag(eval(5))
endif

mode%e_loss = -mode%dpz_damp * ele_ref%value(e_tot$)

if (sigma_mat(5,5) < 0 .or. sigma_mat(6,6) < 0) then
  mode%sig_z = -1
  mode%sigE_E = -1
else
  mode%sig_z = sqrt(sigma_mat(5,5))
  mode%sigE_E = sqrt(sigma_mat(6,6))
endif

! This is Eq (86) from Ohmi, Hirata, & Oide, "From the beam-envelope matrix to synchrotron-radiation integrals".

call mat_eigen (xfer_nodamp_mat, eval, evec, err)
evec = transpose(evec)
call cplx_mat_inverse(evec, cmat)
call mat_inverse(xfer_nodamp_mat, mat6)
cmat = matmul(matmul(matmul(cmat, rmap%damp_dmat), mat6), evec)

mode%a%alpha_damp = -real(cmat(1,1), 8)
mode%b%alpha_damp = -real(cmat(3,3), 8)
mode%z%alpha_damp = -real(cmat(5,5), 8)

! E_loss = 0 can happen in toy lattices without bends.

if (mode%e_loss /= 0) then
  mode%a%j_damp = -2 * mode%a%alpha_damp / mode%dpz_damp
  mode%b%j_damp = -2 * mode%b%alpha_damp / mode%dpz_damp
  mode%z%j_damp = -2 * mode%z%alpha_damp / mode%dpz_damp
endif

if (rf_off) then
  mode%z%emittance = -1
endif

end subroutine emit_6d

!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!+
! Subroutine rad_damp_and_stoc_mats (ele1, ele2, include_opening_angle, rmap, mode, xfer_nodamp_mat, err_flag, closed_orbit, rad_int_branch)
!
! Routine to calculate the damping and stochastic variance matrices from exit end of ele1
! to the exit end of ele2. Use ele1 = ele2 to get 1-turn matrices.
!
! If ele2 is before ele1 the integration range if from ele1 to the branch end plus 
! from the beginning to ele2.
!
! Note: The ele%mat6 matrices will be remade. By convention, these matrices
! do not include damping.
!
! Input:
!   ele1                  -- ele_struct: Start element of integration range.
!   ele2                  -- ele_struct: End element of integration range.
!   include_opening_angle -- logical: If True include the effect of the vertical opening angle of emitted radiation.
!                             Generally use True unless comparing against other codes.
!   closed_orbit(0:)      -- coord_struct, optional: Closed orbit. If not present this routine will calculate it.
!
! Output:
!   rmap                  -- rad_map_struct: Damping and stochastic mats 
!     %stoc_mat               --  stochastic variance matrix.
!   mode                  -- normal_modes_struct:
!     %dpz_damp                 -- Change in pz without RF.
!     %pz_average               -- Average pz due to damping.
!   xfer_nodamp_mat(6,6)  -- real(rp): Transfer matrix without damping.
!   rad_int_branch        -- rad_int_branch_struct, optional: Array of element-by-element radiation integrals.
!   err_flag              -- logical: Set true if there is a problem.
!-

subroutine rad_damp_and_stoc_mats (ele1, ele2, include_opening_angle, rmap, mode, xfer_nodamp_mat, err_flag, closed_orbit, rad_int_branch)

type (ele_struct), allocatable :: eles_save(:)
type (ele_struct), target :: ele1, ele2
type (ele_struct), pointer :: ele1_track, ele2_track, ele3
type (rad_map_struct) rmap
type (normal_modes_struct) mode
type (rad_int_branch_struct), optional :: rad_int_branch
type (branch_struct), pointer :: branch
type (bmad_common_struct) bmad_com_save
type (rad_map_struct), allocatable :: rm1(:)
type (coord_struct), optional, target :: closed_orbit(0:)
type (coord_struct), pointer :: closed_orb(:)
type (coord_struct), allocatable, target :: closed(:)

real(rp) sig_mat(6,6), mt(6,6), xfer_nodamp_mat(6,6), tol, length
integer ie

logical include_opening_angle, err_flag

character(*), parameter :: r_name = 'rad_damp_and_stoc_mats'

!

call find_element_ends(ele1, ele3, ele1_track)
call find_element_ends(ele2, ele3, ele2_track)

branch => ele1_track%branch
allocate (eles_save(0:ubound(branch%ele,1)))
call transfer_eles(branch%ele, eles_save)

allocate (rm1(branch%n_ele_track))

bmad_com_save = bmad_com
bmad_com%radiation_fluctuations_on = .false.
bmad_com%spin_tracking_on = .false.

if (present(closed_orbit)) then
  closed_orb => closed_orbit
  err_flag = .false.
elseif (rf_is_on(branch)) then
  bmad_com%radiation_damping_on = .true.
  call closed_orbit_calc(branch%lat, closed, 6, +1, branch%ix_branch, err_flag)
  closed_orb => closed
else
  bmad_com%radiation_damping_on = .false.
  call closed_orbit_calc(branch%lat, closed, 4, +1, branch%ix_branch, err_flag)
  closed_orb => closed
endif

bmad_com = bmad_com_save

if (err_flag) then
  call out_io (s_error$, r_name, 'Closed orbit calculation failing.', &
                                   'Emittance and other radiation related quantities will not be calculated.')
  goto 9000
endif

if (.not. present(closed_orbit)) then
  call lat_make_mat6 (ele1_track%branch%lat, -1, closed_orb, branch%ix_branch, err_flag)
  if (err_flag) goto 9000  ! Restore and return
endif

! Calculate element-by-element damping and stochastic mats.

tol = 1d-4 / branch%param%total_length
do ie = 1, branch%n_ele_track
  if (closed_orb(ie)%state /= alive$) then
    call out_io (s_error$, r_name, 'Closed orbit calculation failing.', &
                                   'Emittance and other radiation related quantities will not be calculated.')
    goto 9000
  endif

  if (present(rad_int_branch)) then
    call rad1_damp_and_stoc_mats(branch%ele(ie), include_opening_angle, closed_orb(ie-1), closed_orb(ie), rm1(ie), &
                                 tol * branch%param%g2_integral, tol * branch%param%g3_integral, err_flag, branch%ele(ie-1), rad_int_branch%ele(ie))
  else
    call rad1_damp_and_stoc_mats(branch%ele(ie), include_opening_angle, closed_orb(ie-1), closed_orb(ie), rm1(ie), &
                                 tol * branch%param%g2_integral, tol * branch%param%g3_integral, err_flag)
  endif
  if (err_flag) goto 9000
enddo

!

rmap = rad_map_struct(-1, 0, 0, mat6_unit$, 0)
xfer_nodamp_mat = mat6_unit$
mode%dpz_damp = 0
mode%pz_average = 0
length = 0

ie = ele1_track%ix_ele
do 
  ie = ie + 1
  if (ie > branch%n_ele_track) ie = 0
  if (ie /= 0) then
    ele3 => branch%ele(ie)
    mt = rm1(ie)%damp_dmat + ele3%mat6
    rmap%damp_dmat = matmul(mt, rmap%damp_dmat) + matmul(rm1(ie)%damp_dmat, rmap%xfer_damp_mat)
    rmap%xfer_damp_vec = matmul(mt, rmap%xfer_damp_vec) + rm1(ie)%xfer_damp_vec
    rmap%xfer_damp_mat = matmul(mt, rmap%xfer_damp_mat)
    rmap%stoc_mat = matmul(matmul(mt, rmap%stoc_mat), transpose(mt)) + rm1(ie)%stoc_mat
    xfer_nodamp_mat = matmul(ele3%mat6, xfer_nodamp_mat)
    mode%pz_average = mode%pz_average + 0.5_rp * ele3%value(l$) * (closed_orb(ie-1)%vec(6) + closed_orb(ie)%vec(6))
    length = length + ele3%value(l$)
    if (ele3%key /= rfcavity$) then
      mode%dpz_damp = mode%dpz_damp + rm1(ie)%xfer_damp_vec(6)
    endif
  endif
  if (ie == ele2%ix_ele) exit
enddo

if (length /= 0) mode%pz_average = mode%pz_average / length

9000 continue
call transfer_eles(eles_save, branch%ele)

end subroutine rad_damp_and_stoc_mats

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!+
! Subroutine rad1_damp_and_stoc_mats (ele, include_opening_angle, orb_in, orb_out, rad_map, g2_tol, g3_tol, err_flag, ele0, rad_int1)
!
! Routine to calculate the damping and stochastic matrices for a given lattice element.
!
! Input:
!   ele                   -- ele_struct: Element under consideration.
!   include_opening_angle -- logical: If True include the effect of the vertical opening angle of emitted radiation.
!                             Generally use True unless comparing against other codes.
!   orb_in                -- coord_struct: Entrance orbit about which to compute the matrices.
!   orb_out               -- coord_struct: Exit orbit.
!   g2_tol                -- real(rp): Tollerance on g^2 per unit length (damping tolerance).
!   g3_tol                -- real(rp): Tollerance on g^3 per unit length (stocastic tolerance).
!   ele0                  -- ele_struct, optional: Element before `ele`. Needed if and only if rad_int1 is present
!
! Output:
!   rad_map               -- rad_map_struct: Damping and stochastic matrices.
!     %stoc_mat             -- Variance matrix.
!
!   err_flag              -- logical: Set true if there is an error. False otherwise.
!   rad_int1              -- rad_int1_struct, optional: Radiation integrals
!-

subroutine rad1_damp_and_stoc_mats (ele, include_opening_angle, orb_in, orb_out, rad_map, g2_tol, g3_tol, err_flag, ele0, rad_int1)

use super_recipes_mod, only: super_polint

!

type qromb_int_struct
  real(rp) :: h = 0
  real(rp) :: xfer_damp_vec(6) = 0
  real(rp) :: damp_dmat(6,6) = 0
  real(rp) :: stoc_mat(6,6) = 0
  type (rad_int1_struct) :: ri = rad_int1_struct() 
end type

type (ele_struct) ele
type (ele_struct), optional :: ele0
type (coord_struct) orb_in, orb_out, orb0, orb1, orb_end
type (rad_int1_struct), optional :: rad_int1
type (rad_map_struct) :: rad_map
type (fringe_field_info_struct) fringe_info
type (qromb_int_struct) qi(0:4)
type (bmad_common_struct) bmad_com_save
type (private_stash_struct) stash

real(rp) damp_dmat1(6,6), stoc_var_mat1(6,6), g2_tol, g3_tol
real(rp) mat0(6,6), mat1(6,6), ddamp(6,6), dstoc(6,6), damp_dmat_sum(6,6), stoc_var_mat_sum(6,6), mat0_inv(6,6)
real(rp) del_z, l_ref, rel_tol, eps_damp, eps_stoc, z_pos, d_max, mat_end(6,6), damp_vec1(6), damp_vec_sum(6), ddvec(6)
real(rp) kd_coef, kf_coef, radi, length, dum, a, b
real(rp), parameter :: cd = 2.0_rp / 3.0_rp, cf = 55.0_rp * h_bar_planck * c_light / (24.0_rp * sqrt_3)

integer j, j1, i1, i2, n, j_min_test, ll, n_pts
integer :: j_max = 10

logical include_opening_angle, err_flag, err, save_orb_mat

! No radiation cases

rad_map = rad_map_struct(-1, 0, 0, mat6_unit$, 0)
rad_map%ref_orb = orb_out%vec
length = ele%value(l$)
if (length == 0 .or. ele%key == taylor$) return
if (orb_out%vec(2) == orb_in%vec(2) .and. orb_out%vec(4) == orb_in%vec(4) .and. &
       ele%key /= sbend$ .and. ele%key /= rf_bend$ .and. ele%key /= wiggler$ .and. ele%key /= undulator$) return

!

bmad_com_save = bmad_com
bmad_com%radiation_fluctuations_on = .false.
bmad_com%radiation_damping_on = .false.
err_flag = .true.

! Offsets and fringes at upstream end.

orb0 = orb_in

call mat_make_unit(mat0)
call offset_particle (ele, set$, orb0, set_hvkicks = .false., mat6 = mat0, make_matrix = .true.)
call init_fringe_info (fringe_info, ele)
if (fringe_info%has_fringe) then
  fringe_info%particle_at = first_track_edge$
  call apply_element_edge_kick(orb0, fringe_info, ele, ele%branch%param, .false., mat0, .true.)
endif

! Integrate through body.
! This is adapted from qromb and trapzd from Numerical Recipes.

qi(0) = qromb_int_struct()
qi(0)%h = 4

call convert_pc_to (ele%value(p0c$), ele%ref_species, gamma = stash%gamma0)

if (present(rad_int1)) then
  a = ele%x%eta - ele0%x%eta - length * ele0%x%deta_ds;  b = (ele%x%deta_ds - ele0%x%deta_ds) * length
  stash%eta_x_coef(0:3) = [ele0%x%eta, ele0%x%deta_ds, (3.0_rp * a - b) / length**2, (b - 2.0_rp * a) / length**3]

  a = ele%y%eta - ele0%y%eta - length * ele0%y%deta_ds;  b = (ele%y%deta_ds - ele0%y%deta_ds) * length
  stash%eta_y_coef(0:3) = [ele0%y%eta, ele0%y%deta_ds, (3.0_rp * a - b) / length**2, (b - 2.0_rp * a) / length**3]
endif

radi = classical_radius(ele%ref_species)
kd_coef = cd * radi * stash%gamma0**3
kf_coef = cf * radi * stash%gamma0**5 / mass_of(ele%ref_species)

eps_damp = kd_coef * g2_tol / length
eps_stoc = kf_coef * g3_tol / length

j_min_test = 3
if (ele%key == wiggler$ .or. ele%key == undulator$) then
  j_min_test = 5
endif

do j = 1, j_max
  if (j == 1) then
    n_pts = 2
    del_z = length
    l_ref = 0
  else
    n_pts = 2**(j-2)
    del_z = ele%value(l$) / n_pts
    l_ref = del_z / 2
  endif

  damp_dmat_sum = 0
  damp_vec_sum = 0
  stoc_var_mat_sum = 0
  if (present(rad_int1)) rad_int1 = rad_int1_struct()

  do n = 1, n_pts
    z_pos = l_ref + (n-1) * del_z
    save_orb_mat = (j == 1 .and. n == 1)  ! Save if z-position at end of element
    call calc_rad_at_pt(ele, stash, include_opening_angle, orb0, z_pos, damp_dmat1, damp_vec1, stoc_var_mat1, &
                                                                        save_orb_mat, orb_end, mat_end, err, rad_int1)
    if (err) return
    damp_vec_sum = damp_vec_sum + damp_vec1
    damp_dmat_sum = damp_dmat_sum + damp_dmat1
    stoc_var_mat_sum = stoc_var_mat_sum + stoc_var_mat1
  enddo

  j1 = min(j, 4)
  if (j > 4) qi(0:3) = qi(1:4)
  qi(j1)%h = 0.25_rp * qi(j1-1)%h
  qi(j1)%damp_dmat = 0.5_rp * (qi(j1-1)%damp_dmat + del_z * damp_dmat_sum)
  qi(j1)%xfer_damp_vec = 0.5_rp * (qi(j1-1)%xfer_damp_vec + del_z * damp_vec_sum)
  qi(j1)%stoc_mat = 0.5_rp * (qi(j1-1)%stoc_mat + del_z * stoc_var_mat_sum)

  if (present(rad_int1)) then
    qi(j1)%ri%i0 = 0.5_rp * (qi(j1-1)%ri%i0 + del_z * rad_int1%i0)
    qi(j1)%ri%i1 = 0.5_rp * (qi(j1-1)%ri%i1 + del_z * rad_int1%i1)
    qi(j1)%ri%i2 = 0.5_rp * (qi(j1-1)%ri%i2 + del_z * rad_int1%i2)
    qi(j1)%ri%i3 = 0.5_rp * (qi(j1-1)%ri%i3 + del_z * rad_int1%i3)
  endif

  if (j < j_min_test) cycle

  do i1 = 1, 6
    call super_polint(qi(1:j1)%h, qi(1:j1)%xfer_damp_vec(i1), 0.0_rp, damp_vec1(i1), ddvec(i1))
    do i2 = 1, 6
      call super_polint(qi(1:j1)%h, qi(1:j1)%damp_dmat(i1,i2), 0.0_rp, damp_dmat1(i1,i2), ddamp(i1,i2))
      call super_polint(qi(1:j1)%h, qi(1:j1)%stoc_mat(i1,i2), 0.0_rp, stoc_var_mat1(i1,i2), dstoc(i1,i2))
    enddo
  enddo

  d_max = max(eps_damp*maxval(abs(ddamp)), eps_damp*maxval(abs(ddvec)), eps_stoc*maxval(abs(dstoc)))
  if (d_max < 1) exit
enddo

if (present(rad_int1)) then
  call super_polint(qi(1:j1)%h, qi(1:j1)%ri%i0, 0.0_rp, rad_int1%i0, dum)
  call super_polint(qi(1:j1)%h, qi(1:j1)%ri%i1, 0.0_rp, rad_int1%i1, dum)
  call super_polint(qi(1:j1)%h, qi(1:j1)%ri%i2, 0.0_rp, rad_int1%i2, dum)
  call super_polint(qi(1:j1)%h, qi(1:j1)%ri%i3, 0.0_rp, rad_int1%i3, dum)
endif

! Add bend edge if needed

if (ele%key == sbend$ .or. ele%key == rf_bend$) then
  if (ele%orientation == 1) then
    call add_bend_edge_to_damp_mat(damp_dmat1, ele, ele%value(e1$), orb0)
    call add_bend_edge_to_damp_mat(damp_dmat1, ele, ele%value(e2$), orb_end, mat_end)
  else
    call add_bend_edge_to_damp_mat(damp_dmat1, ele, ele%value(e2$), orb0)
    call add_bend_edge_to_damp_mat(damp_dmat1, ele, ele%value(e1$), orb_end, mat_end)
  endif
endif

! Reference position for constructing the matrices in this routine is the upstream end but
! the reference of the output is the downstream end. So need to correct for this.

mat0_inv = mat_symp_conj(mat0)   ! mat0 is transport matrix through the upstream edge

rad_map%xfer_damp_vec = matmul(matmul(ele%mat6, mat0_inv), damp_vec1)
rad_map%damp_dmat     = matmul(matmul(matmul(ele%mat6, mat0_inv), damp_dmat1), mat0)

stoc_var_mat1    = matmul(matmul(mat0_inv, stoc_var_mat1), transpose(mat0_inv))
rad_map%stoc_mat = matmul(matmul(ele%mat6, stoc_var_mat1), transpose(ele%mat6))

bmad_com = bmad_com_save
err_flag = .false.

!---------------------------------------------------------------------------------
contains

subroutine calc_rad_at_pt (ele, stash, include_opening_angle, orb0, z_pos, damp_dmat1, damp_vec1, &
                        stoc_var_mat1, save_orb_mat, orb_save, mat_save, err, rad_int1)

type (ele_struct) ele0, ele, runt
type (rad_int1_struct), optional :: rad_int1
type (coord_struct) orb0, orbz, orb_save  ! Orbit at start, orbit at z.
type (private_stash_struct) stash

real(rp) z_pos, damp_dmat1(6,6), damp_vec1(6), stoc_var_mat1(6,6), g(3), dg(3,3), g1, g2, g3, dg2_dx, dg2_dy
real(rp) mb(6,6), mb_inv(6,6), kf, kv, rel_p, v(6), mat_save(6,6), eta_vec(2)

integer i, j

logical include_opening_angle, save_orb_mat, err

! Note: g from g_bending_strength_from_em_field is g of the particle and not the zero pz particle.
! So g will have a factor of (1 + pz) but the equations for the mats was developed for g of the zero pz particle

call create_element_slice (runt, ele, z_pos, 0.0_rp, ele%branch%param, .false., .false., err, pointer_to_next_ele(ele, -1))
if (err) return
runt%value(dispatch$) = no_misalignment$
runt%old_value(dispatch$) = no_misalignment$   ! So no bookkeeping is done.
call track1(orb0, runt, ele%branch%param, orbz)
call make_mat6(runt, ele%branch%param, orb0, orbz)
if (save_orb_mat) then
  orb_save = orbz
  mat_save = runt%mat6
endif

if (orbz%state /= alive$) then
  err = .true.
  return
endif

mb = mat_symp_conj(runt%mat6)   ! matrix from z_pos back to 0

call g_bending_strength_from_em_field (ele, ele%branch%param, z_pos, orbz, .true., g, dg)
v = orbz%vec
rel_p = 1 + v(6)
g = g * rel_p
dg = dg * rel_p

g1 = norm2(g)
g2 = g1**2
g3 = g1 * g2
dg2_dx = 2 * dot_product(g, dg(:,1))
dg2_dy = 2 * dot_product(g, dg(:,2))

if (ele%key == sbend$ .or. ele%key == rf_bend$) then
  g1 = g1 * (1 + ele%value(g$) * orb0%vec(1))  ! Variation in path length effect
  g2 = g2 * (1 + ele%value(g$) * orb0%vec(1))  ! Variation in path length effect
  g3 = g3 * (1 + ele%value(g$) * orb0%vec(1))  ! Variation in path length effect
  dg2_dx = dg2_dx * (1 + ele%value(g$) * orb0%vec(1)) + g2*ele%value(g$)
endif

! 

if (present(rad_int1)) then
  eta_vec = [poly_eval(stash%eta_x_coef, z_pos), poly_eval(stash%eta_y_coef, z_pos)]
  rad_int1%i0 = rad_int1%i0 + g1 * stash%gamma0
  rad_int1%i1 = rad_int1%i1 + sum(g(1:2) * eta_vec) 
  rad_int1%i2 = rad_int1%i2 + g2
  rad_int1%i3 = rad_int1%i3 + g3
endif

! Damping matrix

damp_dmat1(:,1) = -kd_coef * (mb(:,2)*dg2_dx*v(2)*rel_p + mb(:,4)*dg2_dx*v(4)*rel_p + mb(:,6)*dg2_dx*rel_p**2)
damp_dmat1(:,2) = -kd_coef * (mb(:,2)*g2*rel_p) 
damp_dmat1(:,3) = -kd_coef * (mb(:,2)*dg2_dy*v(2)*rel_p + mb(:,4)*dg2_dy*v(4)*rel_p + mb(:,6)*dg2_dy*rel_p**2)
damp_dmat1(:,4) = -kd_coef * (mb(:,4)*g2*rel_p)
damp_dmat1(:,5) =  0
damp_dmat1(:,6) = -kd_coef * (mb(:,2)*g2*v(2) + mb(:,4)*g2*v(4) + mb(:,6)*2*g2*rel_p)

! Damping vec

damp_vec1 = -kd_coef * matmul(mb, [0.0_rp, v(2)*g2*rel_p, 0.0_rp, v(4)*g2*rel_p, 0.0_rp, g2*rel_p**2])

! Stochastic matrix

kf = kf_coef * rel_p**2 * g3
forall (i = 1:6)
  stoc_var_mat1(i,:) = kf * (rel_p**2 * mb(i,6)*mb(:,6) + v(2)**2 * mb(i,2)*mb(:,2) + v(4)**2 * mb(i,4)*mb(:,4) + &
                                rel_p * v(2) * (mb(i,2) * mb(:,6) + mb(i,6) * mb(:,2)) + &
                                rel_p * v(4) * (mb(i,4) * mb(:,6) + mb(i,6) * mb(:,4)) + &
                                v(2) * v(4) * (mb(i,2) * mb(:,4) + mb(i,4) * mb(:,2)))
end forall

if (include_opening_angle) then
  kv = kf_coef * 13.0_rp * g1 / (55.0_rp * stash%gamma0**2)
  forall (i = 1:6)
    stoc_var_mat1(i,:) = stoc_var_mat1(i,:) + kv * (g(2)**2 * mb(i,2)*mb(:,2) + g(1)**2 * mb(i,4)*mb(:,4) - &
                                  g(1) * g(2) * (mb(i,2) * mb(:,4) + mb(i,4) * mb(:,2)))
  end forall
endif

end subroutine calc_rad_at_pt

!---------------------------------------------------------------------------------
! contains

subroutine add_bend_edge_to_damp_mat(damp_dmat, ele, e_edge, orb, mat)

type (ele_struct) ele
type (coord_struct) orb
real(rp) damp_dmat(6,6), e_edge, dm(6,6), dg2_dx, rel_p
real(rp), optional :: mat(6,6)

!

if (e_edge == 0) return

dg2_dx = -tan(e_edge) * (ele%value(g$) + ele%value(dg$))**2
rel_p = 1 + orb%vec(6)

dm = 0
dm(2,1:3) = -kd_coef * [dg2_dx * orb%vec(2) * rel_p, 0.0_rp, -dg2_dx * orb%vec(2) * rel_p]
dm(4,1:3) = -kd_coef * [dg2_dx * orb%vec(4) * rel_p, 0.0_rp, -dg2_dx * orb%vec(4) * rel_p]
dm(6,1:3) = -kd_coef * [dg2_dx * rel_p**2,           0.0_rp, -dg2_dx * rel_p**2]

if (present(mat)) then
  dm = matmul(matmul(mat_symp_conj(mat), dm), mat)
endif

damp_dmat = damp_dmat + dm

end subroutine add_bend_edge_to_damp_mat

end subroutine rad1_damp_and_stoc_mats

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!+
! Subroutine rad_g_integrals (ele, where, orb_in, orb_out, int_g, int_g2, int_g3, g_tol, g2_tol, g3_tol)
!
! Routine to calculate bending strength integrals (g(s) = 1/trajectory_bending_radius(s)) in
! laboratory coords.
!
! Input:
!   ele                   -- ele_struct: Element under consideration.
!   where                 -- integer: What part of ele to integrate over. 
!                               upstream$ -> 1st half of element, downsteam$ -> 2nd half, all$ -> everything.
!   orb_in                -- coord_struct: Entrance orbit about which to compute the matrices.
!   orb_out               -- coord_struct: Exit orbit.
!   g_tol                 -- real(rp): Tollerance on |g| per unit length.
!   g2_tol                -- real(rp): Tollerance on g^2 per unit length.
!   g3_tol                -- real(rp): Tollerance on g^3 per unit length.
!
! Output:
!   int_g(2)              -- real(rp): Integrals of (gx,gy) vector.
!   gint_g2, int_g3       -- real(rp): integrals of |g|^2 and |g|^3.
!-

subroutine rad_g_integrals (ele, where, orb_in, orb_out, int_g, int_g2, int_g3, g_tol, g2_tol, g3_tol)

use super_recipes_mod, only: super_polint

!

type qromb_int_struct
  real(rp) :: h = 0
  real(rp) :: int_g(2) = 0,  int_g2 = 0,  int_g3 = 0
end type

type (ele_struct) ele
type (coord_struct) orb_in, orb_out, orb0, orb1, orb_end
type (fringe_field_info_struct) fringe_info
type (qromb_int_struct) qi(0:4)
type (bmad_common_struct) bmad_com_save

real(rp) int_g(2), int_g2, int_g3, g_tol, g2_tol, g3_tol, g_vec(3), g2, g3
real(rp) g_sum(2), g2_sum, g3_sum, dgx, dgy, dg2, dg3, len_int, s0, s1, dg(3,3)
real(rp) del_z, l_ref, rel_tol, eps_damp, eps_stoc, z_pos, d_max, mat_end(6,6)
real(rp) tilt, sin_t, cos_t

integer where, j, j1, i1, i2, n, j_min_test, ll, n_pts
integer :: j_max = 10

logical include_opening_angle, save_orb_mat

! No radiation cases

int_g = 0;  int_g2 = 0;  int_g3 = 0 

if (ele%value(l$) == 0 .or. (orb_out%vec(2) == orb_in%vec(2) .and. &
           orb_out%vec(4) == orb_in%vec(4) .and. ele%key /= sbend$ .and. ele%key /= rf_bend$)) return

!

bmad_com_save = bmad_com
bmad_com%radiation_fluctuations_on = .false.
bmad_com%radiation_damping_on = .false.

! Offsets and fringes at upstream end.

orb0 = orb_in

call offset_particle (ele, set$, orb0, set_hvkicks = .false.)
call init_fringe_info (fringe_info, ele)
if (fringe_info%has_fringe) then
  fringe_info%particle_at = first_track_edge$
  call apply_element_edge_kick(orb0, fringe_info, ele, ele%branch%param, .false.)
endif

select case (where)
case (upstream$)
  s0 = 0
  s1 = 0.5_rp * ele%value(l$)
case (downstream$)
  s0 = 0.5_rp * ele%value(l$)
  s1 = ele%value(l$)
case (all$)
  s0 = 0
  s1 = ele%value(l$)
case default
  call err_exit
end select

len_int = s1 - s0

! Integrate through body.
! This is adapted from qromb and trapzd from Numerical Recipes.

qi(0) = qromb_int_struct()
qi(0)%h = 4

j_min_test = 3
if (ele%key == wiggler$ .or. ele%key == undulator$) then
  j_min_test = 5
endif

do j = 1, j_max
  if (j == 1) then
    n_pts = 2
    del_z = len_int
    l_ref = 0
  else
    n_pts = 2**(j-2)
    del_z = len_int / n_pts
    l_ref = del_z / 2
  endif

  g_sum = 0;  g2_sum = 0; g3_sum = 0
  do n = 1, n_pts
    z_pos = s0 + l_ref + (n-1) * del_z
    save_orb_mat = (j == 1 .and. n == 1)  ! Save if z-position at end of element
    call calc_g_at_pt(ele, orb0, z_pos, g_vec, g2, g3, save_orb_mat, orb_end, mat_end)
    g_sum  = g_sum  + g_vec(1:2)
    g2_sum = g2_sum + g2
    g3_sum = g3_sum + g3
  enddo

  j1 = min(j, 4)
  if (j > 4) qi(0:3) = qi(1:4)
  qi(j1)%h = 0.25_rp * qi(j1-1)%h
  qi(j1)%int_g  = 0.5_rp * (qi(j1-1)%int_g  + del_z * g_sum)
  qi(j1)%int_g2 = 0.5_rp * (qi(j1-1)%int_g2 + del_z * g2_sum)
  qi(j1)%int_g3 = 0.5_rp * (qi(j1-1)%int_g3 + del_z * g3_sum)

  if (j < j_min_test) cycle

  call super_polint(qi(1:j1)%h, qi(1:j1)%int_g(1), 0.0_rp, int_g(1), dgx)
  call super_polint(qi(1:j1)%h, qi(1:j1)%int_g(2), 0.0_rp, int_g(2), dgy)
  call super_polint(qi(1:j1)%h, qi(1:j1)%int_g2,   0.0_rp, int_g2,   dg2)
  call super_polint(qi(1:j1)%h, qi(1:j1)%int_g3,   0.0_rp, int_g3,   dg3)
  d_max = max(g_tol*abs(dgx), g_tol*abs(dgy), g2_tol*abs(dg2), g3_tol*abs(dg3)) / len_int
  if (d_max < 1) exit
enddo

! Add bend edge if needed
! And rotate to laboratory coords

if (ele%key == sbend$ .or. ele%key == rf_bend$) then
  if (ele%orientation == 1) then
    if (s0 == 0)             call add_bend_edge_to_ints(int_g, int_g2, int_g3, ele, ele%value(e1$), orb0)
    if (s1 == ele%value(l$)) call add_bend_edge_to_ints(int_g, int_g2, int_g3, ele, ele%value(e2$), orb_end)
  else
    if (s0 == 0)             call add_bend_edge_to_ints(int_g, int_g2, int_g3, ele, ele%value(e2$), orb0)
    if (s1 == ele%value(l$)) call add_bend_edge_to_ints(int_g, int_g2, int_g3, ele, ele%value(e1$), orb_end)
  endif

  tilt = ele%value(tilt$) + ele%value(roll$) ! Should be a good approximation.

else
  tilt = ele%value(tilt$)
endif

! Rotate to laboratory coords

if (tilt /= 0) then
  cos_t = cos(tilt)
  sin_t = sin(tilt)
  int_g(1) = (int_g(1) * cos_t - int_g(2) * sin_t)
  int_g(2) = (int_g(1) * sin_t + int_g(2) * cos_t)
endif

bmad_com = bmad_com_save

!---------------------------------------------------------------------------------
contains

subroutine calc_g_at_pt (ele, orb0, z_pos, g_vec, g2, g3, save_orb_mat, orb_save, mat_save)

type (ele_struct) ele, runt
type (coord_struct) orb0, orbz, orb_save  ! Orbit at start, orbit at z.

real(rp) g2, g3, f
real(rp) z_pos, g_vec(3), mat_save(6,6)

integer i, j

logical save_orb_mat, err_flag

! Note: g from g_bending_strength_from_em_field is g of the actual particle and not the ref particle
! so g has a factor of 1/(1 + pz) in it that is corrected for.

call create_element_slice (runt, ele, z_pos, 0.0_rp, ele%branch%param, .false., .false., err_flag, pointer_to_next_ele(ele, -1))
call track1(orb0, runt, ele%branch%param, orbz)
call make_mat6(runt, ele%branch%param, orb0, orbz)
if (save_orb_mat) then
  orb_save = orbz
  mat_save = runt%mat6
endif

call g_bending_strength_from_em_field (ele, ele%branch%param, z_pos, orbz, .true., g_vec, dg)

g2 = sum(g_vec * g_vec)
g3 = sqrt(g2)**3

if (ele%key == sbend$ .or. ele%key == rf_bend$) then
  f = (1 + ele%value(g$) * orb0%vec(1))  ! Variation in path length effect
  g_vec = g_vec * f
  g2 = g2 * f
  g3 = g3 * f
endif

end subroutine calc_g_at_pt

!---------------------------------------------------------------------------------
! contains

subroutine add_bend_edge_to_ints(int_g, int_g2, int_g3, ele, e_edge, orb)

type (ele_struct) ele
type (coord_struct) orb
real(rp) int_g(2), int_g2, int_g3
real(rp) e_edge, dg(2), dg_mag

!

if (e_edge == 0) return

dg = e_edge * ele%value(g$) * [-orb%vec(1), orb%vec(3)]
dg_mag = norm2(dg)
int_g  = int_g + dg
int_g2 = int_g2 + dg_mag**2
int_g3 = int_g3 + dg_mag**3

end subroutine add_bend_edge_to_ints

end subroutine rad_g_integrals

end module
