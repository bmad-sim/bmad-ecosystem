!+
! Subroutine chrom_calc (lat, delta_e, chrom_x, chrom_y, err_flag, &
!                        pz, low_E_lat, high_E_lat, low_E_orb, high_E_orb, ix_branch, orb0)
!
! Subroutine to calculate the chromaticities by computing the tune change when the energy is changed.
! This will handle open geometry lattices. In this case it is assumed that dbeta/de and dalpha/de are
! zero at the beginning of the lattice.
!
! Input:
!   lat           -- lat_struct: Lat
!   delta_e       -- real(rp): +/- Delta energy used for the calculation. Notice that the energy difference
!                      between high and low is 2 * delta_e. If 0 then default of 1.0d-4 is used.
!   pz            -- real(rp), optional: reference momentum about which to calculate. Default is 0. 
!   ix_branch     -- integer, optional: Index of the lattice branch to use. Default is 0.
!   orb0          -- coord_struct, optional: On-energy orbit at start. Needed if lattice branch has an open geometry.
!
! Output:
!   delta_e       -- real(rp): Set to 1.0d-4 if on input DELTA_E =< 0.
!   chrom_x       -- real(rp): Horizontal chromaticity.
!   chrom_y       -- real(rp): Vertical chromaticity.
!   err_flag      -- logical, optional: Set true if there is an error. False otherwise.
!   low_E_lat     -- lat_struct, optional: Lattice with RF off and matrices computed at E_lat +pz - delta_e
!   high_E_lat    -- lat_struct, optional: Lattice with RF off and matrices computed at E_lat +pz + delta_e
!   low_E_orb(:)  -- coord_struct, allocatable, optional: Orbit computed at E_lat + pz - delta_e.
!   high_E_orb(:) -- coord_struct, allocatable, optional: Orbit computed at E_lat + pz + delta_e.
!-

subroutine chrom_calc (lat, delta_e, chrom_x, chrom_y, err_flag, &
                       pz, low_E_lat, high_E_lat, low_E_orb, high_E_orb, ix_branch, orb0)

use bmad_interface, except_dummy => chrom_calc

implicit none

type (lat_struct), target :: lat
type (lat_struct), optional, target :: low_E_lat, high_E_lat
type (lat_struct), target, save :: this_lat
type (lat_struct), pointer :: lat2
type (coord_struct), allocatable, optional, target :: low_E_orb(:), high_E_orb(:)
type (coord_struct), optional :: orb0
type (coord_struct), allocatable, target :: this_orb(:)
type (coord_struct), pointer :: orb_ptr(:)
type (branch_struct), pointer :: branch, branch2
type (ele_struct), pointer :: ele

real(rp) :: high_tune_x, high_tune_y, low_tune_x, low_tune_y
real(rp) :: pz0, delta_e, chrom_x, chrom_y
real(rp), optional :: pz
real time0, time1

integer, optional :: ix_branch
integer nt, stat, ix_br

logical, optional, intent(out) :: err_flag
logical err, used_this_lat

! Init setup

call cpu_time(time0)

ix_br = integer_option(0, ix_branch)
branch => lat%branch(ix_br)
ele => branch%ele(0)

if (present(err_flag)) err_flag = .true.
if (delta_e <= 0) delta_e = 1.0d-4

! reference momentum
pz0 = real_option(0.0_rp, pz)

nt = branch%n_ele_track

! lower energy tune

if (present(low_E_lat)) then
  lat2 => low_E_lat
  used_this_lat = .false.
else
  lat2 => this_lat
  used_this_lat = .true.
endif

lat2 = lat
branch2 => lat2%branch(ix_br)

call set_on_off (rfcavity$, lat2, off$, ix_branch = ix_br)

if (present(low_E_orb)) then
  call reallocate_coord (low_E_orb, branch%n_ele_max)
  orb_ptr => low_E_orb
else
  call reallocate_coord (this_orb, branch%n_ele_max)
  orb_ptr => this_orb
endif

if (branch%param%geometry == closed$) then
  orb_ptr(0)%vec(6) = pz0-delta_e
  if (present(low_E_orb)) then; call closed_orbit_calc (lat2, low_E_orb, 4, 1, ix_br, err)
  else;                         call closed_orbit_calc (lat2, this_orb, 4, 1, ix_br, err)
  endif
  if (err) return
  call lat_make_mat6 (lat2, -1, orb_ptr, ix_br)
  call twiss_at_start (lat2, stat, ix_br, .false.)
  if (stat /= ok$) return

else
  orb_ptr(0) = orb0
  orb_ptr(0)%vec = orb0%vec + (pz0-orb0%vec(6)-delta_e) * [ele%x%eta, ele%x%etap, ele%y%eta, ele%y%etap, 0.0_rp, 1.0_rp]
  if (present(low_E_orb)) then; call track_all(lat2, low_E_orb, ix_br)
  else;                         call track_all(lat2, this_orb, ix_br)
  endif
  call lat_make_mat6 (lat2, -1, orb_ptr, ix_br)
endif

call twiss_propagate_all (lat2, ix_br)
low_tune_x = branch2%ele(nt)%a%phi / twopi
low_tune_y = branch2%ele(nt)%b%phi / twopi

! higher energy tune

if (present(high_E_lat)) then
  lat2 => high_E_lat
  used_this_lat = .false.
else
  lat2 => this_lat
endif

if (.not. used_this_lat) then
  lat2 = lat
  branch2 => lat2%branch(ix_br)
  call set_on_off (rfcavity$, lat2, off$, ix_branch = ix_br)
endif

if (present(high_E_orb)) then
  call reallocate_coord (high_E_orb, branch%n_ele_max)
  orb_ptr => high_E_orb
else
  call reallocate_coord (this_orb, branch%n_ele_max)
  orb_ptr => this_orb
endif

if (branch%param%geometry == closed$) then
  orb_ptr(0)%vec(6) = pz0+delta_e
  if (present(low_E_orb)) then; call closed_orbit_calc (lat2, high_E_orb, 4, 1, ix_br, err)
  else;                         call closed_orbit_calc (lat2, this_orb, 4, 1, ix_br, err)
  endif
  if (err) return
  call lat_make_mat6 (lat2, -1, orb_ptr, ix_br)
  call twiss_at_start (lat2, stat, ix_br, .false.)
  if (stat /= ok$) return

else
  orb_ptr(0) = orb0
  orb_ptr(0)%vec = orb0%vec + (pz0-orb0%vec(6)+delta_e) * [ele%x%eta, ele%x%etap, ele%y%eta, ele%y%etap, 0.0_rp, 1.0_rp]
  if (present(high_E_orb)) then; call track_all(lat2, high_E_orb, ix_br)
  else;                          call track_all(lat2, this_orb, ix_br)
  endif
  call lat_make_mat6 (lat2, -1, orb_ptr, ix_br)
endif

call twiss_propagate_all (lat2, ix_br)
high_tune_x = branch2%ele(nt)%a%phi / twopi
high_tune_y = branch2%ele(nt)%b%phi / twopi

! compute the chrom

chrom_x = (high_tune_x - low_tune_x) / (2 * delta_e)
chrom_y = (high_tune_y - low_tune_y) / (2 * delta_e)

if (present(err_flag)) err_flag = .false.

call cpu_time(time1)
if (bmad_com%debug) print '(a, f12.2)', 'chrom_calc execution time:', time1 - time0

end subroutine
