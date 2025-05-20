!+
! Subroutine chrom_calc (lat, delta_e, chrom_a, chrom_b, err_flag,
!                        pz, low_E_lat, high_E_lat, low_E_orb, high_E_orb, ix_branch, orb0)
!
! Subroutine to calculate the chromaticities by computing the tune change when the energy (actually pz) is changed.
! This will handle open geometry lattices. In this case, dbeta/dpz and dalpha/dpz are are
! taken to be what is set in the beginning element.
!
! Input:
!   lat           -- lat_struct: Lat
!   delta_e       -- real(rp): +/- Delta energy used for the calculation. Notice that the energy difference
!                      between high and low is 2 * delta_e. If 0 then default of 1.0d-4 is used.
!   pz            -- real(rp), optional: reference momentum about which to calculate. Default is 0. 
!   ix_branch     -- integer, optional: Index of the lattice branch to use. Default is 0.
!   orb0          -- coord_struct, optional: On-energy orbit at start. Only needed if lattice branch has an open geometry.
!
! Output:
!   delta_e       -- real(rp): Set to 1.0d-4 if on input DELTA_E =< 0.
!   chrom_a       -- real(rp): a-mode chromaticity.
!   chrom_b       -- real(rp): b-mode chromaticity.
!   err_flag      -- logical, optional: Set true if there is an error. False otherwise.
!   low_E_lat     -- lat_struct, optional: Lattice with RF off and matrices computed at E_lat +pz - delta_e
!   high_E_lat    -- lat_struct, optional: Lattice with RF off and matrices computed at E_lat +pz + delta_e
!   low_E_orb(:)  -- coord_struct, allocatable, optional: Orbit computed at E_lat + pz - delta_e.
!   high_E_orb(:) -- coord_struct, allocatable, optional: Orbit computed at E_lat + pz + delta_e.
!-

subroutine chrom_calc (lat, delta_e, chrom_a, chrom_b, err_flag, &
                       pz, low_E_lat, high_E_lat, low_E_orb, high_E_orb, ix_branch, orb0)

use bmad_interface, except_dummy => chrom_calc

implicit none

type (lat_struct), target :: lat
type (lat_struct), optional, target :: low_E_lat, high_E_lat
type (lat_struct), target :: this_lat
type (lat_struct), pointer :: lat2
type (coord_struct), allocatable, optional, target :: low_E_orb(:), high_E_orb(:)
type (coord_struct), optional :: orb0
type (coord_struct), allocatable, target :: this_orb(:)
type (coord_struct), pointer :: orb_ptr(:)
type (branch_struct), pointer :: branch, branch2
type (ele_struct), pointer :: ele, ele2

real(rp) :: high_tune_x, high_tune_y, low_tune_x, low_tune_y
real(rp) :: pz0, delta_e, chrom_a, chrom_b
real(rp), optional :: pz
real time0, time1

integer, optional :: ix_branch
integer nt, nm, stat, ix_br, ie, i0

logical, optional, intent(out) :: err_flag
logical err, used_this_lat

character(*), parameter :: r_name = 'chrom_calc'

! Init setup

ix_br = integer_option(0, ix_branch)
branch => lat%branch(ix_br)
ele => branch%ele(0)

if (present(err_flag)) err_flag = .true.
if (delta_e <= 0) delta_e = 1.0d-4

! reference momentum
pz0 = real_option(0.0_rp, pz)

nt = branch%n_ele_track
nm = branch%n_ele_max

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
  i0 = 0
  orb_ptr(0)%vec(6) = pz0-delta_e
  if (present(low_E_orb)) then; call closed_orbit_calc (lat2, low_E_orb, 4, 1, ix_br, err)
  else;                         call closed_orbit_calc (lat2, this_orb, 4, 1, ix_br, err)
  endif
  if (err) then
    call out_io (s_warn$, r_name, 'Closed orbit calc failing for low-energy orbit.')
    return
  endif
  call lat_make_mat6 (lat2, -1, orb_ptr, ix_br)
  call twiss_at_start (lat2, stat, ix_br, .false.)
  if (stat /= ok$) then
    call out_io  (s_warn$, r_name, 'Twiss calc failing for low-energy orbit.')
    return
  endif

else
  i0 = 1
  orb_ptr(0) = orb0
  orb_ptr(0)%vec = orb0%vec + (pz0-orb0%vec(6)-delta_e) * [ele%x%eta, ele%x%etap, ele%y%eta, ele%y%etap, 0.0_rp, 1.0_rp]
  if (present(low_E_orb)) then; call track_all(lat2, low_E_orb, ix_br)
  else;                         call track_all(lat2, this_orb, ix_br)
  endif
  ele2 => lat2%branch(ix_br)%ele(0)
  ele2%a%beta  = ele%a%beta  - delta_e * ele%a%dbeta_dpz
  ele2%b%beta  = ele%b%beta  - delta_e * ele%b%dbeta_dpz
  ele2%a%alpha = ele%a%alpha - delta_e * ele%a%dalpha_dpz
  ele2%b%alpha = ele%b%alpha - delta_e * ele%b%dalpha_dpz

  ele2%a%eta   = ele%a%eta   - delta_e * ele%a%deta_dpz
  ele2%b%eta   = ele%b%eta   - delta_e * ele%b%deta_dpz
  ele2%x%eta   = ele%x%eta   - delta_e * ele%x%deta_dpz
  ele2%y%eta   = ele%y%eta   - delta_e * ele%y%deta_dpz
  ele2%z%eta   = ele%z%eta   - delta_e * ele%z%deta_dpz

  ele2%a%etap  = ele%a%etap  - delta_e * ele%a%detap_dpz
  ele2%b%etap  = ele%b%etap  - delta_e * ele%b%detap_dpz
  ele2%x%etap  = ele%x%etap  - delta_e * ele%x%detap_dpz
  ele2%y%etap  = ele%y%etap  - delta_e * ele%y%detap_dpz
  ele2%z%etap  = ele%z%etap  - delta_e * ele%z%detap_dpz
  call lat_make_mat6 (lat2, -1, orb_ptr, ix_br)
endif

call twiss_propagate_all (lat2, ix_br)
low_tune_x = branch2%ele(nt)%a%phi / twopi
low_tune_y = branch2%ele(nt)%b%phi / twopi

branch%ele(i0:nm)%a%dbeta_dpz  = branch2%ele(i0:nm)%a%beta
branch%ele(i0:nm)%b%dbeta_dpz  = branch2%ele(i0:nm)%b%beta
branch%ele(i0:nm)%a%dalpha_dpz = branch2%ele(i0:nm)%a%alpha
branch%ele(i0:nm)%b%dalpha_dpz = branch2%ele(i0:nm)%b%alpha

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
  if (err) then
    call out_io (s_warn$, r_name, 'Closed orbit calc failing for high-energy orbit.')
    return
  endif
  call lat_make_mat6 (lat2, -1, orb_ptr, ix_br)
  call twiss_at_start (lat2, stat, ix_br, .false.)
  if (stat /= ok$) then
    call out_io  (s_warn$, r_name, 'Twiss calc failing for high-energy orbit.')
    return
  endif

else
  orb_ptr(0) = orb0
  orb_ptr(0)%vec = orb0%vec + (pz0-orb0%vec(6)+delta_e) * [ele%x%eta, ele%x%etap, ele%y%eta, ele%y%etap, 0.0_rp, 1.0_rp]
  if (present(high_E_orb)) then; call track_all(lat2, high_E_orb, ix_br)
  else;                          call track_all(lat2, this_orb, ix_br)
  endif
  ele2 => lat2%branch(ix_br)%ele(0)
  ele2%a%beta  = ele%a%beta  + delta_e * ele%a%dbeta_dpz
  ele2%b%beta  = ele%b%beta  + delta_e * ele%b%dbeta_dpz
  ele2%a%alpha = ele%a%alpha + delta_e * ele%a%dalpha_dpz
  ele2%b%alpha = ele%b%alpha + delta_e * ele%b%dalpha_dpz

  ele2%a%eta   = ele%a%eta   + delta_e * ele%a%deta_dpz
  ele2%b%eta   = ele%b%eta   + delta_e * ele%b%deta_dpz
  ele2%x%eta   = ele%x%eta   + delta_e * ele%x%deta_dpz
  ele2%y%eta   = ele%y%eta   + delta_e * ele%y%deta_dpz
  ele2%z%eta   = ele%z%eta   + delta_e * ele%z%deta_dpz

  ele2%a%etap  = ele%a%etap  + delta_e * ele%a%detap_dpz
  ele2%b%etap  = ele%b%etap  + delta_e * ele%b%detap_dpz
  ele2%x%etap  = ele%x%etap  + delta_e * ele%x%detap_dpz
  ele2%y%etap  = ele%y%etap  + delta_e * ele%y%detap_dpz
  ele2%z%etap  = ele%z%etap  + delta_e * ele%z%detap_dpz
  call lat_make_mat6 (lat2, -1, orb_ptr, ix_br)
endif

call twiss_propagate_all (lat2, ix_br)
high_tune_x = branch2%ele(nt)%a%phi / twopi
high_tune_y = branch2%ele(nt)%b%phi / twopi

! compute the chrom

chrom_a = (high_tune_x - low_tune_x) / (2 * delta_e)
chrom_b = (high_tune_y - low_tune_y) / (2 * delta_e)

branch%ele(i0:nm)%a%dbeta_dpz  = (branch2%ele(i0:nm)%a%beta  - branch%ele(i0:nm)%a%dbeta_dpz) / (2 * delta_e)
branch%ele(i0:nm)%b%dbeta_dpz  = (branch2%ele(i0:nm)%b%beta  - branch%ele(i0:nm)%b%dbeta_dpz) / (2 * delta_e)
branch%ele(i0:nm)%a%dalpha_dpz = (branch2%ele(i0:nm)%a%alpha - branch%ele(i0:nm)%a%dalpha_dpz) / (2 * delta_e)
branch%ele(i0:nm)%b%dalpha_dpz = (branch2%ele(i0:nm)%b%alpha - branch%ele(i0:nm)%b%dalpha_dpz) / (2 * delta_e)

branch%ele(i0:nm)%a%deta_dpz  = (branch2%ele(i0:nm)%a%eta  - branch%ele(i0:nm)%a%eta) / (2 * delta_e)
branch%ele(i0:nm)%b%deta_dpz  = (branch2%ele(i0:nm)%b%eta  - branch%ele(i0:nm)%b%eta) / (2 * delta_e)
branch%ele(i0:nm)%x%deta_dpz  = (branch2%ele(i0:nm)%x%eta  - branch%ele(i0:nm)%x%eta) / (2 * delta_e)
branch%ele(i0:nm)%y%deta_dpz  = (branch2%ele(i0:nm)%y%eta  - branch%ele(i0:nm)%y%eta) / (2 * delta_e)
branch%ele(i0:nm)%z%deta_dpz  = (branch2%ele(i0:nm)%z%eta  - branch%ele(i0:nm)%z%eta) / (2 * delta_e)

branch%ele(i0:nm)%a%detap_dpz = (branch2%ele(i0:nm)%a%etap - branch%ele(i0:nm)%a%etap) / (2 * delta_e)
branch%ele(i0:nm)%b%detap_dpz = (branch2%ele(i0:nm)%b%etap - branch%ele(i0:nm)%b%etap) / (2 * delta_e)
branch%ele(i0:nm)%x%detap_dpz = (branch2%ele(i0:nm)%x%etap - branch%ele(i0:nm)%x%etap) / (2 * delta_e)
branch%ele(i0:nm)%y%detap_dpz = (branch2%ele(i0:nm)%y%etap - branch%ele(i0:nm)%y%etap) / (2 * delta_e)
branch%ele(i0:nm)%z%detap_dpz = (branch2%ele(i0:nm)%z%etap - branch%ele(i0:nm)%z%etap) / (2 * delta_e)

if (present(err_flag)) err_flag = .false.

end subroutine
