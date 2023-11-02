!+
! Subroutine radiation_integrals (lat, orbit, mode, ix_cache, ix_branch, rad_int_by_ele)
!
! Subroutine to calculate the synchrotron radiation integrals along with the
! emittance, and energy spread.
!
! Note: A negative emittance is possible and just means that the beam is
! unstable. That is, you have a negative damping partition number.
!
! The calculation can spend a significant amount of time calculating the integrals.
! To speed up this calculation, the ix_cache argument can be used to stash values for 
! the bending radius, etc. so that repeated calls to radiation_integrals consumes less time.
! To use caching: 
!   1) First call radiation_integrals with ix_cache set to 0. 
!      radiation_integrals will cache values needed to compute the integrals 
!      for the wigglers and assign ix_cache a unique number that is used to 
!      point to the cache. 
!   2) Subsequent calls to radiation_integrals should just pass the value of 
!      ix_cache that has been assigned.
! A new cache, with a unique value for ix_cache, is created each time 
! radiation_integrals is called with ix_cache = 0. This is useful if there are multiple lattices. 
! To release the memory associated with a cache call release_rad_int_cache(ix_cache).
!
! NOTE: The validity of the cache is dependent upon the orbit and parameters like the
! quadrupole strengths not varying too much but is independent of variations in
! the Twiss parameters. 
!
! Recommendation: Do not use caching of non-wiggler elements unless you need it.
!
! Input:
!   lat        -- lat_struct: Lattice to use. The calculation assumes that 
!                    the Twiss parameters have been calculated.
!   orbit(0:)  -- Coord_struct: Closed orbit for the branch.
!   ix_cache   -- Integer, optional: Cache pointer.
!                      = -2 --> No temporary wiggler cache. This is slow so only use as a check.
!                      = -1 --> Use temporary cache for wiggler elements only (default).
!                      =  0 --> Create a new cache for all elements.
!                      >  0 --> Use the corresponding cache. 
!   ix_branch  -- Integer, optional: Lattice branch index. Default is 0.
!
! Output:
!   mode     -- normal_modes_struct: Parameters for the ("horizontal like") a-mode,
!                              ("vertical like") b-mode, and the z-mode
!     %synch_int(0:3) -- Synchrotron integrals. See Bmad manual
!     %sigE_E         -- Sigma_E/E energy spread
!     %sig_z          -- Bunch Length
!     %e_loss         -- Energy loss in eV per turn
!     %a, %b, %z      -- Anormal_mode_struct: Substructure
!       %emittance      -- Emittance. B-mode emit includes photon opening angle (I6) contribution.
!       %synch_int(4:6) -- Synchrotron integrals
!       %j_damp         -- Damping partition factor
!       %alpha_damp     -- Exponential damping coefficient per turn
!     %lin            -- Linac version of the integrals.
!       %i2_E4           -- Integral: g^2 * gamma^4
!       %i3_E7           -- Integral: g^3 * gamma^7
!       %i5a_E6          -- Integral: (g^3 * H_a) * gamma^6
!       %i5b_E6          -- Integral: (g^3 * H_b) * gamma^6
!       %sig_E1          -- Energy spread after 1 pass (eV)
!       %a_emittance_end -- a-mode emittance at end of linac
!       %b_emittance_end -- b-mode emittance at end of linac
!   ix_cache -- Integer, optional: Cache pointer. If ix_cache = 0 at input then
!                   ix_cache is set to a unique number. Otherwise ix_cache 
!                   is not changed.
!   rad_int_by_ele  -- Rad_int_all_ele_struct, optional: Radiation integrals element by element. 
!     %branch(ix_branch)%ele(0:) -- Array of rad_int1_struct structures, one for each element in the branch.
!
! Notes:
!
! 1) %synch_int(1) = momentum_compaction * lat_length
!
! 2) The lin_norm_emit values are running sums from the beginning of the 
!    lattice and include the beginning emittance stored in lat%a%emit and lat%b%emit.
!-

subroutine radiation_integrals (lat, orbit, mode, ix_cache, ix_branch, rad_int_by_ele)

use radiation_mod, except_dummy1 => radiation_integrals
use transfer_map_mod, except_dummy3 => radiation_integrals

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele, slave, field_ele
type (ele_struct) :: ele2, ele_start, ele_end, runt
type (coord_struct), target :: orbit(0:), orb_start, orb_end, orb_here
type (normal_modes_struct) mode
type (bmad_common_struct) bmad_com_save
type (track_struct), target :: track
type (track_point_struct), pointer :: tp
type (rad_int_all_ele_struct), optional, target :: rad_int_by_ele
type (rad_int_branch_struct), target :: rad_int_branch
type (rad_int_branch_struct), pointer :: ri_branch
type (rad_int_info_struct) :: ri_info
type (rad_int_cache_struct), pointer :: cache
type (rad_int_cache1_struct), pointer :: cache_ele ! pointer to cache in use
type (rad_int_track_point_struct), pointer :: c_pt, cp0, cp
type (rad_int_track_point_struct) pt
type (rad_int_track_point_struct), allocatable :: pt_temp(:)
type (rad_int1_struct) int_tot
type (rad_int1_struct), pointer :: rad_int1
type (rad_int_cache_struct), allocatable :: cache_temp(:)
type (spin_orbit_map1_struct) map1

real(rp) i1, i2, i3, i4a, i4b, i4z, i5a, i5b, i6b, m65, G_max, gx_err, gy_err, del_z
real(rp) theta, energy, gamma2_factor, energy_loss, arg, ll, gamma_f, ds_step
real(rp) a_pole(0:n_pole_maxx), b_pole(0:n_pole_maxx), dk(2,2), kx, ky, beta_min
real(rp) v(4,4), v_inv(4,4), z_here, z_start, mc2, gamma, gamma4, gamma6
real(rp) kz, fac, c, s, factor, g2, g_x0, dz_small, z1, const_q, vec0(6), mat6(6,6)
real(rp), parameter :: const_q_factor = 55 * h_bar_planck * c_light / (32 * sqrt_3) ! Cf: Sands Eq 5.46 pg 124.
real time0, time1

integer, optional :: ix_cache, ix_branch
integer i, j, k, n, ix, ixe, ib, ip, ir, key2, n_step, ix_pole_max

character(*), parameter :: r_name = 'radiation_integrals'

logical do_alloc, use_cache, init_cache, cache_only_wig, err, hybrid_warning_given
logical, parameter :: t = .true., f = .false.

!---------------------------------------------------------------------
! Allocate rad_int_by_ele

call cpu_time(time0)

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

! Photon branches do not have associated radiation integrals.

if (present(ix_branch)) then
  branch => lat%branch(ix_branch)
else
  branch => lat%branch(0)
endif

mode = normal_modes_struct()

if (branch%param%particle == photon$) return

if (any(orbit(1:branch%n_ele_track)%species == not_set$)) then
  call out_io (s_error$, r_name, 'ORBIT IN UNINITALIZED STATE! RADIATION INTEGRALS NOT COMPUTED.')
  return
endif

if (any(branch%ele(0:branch%n_ele_track)%a%beta == 0)) then
  n = minloc(branch%ele(0:branch%n_ele_track)%a%beta, 1) - 1
  call out_io (s_error$, r_name, 'TWISS PARAMETERS HAVE NOT BEEN INITIALIZED AT ELEMENT \i0\ OF BRANCH \i0\.', &
                                 '[DUE TO UNSTABLE RING? OR MISSING CALL TO TWISS_PROPAGATE_ALL?]', &
                                 'NO RADIATION INTEGRALS WILL BE COMPUTED.', i_array = [n, ib])
  return
endif

! Init
! To make the calculation go faster turn off radiation fluctuations and damping

allocate (rad_int_branch%ele(0:branch%n_ele_max))

ri_info%branch => branch
ri_info%orbit => orbit

m65 = 0
mode%rf_voltage = 0
int_tot = rad_int1_struct()

bmad_com_save = bmad_com
bmad_com%convert_to_kinetic_momentum = .true.
bmad_com%radiation_fluctuations_on = .false.
bmad_com%radiation_damping_on = .false.
bmad_com%csr_and_space_charge_on = .false.
bmad_com%spin_tracking_on = .false.
bmad_com%rel_tol_adaptive_tracking = 1d-6
bmad_com%abs_tol_adaptive_tracking = 1d-8
bmad_com%high_energy_space_charge_on = .false.

call init_ele (ele2)
call init_ele (ele_start)
call init_ele (ele_end)

!---------------------------------------------------------------------
! Caching

! find a cache

use_cache = .true.
init_cache = .true.
cache_only_wig = .false.

if (.not. allocated(rad_int_cache_common)) allocate (rad_int_cache_common(0:10))

if (.not. present(ix_cache)) then
  cache => rad_int_cache_common(0)
  cache_only_wig = .true.

elseif (ix_cache == -1) then
  cache => rad_int_cache_common(0)
  cache_only_wig = .true.

elseif (ix_cache == 0) then
  n = ubound(rad_int_cache_common, 1)
  if (count(rad_int_cache_common(1:n)%in_use) == n) then
    call move_alloc(rad_int_cache_common, cache_temp)
    allocate (rad_int_cache_common(0:2*n))
    rad_int_cache_common(0:n) = cache_temp
  endif

  do i = 1, ubound(rad_int_cache_common, 1)
    if (rad_int_cache_common(i)%in_use) cycle
    rad_int_cache_common(i)%in_use = .true.
    ix_cache = i
    cache => rad_int_cache_common(i)
    exit
    if (i == size(rad_int_cache_common)) then
      call out_io (s_fatal$, r_name, 'CACHE OUT OF MEMORY!')
      if (global_com%exit_on_error) call err_exit
    endif
  enddo

elseif (ix_cache > 0) then
  if (ix_cache > ubound(rad_int_cache_common, 1)) then
    call out_io (s_fatal$, r_name, 'INVALID IX_CACHE ARGUMENT \i0\ ', ix_cache)
    if (global_com%exit_on_error) call err_exit
  endif
  cache => rad_int_cache_common(ix_cache)
  init_cache = .false.

else  ! ix_cache < -1
  use_cache = .false.
  init_cache = .false.
endif

! allocate cache

if (init_cache) then

  if (allocated(cache%c_ele)) then
    if (size(cache%c_ele) < branch%n_ele_max) deallocate (cache%c_ele)
  endif
  if (.not. allocated(cache%c_ele)) allocate (cache%c_ele(branch%n_ele_max))  ! allocate cache memory

  cache%c_ele%cache_type = no_cache$

  branch%ele%bookkeeping_state%rad_int = stale$
endif

! Now cache the information.
! The cache ri_info is computed with all offsets off.
! Only need to update the cache if ele%bookkeeping_state%rad_int has been altered.

if (use_cache .or. init_cache) then
  do ixe = 1, branch%n_ele_track

    ele => branch%ele(ixe)
    if (ele%value(l$) == 0) cycle
    field_ele => pointer_to_field_ele(ele, 1)

    select case (ele%key)
    case (wiggler$, undulator$, em_field$)
    case (quadrupole$, sol_quad$, sbend$, rf_bend$, sad_mult$, hkicker$, vkicker$)
      if (cache_only_wig) cycle
    case default
      if (cache_only_wig) cycle
      if (attribute_index(ele, 'HKICK') == 0) cycle   ! Has no kick attributes.
      if (ele%value(hkick$) == 0 .or. ele%value(vkick$) == 0) cycle
    end select

    cache_ele => cache%c_ele(ixe)
    cache_ele%cache_type = cache_no_misalign$

    if (branch%ele(ixe)%bookkeeping_state%rad_int /= stale$) cycle
    branch%ele(ixe)%bookkeeping_state%rad_int = ok$

    ele2 = branch%ele(ixe)
    key2 = ele2%key
    if (key2 == undulator$ .or. key2 == em_field$) key2 = wiggler$

    orb_start = orbit(ixe-1)
    call set_tracking_method_for_element_integration(ele2)

    !------------------------------------
    ! For tracking methods that go step-by-step, simply track through the element to get the cache info.
    ! This is faster than evaluating at a set of given cache points.

    if (ele2%tracking_method == runge_kutta$ .or. ele2%tracking_method == time_runge_kutta$ .or. ele2%tracking_method == symp_lie_bmad$) then
      track%n_pt = -1
      call track1 (orb_start, ele2, branch%param, orb_end, track, make_map1 = .true.)
      mat6 = ele2%mat6
      call allocate_cache(cache_ele, track%n_pt)
      do i = 0, track%n_pt
        cache_ele%pt(i)%ref_orb_in  = orb_start
        call cache_fill(ele2, branch, cache_ele%pt(i), track, i, max(0, i-1), min(track%n_pt, i+1))
      enddo

    ! All else...
    else  
      beta_min = min(ele2%a%beta, ele2%b%beta, branch%ele(ixe-1)%a%beta, branch%ele(ixe-1)%b%beta)
      n_step = nint(ele%value(l$) / min(ele2%value(ds_step$), beta_min / 10, abs(ele2%value(l$))/2))
      n_step = min(n_step, nint(1000*ele2%value(l$)))   ! So del_z is at least 1 mm 
      n_step = max(n_step, 3)
      del_z = ele2%value(l$) / n_step
      call allocate_cache(cache_ele, n_step)

      z_start = 0
      ele_start = branch%ele(ixe-1)
      dz_small = min (1e-4_rp, del_z/3)
      call mat_make_unit (mat6)
      vec0 = 0

      orb_here = orb_start
      do k = 0, n_step
        z_here = k * del_z
        z1 = z_here + dz_small
        if (z1 > ele2%value(l$)) z1 = max(0.0_rp, z_here - dz_small)
        call cache_this_point (branch, z_start, z1, z_here, orb_start, orb_here, ele2, mat6, vec0, cache_ele%pt(k), orb_end, ele_end)
        z_start = z_here
        orb_here = orb_end
        ele_start = ele_end
      enddo

      ! Make sure we have enough cache points. If del_z < 10^-3 we have enough points or 
      ! if the difference of the integrated gx and gy averaged over the element using a cubic fit to each interval
      ! versus using a linear fit is within 10^-4.

      outer_loop: do
        if (del_z < 1d-3) exit

        gx_err = 0
        gy_err = 0
        do k = 1, n_step-2  ! Ignore end intervals in error estimate.
          ! gx_err and gy_err are the difference between a linear fit and a cubic fit modulo a factor of 12
          ! which is put in at the end.
          gx_err = gx_err + abs(-cache_ele%pt(k-1)%g_x0 + cache_ele%pt(k)%g_x0 + &
                                          cache_ele%pt(k+1)%g_x0 - cache_ele%pt(k+2)%g_x0)
          gy_err = gy_err + abs(-cache_ele%pt(k-1)%g_y0 + cache_ele%pt(k)%g_y0 + &
                                          cache_ele%pt(k+1)%g_y0 - cache_ele%pt(k+2)%g_y0)
        enddo
        if (gx_err < 12 * 1d-4 * (n_step - 2) .and. gx_err < 12 * 1d-4 * (n_step - 2)) exit

        ! Need more points...
        if (ubound(cache_ele%pt, 1) < 2*n_step) then
          call move_alloc(cache_ele%pt, pt_temp)
          allocate (cache_ele%pt(0:2*n_step))
          cache_ele%pt(0:2*n_step:2) = pt_temp
          deallocate(pt_temp)
        else
          cache_ele%pt(0:2*n_step:2) = cache_ele%pt(0:n_step)
        endif

        n_step = 2*n_step
        cache_ele%n_pt = n_step
        del_z = del_z / 2
        dz_small = min (1e-4_rp, del_z/3)

        orb_here = orb_start
        do k = 1, n_step, 2
          z_here = k * del_z
          z1 = z_here + dz_small
          if (z1 > ele2%value(l$)) z1 = max(0.0_rp, z_here - dz_small)
          cp0 => cache_ele%pt(k-1)
          cp => cache_ele%pt(k)
          z_start = cp0%s_body
          orb_here = cp0%ref_orb_out
          mat6 = cp0%mat6
          vec0 = cp0%vec0
          ele_start%map_ref_orb_in   = orb_start
          ele_start%map_ref_orb_out  = cp0%ref_orb_out
          ele_start%time_ref_orb_in  = orb_start
          ele_start%time_ref_orb_out = cp0%ref_orb_out        
          call cache_this_point (branch, z_start, z1, z_here, orb_start, orb_here, ele2, mat6, vec0, cache_ele%pt(k), orb_end, ele_end)
        enddo
      enddo outer_loop
    endif

  enddo
  
endif ! (init_cache)

!---------------------------------------------------------------------
! Loop over all elements...
! We do the elements that can be integrated quickly to establish a baseline 
! for setting the error tolerance for the elements that take more time to integrate through.

hybrid_warning_given = .false.

do ir = 1, branch%n_ele_track

  ele => branch%ele(ir)
  if (.not. ele%is_on) cycle

  if (ele%key == hybrid$) then
    if (.not. hybrid_warning_given) call out_io (s_error$, r_name, 'CANNOT COMPUTE RADIATION INTEGRALS WHEN THERE IS A HYBRID ELEMENT: ' // ele%name)
    hybrid_warning_given = .true.
    cycle
  endif

  ri_info%ele => ele
  rad_int1 => rad_int_branch%ele(ir)

  nullify (ri_info%cache_ele)
  if (use_cache) then
    if (cache%c_ele(ir)%cache_type /= no_cache$) ri_info%cache_ele => cache%c_ele(ir)
  endif

  pt = rad_int_track_point_struct()  ! zero components

  if (ele%key == rfcavity$) then
    m65 = m65 + ele%mat6(6,5)
    mode%rf_voltage = mode%rf_voltage + ele%value(voltage$)
  endif

  ll = ele%value(l$)
  if (ll == 0) cycle

  ! custom

  if (ele%key == custom$) then
    if (.not. associated(radiation_integrals_custom_ptr)) then
      call out_io (s_error$, r_name, 'RADIATION_INTEGRALS_CUSTOM_PTR HAS NOT BEEN SET IN THIS PROGRAM!', &
                                     'NEEDED FOR CUSTOM ELEMENT: ' // ele%name)
      cycle
    endif
    call radiation_integrals_custom_ptr (lat, ir, orbit, rad_int1, err)
    cycle
  endif

  ! Wigglers and undulators get handled later.

  if (ele%key == wiggler$ .or. ele%key == undulator$ .or. ele%key == em_field$) cycle
  if (ele%key == patch$) cycle
  if (ele%value(hkick$) == 0 .and. ele%value(vkick$) == 0 .and. &
          ele%key /= quadrupole$ .and. ele%key /= sol_quad$ .and. ele%key /= sbend$ .and. &
          ele%key /= rf_bend$ .and. ele%key /= hkicker$ .and. ele%key /= vkicker$) cycle

  ! All other elements

  if (ele%key == sbend$) then
    call transfer_ele(ele, runt)
    theta = ele%value(ref_tilt_tot$)
    pt%g_x0 = cos(theta) * ele%value(g$)
    pt%g_y0 = sin(theta) * ele%value(g$)
    g2 = ele%value(g$)**2
    ! Edge effects for a bend. In this case we ignore any rolls.
    call propagate_part_way (orbit(ir-1), branch%param, pt, ri_info, 0.0_rp, runt)
    rad_int1%i4a = -g2 * tan(ele%value(e1$)) * (ri_info%eta_a(1) * cos(theta) + ri_info%eta_a(3) * sin(theta))
    rad_int1%i4b = -g2 * tan(ele%value(e1$)) * (ri_info%eta_b(1) * cos(theta) + ri_info%eta_b(3) * sin(theta))
    call propagate_part_way (orbit(ir-1), branch%param, pt, ri_info, ll, runt)
    rad_int1%i4a = rad_int1%i4a - g2 * tan(ele%value(e2$)) * (ri_info%eta_a(1) * cos(theta) + ri_info%eta_a(3) * sin(theta))
    rad_int1%i4b = rad_int1%i4b - g2 * tan(ele%value(e2$)) * (ri_info%eta_b(1) * cos(theta) + ri_info%eta_b(3) * sin(theta))
  elseif (ele%key == rf_bend$) then
    theta = ele%value(ref_tilt_tot$)
    pt%g_x0 = cos(theta) * ele%value(g$)
    pt%g_y0 = sin(theta) * ele%value(g$)
    g2 = ele%value(g$)**2
  endif

  ! Integrate for quads, bends and nonzero kicks

  call qromb_rad_int (branch%param, [T, T, T, T, T, T, T, T, T], pt, ri_info, int_tot, rad_int1)

enddo

!----------------------------------------------------------
! Integrate wigglers, undulators.

do ir = 1, branch%n_ele_track

  ele => branch%ele(ir)
  if (.not. ele%is_on) cycle
  field_ele => pointer_to_field_ele(ele, 1)

  ri_info%ele => ele
  rad_int1 => rad_int_branch%ele(ir)

  nullify (ri_info%cache_ele)
  if (use_cache) then
    if (cache%c_ele(ir)%cache_type /= no_cache$) ri_info%cache_ele => cache%c_ele(ir)
  endif

  pt = rad_int_track_point_struct()  ! zero components

  select case (ele%key)
  case (wiggler$, undulator$, em_field$)
    ! For a planar or helical model wiggler we make the approximation that the variation of G is
    ! fast compaired to the variation in eta.
    if (field_ele%field_calc /= fieldmap$ .and. ele%tracking_method == bmad_standard$) then
      if (ele%value(l_period$) == 0) cycle        ! Cannot do calculation
      G_max = ele%value(g_max$)       ! 1/rho at max B
      if (field_ele%field_calc == planar_model$) then
        rad_int1%i0 = (ele%value(e_tot$) / mass_of(branch%param%particle)) * 2 * G_max / 3
        rad_int1%i1 = (G_max * ele%value(l_period$) / twopi)**2 / 2
        rad_int1%i2 = ele%value(l$) * G_max**2 / 2
        rad_int1%i3 = ele%value(l$) * 4 * G_max**3 / (3 * pi)
      else
        rad_int1%i0 = (ele%value(e_tot$) / mass_of(branch%param%particle)) * G_max
        rad_int1%i1 = (G_max * ele%value(l_period$) / twopi)**2
        rad_int1%i2 = ele%value(l$) * G_max**2
        rad_int1%i3 = ele%value(l$) * G_max**3
      endif

      call qromb_rad_int (branch%param, [F, F, F, T, T, T, T, F, T], pt, ri_info, int_tot, rad_int1)
      cycle
    endif

  case default
    cycle
  end select

  !

  ll = ele%value(l$)
  if (ll == 0) cycle

  call qromb_rad_int (branch%param, [T, T, T, T, T, T, T, T, T], pt, ri_info, int_tot, rad_int1)

enddo

!---------------------------------------------------------------------
! Now put everything together...
! Linac radiation integrals:

lat%lord_state%rad_int = ok$

mc2 = mass_of (branch%param%particle)
gamma_f = branch%ele(branch%n_ele_track)%value(e_tot$) / mc2
const_q = const_q_factor / mc2

mode%lin%sig_E1 = 0
mode%lin%i2_E4  = 0
mode%lin%i3_E7  = 0
mode%lin%i5a_E6 = 0
mode%lin%i5b_E6 = 0

factor = 2 * const_q * classical_radius_factor / (3 * mc2)

rad_int_branch%ele%i4z = rad_int_branch%ele%i4a + rad_int_branch%ele%i4b

do i = 0, branch%n_ele_track
  gamma = branch%ele(i)%value(e_tot$) / mc2
  gamma4 = gamma**4
  gamma6 = gamma4 * gamma**2
  rad_int_branch%ele(i)%lin_i2_E4  = rad_int_branch%ele(i)%i2 * gamma4
  rad_int_branch%ele(i)%lin_i3_E7  = rad_int_branch%ele(i)%i3 * gamma6 * gamma
  rad_int_branch%ele(i)%lin_i5a_E6 = rad_int_branch%ele(i)%i5a * gamma6
  rad_int_branch%ele(i)%lin_i5b_E6 = rad_int_branch%ele(i)%i5b * gamma6
  mode%lin%i2_E4  = mode%lin%i2_E4  + rad_int_branch%ele(i)%lin_i2_E4
  mode%lin%i3_E7  = mode%lin%i3_E7  + rad_int_branch%ele(i)%lin_i3_E7
  mode%lin%i5a_E6 = mode%lin%i5a_E6 + rad_int_branch%ele(i)%lin_i5a_E6
  mode%lin%i5b_E6 = mode%lin%i5b_E6 + rad_int_branch%ele(i)%lin_i5b_E6
  rad_int_branch%ele(i)%lin_norm_emit_a = branch%a%emit * gamma + factor * mode%lin%i5a_E6
  rad_int_branch%ele(i)%lin_norm_emit_b = branch%b%emit * gamma + factor * mode%lin%i5b_E6
  rad_int_branch%ele(i)%lin_sig_E = mc2 * sqrt(2 * factor * mode%lin%i3_E7)
enddo

mode%lin%sig_E1 = mc2 * sqrt (2 * factor * mode%lin%i3_E7)
mode%lin%a_emittance_end = factor * mode%lin%i5a_e6 / gamma_f
mode%lin%b_emittance_end = factor * mode%lin%i5b_e6 / gamma_f

! Normal integrals

i1   = int_tot%i1
i2   = int_tot%i2
i3   = int_tot%i3
i4a  = int_tot%i4a
i4b  = int_tot%i4b
i5a  = int_tot%i5a
i5b  = int_tot%i5b
i6b  = int_tot%i6b

i4z = i4a + i4b

energy = branch%ele(0)%value(e_tot$)
gamma2_factor = (energy / mc2)**2
energy_loss = 2 * classical_radius_factor * mode%lin%i2_e4 / 3

mode%synch_int(0) = int_tot%i0
mode%synch_int(1) = i1
mode%synch_int(2) = i2
mode%synch_int(3) = i3

mode%a%synch_int(4) = i4a
mode%b%synch_int(4) = i4b
mode%z%synch_int(4) = i4z

mode%a%synch_int(5) = i5a
mode%b%synch_int(5) = i5b

mode%a%synch_int(6) = 0
mode%b%synch_int(6) = i6b

if (branch%param%geometry == closed$) then
  if (i2 /= 0) then

    mode%a%emittance         = const_q * gamma2_factor * i5a / (i2 - i4a)
    mode%a%emittance_no_vert = mode%a%emittance
    mode%b%emittance         = const_q * (gamma2_factor * i5b + 13 * i6b / 55) / (i2 - i4b)
    mode%b%emittance_no_vert = const_q * (gamma2_factor * i5b) / (i2 - i4b)

    mode%a%j_damp = 1 - i4a / i2
    mode%b%j_damp = 1 - i4b / i2
    mode%z%j_damp = 2 + i4z / i2

    arg = (const_q * i3 * gamma2_factor / (2*i2 + i4z))
    if (arg > 0) then
      mode%sigE_E = sqrt(arg)
    else
      mode%sigE_E = 1d30  ! Something large
    endif

  endif

  mode%a%alpha_damp = energy_loss * mode%a%j_damp / (2 * energy)
  mode%b%alpha_damp = energy_loss * mode%b%j_damp / (2 * energy)
  mode%z%alpha_damp = energy_loss * mode%z%j_damp / (2 * energy)

  mode%e_loss = energy_loss

  if (m65*i1 > 0) then
    mode%sig_z = sqrt(i1/m65) * mode%sigE_E
  else   ! Unstable
    mode%sig_z = 1d30  ! Something large
  endif

  mode%z%emittance = mode%sig_z * mode%sigE_E
  mode%z%emittance_no_vert = mode%z%emittance
endif

bmad_com = bmad_com_save

call deallocate_ele_pointers(ele2)
call deallocate_ele_pointers(ele_start)
call deallocate_ele_pointers(ele_end)

! Fill in rad_int_by_ele

if (present(rad_int_by_ele)) then
  ri_branch => rad_int_by_ele%branch(branch%ix_branch)
  call move_alloc (rad_int_branch%ele, ri_branch%ele)

  do i = branch%n_ele_track+1, branch%n_ele_max
    ele => branch%ele(i)
    if (ele%lord_status /= super_lord$) cycle
    do j = 1, ele%n_slave
      slave => pointer_to_slave(ele, j)
      ib = slave%ix_branch
      k = slave%ix_ele
      ri_branch%ele(i)%i0 = ri_branch%ele(i)%i0 + rad_int_by_ele%branch(ib)%ele(k)%i0
      ri_branch%ele(i)%i1 = ri_branch%ele(i)%i1 + rad_int_by_ele%branch(ib)%ele(k)%i1
      ri_branch%ele(i)%i2 = ri_branch%ele(i)%i2 + rad_int_by_ele%branch(ib)%ele(k)%i2
      ri_branch%ele(i)%i3 = ri_branch%ele(i)%i3 + rad_int_by_ele%branch(ib)%ele(k)%i3
      ri_branch%ele(i)%i4a = ri_branch%ele(i)%i4a + rad_int_by_ele%branch(ib)%ele(k)%i4a
      ri_branch%ele(i)%i4b = ri_branch%ele(i)%i4b + rad_int_by_ele%branch(ib)%ele(k)%i4b
      ri_branch%ele(i)%i4z = ri_branch%ele(i)%i4z + rad_int_by_ele%branch(ib)%ele(k)%i4z
      ri_branch%ele(i)%i5a = ri_branch%ele(i)%i5a + rad_int_by_ele%branch(ib)%ele(k)%i5a
      ri_branch%ele(i)%i5b = ri_branch%ele(i)%i5b + rad_int_by_ele%branch(ib)%ele(k)%i5b
      ri_branch%ele(i)%i6b = ri_branch%ele(i)%i6b + rad_int_by_ele%branch(ib)%ele(k)%i6b
      ri_branch%ele(i)%n_steps = ri_branch%ele(i)%n_steps + rad_int_by_ele%branch(ib)%ele(k)%n_steps
      ri_branch%ele(i)%lin_i2_E4 = ri_branch%ele(i)%lin_i2_E4 + rad_int_by_ele%branch(ib)%ele(k)%lin_i2_E4
      ri_branch%ele(i)%lin_i3_E7 = ri_branch%ele(i)%lin_i3_E7 + rad_int_by_ele%branch(ib)%ele(k)%lin_i3_E7
      ri_branch%ele(i)%lin_i5a_E6 = ri_branch%ele(i)%lin_i5a_E6 + rad_int_by_ele%branch(ib)%ele(k)%lin_i5a_E6
      ri_branch%ele(i)%lin_i5b_E6 = ri_branch%ele(i)%lin_i5b_E6 + rad_int_by_ele%branch(ib)%ele(k)%lin_i5b_E6
    enddo
  enddo
endif

call cpu_time(time1)
if (bmad_com%debug) print '(a, f12.2)', 'radiation_integrals execution time:', time1 - time0

!------------------------------------------------------------------------
contains

subroutine cache_this_point (branch, z_start, z1, z_here, orb_start, orb_here, ele2, mat6, vec0, c_pt, orb_end, ele_end)

type (branch_struct) branch
type (ele_struct) ele2
type (rad_int_track_point_struct) :: c_pt
type (coord_struct) orb_start, orb_here, orb_end, orb_end1
type (ele_struct) ele_end

real(rp) mat6(6,6), vec0(6), z_start, z1, z_here
logical reuse_ele_end

!

c_pt%s_body = z_here
reuse_ele_end = ((ele2%key == wiggler$ .or. ele2%key == undulator$) .and. ele2%tracking_method == bmad_standard$)

call twiss_and_track_intra_ele (ele2, branch%param, z_start, z_here, .true., .false., &
                                                            orb_here, orb_end,  ele_start, ele_end, reuse_ele_end = reuse_ele_end)

call concat_transfer_mat (ele_end%mat6, ele_end%vec0, mat6, vec0, mat6, vec0)
c_pt%mat6 = mat6
c_pt%vec0 = vec0
c_pt%ref_orb_in = orb_start
c_pt%ref_orb_out = orb_end

if ((ele2%key == wiggler$ .or. ele2%key == undulator$) .and. &
            (ele2%field_calc == planar_model$ .or. ele2%field_calc == helical_model$) .and. ele2%tracking_method == bmad_standard$) then
  call calc_wiggler_g_params (ele2, branch%param, z_here, orb_end, c_pt)

else
  call twiss_and_track_intra_ele (ele2, branch%param, z_start, z1, .true., .false., orb_here, orb_end1, ele_start)
  c_pt%g_x0 = -(orb_end1%vec(2) - orb_end%vec(2)) / (z1 - z_here)
  c_pt%g_y0 = -(orb_end1%vec(4) - orb_end%vec(4)) / (z1 - z_here)

  if (ele2%key == sbend$ .or. ele2%key == rf_bend$) then
    c_pt%g_x0 = c_pt%g_x0 + ele2%value(g$) * cos(ele2%value(ref_tilt$))
    c_pt%g_y0 = c_pt%g_y0 + ele2%value(g$) * sin(ele2%value(ref_tilt$))
  endif
endif

call multipole_ele_to_ab (ele2, .true., ix_pole_max, a_pole, b_pole, include_kicks = include_kicks$)

c_pt%dgx_dx = 0
c_pt%dgx_dy = 0
c_pt%dgy_dx = 0
c_pt%dgy_dy = 0

do ip = 0, ix_pole_max
  if (a_pole(ip) == 0 .and. b_pole(ip) == 0) cycle
  call ab_multipole_kick (a_pole(ip), b_pole(ip), ip, branch%param%particle, ele2%orientation, orb_end, kx, ky, dk)
  c_pt%dgx_dx = c_pt%dgx_dx - dk(1,1) / ele2%value(l$)
  c_pt%dgx_dy = c_pt%dgx_dy - dk(1,2) / ele2%value(l$)
  c_pt%dgy_dx = c_pt%dgy_dx - dk(2,1) / ele2%value(l$)
  c_pt%dgy_dy = c_pt%dgy_dy - dk(2,2) / ele2%value(l$)
enddo

end subroutine cache_this_point

!------------------------------------------------------------------------
! contains

subroutine cache_fill (ele2, branch, c_pt, track, ix, ix0, ix1)

type (ele_struct) ele2
type (branch_struct) branch
type (rad_int_track_point_struct) :: c_pt
type (track_struct) track

integer ix, ix0, ix1
real(rp) dorb(6), dmat(6,6), denom

!

c_pt%s_body      = track%pt(ix)%s_body
c_pt%ref_orb_out = track%pt(ix)%orb
c_pt%mat6        = track%pt(ix)%mat6
c_pt%vec0        = track%pt(ix)%vec0
call calc_wiggler_g_params (ele2, branch%param, c_pt%s_body, c_pt%ref_orb_out, c_pt)

end subroutine cache_fill

!------------------------------------------------------------------------
! contains

subroutine allocate_cache (cache_ele, n_pt)

type (rad_int_cache1_struct) :: cache_ele
integer n_pt

!

cache_ele%n_pt = n_pt

if (allocated(cache_ele%pt)) then
  if (ubound(cache_ele%pt, 1) < n_pt) deallocate (cache_ele%pt)
endif

if (.not. allocated (cache_ele%pt)) allocate (cache_ele%pt(0:n_pt))

end subroutine allocate_cache

!------------------------------------------------------------------------
! contains

subroutine set_tracking_method_for_element_integration(ele)

implicit none

type (ele_struct) ele

! The Taylor tracking methods must be changed since it is not possible to partially track with these.
! The symp_lie_ptc tracking method is slow so change this method also.

select case (ele%tracking_method)
case (taylor$, symp_lie_ptc$)
  if (ele%field_calc == fieldmap$ .and. associated(ele%cartesian_map)) then
    ele%tracking_method = symp_lie_bmad$
  else
    select case (ele%key)
    case (wiggler$, undulator$, em_field$)
      ele%tracking_method = runge_kutta$
    case default
      ele%tracking_method = bmad_standard$
    end select
  endif
end select

!

select case (ele%mat6_calc_method)
case (taylor$, symp_lie_ptc$)
  if (ele%field_calc == fieldmap$ .and. associated(ele%cartesian_map)) then
    ele%mat6_calc_method = symp_lie_bmad$
  else
    select case (ele%key)
    case (wiggler$, undulator$, em_field$)
      ele%mat6_calc_method = tracking$
    case default
      ele%mat6_calc_method = bmad_standard$
    end select
  endif
end select

!

!if ((ele%key == wiggler$ .or. ele%key == undulator$) .and. ele%tracking_method == bmad_standard$ .and. &
!    (ele%field_calc == planar_model$ .or. ele%field_calc == helical_model$)) then
!  ele%tracking_method = symp_lie_bmad$
!  ele%mat6_calc_method = tracking$
!endif

end subroutine

end subroutine
