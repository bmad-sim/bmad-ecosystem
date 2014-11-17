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
! To speed up this calculation the ix_cache argument can be used to stash values for 
! the bending radius, etc. so that repeated calls to radiation_integrals consumes less time.
! To use caching: 
!   1) First call radiation_integrals with ix_cache set to 0. 
!      radiation_integrals will cache values needed to compute the integrals 
!      for the wigglers and assign ix_cache a unique number that is used to 
!      point to the cache. 
!   2) Subsequent calls to radiation_integrals should just pass the value of 
!      ix_cache that has been assigned.
! A new cache, with a unique value for ix_cache, is created each time 
! radiation_integrals is called with ix_cache = 0. This is useful if there
! are multiple lattices. To release the memory
! associated with a cache call release_rad_int_cache(ix_cache).
!
! NOTE: The validity of the cache is dependent upon the orbit and parameters like the
! quadrupole strengths not varying too much but is independent of variations in
! the Twiss parameters. 
!
! RECOMMENDATION: Do not use caching of non-wiggler elements unless you need it.
!
! Modules needed:
!   use bmad
!
! Input:
!   lat        -- lat_struct: Lattice to use. The calculation assumes that 
!                    the Twiss parameters have been calculated.
!   orbit(0:)  -- Coord_struct: Closed orbit. for the branch.
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
!       %emittance      -- Emittance
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
!     %ele(:)          -- Array of rad_int1_struct structures, one for each element in the branch.
!       %i0              -- I0 integral for the element. See the Bmad manual.
!       %i1              -- I1 integral for the element.
!       %i2              -- I2 integral for the element.
!       %i3              -- I3 integral for the element.
!       %i4a             -- "A" mode I4 integral for the element.
!       %i4b             -- "B" mode I4 integral for the element.
!       %i5a             -- "A" mode I5 integral for the element.
!       %i5b             -- "B" mode I5 integral for the element.
!       %lin_i2_E4       -- I2 * gamma^4 integral for the element.
!       %lin_i3_E7       -- I3 * gamma^7 integral for the element.
!       %lin_i5a_E6      -- I5a * gamma^6 integral for the element.
!       %lin_i5b_E6      -- I5b * gamma^6 integral for the element.
!       %lin_norm_emit_a -- "A" mode emittance. Running sum from the beginning of the branch.
!       %lin_norm_emit_b -- "B" mode emittance. Running sum from the beginning of the branch.
!
!
! Notes:
!
! 1) %synch_int(1) = momentum_compaction * lat_length
!
! 2) The lin_norm_emit values are running sums from the beginning of the 
!    lattice and include the beginning emittance stored in lat%a%emit and lat%b%emit.
!-       

subroutine radiation_integrals (lat, orbit, mode, ix_cache, ix_branch, rad_int_by_ele)

use radiation_mod, except_dummy => radiation_integrals
use symp_lie_mod, only: symp_lie_bmad
use transfer_map_mod, only: concat_transfer_mat

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele, slave
type (ele_struct) :: ele2, ele_start, ele_end
type (coord_struct), target :: orbit(0:), orb_start, orb_end, orb_end1
type (normal_modes_struct) mode
type (bmad_common_struct) bmad_com_save
type (track_struct) :: track
type (rad_int_all_ele_struct), optional :: rad_int_by_ele
type (rad_int_all_ele_struct), target :: rad_int_all
type (rad_int_info_struct) :: ri_info
type (rad_int_cache_struct), pointer :: cache
type (rad_int_cache1_struct), pointer :: cache_ele ! pointer to cache in use
type (rad_int_track_point_struct), pointer :: c_pt
type (rad_int_track_point_struct) pt
type (rad_int1_struct) int_tot
type (rad_int1_struct), pointer :: rad_int1

real(rp) i1, i2, i3, i4a, i4b, i4z, i5a, i5b, i6b, m65, G_max, g3_ave
real(rp) theta, energy, gamma2_factor, energy_loss, arg, ll, gamma_f
real(rp) a_pole(0:n_pole_maxx), b_pole(0:n_pole_maxx), dk(2,2), kx, ky
real(rp) v(4,4), v_inv(4,4), del_z, z_here, z_start, mc2, gamma, gamma4, gamma6
real(rp) kz, fac, c, s, factor, g2, g_x0, dz, z1, const_q, mat6(6,6), vec0(6)
! Cf: Sands Eq 5.46 pg 124.
real(rp), parameter :: const_q_factor = 55 * h_bar_planck * c_light / (32 * sqrt_3) 


integer, optional :: ix_cache, ix_branch
integer i, j, k, ip, ir, key, key2, n_step

character(20) :: r_name = 'radiation_integrals'

logical do_alloc, use_cache, init_cache, cache_only_wig, err, has_nonzero_pole
logical, parameter :: t = .true., f = .false.

!---------------------------------------------------------------------
! Init
! To make the calculation go faster turn off radiation fluctuations and damping

if (present(ix_branch)) then
  branch => lat%branch(ix_branch)
else
  branch => lat%branch(0)
endif

bmad_com_save = bmad_com
bmad_com%radiation_fluctuations_on = .false.
bmad_com%radiation_damping_on = .false.
bmad_com%space_charge_on = .false.

if (allocated(rad_int_all%ele)) then
  if (ubound(rad_int_all%ele, 1) < branch%n_ele_max) then 
    deallocate (rad_int_all%ele)
    do_alloc = .true.
  else
    do_alloc = .false.
  endif
else
  do_alloc = .true.
endif

if (do_alloc) then
  allocate (rad_int_all%ele(0:branch%n_ele_max))
endif

ri_info%branch => branch
ri_info%orbit => orbit

rad_int_all%ele(:) = rad_int1_struct()

m65 = 0
mode%rf_voltage = 0
int_tot = rad_int1_struct()

call init_ele (ele2)
call init_ele (ele_start)
call init_ele (ele_end)

!---------------------------------------------------------------------
! Caching

! find a cache

use_cache = .true.
init_cache = .true.
cache_only_wig = .false.

if (.not. present(ix_cache)) then
  cache => rad_int_cache_common(0)
  cache_only_wig = .true.
elseif (ix_cache == -1) then
  cache => rad_int_cache_common(0)
  cache_only_wig = .true.
elseif (ix_cache == 0) then
  do i = 1, ubound(rad_int_cache_common, 1)
    if (rad_int_cache_common(i)%set) cycle
    rad_int_cache_common(i)%set = .true.
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

  do i = 1, branch%n_ele_track

    ele => branch%ele(i)
    if (ele%value(l$) == 0) cycle

    select case (ele%key)
    case (wiggler$, undulator$)
    case (quadrupole$, sol_quad$, sbend$, em_field$, sad_mult$, hkicker$, vkicker$)
      if (cache_only_wig) cycle
    case default
      if (cache_only_wig) cycle
      if (attribute_index(ele, 'HKICK') == 0) cycle   ! Has no kick attributes.
      if (ele%value(hkick$) == 0 .or. ele%value(vkick$) == 0) cycle
    end select

    cache_ele => cache%c_ele(i)
    cache_ele%cache_type = cache_no_misalign$

    if (branch%ele(i)%bookkeeping_state%rad_int /= stale$) cycle
    branch%ele(i)%bookkeeping_state%rad_int = ok$

    ! Calculation is effectively done in element reference frame with ele2 having
    ! no offsets.

    ele2 = branch%ele(i)
    key2 = ele2%key
    if (key2 == undulator$) key2 = wiggler$

    call zero_ele_offsets (ele2)
    orb_start = orbit(i-1)
    call offset_particle (branch%ele(i), branch%param, set$, orb_start, set_multipoles = .false., set_hvkicks = .false.)

    if (key2 == wiggler$ .and. ele2%sub_key == periodic_type$) then
      n_step = max(nint(10 * ele2%value(l$) / ele2%value(l_pole$)), 1)
      del_z = ele2%value(l$) / n_step
      ! bmad_standard will not properly do partial track through a periodic_type wiggler so
      ! switch to symp_lie_bmad type tracking.
      ele2%tracking_method = symp_lie_bmad$  
      ele2%mat6_calc_method = symp_lie_bmad$  
    else
      call compute_even_steps (ele2%value(ds_step$), ele2%value(l$), &
                                    bmad_com%default_ds_step, del_z, n_step)
    endif

    cache_ele%del_z = del_z
    if (allocated(cache_ele%pt)) then
      if (ubound(cache_ele%pt, 1) < n_step) deallocate (cache_ele%pt)
    endif
    if (.not. allocated (cache_ele%pt)) allocate (cache_ele%pt(0:n_step))

    ! map_type wiggler

    if (key2 == wiggler$ .and. ele2%sub_key == map_type$) then
      call symp_lie_bmad (ele2, branch%param, orb_start, orb_end, calc_mat6 = .true., track = track)
      do k = 0, track%n_pt
        c_pt => cache_ele%pt(k)
        z_here = track%orb(k)%s - (ele2%s - ele2%value(l$)) 
        orb_end = track%orb(k)
        call calc_wiggler_g_params (ele2, z_here, orb_end, c_pt, ri_info)
        c_pt%mat6 = track%map(k)%mat6
        c_pt%vec0 = track%map(k)%vec0
        c_pt%ref_orb_in  = orb_start
        c_pt%ref_orb_out = orb_end
      enddo

    ! non-wiggler element

    else  

      z_start = 0
      ele_start = branch%ele(i-1)
      dz = min (1e-3_rp, cache_ele%del_z/3)
      cache_ele%pt(:)%ref_orb_in  = orb_start
      call mat_make_unit (mat6)
      vec0 = 0

      do k = 0, n_step

        z_here = k * cache_ele%del_z
        z1 = z_here + dz
        if (z1 > ele2%value(l$)) z1 = max(0.0_rp, z_here - dz)

        c_pt => cache_ele%pt(k)

        call twiss_and_track_intra_ele (ele2, branch%param, z_start, z_here, .true., .false., orb_start, orb_end,  ele_start, ele_end)
        call twiss_and_track_intra_ele (ele2, branch%param, z_start, z1,     .true., .false., orb_start, orb_end1, ele_start)

        z_start = z1
        orb_start = orb_end
        ele_start = ele_end

        call concat_transfer_mat (ele_end%mat6, ele_end%vec0, mat6, vec0, mat6, vec0)
        c_pt%mat6 = mat6
        c_pt%vec0 = vec0

        c_pt%ref_orb_out = orb_end
        c_pt%g_x0 = -(orb_end1%vec(2) - orb_end%vec(2)) / (z1 - z_here)
        c_pt%g_y0 = -(orb_end1%vec(4) - orb_end%vec(4)) / (z1 - z_here)
        c_pt%dgx_dx = 0
        c_pt%dgx_dy = 0
        c_pt%dgy_dx = 0
        c_pt%dgy_dy = 0

        if (key2 == quadrupole$ .or. key2 == sol_quad$ .or. key2 == bend_sol_quad$) then
          c_pt%dgx_dx =  ele2%value(k1$)
          c_pt%dgy_dy = -ele2%value(k1$)

        elseif (key2 == sbend$) then
          c_pt%g_x0   =  c_pt%g_x0 + ele2%value(g$)
          c_pt%dgx_dx =  ele2%value(k1$)
          c_pt%dgy_dy = -ele2%value(k1$)
        endif

        call multipole_ele_to_ab (ele2, .false., has_nonzero_pole, a_pole, b_pole)
        if (has_nonzero_pole) then
          do ip = 0, ubound(a_pole, 1)
            if (a_pole(ip) == 0 .and. b_pole(ip) == 0) cycle
            call ab_multipole_kick (a_pole(ip), b_pole(ip), ip, orb_end, kx, ky, dk)
            c_pt%dgx_dx = c_pt%dgx_dx - dk(1,1) / ele2%value(l$)
            c_pt%dgx_dy = c_pt%dgx_dy - dk(1,2) / ele2%value(l$)
            c_pt%dgy_dx = c_pt%dgy_dx - dk(2,1) / ele2%value(l$)
            c_pt%dgy_dy = c_pt%dgy_dy - dk(2,2) / ele2%value(l$)
          enddo
        endif

      enddo
    endif

  enddo
  
endif ! (init_cache)

!---------------------------------------------------------------------
! Loop over all elements
! We do the elements that can be integrated quickly to establish a baseline for setting 
! the error tolerance for the elements that take more time to integrate through.

do ir = 1, branch%n_ele_track

  ele => branch%ele(ir)
  if (.not. ele%is_on) cycle

  ri_info%ele => ele
  rad_int1 => rad_int_all%ele(ir)

  nullify (ri_info%cache_ele)
  if (use_cache) then
    if (cache%c_ele(ir)%cache_type /= no_cache$) ri_info%cache_ele => cache%c_ele(ir)
  endif

  pt%g_x0 = 0
  pt%g_y0 = 0
  pt%dgx_dx = 0
  pt%dgx_dy = 0
  pt%dgy_dx = 0
  pt%dgy_dy = 0

  key = ele%key
  if (key == undulator$) key = wiggler$

  if (key == rfcavity$) then
    m65 = m65 + ele%mat6(6,5)
    mode%rf_voltage = mode%rf_voltage + ele%value(voltage$)
  endif

  ll = ele%value(l$)
  if (ll == 0) cycle

  ! custom

  if (key == custom$) then
    call radiation_integrals_custom (lat, ir, orbit, err)
    cycle
  endif

  ! map type wigglers get handled later.

  if (key == wiggler$ .and. ele%sub_key == map_type$) cycle

  ! for an periodic type wiggler we make the approximation that the variation of G is
  ! fast compaired to the variation in eta.

  if (key == wiggler$ .and. ele%sub_key == periodic_type$) then
    if (ele%value(l_pole$) == 0) cycle        ! Cannot do calculation
    G_max = sqrt(2*abs(ele%value(k1$)))       ! 1/rho at max B
    g3_ave = 4 * G_max**3 / (3 * pi)
    rad_int1%i0 = (ele%value(e_tot$) / mass_of(branch%param%particle)) * 2 * G_max / 3
    rad_int1%i1 = -ele%value(k1$) * (ele%value(l_pole$) / pi)**2
    rad_int1%i2 = ll * G_max**2 / 2
    rad_int1%i3 = ll * g3_ave

    call qromb_rad_int (branch%param, [F, F, F, T, T, T, T, F, T], pt, ri_info, int_tot, rad_int1)
    cycle
  endif

  ! Only possibilities left are quad, sol_quad and sbend elements, or there
  ! is a non-zero bend angle.

  if (ele%key == patch$) cycle
  if (ele%value(hkick$) == 0 .and. ele%value(vkick$) == 0 .and. &
          key /= quadrupole$ .and. key /= sol_quad$ .and. key /= sbend$ .and. &
          key /= hkicker$ .and. key /= vkicker$) cycle

  if (key == sbend$) then
    theta = ele%value(ref_tilt_tot$) + ele%value(roll_tot$)
    pt%g_x0 =  cos(theta) * ele%value(g$)
    pt%g_y0 = -sin(theta) * ele%value(g$)
    g2 = ele%value(g$)**2
    pt%dgx_dx = ele%value(k1$) * cos(2*theta)
    pt%dgy_dx = ele%value(k1$) * sin(2*theta)
    pt%dgx_dy = pt%dgy_dx
    pt%dgy_dy = -pt%dgx_dx
    ! Edge effects for a bend. In this case we ignore any rolls.
    call propagate_part_way (orbit(ir-1), branch%param, pt, ri_info, 0.0_rp, 1, 1)
    rad_int1%i4a = -ri_info%eta_a(1) * g2 * tan(ele%value(e1$))
    rad_int1%i4b = -ri_info%eta_b(1) * g2 * tan(ele%value(e1$))
    call propagate_part_way (orbit(ir-1), branch%param, pt, ri_info, ll, 1, 2)
    rad_int1%i4a = rad_int1%i4a - ri_info%eta_a(1) * g2 * tan(ele%value(e2$))
    rad_int1%i4b = rad_int1%i4a - ri_info%eta_b(1) * g2 * tan(ele%value(e2$))
  endif

  if (key == quadrupole$ .or. key == sol_quad$) then
    theta = ele%value(tilt_tot$)
    pt%dgx_dx = ele%value(k1$) * cos(2*theta)
    pt%dgy_dx = ele%value(k1$) * sin(2*theta)
    pt%dgx_dy = pt%dgy_dx
    pt%dgy_dy = -pt%dgx_dx
  endif

  ! Integrate for quads, bends and nonzero kicks

  call qromb_rad_int (branch%param, [T, T, T, T, T, T, T, T, T], pt, ri_info, int_tot, rad_int1)

enddo

!----------------------------------------------------------
! For elements that take more time to integrate through.

do ir = 1, branch%n_ele_track

  ele => branch%ele(ir)

  if (.not. ele%is_on) cycle

  select case (ele%key)
  case (wiggler$, undulator$) 
    if (ele%sub_key /= map_type$) cycle
  case (sad_mult$, em_field$)
  case default
    cycle
  end select

  ri_info%ele => ele
  rad_int1 => rad_int_all%ele(ir)

  nullify (ri_info%cache_ele)
  if (use_cache) then
    if (cache%c_ele(ir)%cache_type /= no_cache$) ri_info%cache_ele => cache%c_ele(ir)
  endif

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

rad_int_all%ele%i4z = rad_int_all%ele%i4a + rad_int_all%ele%i4b

do i = 0, branch%n_ele_track
  gamma = branch%ele(i)%value(e_tot$) / mc2
  gamma4 = gamma**4
  gamma6 = gamma4 * gamma**2
  rad_int_all%ele(i)%lin_i2_E4  = rad_int_all%ele(i)%i2 * gamma4
  rad_int_all%ele(i)%lin_i3_E7  = rad_int_all%ele(i)%i3 * gamma6 * gamma
  rad_int_all%ele(i)%lin_i5a_E6 = rad_int_all%ele(i)%i5a * gamma6
  rad_int_all%ele(i)%lin_i5b_E6 = rad_int_all%ele(i)%i5b * gamma6
  mode%lin%i2_E4  = mode%lin%i2_E4  + rad_int_all%ele(i)%lin_i2_E4
  mode%lin%i3_E7  = mode%lin%i3_E7  + rad_int_all%ele(i)%lin_i3_E7
  mode%lin%i5a_E6 = mode%lin%i5a_E6 + rad_int_all%ele(i)%lin_i5a_E6
  mode%lin%i5b_E6 = mode%lin%i5b_E6 + rad_int_all%ele(i)%lin_i5b_E6
  rad_int_all%ele(i)%lin_norm_emit_a = branch%a%emit * gamma + factor * mode%lin%i5a_E6
  rad_int_all%ele(i)%lin_norm_emit_b = branch%b%emit * gamma + factor * mode%lin%i5b_E6
enddo

do i = branch%n_ele_track+1, branch%n_ele_max
  ele => branch%ele(i)
  if (ele%lord_status /= super_lord$) cycle
  do j = 1, ele%n_slave
    slave => pointer_to_slave(ele, j)
    k = slave%ix_ele
    rad_int_all%ele(i)%i0 = rad_int_all%ele(i)%i0 + rad_int_all%ele(k)%i0
    rad_int_all%ele(i)%i1 = rad_int_all%ele(i)%i1 + rad_int_all%ele(k)%i1
    rad_int_all%ele(i)%i2 = rad_int_all%ele(i)%i2 + rad_int_all%ele(k)%i2
    rad_int_all%ele(i)%i3 = rad_int_all%ele(i)%i3 + rad_int_all%ele(k)%i3
    rad_int_all%ele(i)%i4a = rad_int_all%ele(i)%i4a + rad_int_all%ele(k)%i4a
    rad_int_all%ele(i)%i4b = rad_int_all%ele(i)%i4b + rad_int_all%ele(k)%i4b
    rad_int_all%ele(i)%i4z = rad_int_all%ele(i)%i4z + rad_int_all%ele(k)%i4z
    rad_int_all%ele(i)%i5a = rad_int_all%ele(i)%i5a + rad_int_all%ele(k)%i5a
    rad_int_all%ele(i)%i5b = rad_int_all%ele(i)%i5b + rad_int_all%ele(k)%i5b
    rad_int_all%ele(i)%i6b = rad_int_all%ele(i)%i6b + rad_int_all%ele(k)%i6b
    rad_int_all%ele(i)%n_steps = rad_int_all%ele(i)%n_steps + rad_int_all%ele(k)%n_steps
    rad_int_all%ele(i)%lin_i2_E4 = rad_int_all%ele(i)%lin_i2_E4 + rad_int_all%ele(k)%lin_i2_E4
    rad_int_all%ele(i)%lin_i3_E7 = rad_int_all%ele(i)%lin_i3_E7 + rad_int_all%ele(k)%lin_i3_E7
    rad_int_all%ele(i)%lin_i5a_E6 = rad_int_all%ele(i)%lin_i5a_E6 + rad_int_all%ele(k)%lin_i5a_E6
    rad_int_all%ele(i)%lin_i5b_E6 = rad_int_all%ele(i)%lin_i5b_E6 + rad_int_all%ele(k)%lin_i5b_E6
  enddo
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

if (i2 /= 0) then

  mode%a%emittance = const_q * gamma2_factor * i5a / (i2 - i4a)
  mode%b%emittance = const_q * (gamma2_factor * i5b + 13 * i6b / 55) / (i2 - i4b)

  mode%a%j_damp = 1 - i4a / i2
  mode%b%j_damp = 1 - i4b / i2
  mode%z%j_damp = 2 + i4z / i2

  arg = (const_q * i3 * gamma2_factor / (2*i2 + i4z))
  if (arg > 0) then
    mode%sigE_E = sqrt(arg)
  else
    mode%sigE_E = 1e30  ! Something large
  endif

endif

mode%a%alpha_damp = energy_loss * mode%a%j_damp / (2 * energy)
mode%b%alpha_damp = energy_loss * mode%b%j_damp / (2 * energy)
mode%z%alpha_damp = energy_loss * mode%z%j_damp / (2 * energy)

mode%e_loss = energy_loss

if (m65*i1 > 0) then
  mode%sig_z = sqrt(i1/m65) * mode%sigE_E
else   ! Unstable
  mode%sig_z = 1e30  ! Something large
endif

mode%z%emittance = mode%sig_z * mode%sigE_E

bmad_com = bmad_com_save

if (present(rad_int_by_ele)) call move_alloc (rad_int_all%ele, rad_int_by_ele%ele)

call deallocate_ele_pointers(ele2)
call deallocate_ele_pointers(ele_start)
call deallocate_ele_pointers(ele_end)

end subroutine
