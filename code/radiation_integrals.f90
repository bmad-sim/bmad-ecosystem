!+
! Subroutine radiation_integrals (ring, orb_, mode, ix_cache)
!
! Subroutine to calculate the synchrotron radiation integrals along with the
! emittance, and energy spread.
!
! Note: A negative emittance is possible and just means that the beam is
! unstable. That is, you have a negative damping partition number.
!
! The calculation can spend a significant amount of time calculating the
! integrals through any wigglers. To speed up this calculation the ix_cache
! argument can be used to stash values for the bending radius, etc. through a
! wiggler so that repeated calls to radiation_integrals consume less time. 
! [If radiation_integrals is going to be called only once then do not use this
! feature. It will actually slow things down in this case.] 
! To use caching: 
!   1) First call radiation_integrals with ix_cache set to 0. 
!      radiation_integrals will cache values needed to compute the integrals 
!      for the wigglers and assign ix_cache a unique number that is used to 
!      point to the cache. 
!   2) Subsequent calls to radiation_integrals should just pass the value of 
!      ix_cache that has been assigned.
! A new cache, with a unique value for ix_cache, is created each time 
! radiation_integrals is called with ix_cache = 0. To release the memory
! associated with a cache call release_rad_int_cache(ix_cache).
! Note: The validity of the cache is dependent upon the orbit being (more or 
! less) constant but is independent of the Twiss parameters.
!
! Modules needed:
!   use bmad
!
! Input:
!   ring      -- Ring_struct: Ring to use. The calculation assumes that 
!                    the Twiss parameters have been calculated.
!   orb_(0:)  -- Coord_struct: Closed orbit.
!   ix_cache   -- Integer, optional: Cache pointer.
!                              == 0 --> Create a new cache.
!                              /= 0 --> Use the corresponding cache. 
! Output:
!   mode -- Modes_struct: Parameters for the ("horizontal like") a-mode,
!                              ("vertical like") b-mode, and the z-mode
!     %synch_int(1:3) -- Synchrotron integrals.
!     %sig_e          -- Sigma_E/E energy spread
!     %sig_z          -- Bunch Length
!     %e_loss         -- Energy loss in eV per turn
!     %a, %b, %z      -- Amode_struct: Substructure
!       %emittance      -- Emittance
!       %synch_int(4:5) -- Synchrotron integrals
!       %j_damp         -- Damping partition factor
!       %alpha_damp     -- Exponential damping coefficient per turn
!   ix_cache -- Integer, optional: Cache pointer. If ix_cache = 0 at input then
!                   ix_cache is set to a unique number. Otherwise ix_cache 
!                   is not changed.
!                  
!
! Notes:
!   1) %synch_int(1) = momentum_compaction * ring_length
!
!   2) There is a common block where the integrals for the individual elements
!      are saved. To access this common block a use statement is needed:
!         use rad_int_common
!      In the common block:
!         ric%i1_(:)  -- I1 integral for each element.
!         ric%i2_(:)  -- I2 integral for each element.
!         ric%i3_(:)  -- I3 integral for each element.
!         ric%i4a_(:) -- "A" mode I4 integral for each element.
!         ric%i4b_(:) -- "B" mode I4 integral for each element.
!         ric%i5a_(:) -- "A" mode I5 integral for each element.
!         ric%i5b_(:) -- "B" mode I5 integral for each element.
!-       

#include "CESR_platform.inc"

subroutine radiation_integrals (ring, orb_, mode, ix_cache)
                     
  use precision_def
  use nr
  use rad_int_common
  use radiation_mod

  implicit none

  type (ring_struct), target :: ring
  type (coord_struct), target :: orb_(0:), start, end
  type (modes_struct) mode
  type (synch_rad_com) sr_com_save
  type (rad_int_cache_struct), pointer :: cache

  real(rp), parameter :: c_gam = 4.425e-5, c_q = 3.84e-13
  real(rp), save :: i1, i2, i3, i4a, i4b, i4z, i5a, i5b, m65, G_max, g3_ave
  real(rp) theta, energy, gamma2_factor, energy_loss, arg, ll
  real(rp) v(4,4), v_inv(4,4), f0, f1, s

  integer, optional :: ix_cache
  integer i, j, k, ix, ir, key

  logical do_alloc
  logical, parameter :: t = .true., f = .false.

!---------------------------------------------------------------------
! Init
! To make the calculation go faster turn off radiation fluctuations and damping

  sr_com_save = sr_com

  if (allocated(ric%i1_)) then
    if (ubound(ric%i1_, 1) < ring%n_ele_max) then 
      deallocate (ric%i1_)
      deallocate (ric%i2_)
      deallocate (ric%i3_)
      deallocate (ric%i4a_)
      deallocate (ric%i4b_)
      deallocate (ric%i5a_)
      deallocate (ric%i5b_)
      deallocate (ric%n_steps)
      do_alloc = .true.
    else
      do_alloc = .false.
    endif
  else
    do_alloc = .true.
  endif

  if (do_alloc) then
    allocate (ric%i1_(ring%n_ele_max))
    allocate (ric%i2_(ring%n_ele_max))
    allocate (ric%i3_(ring%n_ele_max))
    allocate (ric%i4a_(ring%n_ele_max))
    allocate (ric%i4b_(ring%n_ele_max))
    allocate (ric%i5a_(ring%n_ele_max))
    allocate (ric%i5b_(ring%n_ele_max))
    allocate (ric%n_steps(ring%n_ele_max))
  endif

  ric%ring => ring

  ric%i1_ = 0;   ric%i2_ = 0;  ric%i3_ = 0
  ric%i4a_ = 0;  ric%i4b_ = 0
  ric%i5a_ = 0;  ric%i5b_ = 0
  ric%int_tot = 0; ric%n_steps = 0

  m65 = 0

! Caching

  if (present(ix_cache)) then

    if (ix_cache == 0) then
      do i = 1, size(ric%cache)
        if (ric%cache(i)%set) cycle
        ric%cache(i)%set = .true.
        ix_cache = i
        cache => ric%cache(i)
        exit
      enddo

      if (ix_cache == 0) then
        print *, 'ERROR IN RADIATION_INTEGRALS: CACHE OUT OF MEMORY!'
        call err_exit
      endif

      j = 0
      do i = 1, ring%n_ele_ring
        ric%ele => ring%ele_(i)
        if (ric%ele%key /= wiggler$ .or. ric%ele%sub_key /= map_type$) cycle
        j = j + 1
      enddo
      allocate (cache%ele(j))

      j = 0
      do i = 1, ring%n_ele_ring
        ric%ele => ring%ele_(i)
        if (ric%ele%key /= wiggler$ .or. ric%ele%sub_key /= map_type$) cycle
        j = j + 1
        ric%ele%ixx = j
        cache%ele(j)%ix_ele = i
        cache%ele(j)%ds = ric%ele%value(l$) / ric%ele%num_steps
        allocate (cache%ele(j)%v(0:ric%ele%num_steps))
        do k = 0, ric%ele%num_steps
          s = k * cache%ele(j)%ds
          f0 = (ric%ele%value(l$) - s) / ric%ele%value(l$)
          f1 = s / ric%ele%value(l$)
          start%vec = orb_(i-1)%vec * f0 + orb_(i)%vec * f1
          call calc_g_params (s, start)
          cache%ele(j)%v(k)%g     = ric%g
          cache%ele(j)%v(k)%g2    = ric%g2
          cache%ele(j)%v(k)%g_x   = ric%g_x
          cache%ele(j)%v(k)%g_y   = ric%g_y
          cache%ele(j)%v(k)%dg2_x = ric%dg2_x
          cache%ele(j)%v(k)%dg2_y = ric%dg2_y
        enddo
      enddo
      
    else
      cache => ric%cache(ix_cache)

    endif

    ric%use_cache = .true.

  else
    ric%use_cache = .false.

  endif

!---------------------------------------------------------------------
! Loop over all elements
! We do the non-wiggler elements first since we can do this quickly.

  do ir = 1, ring%n_ele_use

    ric%ele => ring%ele_(ir)
    if (.not. ric%ele%is_on) cycle

    ric%ele0 => ring%ele_(ir-1)
    ric%orb0 => orb_(ir-1)
    ric%orb1 => orb_(ir)

    ll = ric%ele%value(l$)
    if (ll == 0) cycle

    key = ric%ele%key

    ric%g_x0 = -ric%ele%value(hkick$) / ll
    ric%g_y0 = -ric%ele%value(vkick$) / ll

    if (key == rfcavity$) m65 = m65 + ric%ele%mat6(6,5)

! custom

    if (key == custom$) then
      call custom_radiation_integrals (ring, ir, orb_)
      cycle
    endif

! exact calculation involves runge_kutta tracking.

    if (ric%ele%exact_rad_int_calc) then

      ric%d_orb%vec = (/ 1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-5 /)

      rk_com%n_pts = 0
      call track1_runge_kutta (ric%orb0, ric%ele, ric%ring%param, ric%orb1)
      call transfer_rk_track (rk_com, ric%rk_track(0))

      do i = 1, 6
        start = ric%orb0
        start%vec(i) = start%vec(i) + ric%d_orb%vec(i)
        call track1_runge_kutta (start, ric%ele, ric%ring%param, end)
        call transfer_rk_track (rk_com, ric%rk_track(i))
      enddo

! exact integration

      call qromb_rad_int ((/ T, T, T, T, T, T, T /), ir)

      cycle

    endif

! new style wigglers get handled later.

    if (key == wiggler$ .and. ric%ele%sub_key == map_type$) cycle

! for an old style wiggler we make the approximation that the variation of G is
! fast compaired to the variation in eta.

    if (key == wiggler$ .and. ric%ele%sub_key == periodic_type$) then
      G_max = sqrt(2*abs(ric%ele%value(k1$)))       ! 1/rho at max B
      g3_ave = 4 * G_max**3 / (3 * pi)

      ric%i1_(ir) = 0
      ric%i2_(ir) = ll * G_max**2 / 2
      ric%i3_(ir) = ll * g3_ave
      ric%i4a_(ir) = 0   ! Too complicated to calculate. Probably small
      ric%i4b_(ir) = 0   ! Too complicated to calculate. Probably small

      ric%g_x0 = g3_ave**(1.0/3)
      ric%g_y0 = 0
      ric%k1 = 0
      ric%s1 = 0

! integrate for periodic_type wigglers

      call qromb_rad_int ((/ F, F, F, F, F, T, T /), ir)

      cycle

    endif

!

    if (ric%g_x0 == 0 .and. ric%g_y0 == 0 .and. &
          key /= quadrupole$ .and. key /= sol_quad$ .and. key /= sbend$) cycle

    if (key == sbend$) then
      theta = ric%ele%value(tilt$) + ric%ele%value(roll$)
      ric%g_x0 = ric%g_x0 + cos(theta) * ric%ele%value(g$)
      ric%g_y0 = ric%g_y0 - sin(theta) * ric%ele%value(g$)
    endif

    ric%g2 = ric%g_x0**2 + ric%g_y0**2
    ric%g = sqrt(ric%g2)

    ric%i2_(ir)  = ric%g2 * ll
    ric%i3_(ir)  = ric%g2 * ric%g * ll

    if (key == quadrupole$ .or. key == sol_quad$) then
      theta = ric%ele%value(tilt$)
      ric%k1 = ric%ele%value(k1$) * cos(2*theta)
      ric%s1 = ric%ele%value(k1$) * sin(2*theta)
    elseif (key == sbend$) then
      theta = ric%ele%value(tilt$) + ric%ele%value(roll$)
      ric%k1 = ric%ele%value(k1$) * cos(2*theta)
      ric%s1 = ric%ele%value(k1$) * sin(2*theta)
    else
      ric%k1 = 0
      ric%s1 = 0
    endif

! Edge effects for a bend. In this case we ignore any rolls.

    if (key == sbend$) then
      call propagate_part_way(0.0_rp)
      ric%i4a_(ir) = -ric%eta_a(1) * ric%g2 * tan(ric%ele%value(e1$))
      ric%i4b_(ir) = -ric%eta_b(1) * ric%g2 * tan(ric%ele%value(e1$))
      call propagate_part_way(ll)
      ric%i4a_(ir) = ric%i4a_(ir) - &
                           ric%eta_a(1) * ric%g2 * tan(ric%ele%value(e2$))
      ric%i4b_(ir) = ric%i4a_(ir) - &
                           ric%eta_b(1) * ric%g2 * tan(ric%ele%value(e2$))
    endif

! Integrate for quads and bends

    call qromb_rad_int ((/ T, F, F, T, T, T, T /), ir)

  enddo

!----------------------------------------------------------
! For map type wigglers

  do ir = 1, ring%n_ele_use

    ric%ele => ring%ele_(ir)
    if (ric%ele%key /= wiggler$) cycle
    if (ric%ele%sub_key /= map_type$) cycle
    if (.not. ric%ele%is_on) cycle

    ric%ele0 => ring%ele_(ir-1)
    ric%orb0 => orb_(ir-1)
    ric%orb1 => orb_(ir)

    ll = ric%ele%value(l$)
    if (ll == 0) cycle

    ric%g_x0 = -ric%ele%value(hkick$) / ll
    ric%g_y0 = -ric%ele%value(vkick$) / ll

    call make_v_mats (ric%ele0, v, v_inv)
    ric%eta_a0 = &
          matmul(v, (/ ric%ele0%x%eta, ric%ele0%x%etap, 0.0_rp, 0.0_rp /))
    ric%eta_b0 = &
          matmul(v, (/ 0.0_rp, 0.0_rp, ric%ele0%y%eta, ric%ele0%y%etap /))

    call make_v_mats (ric%ele, v, v_inv)
    ric%eta_a1 = &
          matmul(v, (/ ric%ele%x%eta, ric%ele%x%etap, 0.0_rp, 0.0_rp /))
    ric%eta_b1 = &
          matmul(v, (/ 0.0_rp, 0.0_rp, ric%ele%y%eta, ric%ele%y%etap /))

    if (ric%use_cache) ric%cache_ele => cache%ele(ric%ele%ixx)

! radiation integrals calc for the map_type wiggler

    call qromb_rad_int ((/ T, T, T, T, T, T, T /), ir) 

  enddo

!---------------------------------------------------------------------
! now put everything together

  i1   = ric%int_tot(1)
  i2   = ric%int_tot(2)
  i3   = ric%int_tot(3)
  i4a  = ric%int_tot(4)
  i4b  = ric%int_tot(5)
  i5a  = ric%int_tot(6)
  i5b  = ric%int_tot(7)

  i4z = i4a + i4b

  energy = ring%param%beam_energy
  gamma2_factor = (energy * 1956.95e-9)**2
  energy_loss = 1e9 * c_gam * (1e-9 * energy)**4 * i2 / pi

  mode%synch_int(1) = i1
  mode%synch_int(2) = i2
  mode%synch_int(3) = i3

  mode%a%synch_int(4) = i4a
  mode%b%synch_int(4) = i4b
  mode%z%synch_int(4) = i4z

  mode%a%synch_int(5) = i5a
  mode%b%synch_int(5) = i5b

  if (i2 /= 0) then

    mode%a%emittance = c_q * gamma2_factor * i5a / (i2 - i4a)
    mode%b%emittance = c_q * gamma2_factor * i5b / (i2 - i4b)

    mode%a%j_damp = 1 - i4a / i2
    mode%b%j_damp = 1 - i4b / i2
    mode%z%j_damp = 2 + i4z / i2

    arg = (c_q * i3 * gamma2_factor / (2*i2 + i4z))
    mode%sig_e = sqrt(max(0.0_rp, arg))

  endif

  mode%a%alpha_damp = energy_loss * mode%a%j_damp / (2 * energy)
  mode%b%alpha_damp = energy_loss * mode%b%j_damp / (2 * energy)
  mode%z%alpha_damp = energy_loss * mode%z%j_damp / (2 * energy)

  mode%e_loss = energy_loss

  if (abs(m65) > 0) then
    mode%sig_z = sqrt(i1/abs(m65)) * mode%sig_e
  else
    mode%sig_z = 0.
  endif

  mode%z%emittance = mode%sig_z * mode%sig_e

  sr_com = sr_com_save

end subroutine
