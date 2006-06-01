!+
! Subroutine radiation_integrals (ring, orbit, mode, ix_cache)
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
!   ring       -- Ring_struct: Ring to use. The calculation assumes that 
!                    the Twiss parameters have been calculated.
!   orbit(0:)  -- Coord_struct: Closed orbit.
!   ix_cache   -- Integer, optional: Cache pointer.
!                              == 0 --> Create a new cache.
!                              /= 0 --> Use the corresponding cache. 
! Output:
!   mode -- Modes_struct: Parameters for the ("horizontal like") a-mode,
!                              ("vertical like") b-mode, and the z-mode
!     %synch_int(1:3) -- Synchrotron integrals.
!     %sigE_E         -- Sigma_E/E energy spread
!     %sig_z          -- Bunch Length
!     %e_loss         -- Energy loss in eV per turn
!     %a, %b, %z      -- Amode_struct: Substructure
!       %emittance      -- Emittance
!       %synch_int(4:5) -- Synchrotron integrals
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

subroutine radiation_integrals (ring, orbit, mode, ix_cache)

  use nr
  use rad_int_common, except => radiation_integrals
  use radiation_mod
  use symp_lie_mod, only: symp_lie_bmad
  
  implicit none

  type (ring_struct), target :: ring
  type (ele_struct), save :: runt
  type (ele_struct), pointer :: ele, ele0
  type (ele_struct) :: ele2
  type (coord_struct), target :: orbit(0:), start, end, c2
  type (modes_struct) mode
  type (bmad_com_struct) bmad_com_save
  type (rad_int_cache_struct), pointer :: cache
  type (track_struct) track

  real(rp), parameter :: c_gam = 4.425e-5, c_q = 3.84e-13
  real(rp), save :: i1, i2, i3, i4a, i4b, i4z, i5a, i5b, m65, G_max, g3_ave
  real(rp) theta, energy, gamma2_factor, energy_loss, arg, ll, gamma_f
  real(rp) v(4,4), v_inv(4,4), f0, f1, del_z, z_here, mc2, gamma, gamma4, gamma6
  real(rp) kz, fac, c, s

  integer, optional :: ix_cache
  integer i, j, k, ir, key, n_step

  character(20) :: r_name = 'radiation_integrals'

  logical do_alloc
  logical, parameter :: t = .true., f = .false.

!---------------------------------------------------------------------
! Init
! To make the calculation go faster turn off radiation fluctuations and damping

  bmad_com_save = bmad_com
  bmad_com%radiation_fluctuations_on = .false.
  bmad_com%radiation_damping_on = .false.
  bmad_com%trans_space_charge_on = .false.

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

!---------------------------------------------------------------------
! Caching

  ric%use_cache = .false.

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

      ! Count number of elements to cache & allocate memory.

      allocate (cache%ix_ele(ring%n_ele_ring))

      j = 0  ! number of elements to cache
      do i = 1, ring%n_ele_use
        key = ring%ele_(i)%key
        if (key == quadrupole$ .or. key == sol_quad$ .or. key == sbend$ .or. &
                 key == wiggler$ .or. ring%ele_(i)%value(hkick$) /= 0 .or. &
                 ring%ele_(i)%value(vkick$) /= 0) then
          j = j + 1
          cache%ix_ele(i) = j  ! mark element for caching
        else
          cache%ix_ele(i) = -1 ! do not cache this element
        endif          
      enddo
      allocate (cache%ele(j))  ! allocate cache memory

      ! Now cache the information

      j = 0 
      do i = 1, ring%n_ele_use

        if (cache%ix_ele(i) == -1) cycle

        j = j + 1
        cache%ele(j)%ix_ele = i

        ele2 = ring%ele_(i)
        if (.not. ele2%map_with_offsets) call zero_ele_offsets (ele2)

        if (ele2%key == wiggler$ .and. ele2%sub_key == periodic_type$) then
          n_step = 10
          del_z = 2 * ele2%value(l_pole$) / n_step
        else
          call compute_even_steps (ele2%value(ds_step$), ele2%value(l$), &
                                      bmad_com%default_ds_step, del_z, n_step)
        endif

        cache%ele(j)%del_z = del_z
        allocate (cache%ele(j)%pt(0:n_step))

        if (ele2%key == wiggler$) then
          if (ele2%sub_key == map_type$) then
            track%save_track = .true.
            call symp_lie_bmad (ele2, ring%param, orbit(i-1), end, &
                                          calc_mat6 = .true., track = track)
            do k = 0, track%n_pt
              z_here = track%pt(k)%s 
              end = track%pt(k)%orb
              call calc_wiggler_g_params (ele2, z_here, end, cache%ele(j)%pt(k))
              cache%ele(j)%pt(k)%orb  = end
              cache%ele(j)%pt(k)%mat6 = track%pt(k)%mat6
              cache%ele(j)%pt(k)%vec0 = track%pt(k)%vec0
            enddo

          else  ! periodic_type$
            if (ele2%value(l_pole$) == 0) then
              call out_io (s_error$, r_name, &
                              'ERROR: PERIODIC_TYPE WIGGLER: ' // ele2%name, &
                              '       HAS L_POLE = 0')
              cycle
            endif
            kz = pi / ele2%value(l_pole$)
            fac = sqrt(-2 * ele2%value(k1$))
            do k = 0, n_step
              z_here = k * del_z
              c = fac * cos (kz * z_here) 
              s = fac * sin (kz * z_here)
              cache%ele(j)%pt(k)%orb%vec = &
                (/ (c - fac) / kz**2, -s / kz, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp /)
              call mat_make_unit(cache%ele(j)%pt(k)%mat6)
              cache%ele(j)%pt(k)%mat6(1,6) = -cache%ele(j)%pt(k)%orb%vec(1)
              cache%ele(j)%pt(k)%mat6(2,6) = -cache%ele(j)%pt(k)%orb%vec(2)
              cache%ele(j)%pt(k)%vec0 = cache%ele(j)%pt(k)%orb%vec
              cache%ele(j)%pt(k)%g = abs(c)
              cache%ele(j)%pt(k)%g2 = c**2
            enddo
          endif

        else  ! non wiggler elements
          do k = 0, n_step
            z_here = k * cache%ele(j)%del_z
            call twiss_and_track_partial (ring%ele_(i-1), ele2, ring%param, &
                                            z_here, runt, orbit(i-1), end)
            cache%ele(j)%pt(k)%orb = end
            cache%ele(j)%pt(k)%mat6 = runt%mat6
            cache%ele(j)%pt(k)%vec0 = runt%vec0
          enddo
        endif

      enddo
      
    else  ! ix_cache /= 0
      cache => ric%cache(ix_cache)

    endif ! ix_cache /= 0

    ric%use_cache = .true.

  endif ! present(ix_cache)

!---------------------------------------------------------------------
! Loop over all elements
! We do the non-wiggler elements first since we can do this quickly.

  do ir = 1, ring%n_ele_use

    ele => ring%ele_(ir)
    if (.not. ele%is_on) cycle

    if (ric%use_cache) then
      if (cache%ix_ele(ir) > 0) ric%cache_ele => cache%ele(cache%ix_ele(ir))
    endif

    ele0 => ring%ele_(ir-1)
    ric%orb0 => orbit(ir-1)
    ric%orb1 => orbit(ir)

    key = ele%key

    if (key == rfcavity$) m65 = m65 + ele%mat6(6,5)

    ll = ele%value(l$)
    if (ll == 0) cycle

! If there is a non-zero hick$ or vkick$ attribute assume that the kick
! is distributed evenly over the element. The element will then contribute
! to the integrals like a bend.

    if (ele%value(hkick$) /= 0 .or. ele%value(vkick$) /= 0) then
      c2%vec = 0
      call offset_particle (ele, ring%param, c2, set$, &
                          set_canonical = .false., set_multipoles = .false.)
      call offset_particle (ele, ring%param, c2, unset$, &
                          set_canonical = .false., set_multipoles = .false.)
      ric%pt%g_x0 = -c2%vec(2) / ll
      ric%pt%g_y0 = -c2%vec(4) / ll
    else
      ric%pt%g_x0 = 0
      ric%pt%g_y0 = 0
    endif

! custom

    if (key == custom$) then
      call custom_radiation_integrals (ring, ir, orbit)
      cycle
    endif

! new style wigglers get handled later.

    if (key == wiggler$ .and. ele%sub_key == map_type$) cycle

! for an old style wiggler we make the approximation that the variation of G is
! fast compaired to the variation in eta.

    if (key == wiggler$ .and. ele%sub_key == periodic_type$) then
      G_max = sqrt(2*abs(ele%value(k1$)))       ! 1/rho at max B
      g3_ave = 4 * G_max**3 / (3 * pi)
      ric%i1_(ir) = - ele%value(k1$) * (ele%value(l_pole$) / pi)**2
      ric%i2_(ir) = ll * G_max**2 / 2
      ric%i3_(ir) = ll * g3_ave
      ric%i4a_(ir) = 0   ! Too complicated to calculate. Probably small
      ric%i4b_(ir) = 0   ! Too complicated to calculate. Probably small

      ric%pt%g_x0 = g3_ave**(1.0/3)
      ric%pt%g_y0 = 0
      ric%pt%k1   = 0
      ric%pt%s1   = 0

      call qromb_rad_int (ele0, ele, (/ F, F, F, F, F, T, T /), ir)
      cycle

    endif

! Only possibilities left are quad, sol_quad and sbend elements, or there
! is a non-zero bend angle.

    if (ric%pt%g_x0 == 0 .and. ric%pt%g_y0 == 0 .and. &
          key /= quadrupole$ .and. key /= sol_quad$ .and. key /= sbend$) cycle

    if (key == sbend$) then
      theta = ele%value(tilt_tot$) + ele%value(roll$)
      ric%pt%g_x0 = ric%pt%g_x0 + cos(theta) * ele%value(g$)
      ric%pt%g_y0 = ric%pt%g_y0 - sin(theta) * ele%value(g$)
    endif

    ric%pt%g2 = ric%pt%g_x0**2 + ric%pt%g_y0**2
    ric%pt%g = sqrt(ric%pt%g2)

    ric%i2_(ir)  = ric%pt%g2 * ll
    ric%i3_(ir)  = ric%pt%g2 * ric%pt%g * ll

    if (key == quadrupole$ .or. key == sol_quad$) then
      theta = ele%value(tilt_tot$)
      ric%pt%k1 = ele%value(k1$) * cos(2*theta)
      ric%pt%s1 = ele%value(k1$) * sin(2*theta)
    elseif (key == sbend$) then
      theta = ele%value(tilt_tot$) + ele%value(roll$)
      ric%pt%k1 = ele%value(k1$) * cos(2*theta)
      ric%pt%s1 = ele%value(k1$) * sin(2*theta)
    else
      ric%pt%k1 = 0
      ric%pt%s1 = 0
    endif

! Edge effects for a bend. In this case we ignore any rolls.

    if (key == sbend$) then
      call propagate_part_way (ele0, ele, runt, 0.0_rp, 1, 1)
      ric%i4a_(ir) = -ric%eta_a(1) * ric%pt%g2 * tan(ele%value(e1$))
      ric%i4b_(ir) = -ric%eta_b(1) * ric%pt%g2 * tan(ele%value(e1$))
      call propagate_part_way (ele0, ele, runt, ll, 1, 2)
      ric%i4a_(ir) = ric%i4a_(ir) - &
                           ric%eta_a(1) * ric%pt%g2 * tan(ele%value(e2$))
      ric%i4b_(ir) = ric%i4a_(ir) - &
                           ric%eta_b(1) * ric%pt%g2 * tan(ele%value(e2$))
    endif

! Integrate for quads and bends

    if (ric%use_cache) ric%cache_ele => cache%ele(cache%ix_ele(ir))

    call qromb_rad_int (ele0, ele, (/ T, F, F, T, T, T, T /), ir)

  enddo

!----------------------------------------------------------
! For map type wigglers

  do ir = 1, ring%n_ele_use

    ele => ring%ele_(ir)
    if (ele%key /= wiggler$) cycle
    if (ele%sub_key /= map_type$) cycle
    if (.not. ele%is_on) cycle

    ele0 => ring%ele_(ir-1)
    ric%orb0 => orbit(ir-1)
    ric%orb1 => orbit(ir)

    ll = ele%value(l$)
    if (ll == 0) cycle

    ric%pt%g_x0 = -ele%value(hkick$) / ll
    ric%pt%g_y0 = -ele%value(vkick$) / ll

    call make_v_mats (ele0, v, v_inv)
    ric%eta_a0 = &
          matmul(v, (/ ele0%x%eta, ele0%x%etap, 0.0_rp, 0.0_rp /))
    ric%eta_b0 = &
          matmul(v, (/ 0.0_rp, 0.0_rp, ele0%y%eta, ele0%y%etap /))

    call make_v_mats (ele, v, v_inv)
    ric%eta_a1 = &
          matmul(v, (/ ele%x%eta, ele%x%etap, 0.0_rp, 0.0_rp /))
    ric%eta_b1 = &
          matmul(v, (/ 0.0_rp, 0.0_rp, ele%y%eta, ele%y%etap /))

    if (ric%use_cache) ric%cache_ele => cache%ele(cache%ix_ele(ir))

    call qromb_rad_int (ele0, ele, (/ T, T, T, T, T, T, T /), ir) 

  enddo

!---------------------------------------------------------------------
! Now put everything together...
! Linac radiation integrals:

  mc2 = mass_of (ring%param%particle)
  gamma_f = ring%ele_(ring%n_ele_use)%value(beam_energy$) / mc2

  mode%lin%sig_E1 = 0
  mode%lin%i2_E4  = 0
  mode%lin%i3_E7  = 0
  mode%lin%i5a_E6 = 0
  mode%lin%i5b_E6 = 0

  do i = 1, ring%n_ele_use
    gamma = ring%ele_(i)%value(beam_energy$) / mc2
    gamma4 = gamma**4
    gamma6 = gamma4 * gamma**2
    mode%lin%i2_E4  = mode%lin%i2_E4  + ric%i2_(i) * gamma4
    mode%lin%i3_E7  = mode%lin%i3_E7  + ric%i3_(i) * gamma6 * gamma
    mode%lin%i5a_E6 = mode%lin%i5a_E6 + ric%i5a_(i) * gamma6
    mode%lin%i5b_E6 = mode%lin%i5b_E6 + ric%i5b_(i) * gamma6
  enddo

  mode%lin%sig_E1 = mc2 * sqrt (4 * c_q * r_e * mode%lin%i3_E7 / 3)
  mode%lin%a_emittance_end = 2 * c_q * r_e * mode%lin%i5a_e6 / (3 * gamma_f)
  mode%lin%b_emittance_end = 2 * c_q * r_e * mode%lin%i5b_e6 / (3 * gamma_f)

! Normal integrals

  i1   = ric%int_tot(1)
  i2   = ric%int_tot(2)
  i3   = ric%int_tot(3)
  i4a  = ric%int_tot(4)
  i4b  = ric%int_tot(5)
  i5a  = ric%int_tot(6)
  i5b  = ric%int_tot(7)

  i4z = i4a + i4b

  energy = ring%ele_(0)%value(beam_energy$)
  gamma2_factor = (energy * 1956.95e-9)**2
  energy_loss = 1e9 * c_gam * (1e-9 * mc2)**4 * mode%lin%i2_E4 / pi

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
    mode%sigE_E = sqrt(max(0.0_rp, arg))

  endif

  mode%a%alpha_damp = energy_loss * mode%a%j_damp / (2 * energy)
  mode%b%alpha_damp = energy_loss * mode%b%j_damp / (2 * energy)
  mode%z%alpha_damp = energy_loss * mode%z%j_damp / (2 * energy)

  mode%e_loss = energy_loss

  if (abs(m65) > 0) then
    mode%sig_z = sqrt(i1/abs(m65)) * mode%sigE_E
  else
    mode%sig_z = 0.
  endif

  mode%z%emittance = mode%sig_z * mode%sigE_E

  bmad_com = bmad_com_save

end subroutine
