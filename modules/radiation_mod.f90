#include "CESR_platform.inc"

module radiation_mod

  use bmad_struct
  use bmad_interface
  use runge_kutta_mod

  type synch_rad_com
    real(rp) :: scale = 1.0               ! used to scale the radiation
    logical :: damping_on = .false.       ! Radiation damping toggle.
    logical :: fluctuations_on = .false.  ! Radiation fluctuations toggle.
  end type

  type (synch_rad_com), save :: sr_com

contains

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine track1_radiation (start, ele, param, end, edge)
!
! Subroutine to put in radiation dampling and/or fluctuations.
! This routine calculates half the radiation of an element so this routine
! needs to be called before entering an element and just after exiting.
!
! Use the setup_radiation_tracking routine initially prior to tracking.
!
! Modules needed:
!   use bmad
!
! Input:
!   start -- Coord_struct: Starting particle position.
!   ele   -- Ele_struct: Element that causes radiation.
!   edge  -- Integer: Where the particle is: start_edge$ or end_edge$.
!
! Output:
!   end   -- Ele_struct: Particle position after radiation has been applied.
!-

subroutine track1_radiation (start, ele, param, end, edge)

  implicit none

  type (coord_struct), intent(in) :: start
  type (ele_struct), intent(in) :: ele
  type (param_struct), intent(in) :: param
  type (coord_struct), intent(out) :: end

  integer, intent(in) :: edge

!

  type (coord_struct) start2

  real(rp), save :: z_start
  real(rp) s_len, h, h2, h3, h_x, h_y, f, this_ran
  real(rp) x_ave, y_ave, gamma, dE_p, h_bar, fluct_const, fact_d, fact_f

  integer direc

  logical set
  logical :: init_needed = .true.

! Init

  if (init_needed) then
    h_bar = h_bar_planck * c_light**2 / e_charge  ! h_bar in eV*sec
    fluct_const = 55 * r_e * h_bar * c_light / &
                               (24 * sqrt(3.0) * m_electron**6)
    init_needed = .false.
  endif

! If not a magnetic element then nothing to do

  if (.not. any (ele%key == (/ quadrupole$, sextupole$, octupole$, sbend$, &
                                                sol_quad$, wiggler$ /))) then
    end = start
    return
  endif

! The total radiation length is the element length + any change in path length.
! If entering the element then the length over which radiation is generated
! is taken to be 1/2 the element length.
! If leaving the element the radiation length is taken to be 1/2 the element
! length + delta_Z

  if (edge == start_edge$) then
    set = set$
    s_len = ele%value(l$) / 2
    direc = +1
    z_start = start%vec(5)
  elseif (edge == end_edge$) then
    set = unset$
    s_len = ele%value(l$)/2 + (start%vec(5) - z_start)
    direc = -1
  else
    print *, 'ERROR IN RADIATION_MOD: BAD EDGE ARGUMENT:', set
    call err_exit
  endif

! Get the coords in the frame of reference of the element

  start2 = start
  call offset_particle (ele, param, start2, set)

! Calculate the radius of curvature for an on-energy particle

  select case (ele%key)

  case (quadrupole$, sol_quad$)
    x_ave = start2%vec(1) + direc * start2%vec(2) * ele%value(l$) / 4
    y_ave = start2%vec(3) + direc * start2%vec(4) * ele%value(l$) / 4
    h_x =  ele%value(k1$) * x_ave
    h_y = -ele%value(k1$) * y_ave
    h2 = h_x**2 + h_y**2

  case (sextupole$)
    h = ele%value(k2$) * (start2%vec(1)**2 + start2%vec(3)**2)
    h2 = h**2

  case (octupole$)
    h2 = ele%value(k3$)**2 * (start2%vec(1)**2 + start2%vec(3)**2)**3

  case (sbend$)
    if (ele%value(k1$) == 0) then
      h = ele%value(g$) + ele%value(delta_g$)
      h2 = h**2 
    else
      x_ave = start2%vec(1) + direc * start2%vec(2) * ele%value(l$) / 4
      y_ave = start2%vec(3) + direc * start2%vec(4) * ele%value(l$) / 4
      h_x = ele%value(g$) + ele%value(delta_g$) + ele%value(k1$) * x_ave
      h_y = ele%value(k1$) * y_ave
      h2 = h_x**2 + h_y**2
    endif

  case (wiggler$)
    if (ele%sub_key == map_type$) then
      h2 = ele%const(0) + &
                  dot_product(start%vec-ele%const(7:12), ele%const(1:6))
    elseif (ele%sub_key == periodic_type$) then
      h2 = ele%value(k1$)
    endif

  end select

! Apply the radiation kicks

  gamma = param%beam_energy / m_electron

  fact_d = 0
  if (sr_com%damping_on) fact_d = 2 * r_e * gamma**3 * h2 * s_len / 3

  fact_f = 0
  if (sr_com%fluctuations_on) then
    h3 = sqrt(h2)**3
    call gauss_ran (this_ran)
    fact_f = sqrt(fluct_const * s_len * gamma**5 * h3) * this_ran
  endif

  dE_p = (1 + start%vec(6)) * (fact_d + fact_f) * sr_com%scale 

  end = start
  end%vec(2) = end%vec(2) * (1 - dE_p)
  end%vec(4) = end%vec(4) * (1 - dE_p)
  end%vec(6) = end%vec(6)  - dE_p * (1 + end%vec(6))

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine setup_radiation_tracking (ring, closed_orb,
!                                        fluctuations_on, damping_on) 
!
! Subroutine to compute synchrotron radiation parameters prior to tracking.
!
! Note: The fluctuations_on and damping_on switches need to be set for 
! each ring if multiple rings are used. The default if not present is to
! not change the settings for a ring. Initially a ring is set for fluctuations
! and damping off.
!
! Modules needed:
!   use bmad
!
! Input:
!   ring            -- Ring_struct:
!   closed_orb(0:)  -- Coord_struct: Closed_orbit.
!   fluctuations_on -- Logical, optional: If True then radiation fluctuations
!                        will be present. 
!   damping_on      -- Logical, optional: If True then radiation damping
!                        will be present. 
!
! Output:
!   ring           -- Ring_struct: Lattice with radiation parameters computed. 
!-

subroutine setup_radiation_tracking (ring, closed_orb, &
                                        fluctuations_on, damping_on)

  use symp_lie_mod

  implicit none

  type (ring_struct), target :: ring
  type (coord_struct), intent(in) :: closed_orb(0:)
  type (coord_struct) start, end

  real(rp), save :: del_orb(6) = (/ 1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-5 /)

  integer i, j

  logical, optional, intent(in) :: fluctuations_on, damping_on

! Set logicals.

  if (present(fluctuations_on)) sr_com%fluctuations_on = fluctuations_on
  if (present(damping_on))      sr_com%damping_on      = damping_on

! compute wiggler parameters

  do i = 1, ring%n_ele_ring
    if (ring%ele_(i)%key /= wiggler$) cycle
    if (ring%ele_(i)%sub_key == map_type$) then
      sl_com%save_orbit = .true.
      if (.not. associated(ring%ele_(i)%const)) allocate (ring%ele_(i)%const(0:12))
      ring%ele_(i)%const(7:12) = closed_orb(i)%vec
      call symp_lie_bmad (ring%ele_(i), ring%param, closed_orb(i), end, .false.)
      call calc_h2 (ring%ele_(i)%const(0))
      do j = 1, 6
        start = closed_orb(i)
        start%vec(j) = start%vec(j) + del_orb(j)
        call symp_lie_bmad (ring%ele_(i), ring%param, start, end, .false.)
        call calc_h2 (ring%ele_(i)%const(j))
        ring%ele_(i)%const(j) = &
              (ring%ele_(i)%const(j) - ring%ele_(i)%const(0)) / del_orb(j)
      enddo
    endif
  enddo

!-------------------------------------------------------
contains

subroutine calc_h2 (h2)

  real(rp) h2, kick(6)
  integer kk

! h2 is the average over the element.

  h2 = 0
  do kk = 0, ubound(sl_com%s, 1) 
    call derivs_bmad (ring%ele_(i), ring%param, sl_com%s(kk), &
                                            sl_com%orb(kk)%vec, kick)
    if (kk == 0 .or. kk == ubound(sl_com%s, 1)) then
      h2 = h2 + (kick(2)**2 + kick(4)**2) / 2
    else
      h2 = h2 + (kick(2)**2 + kick(4)**2)
    endif
  enddo

  h2 = h2 / ubound(sl_com%s, 1) 

end subroutine

end subroutine

end module
