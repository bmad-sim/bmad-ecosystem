module radiation_mod

use bmad_struct
use bmad_interface
use runge_kutta_mod
use rad_int_common

contains

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!+
! Subroutine release_rad_int_cache (ix_cache)
!
! Subroutine to release the memory associated with caching wiggler values.
! See the radiation_integrals routine for further details.
!
! Modules needed:
!   use radiation_mod
!
! Input:
!   ix_cache -- Integer: Cache number.
!
! Output:
!   ix_cache -- Integer: Cache number set to 0,
!-

subroutine release_rad_int_cache (ix_cache)

implicit none

integer i, ix_cache

!

rad_int_cache_common(ix_cache)%set = .false.
ix_cache = 0

end subroutine release_rad_int_cache 

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine track1_radiation (orb_start, ele, param, orb_end, edge)
!
! Subroutine to put in radiation dampling and/or fluctuations.
! This routine calculates half the radiation of an element so this routine
! needs to be called before entering an element and just after exiting.
!
! Modules needed:
!   use radiation_mod
!
! Input:
!   orb_start -- coord_struct: Particle position before radiation applied.
!   ele       -- ele_struct: Element that causes radiation.
!   edge      -- integer: Where the particle is: start_edge$ or end_edge$.
!
! Output:
!   orb_end   -- coord_struct: Particle position after radiation has been applied.
!-

subroutine track1_radiation (orb_start, ele, param, orb_end, edge)

use random_mod

implicit none

type (coord_struct) :: orb_start, orb_end
type (ele_struct), target :: ele
type (ele_struct), pointer :: ele0
type (lat_param_struct) :: param
type (coord_struct), save :: orb_start2

integer :: edge

real(rp), save :: z_start, g2, g3
real(rp) s_len, g, g_x, g_y, this_ran, mc2
real(rp) x_ave, y_ave, gamma_0, dE_p, fact_d, fact_f
real(rp), parameter :: rad_fluct_const = 55 * classical_radius_factor * &
                                                h_bar_planck * c_light / (24 * sqrt_3)

integer direc

logical set
logical :: init_needed = .true.

character(20) :: r_name = 'track1_radiation'

! If not a magnetic element then nothing to do.
! Also symplectic tracking handles the radiation.

if (ele%tracking_method == symp_lie_bmad$ .or. .not. any (ele%key == &
          [quadrupole$, sextupole$, octupole$, sbend$, sol_quad$, wiggler$, undulator$])) then
  orb_end = orb_start
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
  z_start = orb_start%vec(5)
elseif (edge == end_edge$) then
  set = unset$
  s_len = ele%value(l$)/2 + (orb_start%vec(5) - z_start)
  direc = -1
else
  call out_io (s_fatal$, r_name, 'BAD EDGE ARGUMENT:', set)
  if (global_com%exit_on_error) call err_exit
endif

if (s_len < 0) s_len = 0

! Get the coords in the frame of reference of the element

if (ele%key /= wiggler$ .and. ele%key /= undulator$) then
  orb_start2 = orb_start
  call offset_particle (ele, param, set, orb_start2)
  call canonical_to_angle_coords (orb_start2)
endif

! Calculate the radius of curvature for an on-energy particle

select case (ele%key)

case (quadrupole$, sol_quad$)
  x_ave = orb_start2%vec(1) + direc * orb_start2%vec(2) * ele%value(l$) / 4
  y_ave = orb_start2%vec(3) + direc * orb_start2%vec(4) * ele%value(l$) / 4
  g_x =  ele%value(k1$) * x_ave
  g_y = -ele%value(k1$) * y_ave
  g2 = g_x**2 + g_y**2
  if (bmad_com%radiation_fluctuations_on) g3 = sqrt(g2)**3

case (sextupole$)
  g = ele%value(k2$) * (orb_start2%vec(1)**2 + orb_start2%vec(3)**2)
  g2 = g**2
  if (bmad_com%radiation_fluctuations_on) g3 = g2 * abs(g)

case (octupole$)
  g2 = ele%value(k3$)**2 * (orb_start2%vec(1)**2 + orb_start2%vec(3)**2)**3
  if (bmad_com%radiation_fluctuations_on) g3 = sqrt(g2)**3

case (sbend$)
  if (ele%value(k1$) == 0) then
    g = ele%value(g$) + ele%value(g_err$)
    g2 = g**2 
    if (bmad_com%radiation_fluctuations_on) g3 = g2 * abs(g)
  else
    x_ave = orb_start2%vec(1) + direc * orb_start2%vec(2) * ele%value(l$) / 4
    y_ave = orb_start2%vec(3) + direc * orb_start2%vec(4) * ele%value(l$) / 4
    g_x = ele%value(g$) + ele%value(g_err$) + ele%value(k1$) * x_ave
    g_y = ele%value(k1$) * y_ave
    g2 = g_x**2 + g_y**2
    if (bmad_com%radiation_fluctuations_on) g3 = sqrt(g2)**3
  endif

case (wiggler$, undulator$)
  ! Reuse g2 and g3 values from start_edge
  if (ele%sub_key == map_type$) then
    if (edge == start_edge$) then
      if (.not. associated(ele%rad_int_cache) .or. ele%rad_int_cache%stale) call calc_radiation_tracking_consts(ele, param)
      g2 = ele%rad_int_cache%g2_0 + dot_product(orb_start%vec(1:4)-ele%rad_int_cache%orb0(1:4), ele%rad_int_cache%dg2_dorb(1:4))
      g3 = ele%rad_int_cache%g3_0 + dot_product(orb_start%vec(1:4)-ele%rad_int_cache%orb0(1:4), ele%rad_int_cache%dg3_dorb(1:4))
      if (g3 < 0) g3 = 0
    endif
  elseif (ele%sub_key == periodic_type$) then
    g2 = abs(ele%value(k1$))
    g3 = 4 * sqrt(2*g2)**3 / (3 * pi)  
  endif

end select

! Apply the radiation kicks
! Basic equation is E_radiated = xi * (dE/dt) * sqrt(L) / c_light
! where xi is a random number with sigma = 1.

mc2 = mass_of(param%particle)
gamma_0 = ele%value(e_tot$) / mc2

fact_d = 0
if (bmad_com%radiation_damping_on) then
  fact_d = 2 * classical_radius_factor * gamma_0**3 * g2 * s_len / (3 * mc2)
  if (param%backwards_time_tracking) fact_d = -fact_d
endif

fact_f = 0
if (bmad_com%radiation_fluctuations_on) then
  call ran_gauss (this_ran)
  fact_f = sqrt(rad_fluct_const * s_len * gamma_0**5 * g3) * this_ran / mc2
endif

dE_p = (1 + orb_start%vec(6)) * (fact_d + fact_f) * synch_rad_com%scale 

orb_end = orb_start
orb_end%vec(2) = orb_end%vec(2) * (1 - dE_p)
orb_end%vec(4) = orb_end%vec(4) * (1 - dE_p)
orb_end%vec(6) = orb_end%vec(6)  - dE_p * (1 + orb_end%vec(6))

if (synch_rad_com%i_calc_on) then
  synch_rad_com%i2 = synch_rad_com%i2 + g2 * s_len
  synch_rad_com%i3 = synch_rad_com%i3 + g3 * s_len
  if (associated(ele%branch)) then
    if (edge == start_edge$) then
      ele0 => ele%branch%ele(ele%ix_ele-1)
    else
      ele0 => ele
    endif
    synch_rad_com%i5a = synch_rad_com%i5a + g3 * s_len * (ele0%a%gamma * ele0%a%eta**2 + &
                  2 * ele0%a%alpha * ele0%a%eta * ele0%a%etap + ele0%a%beta * ele0%a%etap**2)
    synch_rad_com%i5b = synch_rad_com%i5b + g3 * s_len * (ele0%b%gamma * ele0%b%eta**2 + &
                  2 * ele0%b%alpha * ele0%b%eta * ele0%b%etap + ele0%b%beta * ele0%b%etap**2)
  endif
endif

end subroutine track1_radiation 

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine calc_radiation_tracking_consts (ele, param)
!
! Routine to compute synchrotron radiation parameters prior to tracking.
! This routine is not meant for general use
!
! Modules needed:
!   use radiation_mod
!
! Input:
!   ele   -- ele_struct: Element 
!   param -- lat_param_struct
!
! Output:
!   ele  -- ele_struct: Element 
!     %rad_int_cache  -- radiation tracking constants.
!-

subroutine calc_radiation_tracking_consts (ele, param)

use symp_lie_mod

implicit none

type (ele_struct) ele
type (lat_param_struct) param
type (coord_struct) start0_orb, start_orb, end_orb
type (track_struct), save :: track

real(rp), parameter :: del_orb = 1e-4
real(rp) g2, g3

integer i, j

! Compute for a map_type wiggler the change in g2 and g3 
!   with respect to transverse orbit for an on-energy particle.

if (ele%key /= wiggler$ .and. ele%key /= undulator$) return
if (ele%sub_key /= map_type$) return

if (.not. associated(ele%rad_int_cache)) allocate (ele%rad_int_cache)
ele%rad_int_cache%orb0 = ele%map_ref_orb_in%vec

track%n_pt = -1
call symp_lie_bmad (ele, param, ele%map_ref_orb_in, end_orb, .false., track)
call calc_g (track, ele%rad_int_cache%g2_0, ele%rad_int_cache%g3_0)

do j = 1, 4
  start_orb = ele%map_ref_orb_in
  start_orb%vec(j) = start_orb%vec(j) + del_orb
  track%n_pt = -1
  call symp_lie_bmad (ele, param, start_orb, end_orb, .false., track)
  call calc_g (track, g2, g3)
  ele%rad_int_cache%dg2_dorb(j) = (g2 - ele%rad_int_cache%g2_0) / del_orb
  ele%rad_int_cache%dg3_dorb(j) = (g3 - ele%rad_int_cache%g3_0) / del_orb
enddo

ele%rad_int_cache%stale = .false.

!-------------------------------------------------------
contains

subroutine calc_g (track, g2, g3)

type (track_struct) track
real(rp) g2, g3, g2_here, g3_here, g(3), f, s0
integer j, n0, n1

! g2 is the average g^2 over the element for an on-energy particle.
! Note: em_field_g_bend assumes orb is lab (not local) coords.

track%orb(:)%vec(6) = 0  ! on-energy

g2 = 0; g3 = 0

n0 = lbound(track%orb, 1)
n1 = track%n_pt
s0 = ele%s + ele%value(z_offset_tot$) - ele%value(l$)

do j = n0, n1

  call em_field_g_bend (ele, param, track%orb(j)%s - s0, 0.0_rp, track%orb(j), g)

  g2_here = g(1)**2 + g(2)**2 ! = g_x^2 + g_y^2
  g3_here = sqrt(g2_here)**3

  if (i == n0 .or. i == n1) then
    g2_here = g2_here / 2
    g3_here = g3_here / 2
  endif

  g2 = g2 + g2_here
  g3 = g3 + g3_here

enddo

g2 = g2 / (n1 - n0 + 1)
g3 = g3 / (n1 - n0 + 1)

end subroutine

end subroutine calc_radiation_tracking_consts

end module
