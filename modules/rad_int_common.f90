!+
! Module rad_int_common
!
! Module needed:
!   use rad_int_common
!-

module rad_int_common               

use ptc_interface_mod
use em_field_mod

! The "cache" is for saving values for g, etc through an element to speed
! up the calculation.

integer, parameter :: no_cache$ = 0, cache_with_misalign$ = 1, cache_no_misalign$ = 2

type rad_int_track_point_struct
  real(rp) mat6(6,6)
  real(rp) vec0(6)
  type (coord_struct) ref_orb_in
  type (coord_struct) ref_orb_out
  real(rp) g_x0, g_y0     ! Additional g factors for bends.
  real(rp) dgx_dx, dgx_dy   ! bending strength gradiant
  real(rp) dgy_dx, dgy_dy   ! bending strength gradiant
end type

type rad_int_cache1_struct
  type (rad_int_track_point_struct), allocatable :: pt(:)
  !! type (ele_struct), pointer :: ele
  real(rp) del_z
  integer :: cache_type = no_cache$
end type

type rad_int_cache_struct
  type (rad_int_cache1_struct), allocatable :: c_ele(:)
  logical set
end type

type (rad_int_cache_struct), target, save :: rad_int_cache_common(0:10)

!

type rad_int_info_struct
  type (branch_struct), pointer :: branch
  type (ele_struct), pointer :: ele
  type (coord_struct), pointer :: orbit(:)
  type (twiss_struct)  a, b
  type (rad_int_cache1_struct), pointer :: cache_ele ! pointer to cache in use
  real(rp) eta_a(4), eta_b(4)
  real(rp) g, g2          ! bending strength (1/bending_radius)
  real(rp) g_x, g_y       ! components in x-y plane
  real(rp) dg2_x, dg2_y
end type

contains

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!+
! Subroutine qromb_rad_int(param, do_int, pt, info, int_tot, rad_int)
!
! Function to do integration using Romberg's method on the 7 radiation 
! integrals.
! This is a modified version of QROMB from Num. Rec.
! See the Num. Rec. book for further details.
!
! This routine is only meant to be called by radiation_integrals and
! is not meant for general use.
!
! There are up to 9 integrals that are calculated:
!          I1, I2, I3, I4a, I4b, I5a, I5b, i0, i6b
! If do_int(1:7) is False for an integral that means that the integral has
! been calculated by the calling routine using a formula and therefore does
! not have to be done by this routine.
!-

subroutine qromb_rad_int (param, do_int, pt, info, int_tot, rad_int1)

use precision_def
use nrtype
use nr, only: polint

implicit none

integer, parameter :: num_int = 9
type ri_array_struct
  real(rp) h(num_int)
  real(rp) sum(num_int)
end type

type (ri_array_struct) ri_array(0:4) ! ri_array(n) holds info for the integrals for a particular step.
type (ele_struct), pointer :: ele
type (coord_struct) start, end
type (rad_int_track_point_struct) pt
type (rad_int_info_struct) info
type (rad_int1_struct) rad_int1, int_tot
type (lat_param_struct) param

integer n,  j, j1, j_min_test, j_max, n_pts

real(rp) :: eps_int, eps_sum, gamma
real(rp) :: ll, del_z, l_ref, z_pos, dint, d0, d_max
real(rp) i_sum(num_int), rad_int_vec(num_int), int_tot_vec(num_int)

logical do_int(num_int), converged

character(*), parameter :: r_name = 'qromb_rad_int'

!

ele => info%ele

eps_int = 1d-4
eps_sum = 1d-6

ri_array(0)%h = 4
ri_array(0)%sum = 0
rad_int_vec = 0
int_tot_vec = [int_tot%i1, int_tot%i2, int_tot%i3, int_tot%i4a, &
               int_tot%i4b, int_tot%i5a, int_tot%i5b, int_tot%i0,int_tot%i6b]

ll = ele%value(l$)

start = info%orbit(ele%ix_ele-1)
end   = info%orbit(ele%ix_ele)
gamma = ele%value(e_tot$) / mass_of(param%particle)

! Go to the local element frame if there has been caching.
if (associated(info%cache_ele)) then
  call offset_particle (ele, param, set$, start, set_multipoles = .false., set_hvkicks = .false.)
  call offset_particle (ele, param, set$, end, set_multipoles = .false., set_hvkicks = .false., ds_pos = ll)
endif

  ! For j >= 3 we test if the integral calculation has converged.
  ! Exception: Since wigglers have a periodic field, the calculation can 
  ! fool itself if we stop before j = 5.

j_min_test = 3
j_max = 14

if (ele%key == wiggler$ .or. ele%key == undulator$) then
  if (ele%sub_key == periodic_type$) then
    j_min_test = 3 + log(max(1.0_rp, ele%value(n_pole$))) / log(2.0_rp)
    j_max = j_min_test + 8
  else
    j_min_test = 5
    j_max = 16
  endif
endif

! Loop until integrals converge.

do j = 1, j_max

  !---------------
  ! This is trapzd from Numerical Recipes

  if (j == 1) then
    n_pts = 2
    del_z = ll
    l_ref = 0
  else
    n_pts = 2**(j-2)
    del_z = ll / n_pts
    l_ref = del_z / 2
  endif

  i_sum = 0

  do n = 1, n_pts
    z_pos = l_ref + (n-1) * del_z
    call propagate_part_way (start, param, pt, info, z_pos)
    i_sum(1) = i_sum(1) + info%g_x * (info%eta_a(1) + info%eta_b(1)) + &
                  info%g_y * (info%eta_a(3) + info%eta_b(3))
    i_sum(2) = i_sum(2) + info%g2
    i_sum(3) = i_sum(3) + info%g2 * info%g
    i_sum(4) = i_sum(4) + &
                  info%g2 * (info%g_x * info%eta_a(1) + info%g_y * info%eta_a(3)) + &
                  info%dg2_x * info%eta_a(1) + info%dg2_y * info%eta_a(3)
    i_sum(5) = i_sum(5) + &
                  info%g2 * (info%g_x * info%eta_b(1) + info%g_y * info%eta_b(3)) + &
                  info%dg2_x * info%eta_b(1) + info%dg2_y * info%eta_b(3)
    i_sum(6) = i_sum(6) + &
                  info%g2 * info%g * (info%a%gamma * info%a%eta**2 + &
                  2 * info%a%alpha * info%a%eta * info%a%etap + &
                  info%a%beta * info%a%etap**2)
    i_sum(7) = i_sum(7) + &
                  info%g2 * info%g * (info%b%gamma * info%b%eta**2 + &
                  2 * info%b%alpha * info%b%eta * info%b%etap + &
                  info%b%beta * info%b%etap**2)
    i_sum(8) = i_sum(8) + gamma * info%g
    i_sum(9) = i_sum(9) + info%g2 * info%g * info%b%beta
  enddo

  if (j <= 4) then
    j1 = j
  else
    ri_array(0:3) = ri_array(1:4)
    j1 = 4
  endif

  ri_array(j1)%h = ri_array(j1-1)%h / 4
  ri_array(j1)%sum = (ri_array(j1-1)%sum + del_z * i_sum) / 2

  !--------------
  ! Back to qromb.

  if (j < j_min_test) cycle

  converged = .true.
  d_max = 0

  do n = 1, num_int
    if (.not. do_int(n)) cycle
    call polint (ri_array(1:j1)%h(n), ri_array(1:j1)%sum(n), 0.0_rp, rad_int_vec(n), dint)
    d0 = eps_int * abs(rad_int_vec(n)) + eps_sum * abs(int_tot_vec(n))
    if (abs(dint) > d0)  converged = .false.
    if (d0 /= 0) d_max = max(d_max, abs(dint) / d0)
  enddo

  ! If we have convergance or we are giving up (when j = j_max) then 
  ! stuff the results in the proper places.

  if (converged .or. j == j_max) then

    rad_int1%n_steps = j

    ! Note that rad_int%i... may already contain a contribution from edge
    ! affects (Eg bend face angles) so add it on to rad_int_vec(i)

    rad_int1%i1  = rad_int1%i1  + rad_int_vec(1)
    rad_int1%i2  = rad_int1%i2  + rad_int_vec(2)
    rad_int1%i3  = rad_int1%i3  + rad_int_vec(3)
    rad_int1%i4a = rad_int1%i4a + rad_int_vec(4)
    rad_int1%i4b = rad_int1%i4b + rad_int_vec(5)
    rad_int1%i5a = rad_int1%i5a + rad_int_vec(6)
    rad_int1%i5b = rad_int1%i5b + rad_int_vec(7)
    rad_int1%i0  = rad_int1%i0  + rad_int_vec(8)
    rad_int1%i6b = rad_int1%i6b + rad_int_vec(9)

    int_tot%i1  = int_tot_vec(1) + rad_int1%i1
    int_tot%i2  = int_tot_vec(2) + rad_int1%i2
    int_tot%i3  = int_tot_vec(3) + rad_int1%i3
    int_tot%i4a = int_tot_vec(4) + rad_int1%i4a
    int_tot%i4b = int_tot_vec(5) + rad_int1%i4b
    int_tot%i5a = int_tot_vec(6) + rad_int1%i5a
    int_tot%i5b = int_tot_vec(7) + rad_int1%i5b
    int_tot%i0  = int_tot_vec(8) + rad_int1%i0
    int_tot%i6b = int_tot_vec(9) + rad_int1%i6b

  endif

  if (converged) return

end do

! We should not be here

call out_io (s_warn$, r_name, 'Note: Radiation Integral is not converging \es12.3\ ', &
                              'For element: ' // ele%name, r_array = [d_max])

end subroutine qromb_rad_int

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

subroutine propagate_part_way (orb_start, param, pt, info, z_here)

implicit none

type (coord_struct) orb, orb_start, orb0, orb1, orb_end, orb_end1
type (ele_struct), pointer :: ele0, ele
type (ele_struct), save :: runt, ele_end
type (twiss_struct) a0, b0, a1, b1
type (rad_int_info_struct) info
type (rad_int_track_point_struct) pt, pt0, pt1
type (lat_param_struct) param

real(rp) z_here, v(4,4), v_inv(4,4), s1, s2, error
real(rp) f0, f1, del_z, c, s, x, y
real(rp) eta_a0(4), eta_a1(4), eta_b0(4), eta_b1(4)
real(rp) dz, z1, tilt
real(rp) a_pole(0:n_pole_maxx), b_pole(0:n_pole_maxx), dk(2,2), kx, ky

integer i0, i1, tm_saved, m6cm_saved
integer i, ix, ip, j_loop, n_pt, n, n1, n2
integer, save :: ix_ele = -1

logical has_nonzero_pole

! Init

ele0 => info%branch%ele(info%ele%ix_ele-1)
ele  => info%ele

if (ix_ele /= info%ele%ix_ele) then
  runt = ele
  ix_ele = info%ele%ix_ele
endif

!--------------------------------------
! With caching
! Remember that here the start and stop coords, etc. are in the local ref frame.

if (associated(info%cache_ele)) then

  ! find cached point info near present z position

  del_z = info%cache_ele%del_z
  i0 = int(z_here/del_z)
  f1 = (z_here - del_z*i0) / del_z 
  f0 = 1 - f1
  i1 = i0 + 1
  if (i1 > ubound(info%cache_ele%pt, 1)) i1 = i0  ! can happen with roundoff
  pt0 = info%cache_ele%pt(i0)
  pt1 = info%cache_ele%pt(i1)

  orb0%vec = matmul(pt0%mat6, orb_start%vec) + pt0%vec0
  orb1%vec = matmul(pt1%mat6, orb_start%vec) + pt1%vec0

  ! Interpolate information

  if ((ele%key == wiggler$ .or. ele%key == undulator$) .and. ele%sub_key == map_type$) then 
    orb = orb_start
    orb%vec = f0 * orb0%vec + f1 * orb1%vec
    call calc_wiggler_g_params (ele, z_here, orb, pt, info)

  else
    pt%dgx_dx = f0 * pt0%dgx_dx + f1 * pt1%dgx_dx
    pt%dgx_dy = f0 * pt0%dgx_dy + f1 * pt1%dgx_dy
    pt%dgy_dx = f0 * pt0%dgy_dx + f1 * pt1%dgy_dx
    pt%dgy_dy = f0 * pt0%dgy_dy + f1 * pt1%dgy_dy

    info%g_x   = f0 * (pt0%g_x0 + pt0%dgx_dx * (orb0%vec(1) - pt0%ref_orb_out%vec(1)) + &
                                  pt0%dgx_dy * (orb0%vec(3) - pt0%ref_orb_out%vec(3))) + &
                 f1 * (pt1%g_x0 + pt1%dgx_dx * (orb1%vec(1) - pt1%ref_orb_out%vec(1)) + &
                                  pt1%dgx_dy * (orb1%vec(3) - pt1%ref_orb_out%vec(3)))
    info%g_y   = f0 * (pt0%g_y0 + pt0%dgy_dx * (orb0%vec(1) - pt0%ref_orb_out%vec(1)) + &
                                  pt0%dgy_dy * (orb0%vec(3) - pt0%ref_orb_out%vec(3))) + &
                 f1 * (pt1%g_y0 + pt1%dgy_dx * (orb1%vec(1) - pt1%ref_orb_out%vec(1)) + &
                                  pt1%dgy_dy * (orb1%vec(3) - pt1%ref_orb_out%vec(3)))
  endif
                 
  info%dg2_x = 2 * (info%g_x * pt%dgx_dx + info%g_y * pt%dgy_dx)
  info%dg2_y = 2 * (info%g_x * pt%dgx_dy + info%g_y * pt%dgy_dy) 
  info%g2 = info%g_x**2 + info%g_y**2
  info%g  = sqrt(info%g2)

  ! Now convert the g calc back to lab coords.
  
  if (ele%key == sbend$) then
    tilt = ele%value(ref_tilt_tot$) + ele%value(roll$)  ! Small angle approx.
  else
    tilt = ele%value(tilt_tot$) 
  endif

  if (tilt /= 0) then
    c = cos(tilt); s = sin(tilt) 
    x = info%g_x; y = info%g_y
    info%g_x = c * x - s * y
    info%g_y = s * x + c * y
    x = info%dg2_x; y = info%dg2_y
    info%dg2_x = c * x - s * y
    info%dg2_y = s * x + c * y
  endif

  ! Interpolate the dispersions

  runt%mat6 = pt0%mat6
  runt%vec0 = pt0%vec0
  runt%map_ref_orb_in  = pt0%ref_orb_in
  runt%map_ref_orb_out = pt0%ref_orb_out

  call mat6_add_offsets (runt, param)  ! back to lab coords
  call twiss_propagate1 (ele0, runt)
  a0 = runt%a; b0 = runt%b
  call make_v_mats (runt, v, v_inv)
  eta_a0 = matmul(v, [runt%a%eta, runt%a%etap, 0.0_rp,   0.0_rp    ])
  eta_b0 = matmul(v, [0.0_rp,   0.0_rp,    runt%b%eta, runt%b%etap ])

  runt%mat6 = pt1%mat6
  runt%vec0 = pt1%vec0
  runt%map_ref_orb_in  = pt1%ref_orb_in
  runt%map_ref_orb_out = pt1%ref_orb_out

  call mat6_add_offsets (runt, param)  ! back to lab coords
  call twiss_propagate1 (ele0, runt)
  a1 = runt%a; b1 = runt%b
  call make_v_mats (runt, v, v_inv)
  eta_a1 = matmul(v, [runt%a%eta, runt%a%etap, 0.0_rp,   0.0_rp    ])
  eta_b1 = matmul(v, [0.0_rp,   0.0_rp,    runt%b%eta, runt%b%etap ])

  info%a%beta  = a0%beta  * f0 + a1%beta  * f1
  info%a%alpha = a0%alpha * f0 + a1%alpha * f1
  info%a%gamma = a0%gamma * f0 + a1%gamma * f1
  info%a%eta   = a0%eta   * f0 + a1%eta   * f1
  info%a%etap  = a0%etap  * f0 + a1%etap  * f1

  info%b%beta  = b0%beta  * f0 + b1%beta  * f1
  info%b%alpha = b0%alpha * f0 + b1%alpha * f1
  info%b%gamma = b0%gamma * f0 + b1%gamma * f1
  info%b%eta   = b0%eta   * f0 + b1%eta   * f1
  info%b%etap  = b0%etap  * f0 + b1%etap  * f1

  info%eta_a = eta_a0 * f0 + eta_a1 * f1
  info%eta_b = eta_b0 * f0 + eta_b1 * f1

  return

endif

!--------------------------------------
! No caching...
! Note: calc_wiggler_g_params assumes that orb is lab (not element) coords.

dz = 1d-3
z1 = z_here + dz
if (z1 > ele%value(l$)) then
  z_here = max(0.0_rp, z_here - dz)
  z1 = min(ele%value(l$), z_here + dz)
endif

! bmad_standard will not properly do partial tracking through a periodic_type wiggler so
! switch to symp_lie_bmad type tracking.

if ((ele%key == wiggler$ .or. ele%key == undulator$) .and. &
                    ele%sub_key == periodic_type$ .and. ele%tracking_method /= custom$) then
  tm_saved = ele%tracking_method  
  m6cm_saved = ele%mat6_calc_method  
  ele%tracking_method = symp_lie_bmad$
  ele%mat6_calc_method = symp_lie_bmad$
endif

call twiss_and_track_intra_ele (ele, info%branch%param, 0.0_rp, z_here, .true., .false., orb_start, orb_end, ele0, ele_end)

info%a = ele_end%a
info%b = ele_end%b

if ((ele%key == wiggler$ .or. ele%key == undulator$) .and. ele%sub_key == periodic_type$) then
  ele%tracking_method  = tm_saved 
  ele%mat6_calc_method = m6cm_saved 
endif

!

call make_v_mats (ele_end, v, v_inv)

info%eta_a = matmul(v, [info%a%eta, info%a%etap, 0.0_rp,   0.0_rp ])
info%eta_b = matmul(v, [0.0_rp,   0.0_rp,    info%b%eta, info%b%etap ])

if ((ele%key == wiggler$ .or. ele%key == undulator$) .and. ele%sub_key == map_type$) then 
  call calc_wiggler_g_params (ele, z_here, orb_end, pt, info)
else
  call twiss_and_track_intra_ele (ele, info%branch%param, z_here, z1,     .false., .false., orb_end, orb_end1)
  info%g_x = pt%g_x0 - (orb_end1%vec(2) - orb_end%vec(2)) / (z1 - z_here)
  info%g_y = pt%g_y0 - (orb_end1%vec(4) - orb_end%vec(4)) / (z1 - z_here)
endif

info%dg2_x = 2 * (info%g_x * pt%dgx_dx + info%g_y * pt%dgy_dx)
info%dg2_y = 2 * (info%g_x * pt%dgx_dy + info%g_y * pt%dgy_dy) 

info%g2 = info%g_x**2 + info%g_y**2
info%g = sqrt(info%g2)

! Add in multipole gradient

if ((ele%key /= wiggler$ .and. ele%key /= undulator$) .or. ele%sub_key /= map_type$) then 
  call multipole_ele_to_ab (ele, .true., has_nonzero_pole, a_pole, b_pole)
  if (has_nonzero_pole) then
    do ip = 0, ubound(a_pole, 1)
      if (a_pole(ip) == 0 .and. b_pole(ip) == 0) cycle
      call ab_multipole_kick (a_pole(ip), b_pole(ip), ip, orb_end, kx, ky, dk)
      info%dg2_x = info%dg2_x - 2 * (info%g_x * dk(1,1) + info%g_y * dk(2,1)) / ele%value(l$)
      info%dg2_y = info%dg2_y - 2 * (info%g_x * dk(1,2) + info%g_y * dk(2,2)) / ele%value(l$)
    enddo
  endif
endif

end subroutine propagate_part_way

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine calc_wiggler_g_params (ele, s_rel, orb, pt, info)

implicit none

type (coord_struct) orb
type (rad_int_track_point_struct) pt
type (rad_int_info_struct) info
type (ele_struct) ele

real(rp) s_rel, g(3), dg(3,3)

! Note: em_field_g_bend assumes orb is lab (not element) coords.

call em_field_g_bend (ele, info%branch%param, s_rel, 0.0_rp, orb, g, dg)

pt%g_x0 = g(1)
pt%g_y0 = g(2)
pt%dgx_dx = dg(1,1)
pt%dgx_dy = dg(1,2)
pt%dgy_dx = dg(2,1)
pt%dgy_dy = dg(2,2)

info%g_x = g(1)
info%g_y = g(2)

end subroutine calc_wiggler_g_params

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine em_field_g_bend (ele, param, s_, orbit, g, dg)
!
! Subroutine to calculate the g bending kick felt by a particle in a element. 
!
! Modules needed:
!   use bmad
!
! Input:
!   ele    -- Ele_struct: Element being tracked thorugh.
!   param  -- lat_param_struct: Lattice parameters.
!   s_rel  -- Real(rp): Distance from the start of the element to the particle.
!   t_rel  -- Real(rp): Time relative to the reference particle.
!   orbit  -- Coord_struct: Particle position in lab (not element) frame.
!
! Output:
!   g(3)    -- Real(rp): (g_x, g_y, g_s) bending radiuses
!   dg(3,3) -- Real(rp), optional: dg(:)/dr gradient. 
!-

subroutine em_field_g_bend (ele, param, s_rel, t_rel, orbit, g, dg)

implicit none

type (ele_struct) ele
type (lat_param_struct) param
type (em_field_struct) field
type (coord_struct) orbit

real(rp), intent(in) :: s_rel, t_rel
real(rp), intent(out) :: g(3)
real(rp), optional :: dg(3,3)
real(rp) vel_unit(3), fact
real(rp) f

! calculate the field

call em_field_calc (ele, param, s_rel, t_rel, orbit, .false., field, present(dg))

!

! vel_unit is the velocity normalized to unit length

vel_unit(1:2) = [orbit%vec(2), orbit%vec(4)] / (1 + orbit%vec(6))
vel_unit(3) = sqrt(1 - vel_unit(1)**2 - vel_unit(2)**2)
fact = 1 / (ele%value(p0c$) * (1 + orbit%vec(6)))
g = g_from_field (field%B, field%E)

! Derivative

if (present(dg)) then
  dg(:,1) = g_from_field (field%dB(:,1), field%dE(:,1))
  dg(:,2) = g_from_field (field%dB(:,2), field%dE(:,2))
  dg(:,3) = g_from_field (field%dB(:,3), field%dE(:,3))
endif

!---------------------------------------------------------------
contains

function g_from_field (B, E) result (g_bend)

real(rp) B(3), E(3), g_bend(3)
real(rp) force(3), force_perp(3)

! force_perp is the perpendicular component of the force.

force = (E + cross_product(vel_unit, B) * orbit%beta * c_light) * charge_of(param%particle)
force_perp = force - vel_unit * (dot_product(force, vel_unit))
g_bend = -force_perp * fact

end function g_from_field

end subroutine em_field_g_bend

end module
