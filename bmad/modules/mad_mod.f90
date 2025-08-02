!+
! Module mad_mod                        
!
! Module containing the routines to calculate 2nd order transport map.
! These routines are adapted from the MAD program.
!
! Note: These routines calculate the trasport map around the zero-orbit.
!-

module mad_mod

use taylor_mod

type mad_energy_struct
  real(rp) total 
  real(rp) beta         ! normalized velocity: v/c
  real(rp) gamma        ! relativistic factor: 1/sqrt(1-beta^2)
  real(rp) kinetic      ! kinetic energy
  real(rp) p0c          ! particle momentum
  integer particle      ! particle species
end type

type mad_map_struct
  real(rp) k(6)         ! 0th order map.
  real(rp) r(6,6)       ! 1st order map.
  real(rp) t(6,6,6)     ! 2nd order map.
end type

!---------------------------------------------------------------------------
contains

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine make_mat6_mad (ele, param, c0, c1)
!
! Subroutine to make the 6x6 transfer matrix for an element from the 
! 2nd order MAD transport map. The map is stored in ele%taylor.
! If the map exists then it is simply used to calculate ele%mat6. 
! If ele%taylor doesn't exist then calculate it.
!
! Input:
!   ele    -- Ele_struct: Element with transfer matrix.
!   param  -- lat_param_struct: Lattice parameters.
!   map    -- mad_map_struct: 2nd order map.
!   c0     -- Coord_struct: Coordinates at the beginning of element.
!
! Output:
!   ele    -- Ele_struct: Element with transfer matrix.
!     %c0(6)      -- 0th order transfer matrix.
!     %mat6(6,6)  -- 6x6 1st order transfer matrix.
!   c1     -- Coord_struct: Coordinates at the end of element.
!-

subroutine make_mat6_mad (ele, param, c0, c1) 

implicit none

type (ele_struct) ele
type (lat_param_struct) param
type (mad_map_struct) map
type (mad_energy_struct) energy
type (coord_struct) c0, c1

! If ele%taylor does not exist then make it.

if (.not. associated(ele%taylor(1)%term)) then
  call make_mad_map (ele, param, energy, map)
  call mad_map_to_taylor (map, energy, ele%taylor)
  ele%taylor_map_includes_offsets = .true.
endif

! make the trasfer map.

call taylor_to_mat6 (ele%taylor, c0%vec, ele%vec0, ele%mat6)
c1%vec = track_taylor (c0%vec, ele%taylor)

end subroutine make_mat6_mad

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine make_mad_map (ele, param, energy, map)
!
! Subroutine to make a 2nd order transport map a la MAD.
!
! Input:
!   ele      -- Ele_struct: Element
!   param    -- lat_param_struct: particle id
!
! Output:
!   energy -- mad_energy_struct: Energy of the particle
!   map    -- Mad_map_struct: Structure holding the transfer map.
!-

subroutine make_mad_map (ele, param, energy, map) 

implicit none

type (ele_struct), target ::  ele
type (mad_energy_struct) energy
type (mad_map_struct) map
type (lat_param_struct) param

real(rp), pointer :: val(:)

character(16), parameter :: r_name = 'make_mad_map'

! energy structure

energy%total = ele%value(E_TOT$)
energy%particle = param%particle
call convert_total_energy_to (energy%total, energy%particle, energy%gamma, energy%kinetic, energy%beta, energy%p0c)

! choose element key

select case (ele%key)

case (drift$, pipe$, ecollimator$, rcollimator$, instrument$, monitor$)
  call mad_drift (ele, energy, map)

case (sbend$)
  call mad_sbend (ele, energy, map)

case (elseparator$)
  call mad_elsep (ele, energy, map)

case (sextupole$)
  call mad_sextupole (ele, energy, map)

case (quadrupole$)
  call mad_quadrupole (ele, energy, map)

case (rfcavity$)
  call mad_rfcavity (ele, energy, map)

case (solenoid$)
  call mad_solenoid (ele, energy, map)

case (marker$, detector$, fixer$, fiducial$)
  call make_unit_mad_map (map)

case default

  call out_io (s_fatal$, r_name, 'ELEMENT NOT IMPLEMENTED: ' // ele%name)
  if (global_com%exit_on_error) call err_exit

end select

! offsets and multipoles

!call mad_add_offsets_and_multipoles (ele, map)

val => ele%value

!if (.not.(val(z_offset_tot$) == 0 .and. val(x_offset_tot$) == 0 .and. &
!    val(x_pitch_tot$) == 0  .and. val(y_offset_tot$) == 0 .and. &
!    val(y_pitch_tot$) == 0 .and. .not. associated(ele%a_pole) .and. &
!    ((val(hkick$) == 0 .and. val(vkick$) == 0) .or. &
!    ele%key == elseparator$)) .and. mad_print_if_misaligned) &
!    call out_io (s_error$, r_name, &
!    'ELEMENT HAS OFFSET, PITCH, OR KICK. MAD DOES NOT SUPPORT THIS. ' // ele%name)

end subroutine make_mad_map

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine mad_add_offsets_and_multipoles (ele, map)
!
! Subroutine to add in the effect of element offsets and/or multipoles
! on the 2nd order transport map for the element.
!
! Input:
!   ele    -- Ele_struct: Drift element.
!   energy -- Mad_energy_struct: particle energy structure.
!
! Output:
!   map -- mad_map_struct: Structure holding the transfer map.
!-

subroutine mad_add_offsets_and_multipoles (ele, map)

implicit none

type (ele_struct), target :: ele
type (mad_map_struct) map
type (mad_map_struct) map2

real(rp), pointer :: val(:)
real(rp) s_here

character(40), parameter :: r_name = 'mad_add_offsets_and_multipoles'

! Setup.

val => ele%value

if (val(z_offset_tot$) == 0 .and. val(x_offset_tot$) == 0 .and. &
    val(x_pitch_tot$) == 0  .and. val(y_offset_tot$) == 0 .and. &
    val(y_pitch_tot$) == 0 .and. .not. associated(ele%a_pole) .and. &
    ((val(hkick$) == 0 .and. val(vkick$) == 0) .or. ele%key == elseparator$)) return

! Front side: Unit map.

call make_unit_mad_map (map2)

! Front side: z_offset

if (val(z_offset_tot$) /= 0) then
  map2%r(1,2) = val(z_offset_tot$)
  map2%r(3,4) = val(z_offset_tot$)
endif

! Front side: Offset and pitch

if (val(x_offset_tot$) /= 0 .or. val(y_offset_tot$) /= 0 .or. &
                val(x_pitch_tot$) /= 0 .or. val(y_pitch_tot$) /= 0) then
  call make_unit_mad_map (map2)
  s_here = -val(l$) / 2
  map2%k(1) = - val(x_offset_tot$) - val(x_pitch_tot$) * s_here
  map2%k(2) = - val(x_pitch_tot$) 
  map2%k(3) = - val(y_offset_tot$) - val(y_pitch_tot$) * s_here
  map2%k(4) = - val(y_pitch_tot$) 
endif

! Front side: HV kicks.
! Note: Separator already has the kicks put in the map

if ((val(hkick$) /= 0 .or. val(vkick$) /= 0) .and. ele%key /= elseparator$) then
  map2%k(2) = map2%k(2) + val(hkick$) / 2
  map2%k(4) = map2%k(4) + val(vkick$) / 2
endif

! Front side: Multipoles.

if (associated(ele%a_pole)) then
  call out_io (s_fatal$, r_name, 'YOUR MOTHER DOES NOT WORK HERE! CLEAN UP THIS MESS!')
  if (global_com%exit_on_error) call err_exit
  map2%k(2) = map2%k(2) - ele%b_pole(0) / 2
  map2%k(4) = map2%k(4) + ele%a_pole(0) / 2
  map2%r(2,1) = -ele%b_pole(1) / 2
  map2%r(2,3) =  ele%a_pole(1) / 2
  map2%r(4,1) =  ele%a_pole(1) / 2
  map2%r(4,3) =  ele%b_pole(1) / 2
  map2%t(2,1,1) = -ele%b_pole(2) / 2
  map2%t(2,3,3) =  ele%b_pole(2) / 2
  map2%t(4,1,3) =  ele%b_pole(2) / 2
  map2%t(2,1,3) =  ele%a_pole(2) / 2
  map2%t(4,1,1) =  ele%a_pole(2) / 2
  map2%t(4,3,3) = -ele%a_pole(2) / 2
  call mad_tmsymm (map2%t)
endif

! concat with map of body

call mad_concat_map2 (map2, map, map)

! Back side: Unit map.

call make_unit_mad_map (map2)

! Back side: Multipoles.

if (associated(ele%a_pole)) then
  call out_io (s_fatal$, r_name, 'YOUR MOTHER DOES NOT WORK HERE! CLEAN UP THIS MESS!')
  if (global_com%exit_on_error) call err_exit
  map2%k(2) = map2%k(2) - ele%b_pole(0) / 2
  map2%k(4) = map2%k(4) + ele%a_pole(0) / 2
  map2%r(2,1) = -ele%b_pole(1) / 2
  map2%r(2,3) =  ele%a_pole(1) / 2
  map2%r(4,1) =  ele%a_pole(1) / 2
  map2%r(4,3) =  ele%b_pole(1) / 2
  map2%t(2,1,1) = -ele%b_pole(2) / 2
  map2%t(2,3,3) =  ele%b_pole(2) / 2
  map2%t(4,1,3) =  ele%b_pole(2) / 2
  map2%t(2,1,3) =  ele%a_pole(2) / 2
  map2%t(4,1,1) =  ele%a_pole(2) / 2
  map2%t(4,3,3) = -ele%a_pole(2) / 2
  call mad_tmsymm (map2%t)
endif

! Back side: HV Kicks.
! Note: Separator already has the kicks put in the map

if ((val(hkick$) /= 0 .or. val(vkick$) /= 0) .and. ele%key /= elseparator$) then
  map2%k(2) = map2%k(2) + val(hkick$) / 2
  map2%k(4) = map2%k(4) + val(vkick$) / 2
endif

! Back side: Offset and pitch

if (val(x_offset_tot$) /= 0 .or. val(y_offset_tot$) /= 0 .or. &
                val(x_pitch_tot$) /= 0 .or. val(y_pitch_tot$) /= 0) then
  call make_unit_mad_map (map2)
  s_here = val(l$) / 2
  map2%k(1) = - val(x_offset_tot$) - val(x_pitch_tot$) * s_here
  map2%k(2) = - val(x_pitch_tot$) 
  map2%k(3) = - val(y_offset_tot$) - val(y_pitch_tot$) * s_here
  map2%k(4) = - val(y_pitch_tot$) 
endif

! Back side: S_offset.

if (val(z_offset_tot$) /= 0) then
  map2%r(1,2) = -val(z_offset_tot$)
  map2%r(3,4) = -val(z_offset_tot$)
endif

! concat with map of body

call mad_concat_map2 (map, map2, map)

end subroutine mad_add_offsets_and_multipoles

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine mad_drift (ele, energy, map)
!
! Subroutine to make a transport map for a drift space.
! The equivalent MAD-8 routine is: TMDRF
!
! Input:
!   ele    -- Ele_struct: Drift element.
!   energy -- Mad_energy_struct: particle energy structure.
!
! Output:
!   map -- Mad_map_struct: Structure holding the transfer map.
!-

Subroutine mad_drift (ele, energy, map)

implicit none

type (ele_struct) ele
type (mad_energy_struct) energy
type (mad_map_struct), target :: map

real(rp) el, beta, gamma, f
real(rp), pointer :: re(:,:), te(:,:,:)

! Init

call make_unit_mad_map (map)

el = ele%value(l$)
beta = energy%beta
gamma = energy%gamma

re => map%r
te => map%t

! First order terms.

re(1,2) = el
re(3,4) = el
re(5,6) = el / (beta * gamma) ** 2

! second-order terms.
 
f = - el / (2.0 * beta)
te(1,2,6) = f
te(1,6,2) = f
te(3,4,6) = f
te(3,6,4) = f
te(5,2,2) = f
te(5,4,4) = f
te(5,6,6) = f * 3.0 / (beta * gamma) ** 2

end subroutine mad_drift

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine mad_elsep (ele, energy, map)
!
! Subroutine to make a transport map for an electric separator. 
! The equivalent MAD-8 routine is: TMSEP
!
! Input:
!   ele    -- Ele_struct: Electric seperator element.
!   energy -- Mad_energy_struct: particle energy structure.
!
! Output:
!   map -- Mad_map_struct: Structure holding the transfer map.
!-

Subroutine mad_elsep (ele, energy, map)

implicit none

type (ele_struct) ele
type (mad_energy_struct) energy            
type (mad_map_struct), target :: map

real(rp) el, beta, gamma, fact
real(rp) ekick, ekl, ch, sh, sy, dy, tilt
     
real(rp), pointer :: ek(:), re(:,:), te(:,:,:)

! check degenerate case

if (ele%value(hkick$) == 0 .and. ele%value(vkick$) == 0) then
  call mad_drift (ele, energy, map)
  return
endif

! Init

call make_unit_mad_map (map)

el = ele%value(l$)
beta = energy%beta
gamma = energy%gamma

ek => map%k
re => map%r
te => map%t

! Prepare linear transformation particle energys.
!    DY = (COSH(K*L) - 1) / K.

ekick = sqrt(ele%value(hkick$)**2 + ele%value(vkick$)**2) * charge_of(energy%particle) / el
ekl   = ekick * el

if (abs(ekl) > 1d-6) then
  ch = cosh(ekl)
  sh = sinh(ekl)
  sy = sh / ekick
  dy = (ch - 1) / ekick**2
else
  ch = (1 + ekl**2 / 2)
  sy = (1 + ekl**2 / 6) * el
  sh = sy * ekick
  dy = (0.5 + ekl**2 / 24) * el**2
endif

! kicks.

ek(3) = dy * (ekick / beta)
ek(4) = sy * (ekick / beta)

! first-order terms.

re(1,2) = el
re(3,3) = ch - ekl * sh / beta**2
re(3,4) = sy
re(3,6) = (dy - el * sy / beta**2) * ekick
re(4,3) = (sh - ekl * ch / beta**2) * ekick
re(4,4) = ch
re(4,6) = (sh - ekl * ch / beta**2)
re(5,3) = - re(4,6)
re(5,4) = - dy * ekick
re(5,6) = - (sy - el * ch / beta**2)

! second-order terms.

fact = el / (2 * beta)
te(1,2,3) = - fact * ekick
te(1,2,6) = - fact

fact = el * (3*sh/gamma**2 + ekl*ch) / (2*beta**3)
te(3,3,3) = fact * ekick**2
te(3,3,6) = fact * ekick
te(3,6,6) = fact

fact = el * (3*ch/gamma**2 + ekl*sh) / (2*beta**3)
te(4,3,3) = fact * ekick**3
te(4,3,6) = fact * ekick**2
te(4,6,6) = fact * ekick
te(5,3,3) = - fact * ekick**2
te(5,3,6) = - fact * ekick
te(5,6,6) = - fact

fact = el * sh / (2 * beta)
te(3,2,2) = fact
te(3,4,4) = fact
te(4,3,4) = - fact * ekick**2
te(4,4,6) = - fact * ekick
te(5,3,4) = fact * ekick
te(5,4,6) = fact

fact = el * ch / (2 * beta)
te(3,3,4) = - fact * ekick
te(3,4,6) = - fact
te(4,2,2) = fact * ekick
te(4,4,4) = fact * ekick
te(5,2,2) = - fact
te(5,4,4) = - fact

call mad_tmsymm(te)

! apply tilt
                  
tilt = -atan2 (ele%value(hkick$), ele%value(vkick$))  
call mad_tmtilt (map, tilt)

end subroutine mad_elsep

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine mad_sextupole (ele, energy, map)
!
! Subroutine to make a transport map for an sextupole.
! The equivalent MAD-8 routine is: TMSEXT
!
! Input:
!   ele    -- Ele_struct: Sextupole element.
!   energy -- Mad_energy_struct: particle energy structure.
!
! Output:
!   map -- Mad_map_struct: Structure holding the transfer map.
!-

subroutine mad_sextupole (ele, energy, map) 

implicit none
                        
type (ele_struct) ele
type (mad_energy_struct) energy
type (mad_map_struct), target :: map

real(rp) el, beta, gamma
real(rp) skl, s1, s2, s3, s4

real(rp), pointer :: re(:,:), te(:,:,:)

! Init

call make_unit_mad_map (map)

el = ele%value(l$)
beta = energy%beta
gamma = energy%gamma

re => map%r
te => map%t

! First-order terms.

re(1,2) = el
re(3,4) = el
re(5,6) = el / (beta * gamma) ** 2

! second-order terms.

skl = ele%value(k2$) * el
if (skl /= 0.0) then
  s1 = skl / 2.0
  s2 = s1 * el / 2.0
  s3 = s2 * el / 3.0
  s4 = s3 * el / 4.0
  te(1,1,1) = - s2
  te(1,1,2) = - s3
  te(1,2,2) = - 2.0 * s4
  te(1,3,3) = + s2
  te(1,3,4) = + s3
  te(1,4,4) = + 2.0 * s4
  te(2,1,1) = - s1
  te(2,1,2) = - s2
  te(2,2,2) = - 2.0 * s3
  te(2,3,3) = + s1
  te(2,3,4) = + s2
  te(2,4,4) = + 2.0 * s3
  te(3,1,3) = + s2
  te(3,1,4) = + s3
  te(3,2,3) = + s3
  te(3,2,4) = + 2.0 * s4
  te(4,1,3) = + s1
  te(4,1,4) = + s2
  te(4,2,3) = + s2
  te(4,2,4) = + 2.0 * s3
endif

te(1,2,6) = - el / (2.0 * beta)
te(3,4,6) = te(1,2,6)
te(5,2,2) = te(1,2,6)
te(5,4,4) = te(1,2,6)
te(5,6,6) = - 3.0 * re(5,6) / (2.0 * beta)

call mad_tmsymm(te)

! Apply tilt.

if (ele%value(tilt_tot$) /= 0.0) call mad_tmtilt(map, ele%value(tilt_tot$))

end subroutine mad_sextupole

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine mad_sbend (ele, energy, map)
!
! Subroutine to make a transport map for a sector bend element.
! The equivalent MAD-8 routine is: TMBEND
!
! Input:
!   ele    -- Ele_struct: Sbend element.
!   energy -- Mad_energy_struct: particle energy structure.
!
! Output:
!   map -- Mad_map_struct: Structure holding the transfer map.
!-

subroutine mad_sbend (ele, energy, map) 

implicit none
                        
type (ele_struct), target :: ele
type (mad_energy_struct) energy
type (mad_map_struct) map2, map_roll
type (mad_map_struct) map

real(rp) angle, roll

! roll

roll = ele%value(roll_tot$) 
if (roll /= 0) then
  angle = ele%value(l$) * ele%value(g$)
  call make_unit_mad_map (map_roll)
  if (abs(roll) < 0.001) then
    map_roll%k(1) = angle * roll**2 / 4
  else
    map_roll%k(1) = angle * (1 - cos(roll)) / 2
  endif
  map_roll%k(3) = -angle * sin(roll) / 2
  call mad_concat_map2 (map_roll, map2, map2)
endif

! Entrance fringe.

call mad_sbend_fringe (ele, energy, .true., map2)
if (roll /= 0) call mad_concat_map2 (map_roll, map2, map2)

call mad_sbend_body (ele, energy, map)
call mad_concat_map2 (map2, map, map)

call mad_sbend_fringe (ele, energy, .false., map2)
call mad_concat_map2 (map, map2, map)

if (roll /= 0) map%k = map%k + map_roll%k

if (ele%value(ref_tilt_tot$) /= 0.0) call mad_tmtilt(map, ele%value(ref_tilt_tot$))

end subroutine mad_sbend

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+                     
! Subroutine mad_sbend_fringe (ele, energy, into, map)
!
! Subroutine to make a transport map for the fringe field of a dipole.
! The equivalent MAD-8 routine is: TMFRNG
!
! Input:
!   ele    -- Ele_struct: Solenoid element.
!   energy -- Mad_energy_struct: particle energy structure.
!   into   -- Logical: If True then map is for particle entering a dipole
!
! Output:
!   map    -- Mad_map_struct: Fringe dipole map.
!     %k(6)     -- 0th order map.
!     %r(6,6)   -- 1st order map.         
!     %t(6,6,6) -- 2nd order map.
!-

subroutine mad_sbend_fringe (ele, energy, into, map) 

implicit none
        
type (ele_struct) ele
type (mad_energy_struct) energy
type (mad_map_struct), target ::  map

real(rp), pointer :: re(:,:), te(:,:,:)
real(rp) h, he, hh, edge, tanedg, secedg, psip, sk1

logical into

! Setup.
! corr is correction factor according to SLAC 75. Not used in this version

call make_unit_mad_map (map)

re => map%r
te => map%t

h = ele%value(g$)
he = 0   ! curvature of the pole face
sk1 = ele%value(k1$)

if (into) then
  edge = ele%value(e1$)
  hh = h / 2
else
  edge = ele%value(e2$)
  hh = -h / 2
endif

tanedg = tan(edge)
secedg = 1.0 / cos(edge)
psip = edge    ! - corr * secedg * (1.0 + sin(edge)**2)

! Linear terms.

re(2,1) = + h * tanedg
re(4,3) = - h * tan(psip)

! Second-order terms.

te(1,1,1) = - hh * tanedg**2
te(1,3,3) = + hh * secedg**2
te(2,1,1) = (h/2) * he * secedg**3 + sk1 * tanedg
te(2,1,2) = - te(1,1,1)
te(2,3,3) = hh * h * tanedg**3 - te(2,1,1)
te(2,3,4) = + te(1,1,1)
te(3,1,3) = - te(1,1,1)
te(4,1,3) = - te(2,1,1)
te(4,1,4) = + te(1,1,1)
te(4,2,3) = - te(1,3,3)

if (into) then
  te(2,3,3) = te(2,3,3) + (h*secedg)**2 * tanedg/2
else
  te(2,1,1) = te(2,1,1) - (h*tanedg)**2 * tanedg/2
  te(4,1,3) = te(4,1,3) + (h*secedg)**2 * tanedg/2
endif

call mad_tmsymm(te)

end subroutine mad_sbend_fringe

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+                     
! Subroutine mad_sbend_body (ele, energy, map)
!
! Subroutine to make a transport map for the body of a sector dipole.
! The equivalent MAD-8 routine is: TMSECT
!
! Input:
!   ele    -- Ele_struct: Solenoid element.
!   energy -- Mad_energy_struct: particle energy structure.
!   into   -- Logical: If True then map is for particle entering a dipole
!
! Output:
!   map -- Mad_map_struct: Structure holding the transfer map.
!-

subroutine mad_sbend_body (ele, energy, map) 

implicit none
        
type (ele_struct) ele
type (mad_energy_struct) energy
type (mad_map_struct), target :: map

real(rp), pointer :: ek(:), re(:,:), te(:,:,:)
real(rp) bi, bi2, bi2gi2, xksq, xk, xkl, xklsq
real(rp) cx, sx, cy, sy, dx, fx, gx, hx, yksq, yk, ykl, yklsq
real(rp) xs6, ys2, h2, t336, t346, t436, t446
real(rp) t116, t126, t166, t216, t226, t266, t516, t526, t566
real(rp) t1, t2, t5, y2ksq, y2klsq                     
real(rp) y0, y1, y2, zc, zs, zd, zf, dd, sumsq, difsq
real(rp) sk1, sk2, h, el, beta, gamma
real(rp) cyy, syy, dyy, fyy, cp, sp, dp, fp, cm, sm, dm, fm

real(rp), parameter :: C1  =   1D0,           C2  =   1D0 / 2D0
real(rp), parameter :: C3  =   1D0 / 24D0,    C4  =   1D0 / 720D0
real(rp), parameter :: S1  =   1D0,           S2  =   1D0 / 6D0
real(rp), parameter :: S3  =   1D0 / 120D0,   S4  =   1D0 / 5040D0
real(rp), parameter :: CG0 =   1D0 / 20D0,    CG1 =   5D0 / 840D0
real(rp), parameter :: CG2 =  21D0 / 60480D0, CH0 =   1D0 / 56D0
real(rp), parameter :: CH1 =  14D0 / 4032D0,  CH2 = 147D0 / 443520D0

! Setup.

call make_unit_mad_map (map)

el = ele%value(l$)
beta = energy%beta
gamma = energy%gamma

sk1 = ele%value(k1$)
sk2 = ele%value(k2$)     ! sextupole strength
h = ele%value(g$)

ek => map%k
re => map%r
te => map%t

bi = 1.0 / beta
bi2 = bi * bi
bi2gi2 = 1.0 / (beta * gamma) ** 2

! Horizontal.

xksq = h**2 + sk1
xk = sqrt(abs(xksq))
xkl = xk * el
xklsq = xksq * el**2

if (abs(xklsq) < 1.0d-2) then
  cx = (c1 - xklsq * (c2 - xklsq*c3))
  sx = (s1 - xklsq * (s2 - xklsq*s3)) * el
  dx = (c2 - xklsq * (c3 - xklsq*c4)) * el**2
  fx = (s2 - xklsq * (s3 - xklsq*s4)) * el**3
  gx = (cg0 - xklsq * (cg1 - xklsq*cg2)) * el**5
  hx = (ch0 - xklsq * (ch1 - xklsq*ch2)) * el**7
else
  if (xklsq > 0.0) then
    cx = cos(xkl)
    sx = sin(xkl) / xk
  else
    cx = cosh(xkl)
    sx = sinh(xkl) / xk
  endif
  dx = (1.0 - cx) / xksq
  fx = (el  - sx) / xksq
  gx = (3.0*el - sx*(4.0-cx)) / (2.0*xksq**2)
  hx = (15.0*el - sx*(22.0-9.0*cx+2.0*cx**2)) / (6.0*xksq**3)
endif

re(1,1) = cx
re(1,2) = sx
re(1,6) = h * dx * bi
re(2,1) = - xksq * sx
re(2,2) = cx
re(2,6) = h * sx * bi
re(5,2) = - re(1,6)
re(5,1) = - re(2,6)
re(5,6) = el * bi2gi2 - h**2 * fx * bi2

! Vertical.

yksq = - sk1
yk = sqrt(abs(yksq))
ykl = yk*el
yklsq = yksq*el**2

if (abs(yklsq) < 1.0d-2) then
  cy = (c1 - yklsq * (c2 - yklsq*c3))
  sy = (s1 - yklsq * (s2 - yklsq*s3)) * el
else if (yklsq > 0.0) then
  cy = cos(ykl)
  sy = sin(ykl) / yk
else
  cy = cosh(ykl)
  sy = sinh(ykl) / yk
endif

re(3,3) = cy
re(3,4) = sy
re(4,3) = - yksq * sy
re(4,4) = cy

ek(3)   = 0.0
ek(4)   = 0.0

! Second-order terms...
! Pure horizontal terms.

xs6 = (sk2 + 2.0*h*sk1) / 6.0
ys2 = (sk2 +   h*sk1) / 2.0
h2 = h / 2.0
t116 = xs6 * (3.0*sx*fx - dx**2) - h * sx**2
t126 = xs6 * (sx*dx**2 - 2.0*cx*gx) - h * sx * dx
t166 = xs6 * (dx**3 - 2.0*sx*gx) - h2 * dx**2
t216 = xs6 * (3.0*cx*fx + sx*dx)
t226 = xs6 * (3.0*sx*fx + dx**2)
t266 = xs6 * (sx*dx**2 - 2.0*cx*gx)
t516 = h * xs6 * (3.0*dx*fx - 4.0*gx) + &
                      (sk1/2.0) * (fx + sx*dx)
t526 = h * xs6 * (dx**3 - 2.0*sx*gx) + (sk1/2.0) * dx**2
t566 = h * xs6 * (3.0*hx - 2.0*dx*gx) + &
                          (sk1/2.0) * gx - fx
t1 = (sk1/2.0) * (dx**2 - sx*fx) - dx
t2 = (sk1/2.0) * (el*dx - fx)
t5 = fx - sk1 * (gx - fx*dx / 2.0)
te(1,1,1) = - xs6 * (sx**2 + dx) - h2*xksq*sx**2
te(1,1,2) = (- xs6*dx + h2*cx) * sx
te(1,2,2) = (- xs6*dx + h2*cx) * dx
te(1,1,6) = (- h2*t116 + (sk1/4.0)*el*sx) * bi
te(1,2,6) = (- h2*t126 + (sk1/4.0) * (el*dx - fx) - sx/2.0) * bi
te(1,6,6) = (- h**2*t166 + h*t1) * bi2 - h2 * dx * bi2gi2
te(2,1,1) = - xs6 * (1.0 + 2.0*cx) * sx
te(2,1,2) = - xs6 * (1.0 + 2.0*cx) * dx
te(2,2,2) = - (2.0*xs6*dx + h2) * sx
te(2,1,6) = (- h2*t216 - (sk1/4.0) * (sx - el*cx)) * bi
te(2,2,6) = (- h2*t226 + (sk1/4.0) * el * sx) * bi
te(2,6,6) = (- h**2*t266 + h*t2) * bi2 - h2 * sx * bi2gi2
te(5,1,1) = (h2*xs6 * (sx*dx + 3.0*fx) - &
                      (sk1/4.0) * (el - cx*sx)) * bi
te(5,1,2) = (h2*xs6*dx**2 + (sk1/4.0)*sx**2) * bi
te(5,2,2) = (h*xs6*gx - sk1 * (fx - sx*dx) / 4.0 - sx/2.0) * bi
te(5,1,6) = h2 * ((t516 - sk1 * (el*dx - fx) / 2.0) * bi2 + &
                              sx * bi2gi2)
te(5,2,6) = h2 * ((t526 - sk1 * (dx**2 - sx*fx) / 2.0) * bi2 + &
                       dx * bi2gi2)
te(5,6,6) = (h**2 * (t566 + t5) * bi2 + &
                 (3.0/2.0) * (h**2*fx - el) * bi2gi2) * bi

! Mixed terms.

y2ksq = 4.0 * yksq
call mad_tmfoc(el, y2ksq, cyy, syy, dyy, fyy)
y2klsq = y2ksq * el**2

if (max(abs(y2klsq),abs(xklsq)) .le. 1.0d-2) then
  y0 = 1.0
  y1 = xklsq + y2klsq
  y2 = xklsq**2 + xklsq*y2klsq + y2klsq**2
  zc = (y0 - (y1 - y2 / 30.0) / 12.0) * el**2 /   2.0
  zs = (y0 - (y1 - y2 / 42.0) / 20.0) * el**3 /   6.0
  zd = (y0 - (y1 - y2 / 56.0) / 30.0) * el**4 /  24.0
  zf = (y0 - (y1 - y2 / 72.0) / 42.0) * el**5 / 120.0
else if (xksq .le. 0.0  .or.  yksq .le. 0.0) then
  dd = xksq - y2ksq
  zc = (cyy - cx) / dd
  zs = (syy - sx) / dd
  zd = (dyy - dx) / dd
  zf = (fyy - fx) / dd
else
  sumsq = (xk/2.0 + yk) ** 2
  difsq = (xk/2.0 - yk) ** 2
  call mad_tmfoc(el, sumsq, cp, sp, dp, fp)
  call mad_tmfoc(el, difsq, cm, sm, dm, fm)
  zc = sp * sm / 2.0
  zs = (sp*cm - cp*sm) / (4.0*xk*yk)
  if (xksq > y2ksq) then
    zd = (dyy - zc) / xksq
    zf = (fyy - zs) / xksq
  else
    zd = (dx - zc) / y2ksq
    zf = (fx - zs) / y2ksq
  endif
endif

t336 = sk2 * (cy*zd - 2.0*sk1*sy*zf) + h * sk1 * fx * sy
t346 = sk2 * (sy*zd - 2.0*cy*zf) + h * fx * cy
t436 = 2.0 * ys2 * fx * cy - sk2 * sk1 * (sy*zd - 2.0*cy*zf)
t446 = 2.0 * ys2 * fx * sy - sk2 * (cy*zd - 2.0*sk1*sy*zf)
te(1,3,3) = + sk2*sk1*zd + ys2*dx
te(1,3,4) = + sk2*zs/2.0
te(1,4,4) = + sk2*zd - h2*dx
te(2,3,3) = + sk2*sk1*zs + ys2*sx
te(2,3,4) = + sk2*zc/2.0
te(2,4,4) = + sk2*zs - h2*sx
te(3,1,3) = + sk2*(cy*zc/2.0 - sk1*sy*zs) + h2*sk1*sx*sy
te(3,1,4) = + sk2*(sy*zc/2.0 - cy*zs) + h2*sx*cy
te(3,2,3) = + sk2*(cy*zs/2.0 - sk1*sy*zd) + h2*sk1*dx*sy
te(3,2,4) = + sk2*(sy*zs/2.0 - cy*zd) + h2*dx*cy
te(3,3,6) = (h2*t336 - sk1*el*sy/4.0) * bi
te(3,4,6) = (h2*t346 - (sy + el*cy) / 4.0) * bi
te(4,1,3) = sk2*sk1*(cy*zs - sy*zc/2.0) + ys2*sx*cy
te(4,1,4) = sk2*(sk1*sy*zs - cy*zc/2.0) + ys2*sx*sy
te(4,2,3) = sk2*sk1*(cy*zd - sy*zs/2.0) + ys2*dx*cy
te(4,2,4) = sk2*(sk1*sy*zd - cy*zs/2.0) + ys2*dx*sy
te(4,3,6) = (h2*t436 + sk1 * (sy - el*cy) / 4.0) * bi
te(4,4,6) = (h2*t446 - sk1*el*sy/4.0) * bi
te(5,3,3) = (- h*sk2*sk1*zf - h*ys2*fx + sk1*(el-cy*sy)/4.0)*bi
te(5,3,4) = (- h*sk2*zd/2.0 - sk1*sy**2/4.0) * bi
te(5,4,4) = (- h*sk2*zf + h*h2*fx - (el + sy*cy)/4.0) * bi

call mad_tmsymm(te)

end subroutine mad_sbend_body

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+                     
! Subroutine mad_tmfoc (el, sk1, c, s, d, f) 
!
! Subroutine to compute the linear focussing functions.  
! The equivalent MAD-8 routine is: TMFOC
!
! Input:
!   el   -- Real(rp): Length.                                 
!   sk1  -- Real(rp): Quadrupole strength.                    
! Output:                                                       
!   c    -- Real(rp): Cosine-like function.             c(k,l)
!   s    -- Real(rp): Sine-like function.               s(k,l)
!   d    -- Real(rp): Dispersion function.              d(k,l)
!   f    -- Real(rp): Integral of dispersion function.  f(k,l)
!-

subroutine mad_tmfoc (el, sk1, c, s, d, f) 

implicit none

real(rp) el, sk1, c, s, d, f
real(rp) qk, qkl, qkl2

!

qk = sqrt(abs(sk1))
qkl = qk * el
qkl2 = sk1 * el**2

if (abs(qkl2) .le. 1.0d-2) then
  c = (1.0 - qkl2 * (1.0 - qkl2 / 12.0) /  2.0)
  s = (1.0 - qkl2 * (1.0 - qkl2 / 20.0) /  6.0) * el
  d = (1.0 - qkl2 * (1.0 - qkl2 / 30.0) / 12.0) * el**2 / 2.0
  f = (1.0 - qkl2 * (1.0 - qkl2 / 42.0) / 20.0) * el**3 / 6.0
else
  if (qkl2 > 0.0) then
    c = cos(qkl)
    s = sin(qkl) / qk
  else
    c = cosh(qkl)
    s = sinh(qkl) / qk
  endif
  d = (1.0 - c) / sk1
  f = (el  - s) / sk1
endif

end subroutine mad_tmfoc

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+                     
! Subroutine mad_quadrupole (ele, energy, map)
!
! Subroutine to make a transport map for an quadrupole element.
! The equivalent MAD-8 routine is: TMSEXT
!
! Input:
!   ele    -- Ele_struct: Quadrupole element.
!   energy -- Mad_energy_struct: particle energy structure.
!
! Output:
!   map -- Mad_map_struct: Structure holding the transfer map.
!-

subroutine mad_quadrupole (ele, energy, map) 

implicit none
                        
type (ele_struct) ele
type (mad_energy_struct) energy
type (mad_map_struct), target :: map

real(rp) el, beta, gamma
real(rp) sk1, qk, qkl, qkl2, cx, sx, cy, sy, biby4

real(rp), pointer :: re(:,:), te(:,:,:)

! Init

call make_unit_mad_map (map)

el = ele%value(l$)
sk1 = ele%value(k1$)

beta = energy%beta
gamma = energy%gamma

re => map%r
te => map%t

! Set up c's and s's.

qk = sqrt(abs(sk1))
qkl = qk * el

if (abs(qkl) < 1.0d-3) then
  qkl2 = sk1 * el**2
  cx = (1.0 - qkl2 / 2.0)
  sx = (1.0 - qkl2 / 6.0) * el
  cy = (1.0 + qkl2 / 2.0)
  sy = (1.0 + qkl2 / 6.0) * el
else if (sk1 > 0.0) then
  cx = cos(qkl)
  sx = sin(qkl) / qk
  cy = cosh(qkl)
  sy = sinh(qkl) / qk
else
  cx = cosh(qkl)
  sx = sinh(qkl) / qk
  cy = cos(qkl)
  sy = sin(qkl) / qk
endif

! first-order terms.

re(1,1) = cx
re(1,2) = sx
re(2,1) = - sk1 * sx
re(2,2) = cx
re(3,3) = cy
re(3,4) = sy
re(4,3) = + sk1 * sy
re(4,4) = cy
re(5,6) = el / (beta * gamma) ** 2

! second-order terms.

biby4 = 1.0 / (4.0 * beta)

te(1,1,6) = + sk1 * el * sx * biby4
te(1,6,1) = te(1,1,6)
te(2,2,6) = te(1,1,6)
te(2,6,2) = te(1,1,6)
te(1,2,6) = - (sx + el*cx) * biby4
te(1,6,2) = te(1,2,6)
te(2,1,6) = - sk1 * (sx - el*cx) * biby4
te(2,6,1) = te(2,1,6)

te(3,3,6) = - sk1 * el * sy * biby4
te(3,6,3) = te(3,3,6)
te(4,4,6) = te(3,3,6)
te(4,6,4) = te(3,3,6)
te(3,4,6) = - (sy + el*cy) * biby4
te(3,6,4) = te(3,4,6)
te(4,3,6) = + sk1 * (sy - el*cy) * biby4
te(4,6,3) = te(4,3,6)

te(5,1,1) = - sk1 * (el - sx*cx) * biby4
te(5,1,2) = + sk1 * sx**2 * biby4
te(5,2,1) = te(5,1,2)
te(5,2,2) = - (el + sx*cx) * biby4
te(5,3,3) = + sk1 * (el - sy*cy) * biby4
te(5,3,4) = - sk1 * sy**2 * biby4
te(5,4,3) = te(5,3,4)
te(5,4,4) = - (el + sy*cy) * biby4
te(5,6,6) = (- 6.0 * re(5,6)) * biby4

! Apply tilt.

if (ele%value(tilt_tot$) /= 0.0) call mad_tmtilt(map, ele%value(tilt_tot$))

end subroutine mad_quadrupole

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine mad_rfcavity (ele, energy, map)
!
! Subroutine to make a transport map for an rfcavity element.
! The equivalent MAD-8 routine is: TMRF
!
! Input:
!   ele    -- Ele_struct: Rfcavity element.
!   energy -- Mad_energy_struct: particle energy structure.
!
! Output:
!   map -- Mad_map_struct: Structure holding the transfer map.
!-

subroutine mad_rfcavity (ele, energy, map) 

implicit none
        
type (ele_struct) ele
type (mad_energy_struct) energy
type (mad_map_struct) map

real(rp) omega, charge, vrf, phirf, c0, c1, c2

! Init.
! Note: orbit is assumed zero.

charge = 1   ! Assume this

omega = charge * twopi * ele%value(rf_frequency$) / c_light
vrf   = ele%value(voltage$) * charge / energy%total
phirf = (ele%value(phi0$) + ele%value(phi0_multipass$)) * twopi   ! - omega * orbit(5)

c0 =  vrf * sin(phirf)
c1 =  vrf * cos(phirf) * omega
c2 = -vrf * sin(phirf) * omega**2 / 2.0

! Map.

call mad_drift (ele, energy, map)

map%k(6) = c0                  ! - c1 * orbit(5) + c2 * orbit(5)**2
map%r(6,5) = c1                ! - 2.0 * c2 * orbit(5)
map%t(6,5,5) = c2

end subroutine mad_rfcavity

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine mad_solenoid (ele, energy, map)
!
! Subroutine to make a transport map for an solenoid.
! The equivalent MAD-8 routine is: TMSEXT
!
! Input:
!   ele    -- Ele_struct: Solenoid element.
!   energy -- Mad_energy_struct: particle energy structure.
!
! Output:
!   map -- Mad_map_struct: Structure holding the transfer map.
!-

subroutine mad_solenoid (ele, energy, map) 

implicit none
        
type (ele_struct) ele
type (mad_energy_struct) energy
type (mad_map_struct), target :: map

real(rp) el, beta, gamma, temp
real(rp) sk, skl, co, si, sibk, sks

real(rp), pointer :: re(:,:), te(:,:,:)

! Init

call make_unit_mad_map (map)

el = ele%value(l$)
beta = energy%beta
gamma = energy%gamma

re => map%r
te => map%t

sks = ele%value(ks$)

! Zeroth order terms.


! Set up C's and S's.

sk = sks / 2.0
skl = sk * el
co = cos(skl)
si = sin(skl)
if (abs(skl) < 1.0d-5) then
  sibk = (1.0 - skl**2/6.0) * el
else
  sibk = si/sk
endif

! First-order terms.

re(1,1) = co**2
re(2,2) = re(1,1)
re(3,3) = re(1,1)
re(4,4) = re(1,1)
re(1,2) = co * sibk
re(3,4) = re(1,2)
re(1,3) = co * si
re(2,4) = re(1,3)
re(3,1) = - re(1,3)
re(4,2) = re(3,1)
re(2,1) = sk * re(3,1)
re(4,3) = re(2,1)
re(1,4) = si * sibk
re(3,2) = - re(1,4)
re(4,1) = sk * si**2
re(2,3) = - re(4,1)
re(5,6) = el / (beta * gamma) ** 2

! Second-order terms.

temp = el * co * si / beta
te(1,4,6) = - temp
te(3,2,6) =   temp
te(1,1,6) =   temp * sk
te(2,2,6) =   temp * sk
te(3,3,6) =   temp * sk
te(4,4,6) =   temp * sk
te(2,3,6) =   temp * sk**2
te(4,1,6) = - temp * sk**2

temp = el * (co**2 - si**2) / (2.0 * beta)
te(1,2,6) = - temp
te(3,4,6) = - temp
te(1,3,6) = - temp * sk
te(2,4,6) = - temp * sk
te(3,1,6) =   temp * sk
te(4,2,6) =   temp * sk
te(2,1,6) =   temp * sk**2
te(4,3,6) =   temp * sk**2

temp = el / (2.0 * beta)
te(5,2,2) = - temp
te(5,4,4) = - temp
te(5,1,4) =   temp * sk
te(5,2,3) = - temp * sk
te(5,1,1) = - temp * sk**2
te(5,3,3) = - temp * sk**2
te(5,6,6) = - 3.0 * re(5,6) / (2.0 * beta)

call mad_tmsymm(te)

end subroutine mad_solenoid

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! subroutine mad_tmsymm (te)
!
! subroutine to symmertrize the 2nd order map t.
! The equivalent MAD-8 routine is: tmsymm
!
! input:
!   te(6,6,6) -- real(rp): array to be symmertrized.
!
! output:
!   te(6,6,6) -- real(rp): symmetrized array.
!-

subroutine mad_tmsymm (te)

implicit none

real(rp) te(6,6,6)
integer k, l

!

do k = 1, 5
  do l = k+1, 6
    te(:,l,k) = te(:,k,l)
  enddo
enddo

end subroutine mad_tmsymm

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine mad_tmtilt (map, tilt)
!
! Subroutine to apply a tilt to a transport map.
! The equivalent MAD-8 routine is: TMTILT
!
! Input:
!   map  -- Mad_map_struct: Unrotated transport map.
!     %k(6)     -- 0th order map.
!     %r(6,6)   -- 1st order map.
!     %t(6,6,6) -- 2nd order map.
!   tilt -- Real(rp): Tilt
!
! Output:
!   map  -- Mad_map_struct: Rotated transport map.
!     %k(6)     -- 0th order map.
!     %r(6,6)   -- 1st order map.
!     %t(6,6,6) -- 2nd order map.
!-

subroutine mad_tmtilt (map, tilt) 

implicit none

type (mad_map_struct), target :: map

real(rp) tilt
real(rp) c, r1j, r2j, ri1, ri2, s, t1jk
real(rp) t2jk, ti1k, ti2k, tij1, tij2, xx

real(rp), pointer :: ek(:), re(:,:), te(:,:,:)   

integer i, j, k

! Setup.

ek => map%k
re => map%r
te => map%t

c =  cos(tilt)
s = sin(tilt)

! Rotate at entrance.

do i = 1, 6

  ri1 = re(i,1)
  re(i,1) = ri1 * c - re(i,3) * s
  re(i,3) = ri1 * s + re(i,3) * c
  ri2 = re(i,2)
  re(i,2) = ri2 * c - re(i,4) * s
  re(i,4) = ri2 * s + re(i,4) * c

  do k = 1, 6
    ti1k = te(i,1,k)
    te(i,1,k) = ti1k * c - te(i,3,k) * s
    te(i,3,k) = ti1k * s + te(i,3,k) * c
    ti2k = te(i,2,k)
    te(i,2,k) = ti2k * c - te(i,4,k) * s
    te(i,4,k) = ti2k * s + te(i,4,k) * c
  enddo

  do j = 1, 6
    tij1 = te(i,j,1)
    te(i,j,1) = tij1 * c - te(i,j,3) * s
    te(i,j,3) = tij1 * s + te(i,j,3) * c
    tij2 = te(i,j,2)
    te(i,j,2) = tij2 * c - te(i,j,4) * s
    te(i,j,4) = tij2 * s + te(i,j,4) * c
  enddo

enddo

! Rotate kick.

xx = ek(1)
ek(1) = xx * c - ek(3) * s
ek(3) = xx * s + ek(3) * c

xx = ek(2)
ek(2) = xx * c - ek(4) * s
ek(4) = xx * s + ek(4) * c

! Rotate at exit.

do j = 1, 6

  r1j = re(1,j)
  re(1,j) = c * r1j - s * re(3,j)
  re(3,j) = s * r1j + c * re(3,j)
  r2j = re(2,j)
  re(2,j) = c * r2j - s * re(4,j)
  re(4,j) = s * r2j + c * re(4,j)

  do k = 1, 6
    t1jk = te(1,j,k)
    te(1,j,k) = c * t1jk - s * te(3,j,k)
    te(3,j,k) = s * t1jk + c * te(3,j,k)
    t2jk = te(2,j,k)
    te(2,j,k) = c * t2jk - s * te(4,j,k)
    te(4,j,k) = s * t2jk + c * te(4,j,k)
  enddo

enddo

end subroutine mad_tmtilt

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine mad_concat_map2 (map1, map2, map3)
!
! Subroutine to concatinate two 2nd order transport maps.
!     map3 = map2(map1)
! The equivalent MAD-8 routine is: TMCAT1
!
! Input:
!   map1 -- Mad_map_struct: First map in the beam line.
!   map2 -- Mad_map_struct: Second map in the beam line.
!
! Output:
!   map3 -- Mad_map_struct: Concatinated map.
!-

subroutine mad_concat_map2 (map1, map2, map3)

type (mad_map_struct) :: map1, map2, map3

real(rp) :: ek(6), re(6,6), te(6,6,6), te2(6,6,6)

integer i, j, k


!-------------------------------------------
! If no constant part in map1

if (all(map1%k == 0)) then

  do i = 1, 6
    te(i,:,:) = matmul(map2%t(i,:,:), map1%r(:,:))
  enddo

  do j = 1, 6
    te2(:,j,:) = matmul(map2%r(:,:), map1%t(:,j,:)) + &
                                 matmul(te(:,:,j), map1%r(:,:))
  enddo

  re(:,:) = matmul (map2%r(:,:), map1%r(:,:))

  map3%k = map2%k
  map3%r = re
  map3%t = te2

  return

endif

!-------------------------------------------
! Here if there is a constant part
! Auxiliary terms.

do k = 1, 6
  re(:,k) = matmul(map2%t(:,:,k), map1%k(:))  ! ek1 * Te2
enddo

do i = 1, 6
  te(i,:,:) = matmul(map2%t(i,:,:), map1%r(:,:))
enddo

! Final values.

ek(:) = map2%k(:) + matmul(map2%r(:,:) + re(:,:), map1%k(:))

do j = 1, 6
  te2(:,j,:) = matmul(map2%r(:,:)+2*re(:,:), map1%t(:,j,:)) + &
                                      matmul(te(:,:,j), map1%r(:,:))
enddo

re(:,:) = matmul (map2%r(:,:) + 2*re(:,:), map1%r(:,:))

!

map3%k = ek
map3%r = re
map3%t = te2

end subroutine mad_concat_map2
      
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine mad_track1 (c0, map, c1)
!
! Subroutine to track through a 2nd order transfer map.
! The equivalent MAD-8 routine is: TMTRAK
!
! Input:
!   c0   -- Coord_struct: Starting coords.
!   map -- Mad_map_struct:  2nd order map.
!
! Output:
!   c1   -- Coord_struct: Ending coords. 
!-

subroutine mad_track1 (c0, map, c1)   

implicit none

type (mad_map_struct) map
type (coord_struct) c0, c1

real(rp) vec0(6), vec1(6), mat(6,6)

integer i

!

vec0 = c0%vec

do i = 1, 6
  mat(i,:) = map%r(i,:) + matmul(vec0, map%t(i,:,:))
enddo

vec1 = map%k + matmul(mat, vec0)

c1%vec = vec1 

end subroutine mad_track1

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine track1_mad (orbit, ele, param)
!
! Subroutine to track through an element using a 2nd order transfer map.
! Note: If map does not exist then one will be created. 
!
! Input:
!   orbit      -- Coord_struct: Starting coords.
!   ele        -- Ele_struct: Element to track through.
!   param      -- lat_param_struct: Lattice parameters.
!
! Output:
!   orbit      -- Coord_struct: Ending coords.
!-

subroutine track1_mad (orbit, ele, param)

implicit none

type (ele_struct) ele
type (lat_param_struct) param
type (coord_struct) orbit, start_orb
type (mad_energy_struct) energy
type (mad_map_struct) map
real(rp) dtime_ref

!

if (.not. associated(ele%taylor(1)%term)) then
  call make_mad_map (ele, param, energy, map)
  call mad_map_to_taylor (map, energy, ele%taylor)
endif

start_orb = orbit
orbit%vec = track_taylor (orbit%vec, ele%taylor)

orbit%s = ele%s
orbit%p0c = ele%value(p0c$)

! If delta_ref_time has not been set then just assume that the particle has constant velocity.

if (ele%value(p0c$) /= ele%value(p0c_start$)) then
  call convert_pc_to (ele%value(p0c$) * (1 + orbit%vec(6)), param%particle, beta = orbit%beta)
endif

dtime_ref = ele%value(delta_ref_time$)
if (dtime_ref == 0) dtime_ref = ele%value(l$) / (orbit%beta * c_light)

orbit%t = start_orb%t + dtime_ref + start_orb%vec(5) / (start_orb%beta * c_light) - orbit%vec(5) / (orbit%beta * c_light)

end subroutine track1_mad

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine mad_map_to_taylor (map, energy, taylor)
!
! Subroutine to convert a MAD order 2 map to a Bmad taylor map.
! The conversion will also convert between MAD's (t, dE) and Bmad's (beta*t, dP) coords.
!
! Input:
!   map       -- Mad_map_struct: Order 2 map.
!   energy    -- mad_energy_struct: Energy numbers.
!
! Output:
!   taylor(6) -- Taylor_struct: Taylor map.
!-

subroutine mad_map_to_taylor (map, energy, taylor)

type (mad_map_struct) map, m
type (mad_energy_struct) energy
type (taylor_struct) taylor(:)

real(rp) beta, dbeta

integer i, j, k, n, nt

! Convert to Bmad phase space coords

taylor%ref = 0
beta = energy%beta
dbeta = 1 / (energy%gamma**2 * energy%total)

m = map

m%k(5) = m%k(5) * beta
m%k(6) = m%k(6) / beta

do i = 1, 6

  do j = 1, 6

    do k = 1, 6
      if (i == 5 .or. j == 6 .or. k == 6) m%t(i,j,k) = m%t(i,j,k) * beta
      if (i == 6 .or. j == 6) m%t(i,j,k) = m%t(i,j,k) / beta
    enddo

    if (i == 5 .or. j == 6) then
      m%r(i,j) = m%r(i,j) * beta
      m%t(i,j,j) = m%t(i,j,j) + m%r(i,j) * dbeta
    endif

    if (i == 6 .or. j == 5) then
      m%r(i,j) = m%r(i,j) / beta
      do k = 1, 6
        m%t(i,j,k) = m%t(i,j,k) * beta
        m%t(i,j,k) = m%t(i,j,k) - m%r(i,j) * m%r(i,k) * dbeta / beta 
      enddo
    endif  

  enddo

enddo

! Count terms and allocate

do i = 1, 6

  nt = 0

  if (m%k(i) /= 0) nt = nt + 1
  do j = 1, 6
    if (m%r(i,j) /= 0) nt = nt + 1
    do k = j, 6
      if (m%t(i,j,k) /= 0) nt = nt + 1
    enddo
  enddo

  if (associated(taylor(i)%term)) then
    if (size(taylor(i)%term) /= nt) then
      deallocate (taylor(i)%term)
      allocate (taylor(i)%term(nt))
    endif
  else
    allocate (taylor(i)%term(nt))
  endif

  do n = 1, nt
    taylor(i)%term(n)%expn = 0
  enddo

enddo

! transfer map to taylor

do i = 1, 6

  nt = 0

  if (m%k(i) /= 0) then
    nt = nt + 1
    taylor(i)%term(nt)%coef = m%k(i)
  endif

  do j = 1, 6

    if (m%r(i,j) /= 0) then
      nt = nt + 1
      taylor(i)%term(nt)%coef = m%r(i,j)
      taylor(i)%term(nt)%expn(j) = 1
    endif

    do k = j, 6

      if (m%t(i,j,k) /= 0) then
        nt = nt + 1
        if (k == j) then
          taylor(i)%term(nt)%coef = m%t(i,j,k)
          taylor(i)%term(nt)%expn(j) = 2
        else
          taylor(i)%term(nt)%coef = 2 * m%t(i,j,k)
          taylor(i)%term(nt)%expn(j) = 1
          taylor(i)%term(nt)%expn(k) = 1
        endif
      endif

    enddo

  enddo

enddo

end subroutine mad_map_to_taylor

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine taylor_to_mad_map (taylor, energy, map)
!
! Subroutine to convert a Taylor map to a mad order 2 map.
! If any of the Taylor terms have order greater than 2 they are ignored.
!
! Input:
!   taylor(6) -- Taylor_struct: Taylor map.
!   energy    -- mad_energy_struct: Energy numbers.
!
! Output:
!   map -- Mad_map_struct: Order 2 map.
!-

subroutine taylor_to_mad_map (taylor, energy, map)

type (mad_map_struct) map
type (mad_energy_struct) energy
type (taylor_struct)  :: taylor(:)
type (taylor_term_Struct) tt

integer i, j, k, n, sm
real(rp) dBmad_dMAD, d2Bmad_dMAD2

character(20), parameter :: r_name = 'taylor_to_mad_map'

!

map%k = 0
map%r = 0
map%t = 0

do i = 1, 6
  n_loop: do n = 1, size(taylor(i)%term)

    tt = taylor(i)%term(n)
    sm = sum(tt%expn)
    select case (sm)

    case (0)
      map%k(i) = tt%coef

    case (1)
      j = maxloc (tt%expn, 1)      
      map%r(i,j) = tt%coef

    case (2)
      j = maxloc (tt%expn, 1)
      if (tt%expn(j) == 2) then
        map%t(i,j,j) = tt%coef
      else
        do k = j+1, 6
          if (tt%expn(k) == 1) then
            map%t(i,j,k) = tt%coef / 2
            map%t(i,k,j) = tt%coef / 2
            cycle n_loop
          endif
          call out_io (s_fatal$, r_name, 'INTERNAL ERROR')
          if (global_com%exit_on_error) call err_exit  ! should not be here
        enddo
      endif

    end select

  enddo n_loop
enddo

! Convert to Bmad phase space coords
! dBmad_dMAD = dE/dP = dct/dz

dBmad_dMAD = 1 / energy%beta
d2Bmad_dMAD2 = -1 / (energy%gamma**2 * energy%p0c)

do i = 1, 6
  do j = 1, 6

    do k = 1, 6
      if (j == 5 .or. j == 6) map%t(i,j,k) = map%t(i,j,k) * dBmad_dMAD
      if (k == 5 .or. k == 6) map%t(i,j,k) = map%t(i,j,k) * dBmad_dMAD
    enddo

    if (j == 5 .or. j == 6) then
      map%r(i,j) = map%r(i,j) * dBmad_dMAD
      map%t(i,j,j) = map%t(i,j,j) + map%r(i,j) * d2Bmad_dMAD2
    endif

    if (i == 5 .or. i == 6) then
      map%r(i,j) = map%r(i,j) / dBmad_dMAD
      do k = 1, 6
        map%t(i,j,k) = map%t(i,j,k) / dBmad_dMAD
        map%t(i,j,k) = map%t(i,j,k) - map%r(i,j) * map%r(i,k) * d2Bmad_dMAD2 / dBmad_dMAD 
      enddo
    endif  
    
  enddo

  if (i == 5 .or. i == 6) map%k(i) = map%k(i) / dBmad_dMAD

enddo

end subroutine taylor_to_mad_map 

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine make_unit_mad_map (map)
!
! Subroutine to initialize a 2nd order transport map to unity.
!
! Input:
!   map -- Mad_map_struct: 2nd order transport map.
!
! Output:
!   map -- Mad_map_struct: Unity 2nd order map.
!-

subroutine make_unit_mad_map (map)

implicit none

type (mad_map_struct) map

integer i

! Make a unit map

map%k = 0
map%r = 0
map%t = 0

forall (i = 1:6) map%r(i,i) = 1

end subroutine make_unit_mad_map

end module
             
