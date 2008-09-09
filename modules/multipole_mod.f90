#include "CESR_platform.inc"

module multipole_mod

use bmad_struct
use matrix_mod

contains

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine multipole_init (ele, zero)
!
! Subroutine to allocate memory for the the ele%a_pole and ele%b_pole multipole 
! vectors.
!
! Modules needed:
!   use bmad
!
! Input:
!   zero -- Logical, optional: If present and True then zero the arrays
!             even if they already exist when this routine is called. 
!             Default is False which means that if the arrays already 
!             exist then this routine will do nothing.
!
! Output:
!   ele -- Ele_struct: Element holding the multipoles.
!     %a_pole(0:n_pole_maxx) -- Multipole An array 
!     %b_pole(0:n_pole_maxx) -- Multipole Bn array
!-

subroutine multipole_init (ele, zero)

implicit none

type (ele_struct) ele
logical, optional :: zero

! If %a_pole and %b_pole already exist then zero them if zero argument present 
! and True.

if (associated (ele%a_pole)) then
  if (logic_option(.false., zero)) then
    ele%a_pole = 0
    ele%b_pole = 0
  endif

! If memory not allocated then allocate and zero.

else
  allocate (ele%a_pole(0:n_pole_maxx), ele%b_pole(0:n_pole_maxx))
  ele%a_pole = 0
  ele%b_pole = 0
endif

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine multipole_ab_to_kt (an, bn, knl, tn)
!
! Subroutine to convert ab type multipoles to kt (MAD standard) multipoles.
!
! Modules needed:
!   use bmad
!
! Input:
!   an(0:n_pole_maxx) -- Real(rp): Skew multipole component.
!   bn(0:n_pole_maxx) -- Real(rp): Normal multipole component.
!
! Output:
!   knl(0:n_pole_maxx) -- Real(rp): Multitude magnatude.
!   tn(0:n_pole_maxx)  -- Real(rp): Multipole angle.
!-

subroutine multipole_ab_to_kt (an, bn, knl, tn)

  implicit none

  real(rp) an(0:), bn(0:)
  real(rp) knl(0:), tn(0:)
  real(rp) n_fact, a, b

  integer n

! use a, b as temp values to avoid problems with a call like:
!   call multipole_ab_to_kt (vec1, vec2, vec1, vec2)

  n_fact = 1

  do n = 0, n_pole_maxx

    if (n /= 0) n_fact = n_fact * n

    a = an(n)
    b = bn(n)

    if (a == 0 .and. b == 0) then
      knl(n) = 0
      tn(n) = 0
    else
      knl(n)  = n_fact * sqrt(a**2 + b**2)
      tn(n) = -atan2(a, b) / (n + 1)
    endif

  enddo

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine multipole_ele_to_kt (ele, particle, knl, tilt, use_ele_tilt)
!
! Subroutine to put the multipole components (strength and tilt)
! into 2 vectors along with the appropriate scaling.
! Note: tilt(:) does includes ele%value(tilt_tot$).
!
! Modules needed:
!   use bmad
!
! Input:
!   ele          -- Ele_struct: Multipole element.
!   particle     -- Integer: Particle species (+1 or -1).
!   use_ele_tilt -- Logical: If True then include ele%value(tilt_tot$) 
!                     in calculations.
!
! Output:
!   knl(0:)  -- Real(rp): Vector of strengths, MAD units.
!   tilt(0:) -- Real(rp): Vector of tilts.
!-

subroutine multipole_ele_to_kt (ele, particle, knl, tilt, use_ele_tilt)

  implicit none

  type (ele_struct)  ele

  real(rp) knl(0:), tilt(0:), a(0:n_pole_maxx), b(0:n_pole_maxx)
  integer particle
  logical use_ele_tilt

!

  if (.not. ele%multipoles_on .or. .not. ele%is_on .or. .not. associated(ele%a_pole)) then
    knl = 0
    tilt = 0
    return
  endif

! multipole
                    
  if (ele%key == multipole$) then
    knl  = ele%a_pole
    tilt = ele%b_pole + ele%value(tilt_tot$)
    return
  endif

! ab_multiple, etc...

  call multipole_ele_to_ab (ele, particle, a, b, .false.)
  call multipole_ab_to_kt (a, b, knl, tilt)
  if (use_ele_tilt) tilt = tilt + ele%value(tilt_tot$)

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine multipole_kt_to_ab (knl, tn, an, bn)
!
! Subroutine to convert kt (MAD standard) multipoles to ab type multipoles.
!
! Modules needed:
!   use bmad
!
! Input:
!   knl(0:) -- Real(rp): Multitude magnatude.
!   tn(0:)  -- Real(rp): Multipole angle.
!
! Output:
!   an(0:) -- Real(rp): Skew multipole component.
!   bn(0:) -- Real(rp): Normal multipole component.
!-

subroutine multipole_kt_to_ab (knl, tn, an, bn)

  implicit none

  real(rp) an(0:), bn(0:)
  real(rp) knl(0:), tn(0:)
  real(rp) n_fact, angle, kl

  integer n

!

  n_fact = 1

  do n = lbound(an, 1), ubound(an, 1)

    if (n /= 0) n_fact = n_fact * n

    kl = knl(n) / n_fact

    if (kl == 0) then
      an(n) = 0
      bn(n) = 0
    else
      angle = -tn(n) * (n + 1)
      an(n) = kl * sin(angle)
      bn(n) = kl * cos(angle)
    endif

  enddo

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine multipole_ele_to_ab (ele, particle, a, b, use_ele_tilt)
!                             
! Subroutine to extract the ab multipole values of an element.
! Note: The ab values will be scalled by the strength of the element.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele          -- Ele_struct: Element.
!     %value()     -- ab_multipole values.
!   particle     -- Integer: Particle species (positron$, etc.).
!                     To be used with electrostatic elements.
!   use_ele_tilt -- Logical: If True then include ele%value(tilt_tot$) 
!                     in calculations.
!
! Output:
!   a(0:n_pole_maxx) -- Real(rp): Array of scalled multipole values.
!   b(0:n_pole_maxx) -- Real(rp): Array of scalled multipole values.
!-

subroutine multipole_ele_to_ab (ele, particle, a, b, use_ele_tilt)

  implicit none

  type (ele_struct) ele

  real(rp) const, radius, factor, a(0:), b(0:)
  real(rp) an, bn, cos_t, sin_t

  integer ref_exp, n, particle

  logical use_ele_tilt

!

  if (.not. ele%multipoles_on .or. .not. ele%is_on .or. .not. associated(ele%a_pole)) then
    a = 0
    b = 0
    return
  endif

! transfer values to a and b vecs

  a = ele%a_pole
  b = ele%b_pole

  if (ele%key == multipole$) call multipole_kt_to_ab (a, b, a, b)

! all zero then we do not need to scale.

  if (all(a == 0) .and. all(b == 0)) return

! flip sign for electrons or antiprotons with a separator.

  if (ele%key == elseparator$ .and. particle < 0) then
    a = -a
    b = -b
  endif

! use tilt?

  if (use_ele_tilt .and. ele%value(tilt_tot$) /= 0) then
    do n = 0, n_pole_maxx
      if (a(n) /= 0 .or. b(n) /= 0) then
        an = a(n); bn = b(n)
        cos_t = cos((n+1)*ele%value(tilt_tot$))
        sin_t = sin((n+1)*ele%value(tilt_tot$))
        b(n) =  bn * cos_t + an * sin_t
        a(n) = -bn * sin_t + an * cos_t
      endif
    enddo
  endif

! radius = 0 defaults to radius = 1

  radius = ele%value(radius$)
  if (radius == 0) radius = 1

! normal case...

  select case (ele%key)

  case (sbend$, rbend$)
    const = ele%value(l$) * (ele%value(g$) + ele%value(g_err$))
    ref_exp = 0

  case (elseparator$, kicker$)
    if (ele%value(hkick$) == 0) then
      const = ele%value(vkick$)
    elseif (ele%value(vkick$) == 0) then
      const = ele%value(hkick$)
    else
      const = sqrt(ele%value(hkick$)**2 + ele%value(vkick$)**2)
    endif
    ref_exp = 0

  case (quadrupole$, sol_quad$)
    const = ele%value(k1$) * ele%value(l$)
    ref_exp = 1

  case (wiggler$)
    const = 2 * ele%value(l$) / &
                    (pi * ele%value(rho$) * ele%value(n_pole$))
    ref_exp = 0

  case (solenoid$)
    const = ele%value(ks$) * ele%value(l$)
    ref_exp = 1

  case (sextupole$)
    const = ele%value(k2$) * ele%value(l$)
    ref_exp = 2
 
  case (octupole$)
    const = ele%value(k3$) * ele%value(l$)
    ref_exp = 3
    
  case (ab_multipole$, multipole$) ! multipoles do not scale
    return

  case default                                  
    print *, 'ERROR IN MULTIPOLE_ELE_TO_AB: ELEMENT NOT A AB_MULTIPOLE, QUAD, ETC.'
    print *, '      ', ele%name
    call err_exit

  end select

! scale multipole values

  do n = 0, n_pole_maxx
    factor = const * radius ** (ref_exp - n)
    a(n) = factor * a(n)
    b(n) = factor * b(n)
  enddo

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine multipole_kicks (knl, tilt, coord, ref_orb_offset)
!
! Subroutine to put in the kick due to a multipole.
!
! Modules Needed:
!   use bmad
!                          
! Input:
!   knl(0:)        -- Real(rp): Multipole strengths (mad units).
!   tilt(0:)       -- Real(rp): Multipole tilts.
!   coord          -- Coord_struct:
!     %vec(1)          -- X position.
!     %vec(3)          -- Y position.
!   ref_orb_offset -- Logical, optional: If present and n = 0 then the
!                       multipole simulates a zero length bend with bending
!                       angle knl.
!
! Output:
!   coord -- Coord_struct: 
!     %vec(2) -- X kick.
!     %vec(4) -- Y kick.
!-

subroutine multipole_kicks (knl, tilt, coord, ref_orb_offset)

  implicit none

  type (coord_struct)  coord

  real(rp) knl(0:), tilt(0:)

  integer n

  logical, optional :: ref_orb_offset

  !
  
  do n = 0, n_pole_maxx
    if (knl(n) == 0) cycle
    call multipole_kick (knl(n), tilt(n), n, coord, ref_orb_offset)
  enddo

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine multipole_kick (knl, tilt, n, coord, ref_orb_offset)
!
! Subroutine to put in the kick due to a multipole.
!
! Modules Needed:
!   use bmad
!                          
! Input:
!   knl   -- Real(rp): Multipole strength (mad units).
!   tilt  -- Real(rp): Multipole tilt.
!   n     -- Real(rp): Multipole order.
!   coord -- Coord_struct:
!     %vec(1) -- X position.
!     %vec(3) -- Y position.
!   ref_orb_offset -- Logical, optional: If present and n = 0 then the
!                       multipole simulates a zero length bend with bending
!                       angle knl.
!
! Output:
!   coord -- Coord_struct: 
!     %vec(2) -- X kick.
!     %vec(4) -- Y kick.
!-

subroutine multipole_kick (knl, tilt, n, coord, ref_orb_offset)

  implicit none

  type (coord_struct)  coord

  real(rp) knl, tilt, x, y, sin_ang, cos_ang
  real(rp) x_vel, y_vel

  integer n, m

  logical, optional :: ref_orb_offset

! simple case

  if (knl == 0) return

! normal case

  if (tilt == 0) then
    sin_ang = 0
    cos_ang = 1
    x = coord%vec(1)
    y = coord%vec(3)
  else
    sin_ang = sin(tilt)
    cos_ang = cos(tilt)
    x =  coord%vec(1) * cos_ang + coord%vec(3) * sin_ang
    y = -coord%vec(1) * sin_ang + coord%vec(3) * cos_ang
  endif

! ref_orb_offset with n = 0 means that we are simulating a zero length dipole.

  if (n == 0 .and. present(ref_orb_offset)) then
    coord%vec(2) = coord%vec(2) + knl * cos_ang * coord%vec(6)
    coord%vec(4) = coord%vec(4) + knl * sin_ang * coord%vec(6)
    coord%vec(5) = coord%vec(5) - knl * &
                    (cos_ang * coord%vec(1) + sin_ang * coord%vec(3))
    return
  endif

! normal case

  x_vel = 0
  y_vel = 0

  do m = 0, n, 2
    x_vel = x_vel + knl * c_multi(n, m) * mexp(x, n-m) * mexp(y, m)
  enddo

  do m = 1, n, 2
    y_vel = y_vel + knl * c_multi(n, m) * mexp(x, n-m) * mexp(y, m)
  enddo

  if (tilt == 0) then
    coord%vec(2) = coord%vec(2) + x_vel
    coord%vec(4) = coord%vec(4) + y_vel
  else
    coord%vec(2) = coord%vec(2) + x_vel * cos_ang - y_vel * sin_ang
    coord%vec(4) = coord%vec(4) + x_vel * sin_ang + y_vel * cos_ang
  endif

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine ab_multipole_kick (a, b, n, coord, kx, ky)
!
! Subroutine to put in the kick due to an ab_multipole.
!
! Modules Needed:
!   use bmad
!                          
! Input:
!   a     -- Real(rp): Multipole skew component.
!   b     -- Real(rp): Multipole normal component.
!   n     -- Real(rp): Multipole order.
!   coord -- Coord_struct:
!
! Output:
!   kx -- Real(rp): X kick.
!   ky -- Real(rp): Y kick.
!-

subroutine ab_multipole_kick (a, b, n, coord, kx, ky)

  implicit none

  type (coord_struct)  coord

  real(rp) a, b, x, y
  real(rp) kx, ky, f

  integer n, m

! simple case

  kx = 0
  ky = 0

  if (a == 0 .and. b == 0) return

! normal case

  x = coord%vec(1)
  y = coord%vec(3)

  do m = 0, n, 2
    f = c_multi(n, m, .true.) * mexp(x, n-m) * mexp(y, m)
    kx = kx + b * f
    ky = ky - a * f
  enddo

  do m = 1, n, 2
    f = c_multi(n, m, .true.) * mexp(x, n-m) * mexp(y, m)
    kx = kx + a * f
    ky = ky + b * f
  enddo

end subroutine

end module
