#include "CESR_platform.inc"

module multipole_mod

  use bmad_struct
  use matrix_mod

contains

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
!   an(0:n_pole_maxx) -- Real(rdef): Skew multipole component.
!   bn(0:n_pole_maxx) -- Real(rdef): Normal multipole component.
!
! Output:
!   knl(0:n_pole_maxx) -- Real(rdef): Multitude magnatude.
!   tn(0:n_pole_maxx)  -- Real(rdef): Multipole angle.
!-

subroutine multipole_ab_to_kt (an, bn, knl, tn)

  implicit none

  real(rdef) an(0:), bn(0:)
  real(rdef) knl(0:), tn(0:)
  real(rdef) n_fact, a, b

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
! Note: tilt(:) does includes ele%value(tilt$).
!
! Modules needed:
!   use bmad
!
! Input:
!   ele          -- Ele_struct: Multipole element.
!   particle     -- Integer: Particle species (+1 or -1).
!   use_ele_tilt -- Logical: If True then include ele%value(tilt$) 
!                     in calculations.
!
! Output:
!   knl(0:)  -- Real(rdef): Vector of strengths, MAD units.
!   tilt(0:) -- Real(rdef): Vector of tilts.
!-

subroutine multipole_ele_to_kt (ele, particle, knl, tilt, use_ele_tilt)

  implicit none

  type (ele_struct)  ele

  real(rdef) knl(0:), tilt(0:), signn, a_n, b_n
  real(rdef) value(n_attrib_maxx), a(0:n_pole_maxx), b(0:n_pole_maxx)

  integer n, particle, n_fact

  logical use_ele_tilt

!

  if (.not. ele%multipoles_on .or. .not. ele%is_on .or. .not. associated(ele%a)) then
    knl = 0
    tilt = 0
    return
  endif

! multipole
                    
  if (ele%key == multipole$) then
    knl  = ele%a
    tilt = ele%b + ele%value(tilt$)
    return
  endif

! ab_multiple, etc...

  call multipole_ele_to_ab (ele, particle, a, b, .false.)
  call multipole_ab_to_kt (a, b, knl, tilt)
  if (use_ele_tilt) tilt = tilt + ele%value(tilt$)

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
!   knl(0:) -- Real(rdef): Multitude magnatude.
!   tn(0:)  -- Real(rdef): Multipole angle.
!
! Output:
!   an(0:) -- Real(rdef): Skew multipole component.
!   bn(0:) -- Real(rdef): Normal multipole component.
!-

subroutine multipole_kt_to_ab (knl, tn, an, bn)

  implicit none

  real(rdef) an(0:), bn(0:)
  real(rdef) knl(0:), tn(0:)
  real(rdef) n_fact, angle, kl

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
!   use_ele_tilt -- Logical: If True then include ele%value(tilt$) 
!                     in calculations.
!
! Output:
!   a(0:n_pole_maxx) -- Real(rdef): Array of scalled multipole values.
!   b(0:n_pole_maxx) -- Real(rdef): Array of scalled multipole values.
!-

subroutine multipole_ele_to_ab (ele, particle, a, b, use_ele_tilt)

  implicit none

  type (ele_struct) ele

  real(rdef) const, radius, factor, a(0:), b(0:)
  real(rdef) an, bn, cos_t, sin_t

  integer ref_exp, n, particle

  logical use_ele_tilt

!

  if (.not. ele%multipoles_on .or. .not. ele%is_on .or. .not. associated(ele%a)) then
    a = 0
    b = 0
    return
  endif

! transfer values to a and b vecs

  a = ele%a
  b = ele%b

  if (ele%key == multipole$) call multipole_kt_to_ab (a, b, a, b)

! all zero then we do not need to scale.

  if (all(a == 0) .and. all(b == 0)) return

! flip sign for electrons or antiprotons with a separator.

  if (ele%key == elseparator$ .and. particle < 0) then
    a = -a
    b = -b
  endif

! use tilt?

  if (use_ele_tilt .and. ele%value(tilt$) /= 0) then
    do n = 0, n_pole_maxx
      if (a(n) /= 0 .or. b(n) /= 0) then
        an = a(n); bn = b(n)
        cos_t = cos((n+1)*ele%value(tilt$))
        sin_t = sin((n+1)*ele%value(tilt$))
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
    const = ele%value(l$) * (ele%value(g$) + ele%value(delta_g$))
    ref_exp = 0

  case (elseparator$, kicker$)
    if (ele%value(hkick$) == 0) then
      const = ele%value(vkick$)
    elseif (ele%value(vkick$) == 0) then
      const = ele%value(hkick$)
    else
      print *, 'ERROR IN MULTIPOLE_ELE_TO_AB: I CANNOT SCALE: ' // ele%name
      print *, '  TO A REFERENCE SEPARATOR WITH HKICK AND VKICK *BOTH* NONZERO'
      call err_exit
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
! Subroutine multipole_kick (knl, tilt, n, coord)
!
! Subroutine to put in the kick due to a multipole.
!
! Modules Needed:
!   use bmad
!                          
! Input:
!   knl   -- Real(rdef): Multipole strength (mad units).
!   tilt  -- Real(rdef): Multipole tilt.
!   n     -- Real(rdef): Multipole order.
!   coord -- Coord_struct:
!     %vec(1) -- X position.
!     %vec(3) -- Y position.
!
! Output:
!   coord -- Coord_struct: 
!     %vec(2) -- X kick.
!     %vec(4) -- Y kick.
!-

subroutine multipole_kick (knl, tilt, n, coord)

  implicit none

  type (coord_struct)  coord

  real(rdef) knl, tilt, x, y, sin_ang, cos_ang
  real(rdef) x_vel, y_vel

  integer n, m

! simple case

  if (knl == 0) return

! normal case

  if (tilt == 0) then
    x = coord%vec(1)
    y = coord%vec(3)
  else
    sin_ang = sin(tilt)
    cos_ang = cos(tilt)
    x =  coord%vec(1) * cos_ang + coord%vec(3) * sin_ang
    y = -coord%vec(1) * sin_ang + coord%vec(3) * cos_ang
  endif

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

end module
