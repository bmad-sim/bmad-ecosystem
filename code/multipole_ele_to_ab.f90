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

!$Id$
!$Log$
!Revision 1.2  2002/02/23 20:32:20  dcs
!Double/Single Real toggle added
!
!Revision 1.1  2002/01/08 21:48:14  dcs
!Align with VMS version
!

#include "CESR_platform.inc"

subroutine multipole_ele_to_ab (ele, particle, a, b, use_ele_tilt)

  use bmad

  implicit none

  type (ele_struct) ele

  real(rdef) const, radius, factor, a(0:n_pole_maxx), b(0:n_pole_maxx)
  real(rdef) an, bn, cos_t, sin_t

  integer ref_exp, n, particle

  logical use_ele_tilt

!

  if (.not. ele%multipoles_on .or. .not. ele%is_on) then
    a = 0
    b = 0
    return
  endif

! transfer values to a and b vecs

  a = ele%value(ix1_m$:ix2_m$-1:2)
  b = ele%value(ix1_m$+1:ix2_m$:2)

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
    const = ele%value(angle$)
    ref_exp = 0

  case (elseparator$, kicker$)
    if (ele%value(hkick$) == 0) then
      const = ele%value(vkick$)
    elseif (ele%value(vkick$) == 0) then
      const = ele%value(hkick$)
    else
      call warning ('I CANNOT SCALE: ' // ele%name, &
        'TO A REFERENCE SEPARATOR WITH HKICK AND VKICK *BOTH* NONZERO:' // &
         ele%name)
      const = 0
    endif
    ref_exp = 0

  case (quadrupole$, sol_quad$)
    const = ele%value(k1$) * ele%value(l$)
    ref_exp = 1

  case (wiggler$)
    const = 2 * ele%value(l$) / (pi * ele%value(rho$) * ele%value(n_pole$))
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
    type *, 'ERROR IN MULTIPOLE_ELE_TO_AB: ELEMENT NOT A AB_MULTIPOLE, QUAD, ETC.'
    type *, '      ', ele%name
    call err_exit

  end select

! scale multipole values

  do n = 0, n_pole_maxx
    factor = const * radius ** (ref_exp - n)
    a(n) = factor * a(n)
    b(n) = factor * b(n)
  enddo

end subroutine
