!+
! Subroutine multipole_ab_scale (ele, particle, a, b)
!
! Subroutine to scale ab_multipole values by the strength of an element.
!
! Modules needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!   ele      -- Ele_struct: Element.
!     %value() -- ab_multipole values.
!   particle -- Integer: Particle species (+1 or -1).
!
! Output:
!   a(0:n_pole_maxx) -- Real: Array of scalled multipole values.
!   b(0:n_pole_maxx) -- Real: Array of scalled multipole values.
!-

!$Id$
!$Log$
!Revision 1.2  2001/09/27 18:31:54  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine multipole_ab_scale (ele, particle, a, b)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ele_struct) ele

  real const, radius, factor, a(0:n_pole_maxx), b(0:n_pole_maxx)

  integer ref_exp, n, particle

! transfer values to a and b vecs

  a = ele%value(ix1_m$:ix2_m$-1:2)
  b = ele%value(ix1_m$+1:ix2_m$:2)

! ab_multipole does not scale

  if (ele%key == ab_multipole$) return

! all zero then we do not need to scale.

  if (all(a == 0) .and. all(b == 0)) return

! flip sign for electrons with a separator.

  if (ele%key == elseparator$ .and. particle == -1) then
    a = -a
    b = -b
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
      return
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
    
  case default                                  
    type *, 'ERROR IN MULTIPOLE_AB_SCALE: ELEMENT NOT A AB_MULTIPOLE, QUAD, ETC.'
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
