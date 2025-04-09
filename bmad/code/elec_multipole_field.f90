!+
! Subroutine elec_multipole_field (a, b, n, coord, Ex, Ey, dE, compute_dE)
!
! Subroutine to put in the field due to an electric_multipole.
!
! Input:
!   a       -- Real(rp): Multipole skew component.
!   b       -- Real(rp): Multipole normal component.
!   n       -- Real(rp): Multipole order.
!   coord   -- Coord_struct:
!
! Output:
!   Ex          -- Real(rp): X field component
!   Ey          -- Real(rp): Y field component.
!   dE(2,2)     -- Real(rp), optional: Field derivatives: dfield(x,y)/d(x,y).
!   compute_dE  -- logical, optional: If False, do not compute the field derivatives 
!                     even if dE is present. Default is True.
!-

subroutine elec_multipole_field (a, b, n, coord, Ex, Ey, dE, compute_dE)

use bmad_routine_interface, dummy => elec_multipole_field

implicit none

type (coord_struct)  coord

real(rp) a, b, x, y
real(rp), optional :: dE(2,2)
real(rp) Ex, Ey, f

integer n, m, n1
logical, optional :: compute_dE
logical compute

! Init

Ex = 0
Ey = 0
compute = (present(dE) .and. logic_option(.true., compute_dE))

if (compute) dE = 0

! simple case

if (a == 0 .and. b == 0) return

! normal case
! Note that c_multi can be + or -

x = coord%vec(1)
y = coord%vec(3)

do m = 0, n, 2
  f = c_multi(n, m, .true.) * mexp(x, n-m) * mexp(y, m)
  Ex = Ex - b * f
  Ey = Ey - a * f
enddo

do m = 1, n, 2
  f = c_multi(n, m, .true.) * mexp(x, n-m) * mexp(y, m)
  Ex = Ex + a * f
  Ey = Ey - b * f
enddo

! dE calc

if (compute) then

  n1 = n - 1
  
  do m = 0, n1, 2
    f = n * c_multi(n1, m, .true.) * mexp(x, n1-m) * mexp(y, m)
    dE(1,1) = dE(1,1) - b * f
    dE(2,1) = dE(2,1) - a * f

    dE(1,2) = dE(1,2) - a * f
    dE(2,2) = dE(2,2) + b * f
  enddo


  do m = 1, n1, 2
    f = n * c_multi(n1, m, .true.) * mexp(x, n1-m) * mexp(y, m)
    dE(1,2) = dE(1,2) - b * f
    dE(2,2) = dE(2,2) - a * f

    dE(1,1) = dE(1,1) + a * f
    dE(2,1) = dE(2,1) - b * f
  enddo

endif

end subroutine elec_multipole_field
