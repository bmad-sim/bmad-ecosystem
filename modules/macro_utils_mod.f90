module macro_utils_mod

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine calc_bunch_emittance (bunch, ele, x_norm_emit, y_norm_emit)
! 
! Subroutine to calculate emittance of a bunch
!
! Modules needed:
!   use macro_utils_mod
!
! Input:
!   bunch   -- bunch_struct: bunch for which emittance should be calculated.
!   ele     -- ele_struct: Element at which to calculate emittance.
!
! Output:   
!   x_norm_emit  -- Real: Normalized Horizontal Emittance of bunch at ele.
!   y_norm_emit  -- Real: Normalized Vertical  Emittance of bunch at ele.
!-

subroutine calc_bunch_emittance (bunch, ele, x_norm_emit, y_norm_emit)

use bmad
use macroparticle_mod

implicit none

type (bunch_struct), target :: bunch
type(ele_struct) ele
type (macro_struct), pointer :: macro(:), mp

integer i, j
real(rp) x_norm_emit, bunch_x, bunch_xp, avg_x2, avg_xxp, avg_xp2
real(rp) y_norm_emit, bunch_y, bunch_yp, avg_y2, avg_yyp, avg_yp2
real(rp) charge, total_charge, avg_energy
real(rp) x, xp, s11, s12, s22, s66, x_eta, x_etap
real(rp) y, yp, s33, s34, s44, y_eta, y_etap
real(rp) zp, s16, s26, s36, s46

! initialize 

bunch_x      = 0
bunch_xp     = 0
bunch_y      = 0
bunch_yp     = 0
total_charge = 0
avg_energy   = 0

x_eta        = ele%x%eta
x_etap       = ele%x%etap
y_eta        = ele%y%eta
y_etap       = ele%y%etap

avg_x2       = 0
avg_xxp      = 0
avg_xp2      = 0
avg_y2       = 0
avg_yyp      = 0
avg_yp2      = 0

! first calculate average x/x'and y/y' of bunch weighted by charge

do i=1, size(bunch%slice)
  macro => bunch%slice(i)%macro(:)
  bunch_x  = bunch_x  + sum (macro%charge * macro%r%vec(1))
  bunch_xp = bunch_xp + sum (macro%charge * macro%r%vec(2))
  bunch_y  = bunch_y  + sum (macro%charge * macro%r%vec(3))
  bunch_yp = bunch_yp + sum (macro%charge * macro%r%vec(4))
  total_charge = total_charge + sum (macro%charge)
enddo

bunch_x  = bunch_x  / total_charge
bunch_xp = bunch_xp / total_charge
bunch_y  = bunch_y  / total_charge
bunch_yp = bunch_yp / total_charge

! now calculate emittance

do i=1, size(bunch%slice)

  do j = 1, size(bunch%slice(i)%macro)

    mp => bunch%slice(i)%macro(j)
    zp = mp%r%vec(6)
    charge = mp%charge

    s11 = mp%sigma(s11$)
    s12 = mp%sigma(s12$)
    s22 = mp%sigma(s22$)
    s33 = mp%sigma(s33$)
    s34 = mp%sigma(s34$)
    s44 = mp%sigma(s44$)
    s66 = mp%sigma(s66$)
    s16 = mp%sigma(s16$)
    s26 = mp%sigma(s26$)
    s36 = mp%sigma(s36$)
    s46 = mp%sigma(s46$)

    x  = mp%r%vec(1) - bunch_x  - x_eta *  zp
    xp = mp%r%vec(2) - bunch_xp - x_etap * zp
    y  = mp%r%vec(3) - bunch_y  - y_eta *  zp
    yp = mp%r%vec(4) - bunch_yp - y_etap * zp

    avg_x2  = avg_x2 + charge * &
             (s11 + x_eta**2 * s66 + x**2 - 2 * s16*x_eta)
    avg_xp2 = avg_xp2+ charge * &
             (s22 + x_etap**2 * s66 + xp**2 - 2 * s26*x_etap)
    avg_xxp = avg_xxp+ charge * &
             (s12 + x_eta*x_etap*s66 + x*xp - s16*x_etap - s26*x_eta)
    avg_y2  = avg_y2 + charge * &
             (s33 + y_eta**2 * s66 + y**2 - 2 * s36*y_eta)
    avg_yp2 = avg_yp2 + charge * &
             (s44 + y_etap**2 * s66 + yp**2 - 2 * s46*y_etap)
    avg_yyp = avg_yyp + charge * &
             (s34 + y_eta*y_etap*s66 + y*yp - s36*y_etap - s46*y_eta)
    avg_energy = avg_energy + charge * (1+mp%r%vec(6)) * ele%value(beam_energy$)

  end do
end do

avg_x2  = avg_x2  / total_charge
avg_xp2 = avg_xp2 / total_charge
avg_xxp = avg_xxp / total_charge
avg_y2  = avg_y2  / total_charge
avg_yp2 = avg_yp2 / total_charge
avg_yyp = avg_yyp / total_charge
avg_energy = avg_energy / total_charge

x_norm_emit  = (avg_energy/m_electron) * sqrt(avg_x2*avg_xp2 - avg_xxp**2)
y_norm_emit  = (avg_energy/m_electron) * sqrt(avg_y2*avg_yp2 - avg_yyp**2)

end subroutine

end module
