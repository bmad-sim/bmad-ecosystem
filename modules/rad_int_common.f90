!+
!-

!$Id$
!$Log$
!Revision 1.2  2001/09/27 18:32:13  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


module rad_int_common               

  use bmad_struct

  type (ring_struct), pointer :: pring
  type (ele_struct), pointer :: ele0, ele
  type (ele_struct) runt
  type (coord_struct), pointer :: orb0, orb1
  type (coord_struct) orb

! i1_(:), etc. are for diagnostics

  type rad_int_common_struct
    real g_x0, g_y0, ll, k1, s1, eta(4), eta_a(4), eta_b(4)
    real g, g2, g_x, g_y, dg2_x, dg2_y 
    real i1_(n_ele_maxx), i2_(n_ele_maxx), i3_(n_ele_maxx)
    real i4a_(n_ele_maxx), i4b_(n_ele_maxx), i5a_(n_ele_maxx), i5b_(n_ele_maxx)
  end type

  type (rad_int_common_struct) rad_com

!

  interface
    function qromb_rad_int (func, a, b, sum)
      use nrtype
      real(sp), intent(in) :: a, b, sum
      real(sp) :: qromb
      interface
        function func(x)
          use nrtype
          real(sp), dimension(:), intent(in) :: x
          real(sp), dimension(size(x)) :: func
        end function func
      end interface
    end function qromb_rad_int
  end interface

end module
