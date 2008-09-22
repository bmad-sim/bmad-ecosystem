!+
!   SUBROUTINE LUMINOSITY_CALC (ELE, coord, param, n_ok, LUM)
!
! routine to compute luminosity for a distribution of particles
! with coordinates COORD, colliding with a gaussian bunch. Assume
! that the weak beam has the same current as the strong beam
!
! INPUT: ELE : ele_struct, beambeam element
!        COORD_ : coord_struct, distribution of weak beam particles
!        param : Param struct
!        n_ok : number of particles
!
! OUTPUT: LUM: real, luminosity (/sec/m^2)
!
! Mod/Commons:
!
! Calls      :
!
! Author     :
!
! Modified   :
!-
!........................................................................
!
! $Id$
!
! $Log$
! Revision 1.2  2007/01/30 16:14:31  dcs
! merged with branch_bmad_1.
!
! Revision 1.1.1.1.2.1  2006/12/22 20:30:42  dcs
! conversion compiles.
!
! Revision 1.1.1.1  2005/06/14 14:59:02  cesrulib
! Beam Simulation Code
!
!
!........................................................................
!
#include "CESR_platform.inc"
  subroutine luminosity_calc (ele, coord, param, n_ok, lum)

  use bmad_struct
  use bmad_interface
  use track1_mod
  use bookkeeper_mod

  implicit none

  type(ele_struct) ele
  type(coord_struct), allocatable :: coord(:)
  type(lat_param_struct) param

  real(rp) z_slice(100), sig_x0, sig_y0, sig_x, sig_y
  real(rp) s_pos, s_pos_old
  real(rp) vec(6)
  real(rp) lum, f
  real(rp) x,y
  real(rp) beta, kx, ky, coef

  integer i, n_slice, j
  integer psize
  integer n_ok

  if (ele%value(charge$) == 0 .or. param%n_part == 0) return
  call attribute_bookkeeper (ele, param)

  lum = 0.
  psize = size(coord) - 1
  do i =1,n_ok

    vec(1:6) = coord(i)%vec(1:6)
    call offset_particle (ele, param, coord(i), set$)

    sig_x0 = ele%value(sig_x$)
    sig_y0 = ele%value(sig_y$)
    n_slice = max(1, nint(ele%value(n_slice$)))
    call bbi_slice_calc (n_slice, ele%value(sig_z$), z_slice)
    s_pos = 0    ! end at the ip
    do j = 1, n_slice
      s_pos_old = s_pos
      s_pos = (vec(5) + z_slice(j)) / 2
      vec(1) = vec(1) + vec(2) * (s_pos - s_pos_old)
      vec(3) = vec(3) + vec(4) * (s_pos - s_pos_old)
      if (ele%a%beta == 0) then
        sig_x = sig_x0
        sig_y = sig_y0
      else
        beta = ele%a%beta - 2 * ele%a%alpha * s_pos + ele%a%gamma * s_pos**2
        sig_x = sig_x0 * sqrt(beta / ele%a%beta)
        beta = ele%b%beta - 2 * ele%b%alpha * s_pos + ele%b%gamma * s_pos**2
        sig_y = sig_y0 * sqrt(beta / ele%b%beta)
      endif

      x = vec(1)/sig_x
      y = vec(3)/sig_y
      f =  exp(-x**2/2.-y**2/2.)
      lum = lum + f/sig_x/sig_y

      call bbi_kick (vec(1)/sig_x, vec(3)/sig_y, sig_y/sig_x, kx, ky)
      coef = ele%value(bbi_const$) / (n_slice * (1 + vec(6)))
      vec(2) = vec(2) + kx * coef
      vec(4) = vec(4) + ky * coef
    enddo
  end do
  lum = lum * (param%n_part)**2/n_slice/psize/twopi *c_light/ param%total_length
 return
end subroutine luminosity_calc


