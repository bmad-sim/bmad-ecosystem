!+
! Subroutine SET_Z_TUNE (RING)
!
! Subroutine to set the longitudinal tune by setting the RF voltages
! in the RF cavities.
!
! Modules Needed:
!   use bmad
!
! Input:
!   RING   -- Ring_struct:
!     %Z%TUNE  -- Longitudinal tune in radians. 
!
! Output:
!   RING
!     %ELE_(I_RF)%VALUE(VOLT$) -- Voltage on the cavity
!
! Notes: 
!   0) The calculation assumes that Q_z << 1.
!   2) If the RF wavelength has not been set (needed for the voltage) then
!      the RF harmonic number is set to 1 and the wavelength set to the ring 
!      circumference.
!   3) By convention a positive tune signifies a clockwise rotation 
!      so that the transverse tunes are positive. This means the
!      longitudinal tune is negative above transition.
!-

!$Id$
!$Log$
!Revision 1.3  2002/02/23 20:32:25  dcs
!Double/Single Real toggle added
!
!Revision 1.2  2001/09/27 18:31:58  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine set_z_tune (ring)

  use bmad

  implicit none

  type (ring_struct), target :: ring
  type (ele_struct), pointer :: ele
  type (modes_struct) mode
  type (coord_struct)  c0

  real(rdef) z_tune, volt_total, delta_volt

  integer i, n_rf, ix_rf(100)


! error detec

  if (ring%z%tune == 0) then
    type *, 'WARNING FROM SET_Z_TUNE: RING%Z%TUNE = 0!'
  endif
      
  if (ring%z%tune .gt. 0) then
    type *, 'WARNING FROM SET_Z_TUNE: RING%Z%TUNE IS POSITIVE!'
    type *, '     I AM ASSUMING THIS IS INCORRECT AND AM SWITCHING THE SIGN'
    ring%z%tune = -ring%z%tune
  endif

  z_tune = ring%z%tune

  call calc_z_tune( ring)


! find rf cavities and calculate the transfer matrices between them

  n_rf = 0
  do i = 1, ring.n_ele_ring
    if (ring%ele_(i)%key == rfcavity$) then
      n_rf = n_rf + 1
      ix_rf(n_rf) = i
      if (ring%ele_(i)%value(rf_wavelength$) == 0) then
        type *, 'ERROR IN SET_Z_TUNE: RF_WAVELENGTH  ATTRIBUTE NOT SET'
        type *, '      FOR: ', ring%ele_(i)%name
        type *, '      WILL SET RF_WAVELENGTH FOR HAMONIC NUMBER = 1'
        ring%ele_(i)%value(harmon$) = 1
        ring%ele_(i)%value(rf_wavelength$) = ring.param.total_length
      endif
    endif
  enddo                         

! now set cavity voltage to get the correct tune

  volt_total = 0.
  do i = 1, n_rf
    ele => ring%ele_(ix_rf(i))
    volt_total = ele%value(volt$) + volt_total
  enddo

  delta_volt = volt_total * ((z_tune / ring%z%tune)**2 -1.)

  do i = 1, n_rf
    ring%ele_(ix_rf(i))%value(volt$) = ring%ele_(ix_rf(i))%value(volt$) &
        + delta_volt/n_rf
    call make_mat6 (ring%ele_(ix_rf(i)), ring%param, c0, c0)
  enddo

  call calc_z_tune( ring)
 
end subroutine
