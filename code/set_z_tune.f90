!+
! Subroutine set_z_tune (ring)
!
! Subroutine to set the longitudinal tune by scalling the RF voltages
! in the RF cavities.
!
! Modules Needed:
!   use bmad
!
! Input:
!   ring   -- Ring_struct:
!     %z%tune  -- Longitudinal tune in radians (must be negative). 
!
! Output:
!   ring
!     %ele_(i_rf)%value(voltage$) -- Voltage on the cavity
!
! Notes: 
!   1) The calculation assumes that Q_z << 1.
!   2) By convention a positive tune signifies a clockwise rotation 
!      in phase space so that the transverse tunes are positive. This means 
!      the longitudinal tune is negative above transition.
!-

#include "CESR_platform.inc"

subroutine set_z_tune (ring)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ring_struct), target :: ring
  type (ele_struct), pointer :: ele
  type (modes_struct) mode
  type (coord_struct)  c0

  real(rp) z_tune_wanted, volt_total, r_volt
  real(rp) coef_tot, volt, E0

  integer i, j, k, ix, n_rf, ix_rf(100), ix_attrib(100)

  logical found_control, rf_is_on

! Error detec and init.

  if (ring%z%tune > 0) then
    print *, 'WARNING FROM SET_Z_TUNE: RING%Z%TUNE IS POSITIVE!'
    print *, '     I AM ASSUMING THIS IS INCORRECT AND AM SWITCHING THE SIGN.'
    ring%z%tune = -ring%z%tune
  endif

  z_tune_wanted = ring%z%tune

! Make a list of controllers for the voltage of the RFcavities.
! The list is:
!   1) RFcavities that are not super_slaves and do not have their voltage
!      controlled by an overlay.
!   2) overlays that control the voltage of an RFcavity


  E0 = ring%param%beam_energy

  n_rf = 0
  coef_tot = 0
  rf_is_on = .false.

  do i = 1, ring%n_ele_max

    ele => ring%ele_(i)

    if (ele%key == rfcavity$) then

      if (ele%control_type == super_slave$) cycle 

      do j = ele%ic1_lord, ele%ic2_lord ! check any overlays.
        ix = ring%ic_(j)
        if (ring%control_(ix)%ix_attrib == voltage$) cycle
      enddo

      if (ele%value(rf_frequency$) /= 0) rf_is_on = .true.

      n_rf = n_rf + 1
      ix_rf(n_rf) = i
      coef_tot = coef_tot + twopi * cos(twopi*ele%value(phi0$)) * &
                                ele%value(rf_frequency$) / (c_light * E0)
      ix_attrib(n_rf) = voltage$

    endif

    if (ele%key == overlay$) then
      found_control = .false.
      do j = ele%ix1_slave, ele%ix2_slave
        ix = ring%control_(j)%ix_slave
        if (ring%ele_(ix)%key == rfcavity$ .and. &
                          ring%control_(j)%ix_attrib == voltage$) then
          if (.not. found_control) n_rf = n_rf + 1
          found_control = .true.
          coef_tot = coef_tot + ring%control_(j)%coef * twopi * &
                   cos(twopi*ring%ele_(ix)%value(phi0$)) * &
                   ring%ele_(ix)%value(rf_frequency$) / (c_light * E0)
          k = ele%ix_value
          ix_attrib(n_rf) = k
        else
          if (found_control) then
            print *, 'WARNING FROM SET_Z_TUNE: FOUND OVERLAY THAT DOES NOT'
            print *, '        PURELY CONTROL RF VOLTAGE: ', ele%name
          endif
        endif
      enddo
    endif

  enddo

!

  if (.not. rf_is_on) then
    print *, 'ERROR IN SET_Z_TUNE: RF_FREQUENCY ATTRIBUTE NOT SET'
    print *, '      FOR ANY RF CAVITIES.'
    print *, '      Z TUNE WILL NOT BE SET.'
    return
  endif


! If the voltage is near zero then start from scratch.
! This is only approximate.

  call calc_z_tune (ring)
  if (abs(ring%z%tune) < 0.001 .or. z_tune_wanted == 0) then
    volt = -z_tune_wanted**2 / (ring%param%t1_mat6(5,6) * coef_tot)
    do i = 1, n_rf
      ring%ele_(ix_rf(i))%value(ix_attrib(i)) = volt
      call ring_make_mat6 (ring, ix_rf(i))
    enddo
  endif

! now set cavity voltage to get the correct tune

  do k = 1, 10

    call calc_z_tune (ring)
    if (abs(ring%z%tune - z_tune_wanted) < 0.001) return

    r_volt = (z_tune_wanted / ring%z%tune)**2 

    do i = 1, n_rf
      ring%ele_(ix_rf(i))%value(ix_attrib(i)) = &
                      ring%ele_(ix_rf(i))%value(ix_attrib(i)) * r_volt
      call ring_make_mat6 (ring, ix_rf(i))
    enddo

    call calc_z_tune (ring)

  enddo

! 

  print *, 'ERROR IN SET_Z_TUNE: I CANNOT SET THE TUNE TO THE CORRECT VALUE.'
  print *, '      VALUE WANTED:   ', z_tune_wanted
  print *, '      VALUE OPTAINED: ', ring%z%tune
  call err_exit

end subroutine

