!+
! Subroutine set_z_tune (lat)
!
! Subroutine to set the longitudinal tune by scalling the RF voltages
! in the RF cavities.
!
! Modules Needed:
!   use bmad
!
! Input:
!   lat   -- lat_struct:
!     %z%tune  -- Longitudinal tune in radians (must be negative). 
!
! Output:
!   lat
!     %ele(i_rf)%value(voltage$) -- Voltage on the cavity
!
! Notes: 
!   1) The calculation assumes that Q_z << 1.
!   2) By convention a positive tune signifies a clockwise rotation 
!      in phase space so that the transverse tunes are positive. This means 
!      the longitudinal tune is negative above transition.
!-

#include "CESR_platform.inc"

subroutine set_z_tune (lat)

  use bmad_struct
  use bmad_interface, except => set_z_tune

  implicit none

  type (lat_struct), target :: lat
  type (ele_struct), pointer :: ele, ele2

  real(rp) z_tune_wanted, r_volt
  real(rp) coef_tot, volt, E0, phase

  integer i, j, k, ix, n_rf, ix_rf(100), ix_attrib(100)

  logical found_control, rf_is_on

! Error detec and init.

  if (lat%z%tune > 0) then
    print *, 'WARNING FROM SET_Z_TUNE: LAT%Z%TUNE IS POSITIVE!'
    print *, '     I AM ASSUMING THIS IS INCORRECT AND AM SWITCHING THE SIGN.'
    lat%z%tune = -lat%z%tune
  endif

  z_tune_wanted = lat%z%tune

! Make a list of controllers for the voltage of the RFcavities.
! The list is:
!   1) RFcavities that are not super_slaves and do not have their voltage
!      controlled by an overlay.
!   2) overlays that control the voltage of an RFcavity


  E0 = lat%ele(0)%value(E_TOT$)

  n_rf = 0
  coef_tot = 0
  rf_is_on = .false.

  do i = 1, lat%n_ele_max

    ele => lat%ele(i)

    if (ele%key == rfcavity$) then

      if (ele%control_type == super_slave$) cycle 

      do j = ele%ic1_lord, ele%ic2_lord ! check any overlays.
        ix = lat%ic(j)
        if (lat%control(ix)%ix_attrib == voltage$) cycle
      enddo

      if (ele%value(rf_frequency$) /= 0) rf_is_on = .true.

      n_rf = n_rf + 1
      ix_rf(n_rf) = i
      phase = twopi * (ele%value(phi0$) + ele%value(dphi0$))
      coef_tot = coef_tot + twopi * cos(phase) * &
                                ele%value(rf_frequency$) / (c_light * E0)
      ix_attrib(n_rf) = voltage$

    endif

    if (ele%key == overlay$) then
      found_control = .false.
      do j = ele%ix1_slave, ele%ix2_slave
        ix = lat%control(j)%ix_slave
        ele2 => lat%ele(ix)
        if (ele2%key == rfcavity$ .and. &
                          lat%control(j)%ix_attrib == voltage$) then
          if (.not. found_control) n_rf = n_rf + 1
          found_control = .true.
          phase = twopi * (ele2%value(phi0$) + ele2%value(dphi0$))
          coef_tot = coef_tot + lat%control(j)%coef * twopi * &
                   cos(phase) * ele2%value(rf_frequency$) / (c_light * E0)
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

  call calc_z_tune (lat)
  if (abs(lat%z%tune) < 0.001 .or. z_tune_wanted == 0) then
    volt = -z_tune_wanted**2 / (lat%param%t1_with_RF(5,6) * coef_tot)
    do i = 1, n_rf
      lat%ele(ix_rf(i))%value(ix_attrib(i)) = volt
      call lat_make_mat6 (lat, ix_rf(i))
    enddo
  endif

! now set cavity voltage to get the correct tune

  do k = 1, 10

    call calc_z_tune (lat)
    if (abs(lat%z%tune - z_tune_wanted) < 0.001) return

    r_volt = (z_tune_wanted / lat%z%tune)**2 

    do i = 1, n_rf
      lat%ele(ix_rf(i))%value(ix_attrib(i)) = &
                      lat%ele(ix_rf(i))%value(ix_attrib(i)) * r_volt
      call lat_make_mat6 (lat, ix_rf(i))
    enddo

    call calc_z_tune (lat)

  enddo

! 

  print *, 'ERROR IN SET_Z_TUNE: I CANNOT SET THE TUNE TO THE CORRECT VALUE.'
  print *, '      VALUE WANTED:   ', z_tune_wanted
  print *, '      VALUE OPTAINED: ', lat%z%tune
  call err_exit

end subroutine

