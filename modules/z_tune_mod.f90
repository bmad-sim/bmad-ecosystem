module z_tune_mod

use bmad_interface

type (lat_struct), private, pointer :: lat_com
type (all_pointer_struct) :: voltage_control(100)
real(rp), private :: volt0(100), z_tune_wanted
integer, private :: ix_rf(100), n_rf

contains

!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------
!+
! Subroutine set_z_tune (lat, z_tune, ok)
!
! Subroutine to set the longitudinal tune by scalling the RF voltages
! in the RF cavities. Note: RF cavity elements that are set OFF will not
! have their voltages varied.
!
! Input:
!   lat    -- lat_struct:
!   z_tune -- Real(rp), optional: Longitudinal tune in radians (must be negative). 
!               If not present, lat%z%tune will be used instead.
!
! Output:
!   lat
!     %ele(i_rf)%value(voltage$) -- Voltage on the cavity.
!     %z%tune                    -- equal to z_tune.
!   ok                           -- logical.  If present, returns true or false if set was
!                                             successful.  If not present, set_z_tune will
!                                             bomb if tune could not be set.
!
! Notes:
!   1) The calculation assumes that Q_z < 1.
!   2) By convention a positive tune signifies a clockwise rotation 
!      in phase space so that the transverse tunes are positive. This means 
!      the longitudinal tune is negative above transition.
!-

subroutine set_z_tune (lat, z_tune, ok)

use expression_mod, only: linear_coef

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele, ele2, lord
type (control_struct), pointer :: ctl

real(rp) dQz_max
real(rp) coef_tot, volt, E0, phase, dz_tune0, coef0, coef, dz_tune
real(rp), optional :: z_tune
logical, optional :: ok

integer i, j, k, ix, is
integer :: loop_max = 10

logical found_control, rf_is_on, err_flag

character(16), parameter :: r_name = 'set_z_tune'

! Error detec and init.

dQz_max = 0.0001

if (present (z_tune)) lat%z%tune = z_tune

if (present (ok)) ok = .true.

if (lat%z%tune > 0) then
  call out_io (s_warn$, r_name, 'LAT%Z%TUNE IS POSITIVE!', &
                  'I AM ASSUMING THIS IS INCORRECT AND AM SWITCHING THE SIGN.')
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

  ! RFcavity element

  if (ele%key == rfcavity$) then

    if (ele%slave_status == super_slave$) cycle 
    if (.not. ele%is_on) cycle
    if (ele%value(rf_frequency$) == 0) cycle

    do j = 1, ele%n_lord ! check any overlays.
      lord => pointer_to_lord (ele, j, ctl)
      if (lord%key == overlay$ .and. ctl%ix_attrib == voltage$) cycle
    enddo

    n_rf = n_rf + 1
    ix_rf(n_rf) = i
    phase = twopi * (ele%value(phi0$) + ele%value(phi0_multipass$))
    coef_tot = coef_tot + twopi * cos(phase) * ele%value(rf_frequency$) / (c_light * E0)
    voltage_control(n_rf)%r => ele%value(voltage$)

  endif

  ! Overlay element

  if (ele%key == overlay$) then
    found_control = .false.
    do is = 1, ele%n_slave
      ele2 => pointer_to_slave(ele, is, ctl)

      if (found_control .and. &
            (ele2%key /= rfcavity$ .or. ctl%ix_attrib /= voltage$)) then
        call out_io (s_error$, r_name, 'FOUND OVERLAY THAT DOES NOT', &
                                       'PURELY CONTROL RF VOLTAGE: ' // ele%name)
        cycle
      endif

      if (ele%value(rf_frequency$) == 0) cycle
      if (.not. ele2%is_on) cycle

      if (.not. found_control) n_rf = n_rf + 1
      found_control = .true.
      phase = twopi * (ele2%value(phi0$) + ele2%value(phi0_multipass$))
      coef_tot = coef_tot + linear_coef(ctl%stack, err_flag) * twopi * &
               cos(phase) * ele2%value(rf_frequency$) / (c_light * E0)
      voltage_control(n_rf)%r => ele%control%var(1)%value
    enddo
  endif

enddo

!

if (coef_tot == 0) then
  call out_io (s_error$, r_name, 'CANNOT FIND ANY RFCAVITY ELEMENTS WHICH ARE:', &
                                 '      1) ON,  AND', &
                                 '      2) HAVE A FINETE RF_FREQUENCY!', &
                                 '      THE Z TUNE WILL NOT BE SET.')
  return
endif


! If the voltage is near zero then start from scratch.
! This is only approximate.

call calc_z_tune (lat)

if (lat%z%tune == 0) then
  call out_io (s_error$, r_name, 'CALCULATED Z TUNE IS ZERO. CANNOT SET THE TUNE.')
  if (global_com%exit_on_error) call err_exit
  return
endif

if (abs(lat%z%tune) < dQz_max .or. z_tune_wanted == 0) then
  volt = -z_tune_wanted**2 / (lat%param%t1_with_RF(5,6) * coef_tot)
  do i = 1, n_rf
    ele => lat%ele(ix_rf(i))
    voltage_control(i)%r = volt
    call set_flags_for_changed_attribute (ele, voltage_control(i)%r)
    call lat_make_mat6 (lat, ix_rf(i))
  enddo
  call calc_z_tune (lat)
endif

! record

do i = 1, n_rf
  ele => lat%ele(ix_rf(i))
  volt0(i) = voltage_control(i)%r
enddo

dz_tune = lat%z%tune - z_tune_wanted
dz_tune0 = dz_tune
lat_com => lat
coef = 1

! now set cavity voltage to get the correct tune

do k = 1, loop_max

  if (abs(dz_tune) < dQz_max) return

  if (dz_tune * dz_tune0 < 0) exit  ! Have bracketed solution

  coef0 = coef
  coef = coef * (z_tune_wanted / (dz_tune + z_tune_wanted))**2 

  dz_tune = dz_tune_func(coef)

  if (k == loop_max) then
    call out_io (s_error$, r_name, 'I CANNOT SET THE TUNE TO THE CORRECT VALUE.', &
                                   '      VALUE WANTED:   \f12.3\ ', '      VALUE OBTAINED: \f12.3\ ', &
                                   r_array = [z_tune_wanted, lat%z%tune])
    if (present(ok)) then
      ok = .false.
    else
      if (global_com%exit_on_error) call err_exit
    endif
    lat%z%stable = .false.
    return
  endif
enddo

lat%z%stable = .true.

! Have bracketed index
! Superfluous? coef = super_zbrent (dz_tune_func, min(coef0, coef), max(coef0, coef), dQz_max)

end subroutine set_z_tune

!-------------------------------------------------------------------------------------

function dz_tune_func (coef) result (dz_tune)

type (ele_struct), pointer :: ele
real(rp), intent(in) :: coef
real(rp) dz_tune
integer i

!

do i = 1, n_rf
  ele => lat_com%ele(ix_rf(i))
  voltage_control(i)%r = volt0(i) * coef
  call set_flags_for_changed_attribute (ele, voltage_control(i)%r)
  call lat_make_mat6 (lat_com, ix_rf(i))
enddo

call calc_z_tune (lat_com)
dz_tune = lat_com%z%tune - z_tune_wanted


end function

end module

