!+
! Subroutine set_z_tune (branch, z_tune, ok, print_err)
!
! Subroutine to set the longitudinal tune by scalling the RF voltages
! in the RF cavities. Note: RF cavity elements that are set OFF will not
! have their voltages varied.
!
! Input:
!   branch      -- branch_struct:
!   z_tune      -- real(rp): Longitudinal tune in radians (must be negative above transition). 
!   print_err   -- logical, optional: Default is True. If False, suppress error messages
!
! Output:
!   branch      -- branch_struct:
!     %ele(i_rf)%value(voltage$) -- Voltage on the cavity.
!     %z%tune                    -- equal to z_tune.
!   ok          -- logical, optional:  If present, returns true or false if set was successful.
!                     If not present, set_z_tune will bomb if tune could not be set.
!
! Notes:
!   1) The calculation assumes that Q_z < 1.
!   2) By convention a positive tune signifies a clockwise rotation 
!      in phase space so that the transverse tunes are positive. This means 
!      the longitudinal tune is negative above transition.
!-

subroutine set_z_tune (branch, z_tune, ok, print_err)

use bmad_interface, dummy => set_z_tune
use expression_mod, only: linear_coef
use super_recipes_mod, only: super_zbrent
use twiss_and_track_mod, only: twiss_and_track

implicit none

type (branch_struct), target :: branch
type (ele_struct), pointer :: ele, ele2, lord
type (control_struct), pointer :: ctl
type (coord_struct), allocatable :: orbit(:)

real(rp) Qz_rel_tol, Qz_abs_tol
real(rp) coef_tot, volt, E0, phase, dz_tune0, coef0, dz_tune, coef
real(rp) :: z_tune
logical, optional :: ok, print_err

integer i, j, k, ix, is, status
integer :: loop_max = 10

logical found_control, err_flag

! common

type (all_pointer_struct) :: voltage_control(100)
real(rp) :: volt0(100)
integer :: ix_rf(100), n_rf

character(16), parameter :: r_name = 'set_z_tune'

! Error detec and init.

Qz_rel_tol = 1d-5
Qz_abs_tol = 1d-7

if (present (ok)) ok = .true.

if (branch%z%tune > 0) then
  call out_io (s_warn$, r_name, 'BRANCH%Z%TUNE IS POSITIVE!', &
                  'I AM ASSUMING THIS IS INCORRECT AND AM SWITCHING THE SIGN.')
  branch%z%tune = -branch%z%tune
endif

! Make a list of controllers for the voltage of the RFcavities.
! The list is:
!   1) RFcavities that are not super_slaves and do not have their voltage
!      controlled by an overlay.
!   2) overlays that control the voltage of an RFcavity

E0 = branch%ele(0)%value(E_TOT$)

n_rf = 0
coef_tot = 0

do i = 1, branch%n_ele_max

  ele => branch%ele(i)

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
    phase = twopi * ele%value(phi0$)
    if (.not. bmad_com%absolute_time_tracking) phase = phase + twopi * ele%value(phi0_multipass$)
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
      phase = twopi * ele%value(phi0$)
      if (.not. bmad_com%absolute_time_tracking) phase = phase + twopi * ele%value(phi0_multipass$)
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

call calc_z_tune (branch)

if (branch%z%tune == 0) then
  call out_io (s_error$, r_name, 'CALCULATED Z TUNE IS ZERO. CANNOT SET THE TUNE.')
  if (global_com%exit_on_error) call err_exit
  return
endif

if (abs(branch%z%tune) < Qz_abs_tol .or. z_tune == 0) then
  volt = -z_tune**2 / (branch%param%t1_with_RF(5,6) * coef_tot)
  do i = 1, n_rf
    ele => branch%ele(ix_rf(i))
    voltage_control(i)%r = volt
    call set_flags_for_changed_attribute (ele, voltage_control(i)%r)
    call lat_make_mat6 (branch%lat, ix_rf(i), ix_branch = branch%ix_branch)
  enddo
  call calc_z_tune (branch)
endif

! record

do i = 1, n_rf
  ele => branch%ele(ix_rf(i))
  volt0(i) = voltage_control(i)%r
enddo

coef = 1
dz_tune = dz_tune_func(coef, status)

coef0 = 1
dz_tune0 = dz_tune

! now set cavity voltage to get the correct tune

do k = 1, loop_max
  if (abs(dz_tune) < Qz_abs_tol) then
    return
  endif

  if (dz_tune * dz_tune0 < 0) exit  ! Have bracketed solution

  coef0 = coef
  dz_tune0 = dz_tune

  coef = coef * (z_tune / (dz_tune + z_tune))**2 
  dz_tune = dz_tune_func(coef, status)

  if (k == loop_max) then
    call out_io (s_error$, r_name, 'I CANNOT SET THE TUNE TO THE CORRECT VALUE.', &
                                   '      VALUE WANTED:   \f12.3\ ', '      VALUE OBTAINED: \f12.3\ ', &
                                   r_array = [z_tune, branch%z%tune])
    if (present(ok)) then
      ok = .false.
    else
      if (global_com%exit_on_error) call err_exit
    endif
    branch%z%stable = .false.
    return
  endif
enddo

! Have bracketed index so now find solution

coef = super_zbrent (dz_tune_func, min(coef0, coef), max(coef0, coef), Qz_rel_tol, Qz_abs_tol, status)
dz_tune = dz_tune_func(coef, status)
branch%z%stable = .true.

!-------------------------------------------------------------------------------------
contains

function dz_tune_func (coef, status) result (dz_tune)

type (ele_struct), pointer :: ele
real(rp), intent(in) :: coef
real(rp) dz_tune
integer status, i

!

do i = 1, n_rf
  ele => branch%ele(ix_rf(i))
  voltage_control(i)%r = volt0(i) * coef
  call set_flags_for_changed_attribute (ele, voltage_control(i)%r)
enddo
call lattice_bookkeeper(branch%lat)

call twiss_and_track(branch%lat, orbit, status, branch%ix_branch, print_err = print_err)
if (status == ok$) status = 0

do i = 1, n_rf
  ele => branch%ele(ix_rf(i))
  call lat_make_mat6 (branch%lat, ix_rf(i), orbit, ix_branch = branch%ix_branch)
enddo

call calc_z_tune (branch)
dz_tune = branch%z%tune - z_tune

end function dz_tune_func

end subroutine set_z_tune


