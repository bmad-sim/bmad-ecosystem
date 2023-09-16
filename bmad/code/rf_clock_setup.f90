!+
! Function rf_clock_setup (branch, n_rf_included, n_rf_excluded) result (ok)
!
! Routine to setup long term tracking to use the "RF clock" when absolute time tracking is used.
! This is to avoid round-off errors when tracking for many turns.
!
! For a cavity to use the RF clock, the cavity's RF frequency must be commensurate with the RF clock period.
! That is: RF clock period * cavity frequency = N integer.
! It may be that the lattice branch has cavities with incommensurate RF frequencies.
! In this case, the RF clock period will be chosen to include the largest number of cavities.
!
! The RF clock period is stored in:
!   bmad_private%rf_clock_period
! The RF clock time for a particle is in stored in coord_struct: 
!   coord%phase(1)  -- Number of RF clock cycles.
!   coord%t         -- Factional time.
!
! Input:
!   branch        -- branch_struct: Lattice branch to setup.
!
! Output:
!   n_rf_included    -- integer: Number of RF cavities with commensurate frequencies using the RF clock.
!   n_rf_excluded    -- integer: Number of RF cavities with incommensurate frequencies not using the RF clock.
!   ok               -- logical: Set True if setup was successful.
!-

function rf_clock_setup (branch, n_rf_included, n_rf_excluded) result (ok)

use bmad_routine_interface, dummy => rf_clock_setup
use attribute_mod, only: has_attribute

implicit none

type rf_ele_info_struct
  type (ele_struct), pointer :: ele
  real(rp) rf_freq
  integer ix_freq    ! Index for ac_kicker where there can be multiple frequencies
  integer ix_basket
end type

type basket_struct
  real(rp) clock_period
  integer :: max_harmonic = 1
  integer :: n_in_basket = 1
end type

type (branch_struct), target :: branch
type (ele_struct), pointer :: ele, ele2
type (rf_ele_info_struct), allocatable :: rf_ele(:)
type (basket_struct), allocatable :: basket(:)

real(rp) ff, clock_period
integer n_rf_included, n_rf_excluded
integer n, ie, ie2, ib, ib0, irf, nx, nn, n_rf, n_basket
logical ok, has_ac_freq

character(*), parameter :: r_name = 'rf_clock_setup'

! Make a list of RF elements

ok = .false.
bmad_private%rf_clock_period = 0
allocate (rf_ele(100))
n_rf = 0

do ie = 1, branch%n_ele_track
  ele => branch%ele(ie)

  if (ele%slave_status == super_slave$) then
    do ie2 = 1, ele%n_lord
      ele2 => pointer_to_lord(ele, ie2)
      call this_ele_setup(ele2, rf_ele, n_rf)
    enddo
  else
    call this_ele_setup(ele, rf_ele, n_rf)
  endif
enddo

! Look for a commensurate frequency to be used as the RF clock frequency.
! A commensurate frequency is such that rf_clock_period * cavity frequency = N
! where N is an integer less than or equal to 100.
! Put groups of RF elements with commensurate frequencies into "baskets".

if (n_rf == 0) return
rf_ele(1)%ix_basket = 1
allocate (basket(n_rf))
n_basket = 1
basket(1)%clock_period = 1.0_rp / rf_ele(1)%rf_freq
basket(1)%max_harmonic = 1

rf_ele_loop: do irf = 2, n_rf
  do ib = 1, n_basket
    do nx = 1, 100/basket(ib)%max_harmonic
      clock_period = nx * basket(ib)%clock_period
      ff = clock_period * rf_ele(irf)%rf_freq
      if (abs(nint(ff) - ff) > 1e-14_rp) cycle
      basket(ib)%n_in_basket = basket(ib)%n_in_basket + 1
      basket(ib)%clock_period = clock_period
      basket(ib)%max_harmonic = max(nx*basket(ib)%max_harmonic, nint(ff))
      rf_ele(irf)%ix_basket = ib
      cycle rf_ele_loop
    enddo

    n_basket = n_basket + 1
    basket(n_basket)%clock_period = 1.0_rp / rf_ele(irf)%rf_freq
    basket(n_basket)%max_harmonic = 1
    rf_ele(irf)%ix_basket = n_basket
  enddo
enddo  rf_ele_loop

! Find basket with most elements. 

ib0 = maxloc(basket(1:n_basket)%n_in_basket, 1)
bmad_private%rf_clock_period = basket(ib0)%clock_period

n_rf_included = basket(ib0)%n_in_basket
n_rf_excluded = n_rf - n_rf_included

if (n_rf_excluded > 0) then
  call out_io (s_warn$, r_name, 'Number of cavities not using RF clock due to incomensurate frequency: ' // &
                                                                                          int_str(n_rf_excluded))
endif

do irf = 1, n_rf
  ele => rf_ele(irf)%ele
  if (rf_ele(irf)%ix_basket == ib0) then
    has_ac_freq = .false.
    if (associated(ele%ac_kick)) has_ac_freq = allocated(ele%ac_kick%frequency)
    if (ele%key == ac_kicker$ .and. has_ac_freq) then
      n = rf_ele(irf)%ix_freq
      ele%ac_kick%frequency(n)%rf_clock_harmonic = nint(rf_ele(irf)%rf_freq * bmad_private%rf_clock_period)
    else
      ele%value(rf_clock_harmonic$) = nint(rf_ele(irf)%rf_freq * bmad_private%rf_clock_period)
    endif
  else
    call out_io (s_warn$, r_name, '  Cavity not using RF clock: ' // ele%name)
  endif
enddo

ok = .true.

!---------------------------------------------------------------------------
contains

subroutine this_ele_setup (ele, rf_ele, n_rf)

type (ele_struct), target :: ele
type (rf_ele_info_struct), allocatable :: rf_ele(:), rf_temp(:)
real(rp) rf_freq
integer i, n_rf, n_new
logical has_ac_freq

!

has_ac_freq = .false.
if (associated(ele%ac_kick)) has_ac_freq = allocated(ele%ac_kick%frequency)

if (ele%key == ac_kicker$ .and. has_ac_freq) then
  n_new = size(ele%ac_kick%frequency)
elseif (has_attribute(ele, 'REPETITION_FREQUENCY')) then
  if (ele%value(repetition_frequency$) == 0) return
  rf_freq = ele%value(repetition_frequency$)
  n_new = 1
elseif (has_attribute(ele, 'RF_FREQUENCY')) then
  rf_freq = ele%value(rf_frequency$)
  n_new = 1
else
  return
endif

if (n_rf+n_new > size(rf_ele)) then
  call move_alloc(rf_ele, rf_temp)
  allocate (rf_ele(2*n_rf + n_new))
  rf_ele(1:n_rf) = rf_temp
endif

if (ele%key == ac_kicker$ .and. has_ac_freq) then
  do i = 1, size(ele%ac_kick%frequency)
    n_rf = n_rf + 1
    rf_ele(n_rf)%ele => ele
    rf_ele(n_rf)%rf_freq = ele%ac_kick%frequency(i)%f
    rf_ele(n_rf)%ix_freq = i
    ele%ac_kick%frequency(i)%rf_clock_harmonic = 0
  enddo
else
  n_rf = n_rf + 1
  rf_ele(n_rf)%ele => ele
  rf_ele(n_rf)%rf_freq = rf_freq
  ele%value(rf_clock_harmonic$) = 0
endif

end subroutine this_ele_setup
end function rf_clock_setup
