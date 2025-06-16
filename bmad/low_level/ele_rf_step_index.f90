!+
! Function ele_rf_step_index(E_ref, s_rel, ele, include_downstream_end) result (ix_step)
!
! Routine to return the step index in the ele%rf%steps(:) array at a particular s-position or referece energy.
!
! E_ref is used in the computation instead of s_rel since s_rel is ambiguous if
! s_rel corresponds to a slice boundary.
!
! Input:
!   E_ref                   -- real(rp): Reference energy of step. If negative, ignore and use s_rel.
!   s_rel                   -- real(rp): S-position relative to the beginning of the element
!   ele                     -- real(rp): RF cavity.
!   include_downstream_end  -- logical, optional: Used to remove ambiguity when creating super and slice slaves and E_ref
!                               is not available. Default is True.
!
! Output:
!   ix_step       -- integer: Corresponding index in the ele%rf%steps(:) array.
!-

function ele_rf_step_index(E_ref, s_rel, ele, include_downstream_end) result (ix_step)

use bmad_routine_interface, dummy => ele_rf_step_index

implicit none

type (ele_struct) :: ele
real(rp) E_ref, dE, dE_rel, s_rel
integer ix_step, n_step
logical, optional :: include_downstream_end
character(*), parameter :: r_name = 'ele_rf_step_index'

!

if (.not. associated(ele%rf)) then
  call out_io(s_error$, r_name, 'MISSING RF STEP BOOKKEEPING PARAMETERS. PLEASE REPORT THIS! FOR ELEMENT ' // ele_full_name(ele))
  return
endif

n_step = ubound(ele%rf%steps, 1) - 1  ! Number of slices
dE = ele%value(E_tot$) - ele%value(E_tot_start$)

! dE == 0 case. Must use s_rel in this case.

if (dE == 0 .or. E_ref < 0) then
  if (s_rel == 0) then
    ix_step = 0
  elseif (s_rel == ele%value(l$) .or. (logic_option(.true.) .and. s_rel > ele%value(l$) - bmad_com%significant_length)) then
    ix_step = n_step + 1
  else
    ix_step = int(n_step * s_rel / ele%value(l$)) + 1
  endif
  return
endif

! dE /= 0 case

ix_step = nint(2.0_rp * n_step * (E_ref - ele%value(E_tot_start$)) / dE)
if (ix_step == 2*n_step) then
  ix_step = n_step + 1
elseif (ix_step > 0) then
  ix_step = (ix_step + 1) / 2
endif

! Sanity check

if (ix_step < 0 .or. ix_step > n_step+1) then
  call out_io(s_error$, r_name, 'RF STEP INDEX BOOKKEEPING FAILURE. PLEASE REPORT THIS!  FOR ELEMENT ' // ele_full_name(ele))
  return
endif

if (abs(E_ref - ele%rf%steps(ix_step)%E_tot0) > 1.0e-10 * (ele%value(E_tot$) + ele%value(E_tot_start$))) then
  call out_io(s_error$, r_name, 'RF STEP ENERGY BOOKKEEPING FAILURE. PLEASE REPORT THIS!  FOR ELEMENT ' // ele_full_name(ele))
  return
endif

end function ele_rf_step_index
