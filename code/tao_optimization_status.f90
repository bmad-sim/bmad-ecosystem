!+
! Function tao_optimization_status (datum) result (why_str)
!
! Routine to return a string saying why a datum's useit_opt component is set to False or
! if the delta_merit calc is using the datum's invalid value.
!
! Input:
!   datum     -- tao_data_struct: Datum to evaluate.
!
! Output:
!   why_str   -- character(60): Optimization status of the datum.
!-

function tao_optimization_status (datum) result (why_str)

use tao_interface, dummy => tao_optimization_status

implicit none

type (tao_data_struct) :: datum
integer ix_uni
character(60) why_str

!

ix_uni = datum%d1%d2%ix_universe

if (datum%useit_opt) then
  if (.not. datum%good_model) then
    why_str = "! good_model = False so delta_merit calc uses invalid value"
  elseif (s%global%opt_with_ref .and. .not. datum%good_design) then
    why_str = "! good_design = False so delta_merit calc uses invalid value"
  elseif (s%global%opt_with_base .and. .not. datum%good_base) then
    why_str = "! good_base = False so delta_merit calc uses invalid value"
  else
    why_str = ''
  endif

elseif (.not. datum%good_opt) then
  why_str = '! useit_opt = False since good_opt = False'

elseif (.not. datum%exists) then
  why_str = '! useit_opt = False since exists = False'

elseif (.not. datum%good_user) then
  why_str = '! useit_opt = False since good_user = False'

elseif (.not. datum%good_meas) then
  why_str = '! useit_opt = False since good_meas = False'

elseif (.not. datum%good_ref .and. s%global%opt_with_ref) then
  why_str = '! useit_opt = False since good_ref = False'

elseif (.not. s%u(ix_uni)%is_on) then
  why_str = "! useit_opt = False since Datum's universe is Off"

else
  why_str = '??? Please report this message!'
endif

end function
