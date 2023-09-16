!+
! Subroutine set_ptc_base_state (component, set_val, old_val)
!
! Routine to set the ptc_private%base_state used for PTC tracking.
!
! Input:
!   component   -- character(*): Name of component. "TOTALPATH", "SPIN", "NOCAVITY", "TIME", etc.
!                   See the PTC internal_state structure for component names. 
!   set_val     -- logical: Value to set to. For TOTALPATH, True => 1, False => 0.
!
! Output:
!   old_val     -- logical, optional: Old value.
!-

subroutine set_ptc_base_state (component, set_val, old_val)

use bmad_struct, only: ptc_private

implicit none

character(*) component
logical set_val
logical, optional :: old_val

!

select case (component)
case ('TOTALPATH')
  if (present(old_val)) old_val = (ptc_private%base_state%totalpath == 1)

  if (set_val) then; ptc_private%base_state%totalpath = 1
  else;              ptc_private%base_state%totalpath = 0
  endif

case("TIME");       call set_it(ptc_private%base_state%time)
case("RADIATION");  call set_it(ptc_private%base_state%radiation)
case("NOCAVITY");   call set_it(ptc_private%base_state%nocavity)
case("FRINGE");     call set_it(ptc_private%base_state%fringe)
case("STOCHASTIC"); call set_it(ptc_private%base_state%stochastic)
case("ENVELOPE");   call set_it(ptc_private%base_state%envelope)
case("PARA_IN");    call set_it(ptc_private%base_state%para_in)
case("ONLY_4D");    call set_it(ptc_private%base_state%only_4d)
case("DELTA");      call set_it(ptc_private%base_state%delta)
case("SPIN");       call set_it(ptc_private%base_state%spin)
case("MODULATION"); call set_it(ptc_private%base_state%modulation)
case("ONLY_2D");    call set_it(ptc_private%base_state%only_2d)
case("FULL_WAY");   call set_it(ptc_private%base_state%full_way)
case default;       call err_exit
end select

!-------------------------------------------------
contains

subroutine set_it (param)
logical param
if (present(old_val)) old_val = param
param = set_val
end subroutine set_it

end subroutine
