!+
! Subroutine tao_ptc_cmd (what, input_str)
!
! Interface to ptc layout.
!
! Input:
!   input_str  -- Character(*): What to show.
!-

subroutine tao_ptc_cmd (what, input_str)

use tao_interface
use ptc_layout_mod
use madx_keywords, only: against_the_method

implicit none

type (lat_struct), pointer :: lat
type (ele_struct), pointer :: ele

integer ix_cmd, iu, ib, ie, order, steps
logical changed

character(*) what, input_str
character(16) :: command(2) = [character(16):: 'init', 'reslice']
character(16) cmd_name
character(*), parameter :: r_name = 'tao_ptc_cmd'

!

call match_word (what, command, ix_cmd, .true., matched_name = cmd_name)
if (ix_cmd == 0) then
  call out_io (s_error$, r_name, 'BAD COMMAND: ' // what)
  return
endif

select case (cmd_name)

! init

case ('init')
  call lat_to_ptc_layout (s%u(1)%model%lat)

! reslice

case ('reslice')
  do iu = 1, size(s%u)
    lat => s%u(iu)%model%lat
    do ib = 0, ubound(lat%branch, 1)
      do ie = 1, lat%branch(ib)%n_ele_max
        ele => lat%branch(ib)%ele(ie)
        if (.not. has_attribute(ele, 'INTEGRATOR_ORDER')) cycle

        select case (ele%key)
        case (wiggler$, undulator$, rfcavity$, lcavity$, crab_cavity$)
        case default
          order = nint(ele%value(integrator_order$))
          steps = nint(ele%value(num_steps$))
          call against_the_method(order, steps, order, steps, 0, changed)
          call out_io (s_blank$, r_name, 'Change for ' // trim(ele%name) // &
                        ' Steps: ' // int_str(nint(ele%value(num_steps$))) // ' -> ' // int_str(steps) // &
                        ', Integrator_order: ' // int_str(nint(ele%value(integrator_order$))) // ' -> ' // int_str(order))
          ele%value(integrator_order$) = order
          ele%value(num_steps$) = steps
          ele%value(ds_step$) = ele%value(l$) / steps
        end select
      enddo
    enddo
  enddo

end select

end subroutine
