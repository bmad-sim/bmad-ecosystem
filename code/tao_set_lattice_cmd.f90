! Subroutine tao_set_lattice_cmd (universe, component, set_value)
!
! Sets a lattice equal to another. This will also recalculate the data in
! the lattice.
!
! Input:
!   universe    -- Character(*): Which universe in the form *;n where
!                  n is the universe number, any word can be before the ';'
!   component   -- Character(*): Which component to set.
!   set_value   -- Character(*): What value to set to.
!
!  Output:
!-

subroutine tao_set_lattice_cmd (universe, component, set_value)

use tao_mod
use quick_plot

implicit none

character(*) universe, component, set_value
character(20) :: r_name = 'tao_set_lattice_cmd'

integer i

logical, automatic :: this_u(size(s%u))
logical err

call tao_pick_universe (universe, universe, this_u, err)
if (err) return

do i = 1, size(s%u)
  if (.not. this_u(i)) cycle
  call set_lat (s%u(i))
  if (err) return
enddo

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
contains

subroutine set_lat (u)

type (tao_universe_struct), target :: u
type (ring_struct), pointer :: u_set_this
type (ring_struct), pointer :: u_to_this

err = .false.

select case (component)
  case ('model')
    u_set_this => u%model
  case ('base')
    u_set_this => u%base
  case default
    call out_io (s_error$, r_name, 'BAD LATTICE: ' // component)
    err = .true.
    return
end select

select case (set_value)
  case ('model')
    u_to_this => u%model
  case ('base')
    u_to_this => u%base
  case ('design')
    u_to_this => u%design
  case default
    call out_io (s_error$, r_name, 'BAD LATTICE: ' // component)
    err = .true.
    return
end select
  
u_set_this = u_to_this

end subroutine set_lat

end subroutine tao_set_lattice_cmd
