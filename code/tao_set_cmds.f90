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
!    su(n)      -- ring_struct: changes specified lattice in specified universe 
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

!------------------------------------------------------------------------------
!+
! Subroutine tao_set_global_cmd (who, set_value)
!
! Routine to set global variables
! 
! Input:
!   who       -- Character(*): which global variable to set
!   set_value -- Character(*): Value to set to.
!
! Output:
!    %global  -- Global variables structure.
!-

subroutine tao_set_global_cmd (who, set_value)

use tao_mod
use quick_plot

implicit none

type (tao_global_struct) global

character(*) who, set_value
character(20) :: r_name = 'tao_set_global_cmd'

integer iu, ios

namelist / params / global

! open a scratch file for a namelist read

iu = lunget()
open (iu, status = 'scratch')
write (iu, *) '&params'
write (iu, *) ' global%' // trim(who) // ' = ' // trim(set_value)
write (iu, *) '/'
rewind (iu)
global = s%global  ! set defaults
read (iu, nml = params, iostat = ios)
close (iu)

if (ios == 0) then
  s%global = global
else
  call out_io (s_error$, r_name, 'BAD COMPONENT OR NUMBER')
endif

end subroutine tao_set_global_cmd

!------------------------------------------------------------------------------
! Subroutine tao_set_plot_cmd (component, set_value)
!
!  Set various aspects of the plotting window
!
! Input:
!   component     -- Character(*): Which component to set.
!   set_value     -- Character(*): What value to set to.
!
!  Output:
!    s%plot_page  -- tao_plot_page_struct:
!-

subroutine tao_set_plot_cmd (component, set_value1, set_value2)

use tao_mod
use quick_plot

implicit none

character(*) component, set_value1
character(*), optional :: set_value2
character(20) :: r_name = 'tao_set_plot_cmd'

real(rp) x, y
integer ix

select case (component)

  case ('title')
    s%plot_page%title(1)%string = trim(set_value1)

  case ('subtitle')
    s%plot_page%title(2)%string = trim(set_value1)
    s%plot_page%title(2)%draw_it = .true.

  case ('subtitle_loc')
    
    if (.not. present(set_value2)) then
      call out_io(s_info$, r_name, "subtitle_loc requires two numbers.")
      return
    endif
    
    read(set_value1, '(f15.10)') x
    read(set_value2, '(f15.10)') y
    s%plot_page%title(2)%x = x
    s%plot_page%title(2)%y = y

end select

end subroutine tao_set_plot_cmd

