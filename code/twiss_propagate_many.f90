!+
! Subroutine twiss_propagate_many (lat, ix_start, ix_end, direction, ix_branch)
!
! Subroutine to propagate the Twiss parameters from one point in the 
! lat to another.
!
! Note: lat%ele()%a Twiss parameters are associated with the "A" mode and
! the lat%ele()%b Twiss parameters are associated with the "B" mode.
!
! Note: Starting and ending points are just after the elements with index
!   IX_START and IX_END. For example, if DIRECTION = +1 then the first element
!   propagateed through is element ix_start+1. If DIRECTION = -1 then the first
!   element propagateed through is element ix_start.
!
! Note: If needed the subroutine will propagate through from the end of the 
!   lat to the beginning (or vice versa) to get to the end point.
!   Also: If IX_START = IX_END then the subroutine will propagate 1 full turn.
!
! Modules Needed:
!   use bmad
!
! Input:
!   lat          -- lat_struct: Lat to propagate through.
!   ix_start     -- Integer: Start index (See Note).
!   ix_end       -- Integer: End index (See Note).
!   direction    -- Integer: Direction to propagate.
!                            = +1  -> Propagate forward
!                            = -1  -> Propagate backward
!   ix_branch    -- Integer, optional: Branch index. Default is 0 (main lattice).
!   bmad_status  -- Common block status structure:
!       %type_out      -- If True then will type a message if the modes are flipped.
!       %exit_on_error -- If True then stop if there is an error.
!
! Output:
!   lat          -- lat_struct:
!     %ele(i)%a       -- X Twiss for i^th element.
!     %ele(i)%b       -- Y Twiss for i^th element.
!     %ele(i)%c_mat   -- C coupling matrix for i^th element.
!     %ele(i)%gamma_c -- Coupling gamma factor for i^th element.
!   bmad_status -- Common block status structure:
!       %ok         -- Set False if an input beta is zero. True otherwise
!-

subroutine twiss_propagate_many (lat, ix_start, ix_end, direction, ix_branch)

use bmad_struct
use bmad_interface, except_dummy => twiss_propagate_many

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch

integer, optional :: ix_branch
integer, intent(in) :: ix_start, ix_end, direction
integer i

!

bmad_status%ok = .true.
branch => lat%branch(integer_option(0, ix_branch))

select case (direction)

case (+1)

  if (ix_start < ix_end) then
    do i = ix_start, ix_end-1
      call twiss_propagate1 (branch%ele(i), branch%ele(i+1))
      if (.not. bmad_status%ok) return
    enddo
    return
  endif

  do i = ix_start, branch%n_ele_track - 1
    call twiss_propagate1 (branch%ele(i), branch%ele(i+1))
    if (.not. bmad_status%ok) return
  enddo

  if (ix_start /= 0) then
    branch%ele(0)%a       = branch%ele(branch%n_ele_track)%a
    branch%ele(0)%b       = branch%ele(branch%n_ele_track)%b
    branch%ele(0)%z       = branch%ele(branch%n_ele_track)%z
    branch%ele(0)%x       = branch%ele(branch%n_ele_track)%x
    branch%ele(0)%y       = branch%ele(branch%n_ele_track)%y
    branch%ele(0)%c_mat   = branch%ele(branch%n_ele_track)%c_mat
    branch%ele(0)%gamma_c = branch%ele(branch%n_ele_track)%gamma_c
  endif

  do i = 0, ix_end-1
    call twiss_propagate1 (branch%ele(i), branch%ele(i+1))
    if (.not. bmad_status%ok) return
  enddo

!

case (-1)

  print *, 'ERROR IN TWISS_PROPAGATE_MANY: BACKWARDS PROPAGATION NOT YET IMPLEMENTED!'
  call err_exit

!

case default

  print *, 'ERROR IN TWISS_PROPAGATE_MANY: BAD DIRECTION:', direction
  call err_exit

end select

end subroutine
