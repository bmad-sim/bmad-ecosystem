!+
! Subroutine twiss_propagate_many (ring, ix_start, ix_end, direction)
!
! Subroutine to propagate the Twiss parameters from one point in the 
! ring to another.
!
! Note: ring%ele_()%x Twiss parameters are associated with the "A" mode and
! the ring%ele_()%y Twiss parameters are associated with the "B" mode.
!
! Note: Starting and ending points are just after the elements with index
!   IX_START and IX_END. For example, if DIRECTION = +1 then the first element
!   propagateed through is element ix_start+1. If DIRECTION = -1 then the first
!   element propagateed through is element ix_start.
!
! Note: If needed the subroutine will propagate through from the end of the 
!   ring to the beginning (or vice versa) to get to the end point.
!   Also: If IX_START = IX_END then the subroutine will propagate 1 full turn.
!
! Modules Needed:
!   use bmad
!
! Input:
!   ring             -- Ring_struct: Ring to propagate through.
!   ix_start         -- Integer: Start index (See Note).
!   ix_end           -- Integer: End index (See Note).
!   direction        -- Integer: Direction to propagate.
!                            = +1  -> Propagate forward
!                            = -1  -> Propagate backward
!   bmad_status -- Common block status structure:
!       %type_out -- If True then will type a message if the modes are flipped.
!       %exit_on_error -- If True then stop if there is an error.
!
! Output:
!   ring          -- Ring_struct:
!     %ele_(i)%x       -- X Twiss for i^th element.
!     %ele_(i)%y       -- Y Twiss for i^th element.
!     %ele_(i)%c_mat   -- C coupling matrix for i^th element.
!     %ele_(i)%gamma_c -- Coupling gamma factor for i^th element.
!   bmad_status -- Common block status structure:
!       %ok         -- Set False if an input beta is zero. True otherwise
!-

subroutine twiss_propagate_many (ring, ix_start, ix_end, direction)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ring_struct) :: ring

  integer, intent(in) :: ix_start, ix_end, direction
  integer i

!

  bmad_status%ok = .true.

  select case (direction)

  case (+1)

    if (ix_start < ix_end) then
      do i = ix_start, ix_end-1
        call twiss_propagate1 (ring%ele_(i), ring%ele_(i+1))
        if (.not. bmad_status%ok) return
      enddo
      return
    endif

    do i = ix_start, ring%n_ele_ring - 1
      call twiss_propagate1 (ring%ele_(i), ring%ele_(i+1))
      if (.not. bmad_status%ok) return
    enddo

    if (ix_start /= 0) then
      ring%ele_(0)%x       = ring%ele_(ring%n_ele_ring)%x
      ring%ele_(0)%y       = ring%ele_(ring%n_ele_ring)%y
      ring%ele_(0)%c_mat   = ring%ele_(ring%n_ele_ring)%c_mat
      ring%ele_(0)%gamma_c = ring%ele_(ring%n_ele_ring)%gamma_c
    endif

    do i = 0, ix_end-1
      call twiss_propagate1 (ring%ele_(i), ring%ele_(i+1))
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
