!+
! Subroutine ring_to_layout (ring, ptc_layout)
!
! Subroutine to create a PTC layout from a BMAD ring.
! Note: If ptc_layout has been already used then you should first do a 
!           call kill(ptc_layout)
! This deallocates the pointers in the layout
!
! Note: Before you call this routine you need to first call:
!    call set_ptc (ring%param, ...)
!
! Modules needed:
!   use accelerator
!
! Input:
!   ring -- Ring_struct: 
!
! Output:
!   ptc_layout -- Layout:
!-

subroutine ring_to_layout (ring, ptc_layout)

  use accelerator

  implicit none

  type (ring_struct), intent(in) :: ring
  type (layout), intent(inout) :: ptc_layout
  type (fibre), pointer :: fib

  real*8 energy, kinetic, beta0, p0c, brho

  integer i

! setup

  call set_up (ptc_layout)

! transfer energy, etc.

!  ptc_layout%energy = ring%param%energy
!  call energy_to_kinetic (ptc_layout%energy, ring%param%particle, &
!                             ptc_layout%kinetic, ptc_layout%beta0, &
!                             ptc_layout%p0c, ptc_layout%brho)
!  ptc_layout%circumference = 0


! transfer elements.

  do i = 1, ring%n_ele_ring
    allocate (fib)
    call ele_to_fibre (ring%ele_(i), fib, ring%param)
    call append (ptc_layout, fib)
    call kill (fib)
  enddo

! circular or not?

  if (ring%param%lattice_type == circular_lattice$) then
    ptc_layout%closed = .true.
    call ring_l (ptc_layout, .true.)
  else
    ptc_layout%closed = .false.
  endif

end subroutine
