!+
! Subroutine type_fibre (fib)
!
! Subroutine to print the global information in a fibre
!
! Modules Needed:
!   use accelerator
!
! Input:
!   fib - fibre: fibre to use.
!+

subroutine type_fibre (fib)

  use accelerator

  implicit none

  type (fibre), intent(in) :: fib

  integer i

!
  
  if (.not. associated (fib%mag)) then
    type *, 'Warning from TYPE_FIBRE: Fibre NOT associated with anything.'
    return
  endif

  type *, 'Name:        ', fib%mag%name
  type *, 'Vorname:     ', fib%mag%vorname
  type *, 'Kind:        ', kind_name(fib%mag%kind)
  type *, 'Knob:        ', fib%magp%knob
    type *, 'L:        ', fib%mag%l

  if (fib%mag%kind == kind4) then
    type *, 'Voltage:  ', fib%mag%volt
    type *, 'Frequency:', fib%mag%freq
    type *, 'Voltage:  ', fib%mag%volt
    type *, 'Phase:    ', fib%mag%phas
    type *, 'Delta_e:  ', fib%mag%delta_e
    type *, 'Thin: ', fib%mag%thin
  endif

  if (fib%mag%kind == kind5) then
    type *, 'KS:       ', fib%mag%b_sol
    type *, 'Thin:     ', fib%mag%thin
  endif

  if (fib%mag%kind == kind2 .and. fib%mag%p%b0 /= 0) then
    type *, 'E1:       ', fib%mag%p%edge(1)
    type *, 'E2:       ', fib%mag%p%edge(2)
    type *, 'Rho:      ', fib%mag%p%b0
    type *, 'L_chord:  ', fib%mag%p%lc
  endif

  type *, 'Integration Order: ', fib%mag%p%method
  type *, 'Integration Steps: ', fib%mag%p%nst


  do i = lbound(fib%mag%bn, 1), ubound(fib%mag%bn, 1)
    if (fib%mag%bn(i) /= 0) type '(a, i2, a, 5x, 1pd12.3)', &
                                  ' BN(', i, '):', fib%mag%bn(i)
  enddo  
  do i = lbound(fib%mag%an, 1), ubound(fib%mag%an, 1)
    if (fib%mag%an(i) /= 0) type '(a, i2, a, 5x, 1pd12.3)', &
                                  ' AN(', i, '):', fib%mag%an(i)
  enddo  


end subroutine
