!+
! Subroutine type2_taylors (bmad_taylor, lines, n_lines)
!
! Subroutine to print the terms in taylor array.
!
! Moudles needed:
!   use bmad
!
! Input:
!   bmad_taylor(6) -- Taylor_struct: Array of taylors.
!
! Output:
!   lines(:)     -- Character*80, allocatable: Character array to hold the 
!                     output. The array size of lines(:) will be set by
!                     this subroutine.
!   n_lines      -- Number of lines in lines(:).
!-

#include "CESR_platform.inc"

subroutine type2_taylors (bmad_taylor, lines, n_lines)

  use bmad

  implicit none

  type (taylor_struct), intent(in), target :: bmad_taylor(6)
  type (taylor_term_struct), pointer :: tt
  type (taylor_struct) tlr

  integer, intent(out) :: n_lines
  integer i, j, k, nl, ix

  character*80, pointer :: lines(:)
  character*40, fmt1, fmt2, fmt

!

  deallocate (lines, stat = ix)
  n_lines = 8 + sum( (/ (size(bmad_taylor(i)%term), i = 1, 6) /) )
  allocate(lines(n_lines))
  

  write (lines(nl+1), *) 'Taylor Terms:'
  write (lines(nl+2), *) &
        'Out     Coef              Exponents           Order        Reference'
  nl = nl+2


  fmt1 = '(i4, a, f20.12, 6i3, i9, f18.9)'
  fmt2 = '(i4, a, 1p, e20.11, 0p, 6i3, i9, f18.9)'

  do i = 1, 6
    write (lines(nl+1), *) &
                      '---------------------------------------------------'
    nl = nl + 1

    nullify (tlr%term)
    call sort_taylor_terms (bmad_taylor(i), tlr)

    do j = 1, size(bmad_taylor(i)%term)

      tt => tlr%term(j)

      if (abs(tt%coef) < 1e5) then
        fmt = fmt1
      else
        fmt = fmt2
      endif

      if (j == 1) then
        write (lines(nl+j), fmt) i, ':', tt%coef, &
                    (tt%exp(k), k = 1, 6), sum(tt%exp), bmad_taylor(i)%ref
      else
        write (lines(nl+j), fmt) i, ':', tt%coef, &
                    (tt%exp(k), k = 1, 6), sum(tt%exp)
      endif
    enddo

    nl = nl + size(bmad_taylor(i)%term)
    deallocate (tlr%term)

  enddo

end subroutine
