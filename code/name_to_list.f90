!+
! Subroutine name_to_list (ring, ele_names, use_ele)
!
! Subroutine to make a list of elements in RING of the elements whose name
! matches the names in ELE_NAMES.
! This subroutine is typiclly used with MAKE_HYBRID_RING
!
! Modules Needed:
!   use bmad
!
! Input:
!   ring         -- Ring_struct: Input ring.
!   ele_names(:) -- Character(*): list of element names. Wild card
!                     characters may be used. The last array element must
!                     be blank.
!
! Output:
!   use_ele(:)   -- Logical array: list elements referenced to the element
!                    list in RING.
!
! Example: The following makes a list of the quads and bends.
!
!   ele_names(1) = 'Q*'      ! quads
!   ele_names(2) = 'B*'      ! bends
!   ele_names(3) = ' '       ! end of ELE_NAMES list
!-

#include "CESR_platform.inc"

subroutine name_to_list (ring, ele_names, use_ele)

  use bmad_struct
  use bmad_interface, except => name_to_list

  implicit none

  type (ring_struct)  ring

  integer n, m, n_names
  integer ic

  logical use_ele(:)

  character(*) ele_names(:)

  logical searching

! find end of lists

  ic = len(ele_names(1))

  n_names = 0
  searching = .true.
  do while (searching)
    n_names = n_names + 1
    if (len_trim(ele_names(n_names)) == 0) searching = .false.
  enddo

! initialize

  do n = 1, ring%n_ele_max
    use_ele(n) = .false.      ! no match yet
  enddo

! match

  do n = 1, ring%n_ele_max

    do m = 1, n_names
      if (match_wild(ring%ele_(n)%name, ele_names(m))) then
        use_ele(n) = .true.
        call update_hybrid_list (ring, n, use_ele)
        goto 1000
      endif
    enddo

1000    continue

  enddo

end subroutine
