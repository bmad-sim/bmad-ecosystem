!+
! Subroutine NAME_TO_LIST (RING, ELE_NAMES, USE_ELE)
!
! Subroutine to make a list of elements in RING of the elements whose name
! matches the names in ELE_NAMES.
! This subroutine is typiclly used with MAKE_HYBRID_RING
!
! Modules Needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!     RING         -- Ring_struct: Input ring.
!     ELE_NAMES(*) -- Character array: list of element names. Wild card
!                     characters may be used. The last array element must
!                     be blank.
!
! Output:
!     USE_ELE(*)   -- Logical array: list elements referenced to the element
!                    list in RING.
!
! Example: The following makes a list of the quads and bends.
!
!     ele_names(1) = 'Q*'      ! quads
!     ele_names(2) = 'B*'      ! bends
!     ele_names(3) = ' '       ! end of ELE_NAMES list
!-

!$Id$
!$Log$
!Revision 1.3  2001/10/02 18:49:12  rwh24
!More compatibility updates; also added many explicit variable declarations.
!
!Revision 1.2  2001/09/27 18:31:55  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"



subroutine name_to_list (ring, ele_names, use_ele)

  use bmad_struct
  use cesr_utils
  implicit none

  type (ring_struct)  ring

  integer n, m, n_names
  integer ic

  logical match_wild, use_ele(*)

  character*(*) ele_names(*)

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

  return
  end
