!+
! Subroutine UPDATE_HYBRID_LIST (RING, N_IN, USE_ELE)
!
! USE_ELE is a list of elements that should not be hyberdized in
! MAKE_HYBRID_RING. If an element is to be used then its controls/controllers
! should also be added to the list. This routine adds the control/controller
! elements for USE_ELE(N_IN).
!
! Input:
!     RING -- Ring_struct: Input ring.
!     N_IN -- Integer: USE_ELE(N_IN) is the element whose controlls/controllers
!             are to be added to USE_ELE.
!
! Output:
!     USE_ELE(:) -- Logical array: list RING elements to be not hyberdized
!                   with MAKE_HYBRID_RING.
!
! Note: If USE_ELE(N_IN) = .false. then no updating is done
!-

!$Id$
!$Log$
!Revision 1.4  2002/07/16 20:44:03  dcs
!*** empty log message ***
!
!Revision 1.3  2002/02/23 20:32:30  dcs
!Double/Single Real toggle added
!
!Revision 1.2  2001/09/27 18:32:01  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"



subroutine update_hybrid_list (ring, n_in, use_ele)

  use bmad
  implicit none      

  type (ring_struct)  ring

  logical use_ele(:)

  integer use_ix(100), ixu
  integer n, ix, n_in, i, j

! see if any work needs to be done

  if (.not. use_ele(n_in)) return
  if (ring%ele_(n_in)%n_slave+ring%ele_(n_in)%n_lord == 0) return

! now go through and put controlled elements in the list and make sure
! all appropriate controllers are on the list

  ixu = 1
  use_ix(1) = n_in

  do while (ixu /= 0)

    n = use_ix(ixu)
    ixu = ixu - 1

    do i = ring%ele_(n)%ix1_slave, ring%ele_(n)%ix2_slave
      ix = ring%control_(i)%ix_slave
      if (.not. use_ele(ix)) then
        use_ele(ix) = .true.
        ixu = ixu + 1
        use_ix(ixu) = ix
      endif
    enddo

    do i = ring%ele_(n)%ic1_lord, ring%ele_(n)%ic2_lord
      j= ring%ic_(i)
      ix = ring%control_(j)%ix_lord
      if (.not. use_ele(ix)) then
        use_ele(ix) = .true.
        ixu = ixu + 1
        use_ix(ixu) = ix
      endif
    enddo
  enddo

  return
  end
