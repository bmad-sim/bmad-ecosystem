!+
! Subroutine RING_MAKE_MAT6 (RING, IX_ELE, COORD_)
!
! Subroutine to make the 6x6 linear transfer matrix for an element
!
! Moudules Needed:
!   use bmad
!
! Input:
!   RING        -- Ring_struct: Ring containing the element.
!   IX_ELE      -- Integer: Index of the element. if < 0 then entire
!                    ring will be made. In this case group elements will
!                    be made up last.
!   COORD_(0:n_ele_maxx) -- [Optional] Coord_struct: Coordinates of the 
!                   center around which the matrix is calculated. If not
!                   present then the orbit is taken to be the orign.
!
! Output:
!   RING%ELE_(I)%MAT6 -- 6x6 transfer matrices.
!-

!$Id$
!$Log$
!Revision 1.4  2002/02/23 20:32:23  dcs
!Double/Single Real toggle added
!
!Revision 1.3  2002/01/16 21:04:18  helms
!Fixed problem with passing optional arguments.
!
!Revision 1.2  2001/09/27 18:31:57  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


recursive subroutine ring_make_mat6 (ring, ix_ele, coord_)

  use bmad_struct
  use bmad_interface, only: control_bookkeeper, make_mat6

  implicit none
                                         
  type (ring_struct)  ring
  type (coord_struct), optional :: coord_(0:n_ele_maxx)
  type (ele_struct)  ele

  integer i, ix_ele, i1, i2, i3, ix1, ix2, ix3

! make entire ring if ix_ele < 0

  if (ix_ele < 0) then         

    do i = ring%n_ele_ring+1, ring%n_ele_max
      if (ring%ele_(i)%control_type /= group_lord$)  &
                                 call control_bookkeeper (ring, i)
    enddo

    do i = ring%n_ele_ring+1, ring%n_ele_max
      if (ring%ele_(i)%control_type == group_lord$)  &
                                 call control_bookkeeper (ring, i)
    enddo

    do i = 1, ring%n_ele_use
      if (ring%ele_(i)%key /= hybrid$) then
         if (present(coord_)) then
            call make_mat6(ring%ele_(i), ring%param, coord_(i-1), coord_(i))
         else
            call make_mat6(ring%ele_(i), ring%param)
         endif
      endif
    enddo

    do i = ring%n_ele_ring+1, ring%n_ele_max
      if (ring%ele_(i)%control_type == component_lord$) then
         if (present(coord_)) then
            call make_mat6(ring%ele_(i), ring%param, coord_(i-1), coord_(i))
         else
            call make_mat6(ring%ele_(i), ring%param)
         endif
      endif
   enddo

    return

  endif

!-----------------------------------------------------------
! otherwise make a single element

  ele = ring%ele_(ix_ele)
  call control_bookkeeper (ring, ix_ele)

! for a regular element

  if (ele%key == hybrid$) return

  if (ix_ele <= ring%n_ele_ring .or.  &
                              ele%control_type == component_lord$) then
     if (present(coord_)) then
        call make_mat6(ring%ele_(ix_ele), ring%param, &
                                         coord_(ix_ele-1), coord_(ix_ele))
     else
        call make_mat6(ring%ele_(ix_ele), ring%param)
     endif

    return
  endif                        

! for a control element

  do i1 = ring%ele_(ix_ele)%ix1_slave, ring%ele_(ix_ele)%ix2_slave
    i = ring%control_(i1)%ix_slave
     if (present(coord_)) then
        call ring_make_mat6 (ring, i, coord_)
     else
        call ring_make_mat6 (ring, i)
     endif
  enddo

end subroutine
