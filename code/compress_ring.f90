!+
! Subroutine COMPRESS_RING (RING, OK)
!
! Subroutine to compress the ele_(*) and control_(*) arrays to remove
! elements no longer used. Note: to mark an element for removal use:
!     ring%ele_(i)%key = -1
!
! Modules Needed:
!   use bmad
!
! Input:
!     RING -- Ring_struct: Ring to compress.
!
! Output:
!     RING -- Ring_struct: Compressed ring.
!     OK   -- Logical: Ring compressed OK.
!-

!$Id$
!$Log$
!Revision 1.4  2003/01/27 14:40:32  dcs
!bmad_version = 56
!
!Revision 1.3  2002/02/23 20:32:13  dcs
!Double/Single Real toggle added
!
!Revision 1.2  2001/09/27 18:31:49  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine compress_ring (ring, ok)

  use bmad_struct
  use bmad_interface

  implicit none
                           
  type (ring_struct), target :: ring
  type (ele_struct), pointer :: ele

  type (control_struct) control_(n_control_maxx)

  integer i, j, ix, i2, ix_(n_control_maxx), ic_(n_control_maxx)
  integer n_ic, n_con, n_lord, n

  logical ok

! remove unwanted ring%ele_() elements

  ok = .true.

  i2 = 0
  do i = 1, ring%n_ele_max
    if (ring%ele_(i)%key == -1) then
      ix_(i) = garbage$
    else
      i2 = i2 + 1
      ix_(i) = i2
      ring%ele_(i2) = ring%ele_(i)
    endif
    if (i == ring%n_ele_ring) ring%n_ele_ring = i2
    if (i == ring%n_ele_use)  ring%n_ele_use = i2
  enddo

  ring%n_ele_max = i2
  if (ring%n_ele_symm > 0) ring%n_ele_symm = ix_(ring%n_ele_symm)

! renumber ring%control_()%ix_ele  

  forall (i = 1:ring%n_control_array) 
    ring%control_(i)%ix_lord = ix_(ring%control_(i)%ix_lord)
    ring%control_(i)%ix_slave = ix_(ring%control_(i)%ix_slave)
  end forall

! compress ring%control_() array

  n_con = 0
  ix_ = 0

  do i = 1, ring%n_ele_max
    ele => ring%ele_(i)
    do j = 1, ele%n_slave
      control_(n_con+j) = ring%control_(ele%ix1_slave+j-1)
      ix_(ele%ix1_slave+j-1) = n_con+j
    enddo
    if (ele%n_slave > 0) then
      ele%ix1_slave = n_con + 1
      ele%ix2_slave = n_con + ele%n_slave
      n_con = n_con + ele%n_slave
    endif
  enddo

! compress ring%ic_() array

  n_ic = 0

  do i = 1, ring%n_ele_max
    ele => ring%ele_(i)
    n_lord = 0
    do i2 = ele%ic1_lord, ele%ic2_lord
      ix = ring%ic_(i2)
      if (ix_(ix) == 0) then
        if (ele%control_type == super_slave$) then
          type *, 'ERROR IN COMPRESS_RING: SUPERPOSITION LORD HAS BEEN REMOVED!'
          ok = .false.
        endif
      else
        n_lord = n_lord + 1                                           
        n_ic = n_ic + 1
        ic_(n_ic) = ix_(ix)
      endif
    enddo

    if (n_lord == 0) then
      ele%ic1_lord = 0
      ele%ic2_lord = -1
      ele%n_lord = 0
      if (i <= ring%n_ele_ring) ele%control_type = free$
    else
      ele%ic1_lord = n_ic - n_lord + 1
      ele%ic2_lord = n_ic
      ele%n_lord = n_lord
    endif

  enddo

  ring%control_ = control_                        
  ring%n_control_array = n_con
  ring%ic_ = ic_
  ring%n_ic_array = n_ic

! do a check

  call check_ring_controls (ring, .true.)

end subroutine
          
