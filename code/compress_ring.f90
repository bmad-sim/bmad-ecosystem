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

#include "CESR_platform.inc"

subroutine compress_ring (ring, ok)

  use bmad_struct
  use bmad_interface, except => compress_ring

  implicit none
                           
  type (ring_struct), target :: ring
  type (ele_struct), pointer :: ele

  type (control_struct), allocatable :: control_(:)

  integer i, j, ix, i2
  integer n_ic, n_con, n_lord, n
  integer, allocatable :: ix_(:), ic_(:)

  logical ok

! allocate

  n = ring%n_control_max
  allocate (control_(n))
  allocate (ix_(n))
  allocate (ic_(ring%n_ic_max))

! remove unwanted ring%ele_() elements

  ok = .true.

  i2 = 0
  do i = 1, ring%n_ele_max
    if (ring%ele_(i)%key == -1) then
      ix_(i) = int_garbage$
    else
      i2 = i2 + 1
      ix_(i) = i2
      ring%ele_(i2) = ring%ele_(i)
    endif
    if (i == ring%n_ele_use) then
       ring%n_ele_use = i2
       ring%n_ele_ring = i2
    endif
  enddo

  ring%n_ele_max = i2

! renumber ring%control_()%ix_ele  

  forall (i = 1:ring%n_control_max) 
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
          print *, 'ERROR IN COMPRESS_RING: SUPERPOSITION LORD HAS BEEN REMOVED!'
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
      if (i <= ring%n_ele_use) ele%control_type = free$
    else
      ele%ic1_lord = n_ic - n_lord + 1
      ele%ic2_lord = n_ic
      ele%n_lord = n_lord
    endif

  enddo

  ring%control_(1:n_con) = control_(1:n_con)                     
  ring%n_control_max = n_con
  ring%ic_(1:n_ic) = ic_(1:n_ic)
  ring%n_ic_max = n_ic

! deallocate and do a check

  deallocate (control_)
  deallocate (ix_)
  deallocate (ic_)

  call check_ring_controls (ring, .true.)

end subroutine
          
