!+
! Subroutine setup_scale_lord (ring, create_new, ref_name, ix_slave, ix_lord)
!
! Subroutine to setup the control information for  a scale_lord of a multipole.
! Note: If ring%ele_(ix_lord)%value(radius$) = 0 then the default value of 1.0
! will be used.
!
! Modules needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!   ring       -- Ring_struct: Ring to be used
!   create_new -- Logical: If True then create a new scale_lord.
!   ref_name   -- Character*16: Name of the reference element.
!   ix_slave   -- Integer: Index of scale_lord element.
!   ix_lord    -- Integer: Index of the lord element if create_new = False.
!
! Output:
!   ix_lord -- Integer: Index of the lord element created if create_new = True.
!-

!$Id$
!$Log$
!Revision 1.2  2001/09/27 18:31:58  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine setup_scale_lord (ring, create_new, ref_name, ix_slave, ix_lord)

  use bmad_struct
  use bmad_interface
                             
  implicit none

  type (ring_struct), target :: ring
  type (ele_struct), pointer :: slave, lord, ref

  integer i, n, nrs, ix_slave, ix_lord, ix_ref

  character*16 ref_name

  logical create_new

! find reference element

  ix_ref = 0
  do i = 1, ring%n_ele_max
    if (ring%ele_(i)%name == ref_name) then
      ix_ref = i
      exit
    endif
  enddo

  if (ix_ref == 0) then
    print *, 'ERROR IN SETUP_SCALE_LORD: CANNOT FIND REFERENCE: ', ref_name
    bmad_status%ok = .false.
    if (bmad_status%exit_on_error) call err_exit
  endif

! set pointers to point to appropriate elements

  ref => ring%ele_(ix_ref)
  slave => ring%ele_(ix_slave)
  slave%control_type = scale_slave$

  if (create_new) call new_control (ring, ix_lord)
  lord => ring%ele_(ix_lord)
  lord%key = slave%key
  lord%control_type = scale_lord$
  lord%value = slave%value
  if (lord%value(radius$) == 0) lord%value(radius$) = 1
  slave%value(radius$) = 0  ! does not get used.
  
! name change

  if (create_new) then
    lord%name = slave%name
  endif

  slave%name = trim(slave%name) // '\' // ref%name

! error check

  if (slave%key /= multipole$ .and. slave%key /= ab_multipole$) then
    type *, 'ERROR IN SETUP_SCALE_LORD: SLAVE IS NOT A MULTIPOLE: ', &
                                                                slave%name
    bmad_status%ok = .false.
    if (bmad_status%exit_on_error) call err_exit
  endif

! make room in the control structure for slave/lord info

  ref%n_slave = ref%n_slave + 1
  call adjust_control_struct (ring, ix_ref)

  lord%n_slave = lord%n_slave + 1
  call adjust_control_struct (ring, ix_lord)

  slave%n_lord = 2
  call adjust_control_struct (ring, ix_slave)

! put ref/lord info in the control structure

  nrs = ref%ix2_slave
  ring%control_(nrs)%ix_lord   = ix_ref
  ring%control_(nrs)%ix_slave  = ix_slave

  n = slave%ic1_lord
  ring%ic_(n) = nrs

! put lord/slave info in the control structure

  nrs = lord%ix2_slave
  ring%control_(nrs)%ix_lord   = ix_lord
  ring%control_(nrs)%ix_slave  = ix_slave

  n = slave%ic2_lord
  ring%ic_(n) = nrs

! make multipole electric if ref is a sep

  if (ref%key == elseparator$) then
    lord%value(electric$) = 1
    slave%value(electric$) = 1
  endif

end subroutine
