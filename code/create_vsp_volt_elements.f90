!+
! Subroutine CREATE_VSP_VOLT_ELEMENTS (RING, ELE_TYPE)
!
! Subroutine to create elements corresponding to the 6 data base elements
! in CSR VSP VOLT. For each vertical separator 3 controller (lord) elements are
! created corresponding to CSR VSP VOLT 1, 3, and 4 for the west vsep and
! 2, 5, and 6 for the east vsep. The controllers control the VKICK attribute
! with coefficients of 1.0. If ELE_TYPE = OVERLAY$ then any VKICK$ in the
! vertical seps will be divided between elemtns 3 and 4 for the west and 5 and
! 6 for the east.
!
! Use BMAD_TO_DB or BMAD_TO_CESR to find where the elements are located.
!
! Modules Needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!     RING     -- Ring_struct: Ring to be modified
!     ELE_TYPE -- Integer: Type of elements to make
!                   = group$      ! make group controller elements
!                   = overlay$    ! make overlay controller elements
!
! Output:
!     RING -- Ring_struct: Modified ring.
!-

!$Id$
!$Log$
!Revision 1.3  2001/12/04 20:28:57  helms
!Changes from DCS
!
!Revision 1.2  2001/09/27 18:31:50  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"



subroutine create_vsp_volt_elements (ring, ele_type)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ring_struct)  ring

  integer ele_type, ix_west(3) / 1, 3, 4 /, ix_east(3) / 2, 5, 6 /
  integer i

  character*16 vsep_west / 'V_SEP_48W' /, vsep_east / 'V_SEP_48E' /

  logical found_west

! find vseps

  found_west = .false.

  do i = 1, ring%n_ele_ring

    if (ring%ele_(i)%name == vsep_west) then
      found_west = .true.
      call do_vsp_eles (ring, i, ix_west, ele_type)
                 
    elseif (ring%ele_(i)%name == vsep_east) then
      call do_vsp_eles (ring, i, ix_east, ele_type)
      if (.not. found_west) then
        type *, 'ERROR IN CREATE_VSP_VOLT_ELEMENTS: CANNOT FIND WEST VSEP!'
        if (bmad_status%exit_on_error) call err_exit
        bmad_status%ok = .false.
      endif
      return
    endif

  enddo

  type *, 'ERROR IN CREATE_VSP_VOLT_ELEMENTS: CANNOT FIND EAST VSEP!'
  if (bmad_status%exit_on_error) call err_exit
  bmad_status%ok = .false.

  return

  end


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------


subroutine do_vsp_eles (ring, i_vsep, ix_, ele_type)

  use bmad_struct
  use bmad_interface
  implicit none

  type (ring_struct)  ring
  type (control_struct)  con_(1)

  integer i_vsep, ix_(3), ele_type, i, i_con

  real vkick

!
                 
  if (ring%ele_(i_vsep)%control_type /= free$) then
    type *, 'ERROR IN CREATE_VSP_VOLT_ELEMENTS: VSEP NOT FREE!', i_vsep
    return
  endif

  ring%ele_(i_vsep)%type = ' '

  con_(1)%ix_attrib = vkick$
  con_(1)%coef = 1.0
  con_(1)%ix_slave = i_vsep
  vkick = ring%ele_(i_vsep)%value(vkick$)

  do i = 1, 3

    call new_control (ring, i_con)
    write (ring%ele_(i_con)%name, '(a, i1)') 'VSP_VOLT_', ix_(i)
    write (ring%ele_(i_con)%type, '(a, i4)') 'CSR VSP VOLT', ix_(i)

    if (ele_type == group$) then
      call create_group (ring, i_con, 1, con_)
    elseif (ele_type == overlay$) then
      call create_overlay (ring, i_con, vkick$, 1, con_)
      if (i == 2 .or. i == 3) ring%ele_(i_con)%value(vkick$) = vkick / 2
    else
      type *, 'ERROR IN CREATE_VSP_VOLT_ELEMENTS: BAD ELE_TYPE: ', ele_type
      call err_exit
    endif

  enddo

  return
  end











