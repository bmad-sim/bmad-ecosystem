!+
! Subroutine db_group_to_bmad (ing_name, ing_num, biggrp_set, ring, db, &
!                                                con_, n_con, ok, type_err)
!
! Subroutine to take a data base group element and find the elements
! controlled along with the coefficients.
!
! Modules used:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!   ing_name   -- Character*12: DB node name (e.g. 'CSR PRETZING')
!   ing_num    -- Integer: DB element number (e.g. 13)
!   biggrp_set -- Integer: Biggrp set number. 0 => read from the data base
!   ring       -- Ring_struct: BMAD ring to use.
!   db         -- Db_struct: Info about ring obtained with a 
!                            "call bmad_to_db (ring, db)"
!   type_err   -- Logical: If true then error messages are typed.
!
! Output:
!   con_(*) -- Control_struct: Control structure.
!     %ix_slave    -- Index to ring element controlled.
!     %ix_attrib -- Index of attribute controlled.
!     %coef      -- Control coefficient
!   n_con   -- Integer: number of con_(*) elements used by the group.
!   ok      -- Logical: True if everything was OK.
!
! Note: See also DB_GROUP_TO_BMAD_GROUP which uses the results of this
! subroutine to create a BMAD group in the ring structure.
!-

!$Id$
!$Log$
!Revision 1.2  2001/09/27 18:31:50  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine db_group_to_bmad (ing_name, ing_num, biggrp_set, ring, db, &
                                              con_, n_con, ok, type_err)

  use bmad_struct
  use bmad_interface

  implicit none

  #include "group.inc"

  type (ring_struct) ring
  type (db_struct) db
  type (group_info_struct) grp
  type (control_struct) con_(*)
  type (db_element_struct), pointer :: db_ptr(:)

  integer k, n, nn, n_con, ing_num, biggrp_set

  character*12 ing_name

  logical ok, type_err

!

  call setup_group (ing_name, ing_num, 0, grp, ok, .true.)

  if (.not. ok) return
  if (grp%ratio_mode) then
    if (type_err) type *, 'ERROR: RATIO MODE NOT SUPPORTED.'
    ok = .false.
    return
  endif

  n_con = 0

  do k = 1, grp%n_node

    call identify_db_node (grp%node(k)%name, db, db_ptr, ok, .true.)

    if (.not. ok) then
      if (type_err) type *, 'ERROR IN DB_GROUP_TO_BMAD:', &
                                   ' CANNOT FIND NODE: ', grp%node(k)%name
      return
    endif

    do n = 1, grp%node(k)%n_ele
      nn = n + grp%node(k)%offset
      if (grp%ele(nn)%coef == 0) cycle
      if (db_ptr(n)%ix_ring == 0) then
        if (type_err) then
          type *, 'ERROR IN DB_GROUP_TO_BMAD: CANNOT FIND IN RING ELEMENT'
          type *, '      FOR: ', db_ptr(n)%db_node_name, db_ptr(n)%ix_db
        endif
        ok = .false.               
        return
      endif
      n_con = n_con + 1
      con_(n_con)%ix_slave = db_ptr(n)%ix_ring
      con_(n_con)%ix_attrib = db_ptr(n)%ix_attrib
      con_(n_con)%coef = grp%ele(nn)%coef * db_ptr(n)%dvar_dcu  
    enddo

  enddo

end subroutine          
