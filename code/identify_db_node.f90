!+
! Subroutine IDENTIFY_DB_NODE (DB_NAME, DB, DP_PTR, OK, TYPE_ERR)
!
! Subroutine to find which array in DB is associated with DB_NAME.
!
! Modules Needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!   DB_NAME -- Character*12: Data base name (eg. "CSR HSP VOLT")
!   DB      -- Db_node_struct: Data base structure.
!  
! Output:
!   DB_PTR   -- Db_element_struct, pointer(:) : Pointer to, eg., DB%HSP_VOLT(:)
!   OK       -- Logical: Set True if everything ok
!   TYPE_ERR -- Logical: If True then error message will be typed if needed.
!-

!$Id$
!$Log$
!Revision 1.2  2001/09/27 18:31:52  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine identify_db_node (db_name, db, db_ptr, ok, type_err)

  use bmad_struct
  use bmad_interface

  implicit none

  type (db_struct), target :: db
  type (db_element_struct), pointer :: db_ptr(:)

  integer k

  character*(*) db_name

  logical ok, type_err

!

  do k = 1, size(db%node)
    if (db%node(k)%ptr(1)%db_node_name /= db_name) cycle
    db_ptr => db%node(k)%ptr
    ok = .true.
    return
  enddo

  if (type_err) type *, &
        'ERROR IN IDENTIFY_DB_NODE: CANNOT FIND DATABASE NODE: ', db_name
  ok = .false.
  return

end subroutine
