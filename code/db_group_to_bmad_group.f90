!+     
! Subroutine db_group_to_bmad_group (group_name, group_num, i_biggrp, 
!                                               ring, db, ix_ele, ok, type_err)
!
! Subroutine to set up a data base group knob in a bmad ring structure.
!
! Modules used:
!   use bmad
!
! Input:
!   group_name -- Character*12: Group node name (eg. "CSR PRETZING")
!   group_num  -- Integer: Group node number
!   i_biggrp   -- Integer: Biggrp number. 0 => read from CESR database
!   ring       -- Ring_struct: Ring to modify
!   db         -- Db_struct: Info about ring obtained with a 
!                            "call bmad_to_db (ring, db)"
!   type_err   -- Logical: If true then error messages are typed.
!
! Output:
!   ring    -- Ring_struct: Modified ring.
!   ix_ele  -- Integer: RING%ELE_(IX_ELE) is the group element
!   ok      -- Logical: True if everything was OK.
!
! Note: See also DB_GROUP_TO_BMAD.
!-

!$Id$
!$Log$
!Revision 1.4  2003/01/27 14:40:33  dcs
!bmad_version = 56
!
!Revision 1.3  2002/02/23 20:32:14  dcs
!Double/Single Real toggle added
!
!Revision 1.2  2001/09/27 18:31:51  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine db_group_to_bmad_group (group_name, group_num, i_biggrp, &
                                       ring, db, ix_ele, ok, type_err)

  use bmad_struct
  use bmad_interface

  implicit none              

  type (ring_struct) ring                                       
  type (control_struct) con_(100)
  type (db_struct) db

  integer n_con, k, j, n, nn, ix
  integer group_num, ix_ele, i_biggrp

  character*12 group_name
                                                
  logical ok, type_err

!                                          

  ix_ele = -1

  call db_group_to_bmad (group_name, group_num, i_biggrp, ring, db, &
                                                    con_, n_con, ok, type_err)
  if (.not. ok) return
  call new_control (ring, ix_ele)

  call vstget (group_name, group_num, group_num, ring%ele_(ix_ele)%name, ix)

  write (ring%ele_(ix_ele)%type, '(a, i4)') group_name, group_num

  ring%ele_(ix_ele)%name = ring%ele_(ix_ele)%type
  do
    ix = index(ring%ele_(ix_ele)%name, ' ')
    if (ix == 0) exit
    ring%ele_(ix_ele)%name(ix:ix) = '_'
  enddo

  call create_group (ring, ix_ele, n_con, con_)

end subroutine
                                                             
