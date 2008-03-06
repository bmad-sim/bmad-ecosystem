module cesr_group_mod

use group_struct
use cesr_db_mod

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine db_group_to_bmad (ing_name, ing_num, biggrp_set, csr_set, &
!                                      lat, db, con, n_con, ok, type_err)
!
! Subroutine to take a data base group element and find the elements
! controlled along with the coefficients.
!
! Modules used:
!   use cesr_mod
!
! Input:
!   ing_name   -- Character(12): DB node name (e.g. 'CSR PRETZING')
!   ing_num    -- Integer: DB element number (e.g. 13)
!   biggrp_set -- Integer: Biggrp set number. 0 => read from the data base
!   csr_set    -- Integer: CSR set number. 0 => read from the data base
!   lat       -- lat_struct: BMAD lat to use.
!   db         -- Db_struct: Info about lat obtained with a 
!                            "call bmad_to_db (lat, db)"
!   type_err   -- Logical: If true then error messages are typed.
!
! Output:
!   con(:) -- Control_struct: Control structure.
!     %ix_slave    -- Index to lat element controlled.
!     %ix_attrib -- Index of attribute controlled.
!     %coef      -- Control coefficient
!   n_con   -- Integer: number of con(:) elements used by the group.
!   ok      -- Logical: True if everything was OK.
!
! Note: See also DB_GROUP_TO_BMAD_GROUP which uses the results of this
! subroutine to create a BMAD group in the lat structure.
!-

subroutine db_group_to_bmad (ing_name, ing_num, biggrp_set, csr_set, &
                                      lat, db, con, n_con, ok, type_err)


  implicit none

  type (lat_struct) lat
  type (db_struct) db
  type (group_info_struct) grp
  type (control_struct) con(:)
  type (db_element_struct), pointer :: db_ptr(:)

  integer k, n, nn, n_con, ing_num, biggrp_set, csr_set

  character(12) ing_name

  logical ok, type_err

!

  call setup_group (ing_name, ing_num, biggrp_set, csr_set, grp, ok, .true.)

  if (.not. ok) return
  if (grp%ratio_mode) then
    if (type_err) print *, 'ERROR: RATIO MODE NOT SUPPORTED.'
    ok = .false.
    return
  endif

  n_con = 0

  do k = 1, grp%n_node

    call identify_db_node (grp%node(k)%name, db, db_ptr, ok, .true.)

    if (.not. ok) then
      if (type_err) print *, 'ERROR IN DB_GROUP_TO_BMAD:', &
                                   ' CANNOT FIND NODE: ', grp%node(k)%name
      return
    endif

    do n = 1, grp%node(k)%n_ele
      nn = n + grp%node(k)%offset
      if (grp%ele(nn)%coef == 0) cycle
      if (db_ptr(n)%ix_lat == 0) then
        if (type_err) then
          print *, 'ERROR IN DB_GROUP_TO_BMAD: CANNOT FIND A LATTICE ELEMENT'
          print *, '      CORRESPONDING TO: ', db_ptr(n)%db_node_name, db_ptr(n)%ix_db
          print *, '      FOR GROUP:        ', ing_name, ing_num
        endif
        ok = .false.               
        return
      endif
      n_con = n_con + 1
      con(n_con)%ix_slave = db_ptr(n)%ix_lat
      con(n_con)%ix_attrib = db_ptr(n)%ix_attrib
      con(n_con)%coef = grp%ele(nn)%coef * db_ptr(n)%dvar_dcu  
    enddo

  enddo

end subroutine 

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+     
! Subroutine db_group_to_bmad_group (group_name, group_num, 
!                   biggrp_set, csr_set, lat, db, ix_ele, ok, type_err)
!
! Subroutine to set up a data base group knob in a bmad lat structure.
!
! Modules used:
!   use cesr_mod
!
! Input:
!   group_name -- Character(12): Group node name (eg. "CSR PRETZING")
!   group_num  -- Integer: Group node number
!   biggrp_set -- Integer: Biggrp number. 0 => read from CESR database
!   csr_set    -- Integer: CSR set number. 0 => read from the data base
!   lat       -- lat_struct: Lat to modify
!   db         -- Db_struct: Info about lat obtained with a 
!                            "call bmad_to_db (lat, db)"
!   type_err   -- Logical: If true then error messages are typed.
!
! Output:
!   lat    -- lat_struct: Modified lat.
!   ix_ele  -- Integer: LAT%ele(IX_ELE) is the group element
!   ok      -- Logical: True if everything was OK.
!
! Note: See also DB_GROUP_TO_BMAD.
!-

subroutine db_group_to_bmad_group (group_name, group_num, &
                   biggrp_set, csr_set, lat, db, ix_ele, ok, type_err)

  implicit none              

  type (lat_struct) lat                                       
  type (control_struct) con(150)
  type (db_struct) db

  integer n_con, group_num, ix_ele, biggrp_set, ix, ixs(1), csr_set
  integer j, endj
  character(12) group_name
  logical ok, type_err, err

!                                          

  ix_ele = -1

  call db_group_to_bmad (group_name, group_num, biggrp_set, csr_set, &
                                         lat, db, con, n_con, ok, type_err)
  if (.not. ok) return
  call new_control (lat, ix_ele)

  call vstget (group_name, group_num, group_num, lat%ele(ix_ele:ix_ele)%name, ixs)

  write (lat%ele(ix_ele)%type, '(a, i4)') group_name, group_num

  lat%ele(ix_ele)%name = lat%ele(ix_ele)%type
  endj = len(lat%ele(ix_ele)%name)
  if (endj > 16) endj = 16
  do j = 1, endj
     if (lat%ele(ix_ele)%name(j:j) == ' ') lat%ele(ix_ele)%name(j:j) = '_'
  enddo

  call create_group (lat, ix_ele, con(1:n_con), err)
  if (err) ok = .false.

end subroutine
                                                             
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine identify_db_node (db_name, db, dp_ptr, ok, type_err)
!
! Subroutine to find which array in DB is associated with DB_NAME.
!
! Modules Needed:
!   use cesr_mod
!
! Input:
!   db_name -- Character(12): Data base name (eg. "CSR HSP VOLT")
!   db      -- Db_node_struct: Data base structure.
!  
! Output:
!   db_ptr   -- Db_element_struct, pointer(:) : Pointer to, eg., DB%HSP_VOLT(:)
!   ok       -- Logical: Set True if everything ok
!   type_err -- Logical: If True then error message will be typed if needed.
!-

subroutine identify_db_node (db_name, db, db_ptr, ok, type_err)

  implicit none

  type (db_struct), target :: db
  type (db_element_struct), pointer :: db_ptr(:)

  integer k

  character(*) db_name

  logical ok, type_err

!

  do k = 1, size(db%node)
    if (db%node(k)%ptr(1)%db_node_name /= db_name) cycle
    db_ptr => db%node(k)%ptr
    ok = .true.
    return
  enddo

  if (type_err) print *, &
        'ERROR IN IDENTIFY_DB_NODE: CANNOT FIND DATABASE NODE: ', db_name
  ok = .false.
  return

end subroutine

end module
