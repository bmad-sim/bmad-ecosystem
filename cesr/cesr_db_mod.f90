module cesr_db_mod

  use cesr_basic_mod
  use group_struct
  use cesr_quad_mod

!---------------------------------------------------------------------------
! DB_STRUCT:                       
! This structrue holds info on the correspondence between CESR data base
! elements and a BMAD ring. Use BMAD_TO_DB to initialize this structure.
! For completeness there are, in addition, arrays for the detectors and
! wigglers, etc.

! for an individual element

  type db_element_struct
    character(16) bmad_name    ! bmad name of element
    real dvar_dcu             ! calibration factor
    real var_theory           ! theory var value
    real var_0                ! extrapolated var value at CU = 0
    integer ix_ring           ! index to element array in ring struct
    integer ix_attrib         ! index to element attribute
    integer ix_cesrv          ! index to cesr_struct arrays
    character(12) db_node_name ! node name ("CSR QUAD CUR")
    integer ix_db             ! element index for data base node (5 for Q05W)
    character(16) db_ele_name  ! element name
    integer cu_high_lim       ! high limit
    integer cu_low_lim        ! low limit
    integer cu_now            ! current CU
    integer ix                ! Integer for general use.
  end type              

! for pointing to a db_struct array

  type db_node_struct
    type (db_element_struct), pointer :: ptr(:) => null()
  end type
                                                       
! db_struct def

  integer, parameter :: n_csr_sqewsext_maxx = 4

  type db_struct
    type (db_element_struct) :: csr_quad_cur(98)
    type (db_element_struct) :: csr_qadd_cur(n_qadd_maxx)
    type (db_element_struct) :: csr_horz_cur(98)
    type (db_element_struct) :: csr_hbnd_cur(n_hbnd_maxx)
    type (db_element_struct) :: csr_vert_cur(98)
    type (db_element_struct) :: csr_hsp_volt(n_sep_maxx)
    type (db_element_struct) :: csr_vsp_volt(n_sep_maxx)
    type (db_element_struct) :: csr_sext_cur(98)
    type (db_element_struct) :: csr_octu_cur(n_oct_maxx)
    type (db_element_struct) :: csr_sqewquad(98)
    type (db_element_struct) :: csr_sqewsext(n_csr_sqewsext_maxx)
    type (db_element_struct) :: scir_quadcur(n_scir_quad_maxx)
    type (db_element_struct) :: scir_skqucur(n_scir_quad_maxx)
    type (db_element_struct) :: scir_vertcur(n_scir_quad_maxx)
    type (db_element_struct) :: scir_sksxcur(1)
    type (db_node_struct) :: node(15)  ! does not include stuff below
! in db but without corresponding BMAD element
    type (db_element_struct) :: scir_pos_stp(n_scir_cam_maxx)
    type (db_element_struct) :: scir_enc_cnt(n_scir_cam_maxx)
    type (db_element_struct) :: scir_pos_rd(3*n_scir_cam_maxx)
! non data base stuff
    type (db_element_struct) :: quad(0:120) ! combinded csr_quad_cur, etc.
    type (db_element_struct) :: detector(0:120)
    type (db_element_struct) :: wiggler(1:n_wig_maxx)
    type (db_element_struct) :: scir_cam_rho(n_scir_cam_maxx)
    type (db_element_struct) :: scir_tilt(n_scir_tilt_maxx)
  end type

! for the synchrotron

  type synch_db_struct
    type (db_element_struct) :: bend(1:100)
    type (db_element_struct) :: detector(1:200)
    type (db_element_struct) :: hkicker(1:200)
    type (db_element_struct) :: vkicker(1:200)
    type (db_element_struct) :: htable(1:200)
    type (db_element_struct) :: vtable(1:200)
  end type

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine bmad_to_db (ring, db, calib_date)
!
! Subroutine to return information on the data base that pertains to
! CESR elements.
!
! Note: This subroutine uses database routines. You must setup communication
! with the database before you use this routine.
!
! Modules needed:
!   use cesr_mod
!
! Input:
!   ring       -- Ring_struct: Ring structure
!   calib_date -- Character(10), optional: Date to use when looking up the 
!                   calibration constants. 
!                   calib_date is in the form: "YYYY-MM-DD".
!                   the present calibrations will be used if calib_date is 
!                   not present or is "NOW"
!
! Output:
!   db   -- DB_struct: Data base structure. See the module file 
!                      BMAD_MODULES.F90 for more details. 
!
! Notes:
!
! %CU_NOW is *NOT* filled in by this subroutine.
!
! To check if an element is in the data base check that %DVAR_DCU is not zero.
! To check if an element has a corresponding ring element check that %IX_RING
!   is not zero.
!
! For completeness the REQ's are included in DB%CSR_QUAD_CUR as 
! elements 0 and 99. 
! These elements have %IX_RING /= 0 but %DVAR_DCU = 0. 
!
! The arrays DB%QUAD_TILT and DB%QADD_TILT are to record the theory rotations
! and do not have corresponding data base elements
!
! The %DVAR_DCU calibrations for the horizontal steerings take +kick as
! radially outward. var units use the BMAD standard so, for example, 
! dvar_dcu for the steerings are radians/CU.
!
! Use the DB%NODE(i)%PTR(1)%DB_NAME array to search for a particular node name
!-

subroutine bmad_to_db (ring, db, calib_date)

  implicit none

  type (ring_struct) ring
  type (db_struct) db
  type (cesr_struct) cesr

  real h_stren(120), v_stren(120), gev, cu_per_k_gev(0:120)
  real k_theory(0:120), k_base(0:120), len_quad(0:120)
  real quad_tilt(0:120), dk_gev_dcu(0:120)

  integer i, ix, ios, cu_theory(0:120), nq100

  character(*), optional :: calib_date
  character(16) type_str
  character(10) this_date

! general init

  nq100 = 100+n_qadd_maxx

  call bmad_to_cesr (ring, cesr)

  call db_init_it (db%csr_quad_cur, 0, 'INIT_IT', 0, db%node)  ! init call

  call db_init_it (db%csr_quad_cur, lbound(db%csr_quad_cur, 1), &
          'CSR QUAD CUR',    k1$,    db%node, cesr%quad_(1:98), 1)

  call db_init_it (db%csr_qadd_cur, lbound(db%csr_qadd_cur, 1), &
          'CSR QADD CUR',    k1$,    db%node, cesr%quad_(101:nq100), 101)

  call db_init_it (db%scir_quadcur, lbound(db%scir_quadcur, 1), &
          'SCIR QUADCUR',    k1$,    db%node, cesr%quad_(111:114), 111)

  call db_init_it (db%csr_horz_cur, lbound(db%csr_horz_cur, 1), &
          'CSR HORZ CUR', hkick$, db%node, cesr%h_steer_(1:98), 1)

  call db_init_it (db%csr_hbnd_cur, lbound(db%csr_hbnd_cur, 1), &
          'CSR HBND CUR', hkick$, db%node, cesr%h_steer_(101:106), 101)

  call db_init_it (db%csr_vert_cur, lbound(db%csr_vert_cur, 1), &
          'CSR VERT CUR', vkick$, db%node, cesr%v_steer_(1:98), 1)

  call db_init_it (db%csr_hsp_volt, lbound(db%csr_hsp_volt, 1), &
          'CSR HSP VOLT', hkick$, db%node)

  call db_init_it (db%csr_vsp_volt, lbound(db%csr_vsp_volt, 1), &
          'CSR VSP VOLT', vkick$, db%node)

  call db_init_it (db%csr_sext_cur, lbound(db%csr_sext_cur, 1), &
          'CSR SEXT CUR',    k2$,    db%node, cesr%sex_(1:98), 1)

  call db_init_it (db%csr_octu_cur, lbound(db%csr_octu_cur, 1), &
          'CSR OCTU CUR',    k3$,    db%node, cesr%oct_, 1)

  call db_init_it (db%csr_sqewquad, lbound(db%csr_sqewquad, 1), &
          'CSR SQEWQUAD',    k1$,    db%node, cesr%skew_quad_(1:98), 1)
                           
  call db_init_it (db%scir_skqucur, lbound(db%scir_skqucur, 1), &
          'SCIR SKQUCUR',    k1$,    db%node, cesr%skew_quad_(111:114), 111)
                           
  call db_init_it (db%csr_sqewsext, lbound(db%csr_sqewsext, 1), &
          'CSR SQEWSEXT',    k2$,    db%node, cesr%skew_sex_(1:4), 1)

  call db_init_it (db%scir_vertcur, lbound(db%scir_vertcur, 1), &
          'SCIR VERTCUR', vkick$, db%node, cesr%v_steer_(111:114), 111)

  call db_init_it (db%scir_sksxcur, lbound(db%scir_sksxcur, 1), &
          'SCIR SKSXCUR',   k1$,    db%node, cesr%skew_sex_(11:11), 11)

!-----------------------------------------------------------------------------
! find the separators

  do i = 1, ring%n_ele_max

    type_str = ring%ele_(i)%type

    if (type_str(:12) == db%csr_vsp_volt(1)%db_node_name) then 
      if (type_str(13:) /= '    ') then
        read (type_str(13:), *, iostat = ios) ix
        if (ios .ne. 0) then
          print *, 'ERROR IN BAMD_TO_DB: READ ERROR FOR NODE INDEX: ', type_str
        endif
        db%csr_vsp_volt(ix)%bmad_name = ring%ele_(i)%name
        db%csr_vsp_volt(ix)%ix_ring = i
      endif
    endif

    if (type_str(:12) == db%csr_hsp_volt(1)%db_node_name) then 
      if (type_str(13:) /= '    ') then
        read (type_str(13:), *, iostat = ios) ix
        if (ios .ne. 0) then
          print *, 'ERROR IN BAMD_TO_DB: READ ERROR FOR NODE INDEX: ', type_str
        endif
        db%csr_hsp_volt(ix)%bmad_name = ring%ele_(i)%name
        db%csr_hsp_volt(ix)%ix_ring = i
      endif
    endif

  enddo

!--------------------------------------------------------------------
! get calibrations...
! steerings
!!  call getcs (v_stren, h_stren)

  this_date = 'NOW'
  if (present(calib_date)) this_date = calib_date

  gev = 1e-9 * ring%ele_(0)%value(beam_energy$)
  if (gev == 0) then
    print *, 'ERROR IN BMAD_TO_DB: BEAM ENERGY IS ZERO!'
    call err_exit
  endif

  call get_steering_calibrations (h_stren, v_stren, this_date)
  db%csr_horz_cur(:)%dvar_dcu = 1.0e-6 * h_stren(1:98) / gev  
  db%csr_hbnd_cur(:)%dvar_dcu = 1.0e-6 * h_stren(101:106) / gev  
  db%csr_vert_cur(:)%dvar_dcu = 1.0e-6 * v_stren(1:98) / gev  
  db%scir_vertcur(:)%dvar_dcu = 1.0e-6 * v_stren(111:114) / gev  

! separator calibration

  call get_sep_strength (db%csr_hsp_volt(:)%dvar_dcu, &
                                       db%csr_vsp_volt(:)%dvar_dcu)
  db%csr_hsp_volt(:)%dvar_dcu = 1.0e-6 * db%csr_hsp_volt(:)%dvar_dcu / gev
  db%csr_vsp_volt(:)%dvar_dcu = 1.0e-6 * db%csr_vsp_volt(:)%dvar_dcu / gev

  call get_ele_theory (ring, db%csr_hsp_volt)
  call get_ele_theory (ring, db%csr_vsp_volt)

! sextupole calibration

  call sext_calib (gev, db%csr_sext_cur(:)%dvar_dcu)
  call get_ele_theory (ring, db%csr_sext_cur)

! octupole calibration

  db%csr_octu_cur(:)%dvar_dcu = 0.213 / gev
  call get_ele_theory (ring, db%csr_octu_cur)

! skew quad calibration

  call skew_quad_calib (db%csr_sqewquad(:)%dvar_dcu, &
                                             db%scir_skqucur(:)%dvar_dcu)
  db%csr_sqewquad(:)%dvar_dcu = db%csr_sqewquad(:)%dvar_dcu / gev
  db%scir_skqucur(:)%dvar_dcu = db%scir_skqucur(:)%dvar_dcu / gev

  call get_ele_theory (ring, db%csr_sqewquad) 
  call get_ele_theory (ring, db%scir_skqucur) 

! skew sextupole calibration

  call skew_sex_calib (db%csr_sqewsext(:)%dvar_dcu, &
                                            db%scir_sksxcur(:)%dvar_dcu)
  db%csr_sqewsext(:)%dvar_dcu = db%csr_sqewsext(:)%dvar_dcu / gev
  db%scir_sksxcur(:)%dvar_dcu = db%scir_sksxcur(:)%dvar_dcu / gev

  call get_ele_theory (ring, db%csr_sqewsext) 
  call get_ele_theory (ring, db%scir_sksxcur) 

! quad calibrations 

  call ring_to_quad_calib (ring, cesr, k_theory, k_base, len_quad, &
                              cu_per_k_gev, quad_tilt, dk_gev_dcu, cu_theory)

  call non_db_set (db%quad, cesr%quad_, k1$, 0)
  db%quad(1:98) = db%csr_quad_cur 
  db%quad(101:nq100) = db%csr_qadd_cur
  db%quad(111:114) = db%scir_quadcur(:)

  db%quad(:)%var_theory = k_theory
  db%quad(:)%dvar_dcu = dk_gev_dcu / gev
  db%quad(:)%var_0 = db%quad(:)%var_theory - &
                                       cu_theory * db%quad(:)%dvar_dcu

  db%csr_quad_cur = db%quad(1:98)
  db%csr_qadd_cur = db%quad(101:nq100)
  db%scir_quadcur = db%quad(111:114)

! Things not part of the CESR DB

  call non_db_set (db%detector, cesr%det_, 0, 0)
  call non_db_set (db%wiggler, cesr%wig_, b_max$, 1)
  call non_db_set (db%scir_cam_rho, cesr%scir_cam_rho_, rho$, 1)
  call non_db_set (db%scir_tilt, cesr%scir_tilt_, tilt$, 1)


  deallocate( cesr%ix_cesr )

!---------------------------------------------------------------------
contains
                        
subroutine db_init_it (node, n1, node_name, ix_attrib, node_array, &
                                                      cesr_ele, ix0_cesrv)

  implicit none

  integer ix_attrib, n1, i, vnumbr, iv, ix_(200)

  type (db_element_struct), target :: node(n1:)
  type (db_node_struct) :: node_array(:)
  type (cesr_element_struct), optional :: cesr_ele(:)

  character(*) node_name

  integer, optional :: ix0_cesrv
  integer, save :: i_array

! init array index at the start

  if (node_name == 'INIT_IT') then
    i_array = 0
    return
  endif

!

  i_array = i_array + 1
  node_array(i_array)%ptr => node

  do i = lbound(node, 1), ubound(node, 1)
    node(i)%bmad_name = ' '
    node(i)%ix_ring = 0
    node(i)%ix_attrib = ix_attrib
    node(i)%db_node_name = node_name       
    node(i)%ix_db = i
    node(i)%db_ele_name = ' '           
    node(i)%dvar_dcu = 0
    node(i)%var_theory = 0
    node(i)%cu_high_lim = 0
    node(i)%cu_low_lim = 0
  enddo

  iv = vnumbr(node_name)
  if (size(node) < iv) then
    print *, 'WARNING IN BMAD_TO_DB: ARRAY IS TOO SMALL TO HOLD DATA BASE'
    print *, '        NODE ELEMENTS: ', node_name
    iv = size(node)
  endif

  if (iv .gt. 0) then
    call vmglim (node_name, &
                   1, iv, node(1:iv)%cu_low_lim, node(1:iv)%cu_high_lim)
    call vstget (node_name, 1, iv, node(1:iv)%db_ele_name, ix_)
  endif

  if (present(cesr_ele)) then
    if (ubound(node, 1)-n1 .ne. ubound(cesr_ele, 1)-1) then
      print *, 'ERROR IN DB_INIT_IT: ARRAYS DO NOT MATCH FOR: ', node_name
      call err_exit
    endif
    node(:)%ix_ring = cesr_ele(:)%ix_ring
    node(:)%bmad_name = cesr_ele(:)%name
    do i = lbound(node, 1), ubound(node, 1)
      node(i)%ix_cesrv = i - lbound(node, 1) + ix0_cesrv 
    enddo
  else
    do i = lbound(node, 1), ubound(node, 1)
      node(i)%ix_cesrv = i - lbound(node, 1) + 1
    enddo
  endif

  return

9000 print *
  print *, 'ERROR IN DB_INIT_IT: CANNOT READ INDEX: ', cesr_ele(i)%db_node_name
  call err_exit

end subroutine db_init_it
                         
!---------------------------------------------------------------------
! contains
                            
subroutine get_ele_theory (ring, db_ele)
                                                               
  implicit none

  type (db_element_struct) :: db_ele(:)
  type (ring_struct) :: ring

  integer i, ie, iv

!

  do i = 1, ubound(db_ele, 1)
    ie = db_ele(i)%ix_ring
    if (ie /= 0) then
      iv = db_ele(i)%ix_attrib
      db_ele(i)%var_theory = ring%ele_(ie)%value(iv)
    endif
  enddo

end subroutine get_ele_theory

!---------------------------------------------------------------------
! contains

subroutine non_db_set (db_ele, cesr_ele, ix_attrib, ix0_cesrv)

  implicit none

  type (db_element_struct) :: db_ele(:)
  type (cesr_element_struct) :: cesr_ele(:)

  integer i, ix0_cesrv, ix_attrib

  db_ele(:)%bmad_name = cesr_ele(:)%name
  db_ele(:)%ix_ring = cesr_ele(:)%ix_ring
  db_ele(:)%ix_db = 0
  db_ele(:)%ix_attrib = ix_attrib
  db_ele(:)%db_node_name = '----------'
  db_ele(:)%db_ele_name  = '----------'
  db_ele(:)%dvar_dcu = 0
  db_ele(:)%var_theory = 0
  db_ele(:)%var_0 = 0
  db_ele(:)%cu_high_lim = 0
  db_ele(:)%cu_low_lim = 0
  forall (i = 1:size(db_ele, 1)) db_ele(i)%ix_cesrv = i + ix0_cesrv - 1

end subroutine non_db_set

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine db_group_to_bmad (ing_name, ing_num, biggrp_set, csr_set, &
!                                      ring, db, con_, n_con, ok, type_err)
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
!   ring       -- Ring_struct: BMAD ring to use.
!   db         -- Db_struct: Info about ring obtained with a 
!                            "call bmad_to_db (ring, db)"
!   type_err   -- Logical: If true then error messages are typed.
!
! Output:
!   con_(:) -- Control_struct: Control structure.
!     %ix_slave    -- Index to ring element controlled.
!     %ix_attrib -- Index of attribute controlled.
!     %coef      -- Control coefficient
!   n_con   -- Integer: number of con_(:) elements used by the group.
!   ok      -- Logical: True if everything was OK.
!
! Note: See also DB_GROUP_TO_BMAD_GROUP which uses the results of this
! subroutine to create a BMAD group in the ring structure.
!-

subroutine db_group_to_bmad (ing_name, ing_num, biggrp_set, csr_set, &
                                      ring, db, con_, n_con, ok, type_err)


  implicit none

  type (ring_struct) ring
  type (db_struct) db
  type (group_info_struct) grp
  type (control_struct) con_(:)
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
      if (db_ptr(n)%ix_ring == 0) then
        if (type_err) then
          print *, 'ERROR IN DB_GROUP_TO_BMAD: CANNOT FIND A LATTICE ELEMENT'
          print *, '      CORRESPONDING TO: ', db_ptr(n)%db_node_name, db_ptr(n)%ix_db
          print *, '      FOR GROUP:        ', ing_name, ing_num
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

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+     
! Subroutine db_group_to_bmad_group (group_name, group_num, 
!                   biggrp_set, csr_set, ring, db, ix_ele, ok, type_err)
!
! Subroutine to set up a data base group knob in a bmad ring structure.
!
! Modules used:
!   use cesr_mod
!
! Input:
!   group_name -- Character(12): Group node name (eg. "CSR PRETZING")
!   group_num  -- Integer: Group node number
!   biggrp_set -- Integer: Biggrp number. 0 => read from CESR database
!   csr_set    -- Integer: CSR set number. 0 => read from the data base
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

subroutine db_group_to_bmad_group (group_name, group_num, &
                   biggrp_set, csr_set, ring, db, ix_ele, ok, type_err)

  implicit none              

  type (ring_struct) ring                                       
  type (control_struct) con_(100)
  type (db_struct) db

  integer n_con, group_num, ix_ele, biggrp_set, ix, ixs(1), csr_set
  character(12) group_name
  logical ok, type_err

!                                          

  ix_ele = -1

  call db_group_to_bmad (group_name, group_num, biggrp_set, csr_set, &
                                         ring, db, con_, n_con, ok, type_err)
  if (.not. ok) return
  call new_control (ring, ix_ele)

  call vstget (group_name, group_num, group_num, ring%ele_(ix_ele:ix_ele)%name, ixs)

  write (ring%ele_(ix_ele)%type, '(a, i4)') group_name, group_num

  ring%ele_(ix_ele)%name = ring%ele_(ix_ele)%type
  do
    ix = index(ring%ele_(ix_ele)%name, ' ')
    if (ix == 0) exit
    ring%ele_(ix_ele)%name(ix:ix) = '_'
  enddo

  call create_group (ring, ix_ele, con_(1:n_con))

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

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine READ_BUTNS_FILE (butns_num, nonlinear_calc, butns, db, 
!                                                         read_ok, type_err)
!
! Subroutine to read the raw information from a BUTNS.nnnnn file and convert
! it to orbit x-y positions.
!
! Modules Needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!   butns_num      -- Integer: Number in BUTNS.nnnnn file name.
!   nonlinear_calc -- Logical: Calculate orbit using Rich Helms nonlinear 
!                       routines? Recomendation: True.
!   type_err       -- Logical: If True then open error message will be printed.
!                         If False then open error message will be supressed.
!
! Output:
!   butns -- Butns_struct: Orbit information.
!       %lattice    -- Character(40): Lattice name.
!       %save_set   -- Integer: Save set number.
!       %date       -- Character(20): Date orbit was taken
!       %file_num   -- Integer: Equal to butns_num.
!       %turn       -- Integer: Turn number for injection orbits. 0 otherwise.
!       %comment(5) -- Character(72): Comment.
!       %det(0:99)%amp(4) -- Integer: raw button numbers.
!       %det(0:99)%x_orb  -- Real: Horizontal orbit in meters.
!       %det(0:99)%y_orb  -- Real: Horizontal orbit in meters.
!   db    -- Db_struct: Structure holding the steering settings.
!     %csr_horz_cur(i)%cu_now -- CU settings for CSR HORZ CUR
!     %csr_hbnd_cur(i)%cu_now -- CU settings for CSR HBND CUR
!     %csr_vert_cur(i)%cu_now -- CU settings for CSR VERT CUR
!     %csr_hsp_volt(i)%cu_now -- CU settings for CSR HSP VVAL
!     %csr_vsp_volt(i)%cu_now -- CU settings for CSR VSP VOLT
!     %scir_vertcur(i)%cu_now -- CU settings for SCIR VERTCUR
!     %scir_pos_stp(i)%cu_now -- CU settings for SCIR POS STP
!     %scir_enc_cnt(i)%cu_now -- CU settings for SCIR ENC CNT
!     %scir_pos_rd(i)%cu_now  -- CU settings for SCIR POS RD
!   read_ok -- Logical: Set True if butns file was successfuly parsed.
!
! Note: orbit numbers are from 0 to 99
!
! Note: db%csr_hsp_volt is actually obtained from the node CSR HSP VVAL which
! records the actual voltage as opposed to the command. Since CSR HSP VVAL
! is a readback node the values put in db%csr_hsp_volt will not be exactly the
! actual commands (and if the separators have been turned off they will not
! even be approximately the same).
!-

subroutine read_butns_file (butns_num, nonlinear_calc, butns, db, &
                                                          read_ok, type_err)

  implicit none

  type (db_struct) db
  type (butns_struct) butns

  integer vec(120), det_type(120)
  integer i, ix, j, butns_num, iu, ios, raw(4, 120)
  integer n_node, n_ele

  real x_orbit(120), y_orbit(120), rdummy

  character line_in*130, file_in*40

  logical nonlinear_calc, read_ok, type_err, err_flag, is_ok(120)

! init comments

  butns%comment = ' '

! compute filename

  read_ok = .false.

  call form_file_name_with_number ('ORBIT', butns_num, file_in, err_flag)
  if (err_flag) return

  butns%file_num = butns_num

! read header line in the data file (to get lattice name)
! and sterring strengths

  iu = lunget()
  open(unit = iu, file = file_in, status = 'old', action = 'READ', iostat = ios)
  if (ios /= 0) then
    if (type_err) print *, &
          'ERROR IN READ_BUTNS_FILE: ERROR OPENING: ', trim(file_in)
    return
  endif

  read (iu, '(a)') line_in
  butns%lattice = line_in(61:)
  butns%date = line_in(30:)                     ! get date
  read (line_in(54:), *) butns%save_set

  do
    read (iu, '(a)', iostat = ios) line_in
    if (ios /= 0) then
      print *, 'ERROR IN ORBIT_READ: ERROR READING STEERINGS IN ORBIT FILE'
      goto 1000
    endif
    ix = index(line_in, 'END BUTNS')
    if (ix /= 0) exit
    ix = index(line_in, 'TURN=')
    if (ix /= 0) then
      read (line_in(ix+5:), *, iostat = ios) butns%turn
      if (ios /= 0) then
        print *, 'ERROR IN READ_BUTNS_FILE: ERROR READING TURN NUMBER.'
        print *, '      IN FILE: ', trim(file_in)
      endif
    endif
  enddo

  read (line_in(ix+10:), *, iostat = ios) n_node
  if (ios /= 0) then
    print *, 'ERROR IN ORBIT_READ: ERROR READING STEERING NODE NUMBER IN ORBIT FILE'
    goto 1000
  endif

! read data base valuse stored in the orbit file

  do i = 1, n_node

    read (iu, '(a)', iostat = ios) line_in
    if (ios /= 0) then
      print *, 'ERROR IN ORBIT_READ: ERROR READING STEERING NODE IN ORBIT FILE'
      exit
    endif

    read (line_in(14:), *, iostat = ios) n_ele

    if (line_in(2:13) == 'CSR COMMENTS') then
      do j = 1, n_ele
        read (iu, '(1x, a)', iostat = ios) butns%comment(j)
        if (ios /= 0) then
          print *, 'ERROR IN ORBIT_READ: ERROR READING COMMENT #', j
          exit
        endif
      enddo
      cycle
    endif

    read (iu, '(12x, 10i6)') vec(1:n_ele)

    if (line_in(2:13) == 'CSR HORZ CUR') then
      db%csr_horz_cur(1:n_ele)%cu_now = vec(1:n_ele)
    elseif (line_in(2:13) == 'CSR VERT CUR') then
      db%csr_vert_cur(1:n_ele)%cu_now = vec(1:n_ele)
    elseif (line_in(2:13) == 'CSR HBND CUR') then
      db%csr_hbnd_cur(1:n_ele)%cu_now = vec(1:n_ele)
    elseif (line_in(2:13) == 'CSR HSP VVAL') then
      call hsp_vval_to_volt (vec, vec)
      n_ele = size(db%csr_hsp_volt)
      db%csr_hsp_volt(1:n_ele)%cu_now = vec(1:n_ele)
    elseif (line_in(2:13) == 'CSR VSP VOLT') then
      db%csr_vsp_volt(1:n_ele)%cu_now = vec(1:n_ele)
    elseif (line_in(2:13) == 'SCIR VERTCUR') then
      db%scir_vertcur(1:n_ele)%cu_now = vec(1:n_ele)
    elseif (line_in(2:13) == 'SCIR POS STP') then
      db%scir_pos_stp(1:n_ele)%cu_now = vec(1:n_ele)
    elseif (line_in(2:13) == 'SCIR ENC CNT') then
      db%scir_enc_cnt(1:n_ele)%cu_now = vec(1:n_ele)
    elseif (line_in(2:13) == 'SCIR POS RD ') then
      db%scir_pos_rd(1:n_ele)%cu_now = vec(1:n_ele)
    else
      print *, 'ERROR IN ORBIT_READ: UNKNOWN NODE IN ORBIT FILE: ', line_in(2:13)
      goto 1000
    endif

  enddo

! close

  1000 continue
  close (unit = iu)

! read in the raw data

  call butfilget (raw, file_in, rdummy, det_type)      ! read in raw data
  if (nonlinear_calc) then
    call nonlin_butcon (raw, 1, 100, y_orbit, x_orbit)
  else
    call butcon (raw, 1, 100, y_orbit, x_orbit)
  endif

  is_ok = .false.
  call det_ok (raw, 1, 100, det_type, is_ok)

  do i = 1, 100
    j = i
    if (i == 100) j = 0
    butns%det(j)%amp = raw(1:4, i)
    butns%det(j)%x_orb = x_orbit(i) / 1000.0   ! convert to m
    butns%det(j)%y_orb = y_orbit(i) / 1000.0   ! convert to m
    butns%det(j)%ok = is_ok(i)
  end do

  read_ok = .true.

end subroutine

end module
