!+
! Module: cesr_mod
!
! This module holds all the CESR specific stuff in BMAD.
!-

#include "CESR_platform.inc"

module cesr_mod

  use bmad_struct
  use bmad_interface
  use group_struct

!-------------------------------------------------------------------------
! CESR logical names

  integer, parameter :: q49aw$ = 101, q47aw$ = 102
  integer, parameter :: q47ae$ = 103, q49ae$ = 104
  integer, parameter :: q43aw$ = 105, q43ae$ = 106
  integer, parameter :: q08aw$ = 107

  integer, parameter :: h_sep_08w$ = 1, h_sep_45w$ = 2
  integer, parameter :: h_sep_45e$ = 3, h_sep_08e$ = 4
  integer, parameter :: v_sep_48w$ = 5, v_sep_48e$ = 6

  integer, parameter :: wig_w$ = 1, wig_e$ = 2

  integer, parameter :: rf_w1$ = 1, rf_w2$ = 2, rf_e1$ = 3, rf_e2$ = 4

  integer, parameter :: scir_tilt_w$ = 1, scir_tilt_sk_w$ = 2
  integer, parameter :: scir_tilt_e$ = 3, scir_tilt_sk_e$ = 4

!---------------------------------------------------------------
! the cesr_struct is used by bmad_to_cesr.
! the element ordering here is similar to what is used by cesrv

  type b_struct
    real(rp) beta_mid             ! beta at the midpoint
    real(rp) beta_ave             ! beta averaged over the element
    real(rp) pos_mid              ! position at midpoint
    real(rp) pos_inj              ! position in injection lattice
  end type

  type cesr_element_struct 
    character*16 name              ! bmad name
    type (b_struct)  x, y          ! beta's and positions
    character*12 db_node_name
    integer ix_db                  ! element index for data base node
    integer ix_ring                ! index to element in ring structure
  end type

  integer, parameter :: n_quad_maxx = 120
  integer, parameter :: n_oct_maxx = 4
  integer, parameter :: n_hbnd_maxx = 6
  integer, parameter :: n_rf_maxx = 4
  integer, parameter :: n_skew_sex_maxx = 20
  integer, parameter :: n_wig_maxx = 2
  integer, parameter :: n_sep_maxx = 6
  integer, parameter :: n_qadd_maxx = 7
  integer, parameter :: n_scir_cam_maxx = 10
  integer, parameter :: n_scir_quad_maxx = 4
  integer, parameter :: n_scir_tilt_maxx = 8

! %ix_cesr is a pointer to cesr_struct for given type

  type cesr_struct
    type (cesr_element_struct) quad_(0:120)
    type (cesr_element_struct) skew_quad_(0:120)
    type (cesr_element_struct) sex_(0:120)
    type (cesr_element_struct) det_(0:120)
    type (cesr_element_struct) skew_sex_(n_skew_sex_maxx)
    type (cesr_element_struct) oct_(n_oct_maxx)
    type (cesr_element_struct) rf_(n_rf_maxx)
    type (cesr_element_struct) wig_(n_wig_maxx)
    type (cesr_element_struct) sep_(n_sep_maxx)
    type (cesr_element_struct) h_steer_(0:120)
    type (cesr_element_struct) v_steer_(0:120)
    type (cesr_element_struct) solenoid       ! solenoid struct
    type (cesr_element_struct) scir_cam_rho_(n_scir_cam_maxx)
    type (cesr_element_struct) scir_tilt_(n_scir_tilt_maxx)
    integer ix_ip_l3                     ! pointer to IP_L3
    integer, pointer :: ix_cesr(:) => null() 
  end type

!---------------------------------------------------------------------------
! DB_STRUCT:                       
! This structrue holds info on the correspondence between CESR data base
! elements and a BMAD ring. Use BMAD_TO_DB to initialize this structure.
! For completeness there are, in addition, arrays for the detectors and
! wigglers, etc.

! for an individual element

  type db_element_struct
    character*16 bmad_name    ! bmad name of element
    real(rp) dvar_dcu             ! calibration factor
    real(rp) var_theory           ! theory var value
    real(rp) var_0                ! extrapolated var value at CU = 0
    integer ix_ring           ! index to element array in ring struct
    integer ix_attrib         ! index to element attribute
    integer ix_cesrv          ! index to cesr_struct arrays
    character*12 db_node_name ! node name ("CSR QUAD CUR")
    integer ix_db             ! element index for data base node (5 for Q05W)
    character*16 db_ele_name  ! element name
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

!-----------------------------------------------------------------------------
! For butns.nnnnn files

  type detector_struct
    real(rp) x_orb, y_orb
    integer amp(4)
    integer type
    logical ok
  endtype

  type butns_struct
    character*40 lattice
    character*20 date
    character*72 comment(5)
    type (detector_struct) det(0:120)
    integer save_set
    integer file_num
    integer turn    ! turn number for injection data
  end type

contains

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine IDENTIFY_DB_NODE (DB_NAME, DB, DP_PTR, OK, TYPE_ERR)
!
! Subroutine to find which array in DB is associated with DB_NAME.
!
! Modules Needed:
!   use bmad
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

subroutine identify_db_node (db_name, db, db_ptr, ok, type_err)

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

  if (type_err) print *, &
        'ERROR IN IDENTIFY_DB_NODE: CANNOT FIND DATABASE NODE: ', db_name
  ok = .false.
  return

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine db_group_to_bmad (ing_name, ing_num, biggrp_set, ring, db, &
!                                                con_, n_con, ok, type_err)
!
! Subroutine to take a data base group element and find the elements
! controlled along with the coefficients.
!
! Modules used:
!   use bmad
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

subroutine db_group_to_bmad (ing_name, ing_num, biggrp_set, ring, db, &
                                              con_, n_con, ok, type_err)


  implicit none

  type (ring_struct) ring
  type (db_struct) db
  type (group_info_struct) grp
  type (control_struct) con_(:)
  type (db_element_struct), pointer :: db_ptr(:)

  integer k, n, nn, n_con, ing_num, biggrp_set

  character*12 ing_name

  logical ok, type_err

!

  call setup_group (ing_name, ing_num, 0, grp, ok, .true.)

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
          print *, 'ERROR IN DB_GROUP_TO_BMAD: CANNOT FIND IN RING ELEMENT'
          print *, '      FOR: ', db_ptr(n)%db_node_name, db_ptr(n)%ix_db
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

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
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

subroutine db_group_to_bmad_group (group_name, group_num, i_biggrp, &
                                       ring, db, ix_ele, ok, type_err)

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
                                                             
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine READ_BUTNS_FILE (butns_num, butns, db, read_ok, type_err)
!
! Subroutine to read in the information in a BUTNS.nnnnn file.
!
! Modules Needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!   butns_num -- Integer: Number in BUTNS.nnnnn file name.
!   type_err  -- Logical: If True then open error message will be printed.
!                         If False then open error message will be supressed.
!
! Output:
!   butns -- Butns_struct: Orbit information.
!       %lattice    -- Character*40: Lattice name.
!       %save_set   -- Integer: Save set number.
!       %date       -- Character*20: Date orbit was taken
!       %file_num   -- Integer: Equal to butns_num.
!       %turn       -- Integer: Turn number for injection orbits. 0 otherwise.
!       %comment(5) -- Character*72: Comment.
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

subroutine read_butns_file (butns_num, butns, db, read_ok, type_err)

  implicit none

  type (db_struct) db
  type (butns_struct) butns

  integer vec(120), det_type(120)
  integer i, ix, j, butns_num, iu, lunget, ios, raw(4, 120)
  integer n_node, n_ele

  real(rp) x_orbit(120), y_orbit(120), rdummy

  character line_in*130, file_in*40

  logical read_ok, type_err, err_flag, is_ok(120)

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
      read (line_in(ix+5:), '(i)', iostat = ios) butns%turn
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
  call butcon (raw, 1, 100, y_orbit, x_orbit)
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

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine BMAD_TO_CESR (RING, CESR)
!
! Subroutine to transfer information from the RING structure returned from
! BMAD_PARSER to a structure for the CESR ring
!
! WARNING! cesr%ix_cesr is allocated in this subroutine.  It should be
! deallocated in the calling function!

!
! Modules Needed:
!   use bmad
!
! Input:
!     RING      -- Ring_struct: Ring to parse.
!     STATUS    -- Common block status structure
!       .TYPE_OUT    -- If .true. then type error messages if all the elements
!                       are not found
!
! Output:
!     CESR        -- Cesr_struct:
!     STATUS    -- Common block status structure
!       .OK         -- Set .false. if failure to find all the elements.
!
! Notes:
!
! Hardbend steerings are put in CESR.H_STEER(101) through CESR.H_STEER(106)
!-

subroutine bmad_to_cesr (ring, cesr)

  implicit none

  type (ring_struct)  ring
  type (cesr_struct)  cesr
  type (ele_struct)  ele

  character cc1*1, cc2*2, cc4*4
  character*16 hsteer_name(0:120), vsteer_name(0:99)

  integer j, i, vnumbr

! load names

  if (associated(cesr%ix_cesr)) deallocate (cesr%ix_cesr)
  allocate (cesr%ix_cesr(ring%n_ele_max))

  cesr%quad_(:)%ix_ring         = 0
  cesr%skew_quad_(:)%ix_ring    = 0
  cesr%skew_sex_(:)%ix_ring     = 0
  cesr%sep_(:)%ix_ring          = 0
  cesr%sex_(:)%ix_ring          = 0
  cesr%det_(:)%ix_ring          = 0
  cesr%oct_(:)%ix_ring          = 0
  cesr%wig_(:)%ix_ring          = 0
  cesr%rf_(:)%ix_ring           = 0
  cesr%h_steer_(:)%ix_ring      = 0
  cesr%v_steer_(:)%ix_ring      = 0
  cesr%scir_cam_rho_(:)%ix_ring = 0
  cesr%scir_tilt_(:)%ix_ring    = 0

  cesr%quad_(:)%name         = 'DUMMY'            ! assume nothing here
  cesr%skew_quad_(:)%name    = 'DUMMY'
  cesr%skew_sex_(:)%name     = 'DUMMY'
  cesr%sep_(:)%name          = 'DUMMY'
  cesr%sex_(:)%name          = 'DUMMY'
  cesr%det_(:)%name          = 'DUMMY'
  cesr%oct_(:)%name          = 'DUMMY'
  cesr%wig_(:)%name          = 'DUMMY'
  cesr%rf_(:)%name           = 'DUMMY'
  cesr%h_steer_(:)%name      = 'DUMMY'
  cesr%v_steer_(:)%name      = 'DUMMY'
  cesr%scir_cam_rho_(:)%name = 'DUMMY'            
  cesr%scir_tilt_(:)%name    = 'DUMMY'            

!

  do i = 0, 49
    write (cc2, '(i2.2)') i
    cesr%quad_(i)%name    = 'Q' // cc2 // 'W'
    cesr%quad_(99-i)%name = 'Q' // cc2 // 'E'
    cesr%skew_quad_(i)%name    = 'SK_Q' // cc2 // 'W'
    cesr%skew_quad_(99-i)%name = 'SK_Q' // cc2 // 'E'
    cesr%sex_(i)%name    = 'SEX_' // cc2 // 'W'
    cesr%sex_(99-i)%name = 'SEX_' // cc2 // 'E'
    cesr%det_(i)%name    = 'DET_' // cc2 // 'W'
    cesr%det_(99-i)%name = 'DET_' // cc2 // 'E'
  enddo

  do i = 0, 99
    write (cc4, '(i4)') i
    hsteer_name(i) = 'CSR HORZ CUR' // cc4
    vsteer_name(i) = 'CSR VERT CUR' // cc4
  enddo

  do i = 1, n_hbnd_maxx
    write (cc4, '(i4)') i
    hsteer_name(i+100) = 'CSR HBND CUR' // cc4
  enddo

  cesr%quad_(q49aw$)%name = 'Q49AW'
  cesr%quad_(q47aw$)%name = 'Q47AW'
  cesr%quad_(q47ae$)%name = 'Q47AE'
  cesr%quad_(q49ae$)%name = 'Q49AE'
  cesr%quad_(q43aw$)%name = 'Q43AW'
  cesr%quad_(q43ae$)%name = 'Q43AE'
  cesr%quad_(q08aw$)%name = 'Q08AW'

  cesr%sep_(h_sep_08w$)%name = 'H_SEP_08W'
  cesr%sep_(h_sep_08e$)%name = 'H_SEP_08E'
  cesr%sep_(h_sep_45w$)%name = 'H_SEP_45W'
  cesr%sep_(h_sep_45e$)%name = 'H_SEP_45E'
  cesr%sep_(v_sep_48w$)%name = 'V_SEP_48W'
  cesr%sep_(v_sep_48e$)%name = 'V_SEP_48E'

  cesr%oct_(1)%name = 'OCT_45W'
  cesr%oct_(2)%name = 'OCT_48W'
  cesr%oct_(3)%name = 'OCT_48E'
  cesr%oct_(4)%name = 'OCT_45E'

  cesr%rf_(rf_w1$)%name = 'RF_W1'
  cesr%rf_(rf_w2$)%name = 'RF_W2'
  cesr%rf_(rf_e1$)%name = 'RF_E1'
  cesr%rf_(rf_e2$)%name = 'RF_E2'


  cesr%skew_sex_(1)%name = 'DUMMY'
  cesr%skew_sex_(2)%name = 'SK_SEX_23W'
  cesr%skew_sex_(3)%name = 'SK_SEX_23E'
  cesr%skew_sex_(4)%name = 'SK_SEX_07E'

  cesr%wig_(wig_w$)%name = 'WIG_W'
  cesr%wig_(wig_e$)%name = 'WIG_E'
 
  cesr%solenoid%name = 'CLEO_SOL'

! phase_iii

  cesr%v_steer_(111)%name = 'SC_V01W'
  cesr%v_steer_(112)%name = 'SC_V02W'
  cesr%v_steer_(113)%name = 'SC_V02E'
  cesr%v_steer_(114)%name = 'SC_V01E'

  cesr%skew_quad_(111)%name = 'SC_SK_Q01W'
  cesr%skew_quad_(112)%name = 'SC_SK_Q02W'
  cesr%skew_quad_(113)%name = 'SC_SK_Q02E'
  cesr%skew_quad_(114)%name = 'SC_SK_Q01E'

  cesr%quad_(111)%name = 'SC_Q01W'
  cesr%quad_(112)%name = 'SC_Q02W'
  cesr%quad_(113)%name = 'SC_Q02E'
  cesr%quad_(114)%name = 'SC_Q01E'

  do i = 1, 5                          
    write (cesr%scir_cam_rho_(i)%name,   '(a, i1, a)') 'SC_CAM_', i, 'W'
    write (cesr%scir_cam_rho_(i+5)%name, '(a, i1, a)') 'SC_CAM_', i, 'E'
  enddo

  cesr%scir_tilt_(scir_tilt_w$)%name    = 'SC_TILT_W'
  cesr%scir_tilt_(scir_tilt_e$)%name    = 'SC_TILT_E'
  cesr%scir_tilt_(scir_tilt_sk_w$)%name = 'SC_TILT_SK_W'
  cesr%scir_tilt_(scir_tilt_sk_e$)%name = 'SC_TILT_SK_E'

  cesr%skew_sex_(11)%name = 'SK_SEX_02E'

!-------------------------------------------------------------
! Load elements from RING to CESR

  ele_loop: do i = 1, ring%n_ele_max

    ele = ring%ele_(i)

! quads and skew quads

    if (ele%key == quadrupole$) then

      if (ele%name(1:1) == 'Q' .or. ele%name(1:4) == 'SC_Q') then
        do j = 0, 120
          if (ele%name == cesr%quad_(j)%name) then
            cesr%ix_cesr(i) = j
            call insert_info (cesr%quad_(j), ele, i)
            cycle ele_loop
          endif
        enddo

      elseif (ele%name(:2) == 'SK' .or. ele%name(1:5) == 'SC_SK') then
        do j = 0, 120
          if (ele%name == cesr%skew_quad_(j)%name) then
            cesr%ix_cesr(i) = j
            call insert_info (cesr%skew_quad_(j), ele, i)
            cycle ele_loop
          endif
        enddo
      endif

    endif

! sex and skew sex

    if (ele%key == sextupole$) then

      if (ele%name(:3) == 'SEX') then
        do j = 0, 120
          if (ele%name == cesr%sex_(j)%name) then
            cesr%ix_cesr(i) = j
            call insert_info (cesr%sex_(j), ele, i)
            cycle ele_loop
          endif
        enddo

      elseif (ele%name(:6) == 'SK_SEX') then
        do j = 1, n_skew_sex_maxx
          if (ele%name == cesr%skew_sex_(j)%name) then
            cesr%ix_cesr(i) = j
            call insert_info (cesr%skew_sex_(j), ele, i)
            cycle ele_loop
          endif
        enddo
      endif

    endif

! octupoles

    if (ele%key == octupole$) then

      do j = 1, n_oct_maxx
        if (ele%name == cesr%oct_(j)%name) then
          cesr%oct_(j)%ix_ring = i
          cesr%ix_cesr(i) = j
          call insert_info (cesr%oct_(j), ele, i)
          cycle ele_loop
        endif
      enddo

    endif

! separators

    if (ele%key == elseparator$) then

      do j = 1, n_sep_maxx
        if (ele%name == cesr%sep_(j)%name) then
          cesr%ix_cesr(i) = j
          call insert_info (cesr%sep_(j), ele, i)
          cycle ele_loop
        endif
      enddo

    endif

! markers... Detectors or IP_L3
! detector markers

    if (ele%key == marker$) then

      if (ele%name(:3) == 'DET') then
        do j = 0, 120
          if (ele%name == cesr%det_(j)%name) then
            cesr%ix_cesr(i) = j
            call insert_info (cesr%det_(j), ele, i)
            cycle ele_loop
          endif
        enddo

      elseif (ele%name(:5) == 'IP_L3') then
        cesr%ix_ip_l3 = i
      endif

    endif


! horz and vert steering overlays
! scir cam and tilts

    if (ele%key == overlay_lord$) then

      if (ele%name(:1) == 'H') then
        do j = 0, 120
          if (ele%type == hsteer_name(j)) then
            call insert_info (cesr%h_steer_(j), ele, i)
            cycle ele_loop
          endif
        enddo

      elseif (ele%name(1:1) == 'V') then
        do j = 0, 99
          if (ele%type == vsteer_name(j)) then
            call insert_info (cesr%v_steer_(j), ele, i)
            cycle ele_loop
          endif
        enddo

      elseif (ele%name(1:4) == 'SC_V') then
        do j = 101, size(cesr%v_steer_)
          if (ele%name == cesr%v_steer_(j)%name) then
            call insert_info (cesr%v_steer_(j), ele, i)
            cycle ele_loop
          endif
        enddo     

      elseif (ele%name(1:6) == 'SC_CAM') then
        do j = 1, size(cesr%scir_cam_rho_)
          if (ele%name == cesr%scir_cam_rho_(j)%name) then
            call insert_info (cesr%scir_cam_rho_(j), ele, i)
            cycle ele_loop
          endif
        enddo

      elseif (ele%name(1:7) == 'SC_TILT') then
        do j = 1, size(cesr%scir_tilt_)
          if (ele%name == cesr%scir_tilt_(j)%name) then
            call insert_info (cesr%scir_tilt_(j), ele, i)
            cycle ele_loop
          endif
        enddo
      endif

    endif

! rf

    if (ele%key == rfcavity$) then
      do j = 1, n_rf_maxx
        if (ele%name == cesr%rf_(j)%name) then
          call insert_info (cesr%rf_(j), ele, i)
          cycle ele_loop
        endif
      enddo
    endif

! wiggler

    if (ele%key == wiggler$) then
      do j = 1, n_wig_maxx
        if (ele%name == cesr%wig_(j)%name) then
          call insert_info (cesr%wig_(j), ele, i)
          cycle ele_loop
        endif
      enddo
    endif

! solenoid

    if (ele%name == cesr%solenoid%name) then
      cesr%solenoid%ix_ring = i
      cycle ele_loop
    endif

  enddo ele_loop
           
!-------------------------------------------------------------------
! check that we have loaded everything...
! do not check Q01 and Q02's

  if (cesr%quad_( 1)%ix_ring == 0) cesr%quad_( 1)%name = 'DUMMY'
  if (cesr%quad_( 2)%ix_ring == 0) cesr%quad_( 2)%name = 'DUMMY'
  if (cesr%quad_(97)%ix_ring == 0) cesr%quad_(97)%name = 'DUMMY'
  if (cesr%quad_(98)%ix_ring == 0) cesr%quad_(98)%name = 'DUMMY'

  call bmad_to_cesr_err_type (cesr%quad_,        'QUADRUPOLE')
  call bmad_to_cesr_err_type (cesr%sep_,         'SEPARATOR')
  call bmad_to_cesr_err_type (cesr%skew_sex_,    'SKEW SEXTUPOLE')
  call bmad_to_cesr_err_type (cesr%oct_,         'OCTUPOLE')
  call bmad_to_cesr_err_type (cesr%wig_,         'WIGGLER')
  call bmad_to_cesr_err_type (cesr%rf_,          'RF CAVITY')
  call bmad_to_cesr_err_type (cesr%scir_cam_rho_,    'SCIR CAM')
  call bmad_to_cesr_err_type (cesr%scir_tilt_,   'SCIR TILT')

!----------------------------------------------------------------
contains

subroutine bmad_to_cesr_err_type (cesr_ele, str)

  type (cesr_element_struct) :: cesr_ele(:)
  integer i
  character*(*) str

!

  do i = lbound(cesr_ele, 1), ubound(cesr_ele, 1)
    if (cesr_ele(i)%ix_ring == 0 .and.  &
                                  cesr_ele(i)%name(:5) /= 'DUMMY') then
      bmad_status%ok = .false.
      if (bmad_status%type_out) then
        print *, 'WARNING FROM BMAD_TO_CESR. ELEMENT: ', cesr_ele(i)%name
        print *, '        NOT LOADED INTO CESR STRUCT: ', str
      endif
    endif
  enddo

end subroutine

end subroutine

!------------------------------------------------------------------------

subroutine insert_info (cesr_ele, ele, i_ele)

  use bmad_struct
  use bmad_interface
  
  implicit none

  type (cesr_element_struct)  cesr_ele
  type (ele_struct)  ele
  integer i_ele, ios

!

  cesr_ele%ix_ring = i_ele
  cesr_ele%db_node_name = ele%type(:12)
  cesr_ele%name = ele%name
  if (ele%type(13:) /= '    ') then
    read (ele%type(13:), *, iostat = ios) cesr_ele%ix_db
    if (ios /= 0) then
      print *, 'ERROR IN INSERT_INFO: READ ERROR FOR NODE INDEX: ', ele%type
    endif
  endif

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine BMAD_TO_DB (RING, DB)
!
! Subroutine to return information on the data base that pertains to
! CESR elements.
!
! Note: This subroutine uses database routines. You must setup communication
! with the database before you use this routine.
!
! Modules needed:
!   use bmad
!   USE BMAD_INTERFACE
!
! Input:
!   RING -- Ring_struct: Ring structure
!
! Output:
!   DB   -- DB_struct: Data base structure. See the module file 
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

subroutine bmad_to_db (ring, db)

  implicit none

  type (ring_struct) ring
  type (db_struct) db
  type (cesr_element_struct) :: cesr_ele(4)

  call bmad_to_db_main (ring, db)

!-------------------------------------------------------------------------
contains

subroutine bmad_to_db_main (ring, db)

  implicit none

  type (ring_struct) ring
  type (cesr_struct) cesr
  type (db_struct) db

  real(rp) h_stren(120), v_stren(120), gev
  real(rp) k_theory(0:120), k_base(0:120), len_quad(0:120), cu_per_k_gev(0:120)
  real(rp) quad_tilt(0:120), dk_gev_dcu(0:120)

  integer n1, i, ix, ios, cu_theory(0:120), nq100

  character*16 type_str

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

  gev = 1e-9 * ring%param%beam_energy
  if (gev == 0) then
    print *, 'ERROR IN BMAD_TO_DB: BEAM ENERGY IS ZERO!'
    call err_exit
  endif

  call getcs (v_stren, h_stren)
  db%csr_horz_cur(:)%dvar_dcu = -1.0e-6 * h_stren(1:98) / gev  
  db%csr_hbnd_cur(:)%dvar_dcu = -1.0e-6 * h_stren(101:106) / gev  
  db%csr_vert_cur(:)%dvar_dcu =  1.0e-6 * v_stren(1:98) / gev  
  db%scir_vertcur(:)%dvar_dcu =  1.0e-6 * v_stren(111:114) / gev  

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

end subroutine

!---------------------------------------------------------------------
! contains
                        
subroutine db_init_it (node, n1, node_name, ix_attrib, node_array, &
                                                      cesr_ele, ix0_cesrv)

  implicit none

  integer ix_attrib, n1, i, vnumbr, ix, iv, ix_(200)

  type (db_element_struct), target :: node(n1:)
  type (db_node_struct) :: node_array(:)
  type (cesr_element_struct), optional :: cesr_ele(:)

  character*(*) node_name

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

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine ring_to_quad_calib (ring, cesr, k_theory, k_base, len_quad,
!                               cu_per_k_gev, quad_rot, dk_gev_dcu, cu_theory)
!
! This subroutine returns the calibration constants for the CESR quads.
!
! Modules Needed:
!   use bmad
!
! Input:
!   ring  -- Ring_struct: Ring with lattice loaded.
!   cesr  -- Cesr_struct: Locations of the quads. 
!              Need previous call to bmad_to_cesr.
!
! Output:
! k_theory(0:120)     -- Real(rp): Theory K of quad_i,
! k_base(0:120)       -- Real(rp): Extrapolated K for zero cu command
! cu_per_k_gev(0:120) -- Real(rp): CU to K*GEV calibration
! len_quad(0:120)     -- Real(rp): Length of quad_i
! quad_rot(0:120)     -- Real(rp): Quad rotation angle in degrees
! dk_gev_dcu(0:120)   -- Real(rp): Derivative of K*GEV vs CU curve.
! cu_theory(0:120)    -- Integer: Scaler needed to get K_THEORY.
!
! Notes:
!
! 0) Without corrections
!         CU_THEORY = (K_THEORY - K_BASE) * GEV * CU_PER_K_GEV
!    In actuality there are corrections due to the calibration of CU with
!    current so the above equation will be off slightly.
!
! 1)      Index     Quad
!             0     REQ W
!            99     REQ E
!           101    QADD 50 W
!           102    QADD 47 W
!           103    QADD 47 E
!           104    QADD 50 E
!
! 2) Nonzero K_BASE is due to using current from the dipoles in a quadrupole.
!
! 3) DK_GEV_DCU = 1 / CU_PER_K_GEV except if there are iron saturation effects:
!
!    DK_GEV_DCU is the derivative when the actual k = K_THEORY.
!    I.e. DK_GEV_DCU is the tangent of the curve.
!
!    CU_PER_K_GEV is the average slope from K_BASE to K_THEORY.
!    I.e. CU_PER_K_GEV is the chord between 2 points on the curve.
!-

subroutine ring_to_quad_calib (ring, cesr, k_theory, k_base,  &
                 len_quad, cu_per_k_gev, quad_rot, dk_gev_dcu, cu_theory)

  implicit none

  record /cesr_struct/ cesr
  record /ring_struct/ ring
  real(rp) gev, k_theory(0:*), k_base(0:*), len_quad(0:*)
  real(rp) cu_per_k_gev(0:*), dk_gev_dcu(0:*), quad_rot(0:*)
  integer cindex, rindex, cu_theory(0:*)

! init  &
  do cindex = 0, 120
    quad_rot(cindex) = 0.0
  enddo

! read lattice file
  gev = 1e-9 * ring%param%beam_energy

  do cindex = 0, 120
    rindex = cesr%quad_(cindex)%ix_ring
    if(rindex/=0) then
      k_theory(cindex) = ring%ele_(rindex)%value(k1$)
      len_quad(cindex) = ring%ele_(rindex)%value(l$)
      quad_rot(cindex) = ring%ele_(rindex)%value(tilt$)*(180./pi)
    endif
  enddo
      
! convert k_theory to scalar computer units given the specified design energy

  call k_to_quad_calib (k_theory, gev, cu_theory, k_base, dk_gev_dcu,  &
                                                                cu_per_k_gev)

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine QUAD_CALIB (LATTICE, K_THEORY, K_BASE, LEN_QUAD,
!                         CU_PER_K_GEV, QUAD_ROT, DK_GEV_DCU, CU_THEORY)
!
! This subroutine returns the calibration constants for the CESR quads.
!
! Modules Needed:
!   use bmad
!
! Input:
!     char*(*)  LATTICE        -- lattice name
!
! Output:
!     real(rp)  K_THEORY(0:120)     -- Theory K of quad_i,
!     real(rp)  K_BASE(0:120)       -- Extrapolated K for zero cu command
!     real(rp)  CU_PER_K_GEV(0:120) -- CU to K*GEV calibration
!     real(rp)  LEN_QUAD(0:120)     -- Length of quad_i
!     real(rp)  QUAD_ROT(0:120)     -- Quad rotation angle in degrees
!     real(rp)  DK_GEV_DCU(0:120)   -- Derivative of K*GEV vs CU curve.
!     integer CU_THEORY(0:120)  -- Scaler needed to get K_THEORY.
!
! Notes:
!
! 0) Without corrections
!         CU_THEORY = (K_THEORY - K_BASE) * GEV * CU_PER_K_GEV
!    In actuality there are corrections due to the calibration of CU with
!    current so the above equation will be off slightly.
!
! 1)      Index     Quad
!             0     REQ W
!            99     REQ E
!           101    QADD 50 W
!           102    QADD 47 W
!           103    QADD 47 E
!           104    QADD 50 E
!
! 2) Nonzero K_BASE is due to using current from the dipoles in a quadrupole.
!
! 3) DK_GEV_DCU = 1 / CU_PER_K_GEV except if there are iron saturation effects:
!
!    DK_GEV_DCU is the derivative when the actual k = K_THEORY.
!    I.e. DK_GEV_DCU is the tangent of the curve.
!
!    CU_PER_K_GEV is the average slope from K_BASE to K_THEORY.
!    I.e. CU_PER_K_GEV is the chord between 2 points on the curve.
!
! 4) Positive QUAD_ROT => CCW rotation from E+ viewpoint
!-

subroutine quad_calib (lattice, k_theory, k_base,  &
                 len_quad, cu_per_k_gev, quad_rot, dk_gev_dcu, cu_theory)

  implicit none

  type (cesr_struct)  cesr
  type (ring_struct)  ring
  character lattice*(*), latfil*50, bmad_lat*40
  real(rp) energy, k_theory(0:*), k_base(0:*), len_quad(0:*)
  real(rp) cu_per_k_gev(0:*), dk_gev_dcu(0:*), quad_rot(0:*)
  integer ix, rindex, cu_theory(0:*)

! read lattice file

  call lattice_to_bmad_file_name (lattice, latfil)
  call bmad_parser(latfil, ring)
  call bmad_to_cesr(ring, CESR)
  energy = 1e-9 * ring%param%beam_energy


  k_theory(0:120) = 0
  len_quad(0:120) = 0
  quad_rot(0:120) = 0

  do ix = 0, 120
    rindex = cesr%quad_(ix)%ix_ring
    if(rindex==0) cycle
    k_theory(ix) = ring%ele_(rindex)%value(k1$)
    len_quad(ix) = ring%ele_(rindex)%value(l$)
    quad_rot(ix) = ring%ele_(rindex)%value(tilt$)*(180./pi)
  enddo

! convert k_theory to scalar computer units given the specified design energy

  call k_to_quad_calib (k_theory, energy, cu_theory, k_base, dk_gev_dcu,  &
                                                                cu_per_k_gev)

end subroutine

end module
