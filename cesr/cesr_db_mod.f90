module cesr_db_mod

  use cesr_basic_mod
  use cesr_quad_mod

  private db_init_it, get_ele_theory, non_db_set

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine bmad_to_db (lat, db, calib_date, use_mpm)
!
! Subroutine to return information on the data base that pertains to
! CESR elements.
!
! Note: This subroutine uses database routines. You must setup communication
! with the database before you use this routine.
!
! Modules needed:
!   use cesr_db_mod
!
! Input:
!   lat        -- lat_struct: Lat structure
!   calib_date -- Character(10), optional: Date to use when looking up the 
!                   calibration constants. 
!                   calib_date is in the form: "YYYY-MM-DD".
!                   the present calibrations will be used if calib_date is 
!                   not present or is "NOW"
!   use_mpm    -- Logical, optional: If present and False then this routine 
!                   will assume that there is no mpm connection. 
!                   In this case, the element name and software limits will
!                   not be loaded.
!
! Output:
!   db   -- DB_struct: Data base structure. See the module file 
!                      cesr_db_mod.F90 for more details. 
!
! Notes:
!
! %CU_NOW is *NOT* filled in by this subroutine.
!
! To check if an element is in the data base check that %DVAR_DCU is not zero.
! To check if an element has a corresponding lat element check that %IX_LAT
!   is not zero.
!
! For completeness the REQ's are included in DB%CSR_QUAD_CUR as 
! elements 0 and 99. 
! These elements have %IX_LAT /= 0 but %DVAR_DCU = 0. 
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

subroutine bmad_to_db (lat, db, calib_date, use_mpm)

  implicit none

  type (lat_struct) lat
  type (db_struct) db
  type (cesr_struct) cesr

  real(rp) gev, cu_per_k_gev(0:120)
  real(rp) k_theory(0:120), k_base(0:120), len_quad(0:120)
  real(rp) quad_tilt(0:120), dk_gev_dcu(0:120)

  real h_stren(120), v_stren(120), gev_sp
  real hsp_dvar_dcu(n_sep_maxx), vsp_dvar_dcu(n_sep_maxx)
  real sqewquad_dvar_dcu(98), scir_skqucur_dvar_dcu(n_scir_quad_maxx)
  real sext_dvar_dcu(98), scsolcur_dvar_dcu(5*n_sc_sol_maxx)
  real sqewsext_dvar_dcu(n_csr_sqewsext_maxx), ir_sksxcur_dvar_dcu(1)

  integer i, ix, ios, cu_theory(0:120), nq100

  character(*), optional :: calib_date
  character(16) type_str
  character(10) this_date

  logical, optional :: use_mpm

! general init

  nq100 = 100+n_qadd_maxx

  call bmad_to_cesr (lat, cesr)

  call db_init_it (db%csr_quad_cur, 0, 'INIT_IT', 0, db%node, use_mpm = use_mpm)  

  call db_init_it (db%csr_quad_cur, lbound(db%csr_quad_cur, 1), &
          'CSR QUAD CUR',    k1$,    db%node, cesr%quad(1:98), 1, use_mpm)

  call db_init_it (db%nir_shuntcur, lbound(db%nir_shuntcur, 1), &
          'NIR SHUNTCUR',    k1$,    db%node, cesr%nir_shuntcur, 1, use_mpm)

  call db_init_it (db%csr_qadd_cur, lbound(db%csr_qadd_cur, 1), &
          'CSR QADD CUR',    k1$,    db%node, cesr%quad(101:nq100), 101, use_mpm)

  call db_init_it (db%scir_quadcur, lbound(db%scir_quadcur, 1), &
          'SCIR QUADCUR',    k1$,    db%node, cesr%quad(111:114), 111, use_mpm)

  call db_init_it (db%csr_horz_cur, lbound(db%csr_horz_cur, 1), &
          'CSR HORZ CUR', hkick$, db%node, cesr%h_steer(1:98), 1, use_mpm)

  call db_init_it (db%csr_hbnd_cur, lbound(db%csr_hbnd_cur, 1), &
          'CSR HBND CUR', hkick$, db%node, cesr%h_steer(101:106), 101, use_mpm)

  call db_init_it (db%csr_vert_cur, lbound(db%csr_vert_cur, 1), &
          'CSR VERT CUR', vkick$, db%node, cesr%v_steer(1:98), 1, use_mpm)

  call db_init_it (db%csr_hsp_volt, lbound(db%csr_hsp_volt, 1), &
          'CSR HSP VOLT', hkick$, db%node, use_mpm = use_mpm)

  call db_init_it (db%csr_vsp_volt, lbound(db%csr_vsp_volt, 1), &
          'CSR VSP VOLT', vkick$, db%node, use_mpm = use_mpm)

  call db_init_it (db%csr_sext_cur, lbound(db%csr_sext_cur, 1), &
          'CSR SEXT CUR',    k2$,    db%node, cesr%sex(1:98), 1, use_mpm)

  call db_init_it (db%csr_octu_cur, lbound(db%csr_octu_cur, 1), &
          'CSR OCTU CUR',    k3$,    db%node, cesr%oct, 1, use_mpm)

  call db_init_it (db%csr_sqewquad, lbound(db%csr_sqewquad, 1), &
          'CSR SQEWQUAD',    k1$,    db%node, cesr%skew_quad(1:98), 1, use_mpm)
                           
  call db_init_it (db%csr_scsolcur, lbound(db%csr_scsolcur, 1), &
          'SCSOL CONTRL',    ks$,    db%node, cesr%scsol_cur(1:10), 1, use_mpm)
                           
  call db_init_it (db%scir_skqucur, lbound(db%scir_skqucur, 1), &
          'SCIR SKQUCUR',    k1$,    db%node, cesr%skew_quad(111:114), 111, use_mpm)
                           
  call db_init_it (db%csr_sqewsext, lbound(db%csr_sqewsext, 1), &
          'CSR SQEWSEXT',    k2$,    db%node, &
           cesr%skew_sex(1:n_csr_sqewsext_maxx), 1, use_mpm)

  call db_init_it (db%scir_vertcur, lbound(db%scir_vertcur, 1), &
          'SCIR VERTCUR', vkick$, db%node, cesr%v_steer(111:114), 111, use_mpm)

  call db_init_it (db%ir_sksxcur, lbound(db%ir_sksxcur, 1), &
          'IR SKSXCUR  ',   k2$,    db%node, cesr%skew_sex(11:11), 11, use_mpm)

!-----------------------------------------------------------------------------
! find the separators

  do i = 1, lat%n_ele_max

    type_str = lat%ele(i)%type

    if (type_str(:12) == db%csr_vsp_volt(1)%db_node_name) then 
      if (type_str(13:) /= '    ') then
        read (type_str(13:), *, iostat = ios) ix
        if (ios .ne. 0) then
          print *, 'ERROR IN BAMD_TO_DB: READ ERROR FOR NODE INDEX: ', type_str
        endif
        db%csr_vsp_volt(ix)%bmad_name = lat%ele(i)%name
        db%csr_vsp_volt(ix)%ix_lat = i
      endif
    endif

    if (type_str(:12) == db%csr_hsp_volt(1)%db_node_name) then 
      if (type_str(13:) /= '    ') then
        read (type_str(13:), *, iostat = ios) ix
        if (ios .ne. 0) then
          print *, 'ERROR IN BAMD_TO_DB: READ ERROR FOR NODE INDEX: ', type_str
        endif
        db%csr_hsp_volt(ix)%bmad_name = lat%ele(i)%name
        db%csr_hsp_volt(ix)%ix_lat = i
      endif
    endif

  enddo

!--------------------------------------------------------------------
! get calibrations...
! steerings
!!  call getcs (v_stren, h_stren)

  this_date = 'NOW'
  if (present(calib_date)) this_date = calib_date

  gev = 1e-9 * lat%ele(0)%value(E_TOT$)
  gev_sp = gev
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

  call get_sep_strength (hsp_dvar_dcu, vsp_dvar_dcu)
  db%csr_hsp_volt(:)%dvar_dcu = 1.0e-6 * hsp_dvar_dcu / gev
  db%csr_vsp_volt(:)%dvar_dcu = 1.0e-6 * vsp_dvar_dcu / gev

  call get_ele_theory (lat, db%csr_hsp_volt)
  call get_ele_theory (lat, db%csr_vsp_volt)

! sextupole calibration

  call sext_calib (gev_sp, sext_dvar_dcu)
  db%csr_sext_cur%dvar_dcu = sext_dvar_dcu
  call get_ele_theory (lat, db%csr_sext_cur)

! octupole calibration

  db%csr_octu_cur(:)%dvar_dcu = 0.213 / gev
  call get_ele_theory (lat, db%csr_octu_cur)

! skew quad calibration

  call skew_quad_calib (sqewquad_dvar_dcu, scir_skqucur_dvar_dcu)
  db%csr_sqewquad(:)%dvar_dcu = sqewquad_dvar_dcu / gev
  db%scir_skqucur(:)%dvar_dcu = scir_skqucur_dvar_dcu / gev

  call get_ele_theory (lat, db%csr_sqewquad) 
  call get_ele_theory (lat, db%scir_skqucur) 

! skew sextupole calibration

  call skew_sex_calib (sqewsext_dvar_dcu, ir_sksxcur_dvar_dcu)
  db%csr_sqewsext(:)%dvar_dcu = sqewsext_dvar_dcu / gev
  db%ir_sksxcur(:)%dvar_dcu   = ir_sksxcur_dvar_dcu / gev

  call get_ele_theory (lat, db%csr_sqewsext) 
  call get_ele_theory (lat, db%ir_sksxcur) 

! quad calibrations 

  call lat_to_quad_calib (lat, cesr, k_theory, k_base, len_quad, &
                              cu_per_k_gev, quad_tilt, dk_gev_dcu, cu_theory)

  call non_db_set (db%quad, cesr%quad, k1$, 0)
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

  call nir_shuntcur_calib (db%nir_shuntcur%dvar_dcu)
  db%nir_shuntcur%dvar_dcu = &
                    db%nir_shuntcur%dvar_dcu * db%quad(48:51)%dvar_dcu
  db%nir_shuntcur%ix_lat = db%quad(48:51)%ix_lat

! sc solenoid calibration

  call sc_sol_calib (scsolcur_dvar_dcu)
  db%csr_scsolcur(:)%dvar_dcu = scsolcur_dvar_dcu / gev

  call get_ele_theory (lat, db%csr_scsolcur) 


! Things not part of the CESR DB

  call non_db_set (db%detector, cesr%det, 0, 0)
  call non_db_set (db%wiggler, cesr%wig, b_max$, 1)
  call non_db_set (db%scir_cam_rho, cesr%scir_cam_rho, rho$, 1)
  call non_db_set (db%scir_tilt, cesr%scir_tilt, tilt$, 1)


  deallocate( cesr%ix_cesr )

end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
                        
subroutine db_init_it (node, n1, node_name, ix_attrib, node_array, &
                                   cesr_ele, ix0_cesrv, use_mpm)

  implicit none

  integer ix_attrib, n1, i, vnumbr, iv, ix(200)

  type (db_element_struct), target :: node(n1:)
  type (db_node_struct) :: node_array(:)
  type (cesr_element_struct), optional :: cesr_ele(:)

  character(*) node_name

  integer, optional :: ix0_cesrv
  integer, save :: i_array

  logical, optional :: use_mpm

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
    node(i)%ix_lat = 0
    node(i)%ix_attrib = ix_attrib
    node(i)%db_node_name = node_name       
    node(i)%ix_db = i
    node(i)%db_ele_name = ' '           
    node(i)%dvar_dcu = 0
    node(i)%var_theory = 0
    node(i)%var_0 = 0
    node(i)%cu_high_lim = 0
    node(i)%cu_low_lim = 0
    node(i)%ix = 0
  enddo

  if (logic_option (.true., use_mpm)) then
    iv = vnumbr(node_name)
    if (size(node) < iv) then
      print *, 'WARNING IN BMAD_TO_DB: ARRAY IS TOO SMALL TO HOLD DATA BASE'
      print *, '        NODE ELEMENTS: ', node_name
      call err_exit
      iv = size(node)
    endif

    if (iv .gt. 0) then
      call vmglim (node_name, &
                     1, iv, node(1:iv)%cu_low_lim, node(1:iv)%cu_high_lim)
      call vstget (node_name, 1, iv, node(1:iv)%db_ele_name, ix)
    endif
  endif

  if (present(cesr_ele)) then
    if (ubound(node, 1)-n1 .ne. ubound(cesr_ele, 1)-1) then
      print *, 'ERROR IN DB_INIT_IT: ARRAYS DO NOT MATCH FOR: ', node_name
      call err_exit
    endif
    node(:)%ix_lat = cesr_ele(:)%ix_lat
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
!---------------------------------------------------------------------
!---------------------------------------------------------------------
                            
subroutine get_ele_theory (lat, db_ele)
                                                               
  implicit none

  type (db_element_struct) :: db_ele(:)
  type (lat_struct) :: lat

  integer i, ie, iv

!

  do i = 1, ubound(db_ele, 1)
    ie = db_ele(i)%ix_lat
    if (ie /= 0) then
      iv = db_ele(i)%ix_attrib
      db_ele(i)%var_theory = lat%ele(ie)%value(iv)
    endif
  enddo

end subroutine get_ele_theory

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

subroutine non_db_set (db_ele, cesr_ele, ix_attrib, ix0_cesrv)

  implicit none

  type (db_element_struct) :: db_ele(:)
  type (cesr_element_struct) :: cesr_ele(:)

  integer i, ix0_cesrv, ix_attrib

  db_ele(:)%bmad_name = cesr_ele(:)%name
  db_ele(:)%ix_lat = cesr_ele(:)%ix_lat
  db_ele(:)%ix_db = 0
  db_ele(:)%ix_attrib = ix_attrib
  db_ele(:)%db_node_name = '----------'
  db_ele(:)%db_ele_name  = '----------'
  db_ele(:)%dvar_dcu = 0
  db_ele(:)%var_theory = 0
  db_ele(:)%var_0 = 0
  db_ele(:)%cu_high_lim = 0
  db_ele(:)%cu_low_lim = 0
  db_ele(:)%ix = 0
  forall (i = 1:size(db_ele, 1)) db_ele(i)%ix_cesrv = i + ix0_cesrv - 1

end subroutine non_db_set

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine read_butns_file (butns_num, nonlinear_calc, butns, db, 
!                                                         read_ok, type_err)
!
! Subroutine to read the raw information from a BUTNS.nnnnn file and convert
! it to orbit x-y positions.
!
! Modules Needed:
!   use cesr_mod
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
!       %det(0:99)%ok     -- Logical: Was there a valid orbit reading?
!       %det(0:99)%amp(4) -- Integer: raw button numbers.
!       %det(0:99)%x_orb  -- Real(rp): Horizontal orbit in meters.
!       %det(0:99)%y_orb  -- Real(rp): Horizontal orbit in meters.
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

  character(130) line_in, file_in

  logical nonlinear_calc, read_ok, type_err, err_flag, is_ok(120)

! init

  butns%comment = ' '
  call db_struct_init (db)

! compute filename

  read_ok = .false.

  call form_file_name_with_number ('ORBIT', butns_num, file_in, err_flag)
  if (err_flag) return

  butns%file_num = butns_num

! read header line in the data file (to get lattice name)
! and sterlat strengths

  iu = lunget()
  open(unit = iu, file = file_in, status = 'old', action = 'READ', iostat = ios)
  if (ios /= 0) then
    if (type_err) print *, &
          'ERROR IN READ_BUTNS_FILE: ERROR OPENING: ', trim(file_in)
    return
  endif

  read (iu, '(a)') line_in
  if (butns_num >= 100000) then
     butns%lattice = line_in(62:)
  else
     butns%lattice = line_in(61:)
  endif
  call string_trim (butns%lattice, butns%lattice, ix)

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
