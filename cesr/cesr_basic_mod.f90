module cesr_basic_mod

use bmad_struct
use bmad_interface

!---------------------------------------------------------------------------
! DB_STRUCT:                       
! This structrue holds info on the correspondence between CESR data base
! elements and a BMAD lat. Use BMAD_TO_DB to initialize this structure.
! For completeness there are, in addition, arrays for the detectors and
! wigglers, etc.

integer, parameter :: n_quad_maxx = 120
integer, parameter :: n_det_maxx = 120
integer, parameter :: n_sex_maxx = 120
integer, parameter :: n_steer_maxx = 120
integer, parameter :: n_oct_maxx = 4
integer, parameter :: n_hbnd_maxx = 6
integer, parameter :: n_rf_maxx = 4
integer, parameter :: n_skew_sex_maxx = 20
integer, parameter :: n_wig_maxx = 2
integer, parameter :: n_sep_maxx = 6
integer, parameter :: n_qadd_maxx = 20
integer, parameter :: n_scir_cam_maxx = 10
integer, parameter :: n_scir_quad_maxx = 4
integer, parameter :: n_scir_tilt_maxx = 8
integer, parameter :: n_nir_shuntcur_maxx = 4
integer, parameter :: n_sc_sol_maxx = 2

! for an individual element

type db_element_struct
  character(16) bmad_name     ! bmad name of element
  real(rp) dvar_dcu           ! calibration factor
  real(rp) var_theory         ! theory var value
  real(rp) var_0              ! extrapolated var value at CU = 0
  integer ix_lat              ! index to element array in lat struct
  integer ix_attrib           ! index to element attribute
  integer ix_cesrv            ! index to cesr_struct arrays
  character(12) db_node_name  ! node name ("CSR QUAD CUR")
  integer ix_db               ! element index for data base node (5 for Q05W)
  character(16) db_ele_name   ! element name
  integer cu_high_lim         ! high limit
  integer cu_low_lim          ! low limit
  integer cu_now              ! current CU
  integer ix                  ! Integer for general use.
  logical valid_cu_now        ! Valid cu_now value.
  integer cu_meas             ! Value when measured data is taken.
  integer cu_ref              ! Value when reference data is taken
  integer cu_target           ! Target to load into the data base.
  integer cu_design           ! Design value.
end type              

! for pointing to a db_struct array

type db_node_struct
  type (db_element_struct), pointer :: ptr(:) => null()
end type
                                                     
! db_struct def

integer, parameter :: n_csr_sqewsext_maxx = 8

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
  type (db_element_struct) :: csr_scsolcur(5*n_sc_sol_maxx)
  type (db_element_struct) :: csr_sqewsext(n_csr_sqewsext_maxx)
  type (db_element_struct) :: scir_quadcur(n_scir_quad_maxx)
  type (db_element_struct) :: scir_skqucur(n_scir_quad_maxx)
  type (db_element_struct) :: scir_vertcur(n_scir_quad_maxx)
  type (db_element_struct) :: nir_shuntcur(n_nir_shuntcur_maxx)
  type (db_element_struct) :: ir_sksxcur(1)
  type (db_node_struct) :: node(17)  ! does not include stuff below
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
  character(16) name              ! bmad name
  type (b_struct)  x, y          ! beta's and positions
  character(12) db_node_name
  integer ix_db                  ! element index for data base node
  integer ix_lat                 ! index to element in lat structure
end type

! %ix_cesr is a pointer to cesr_struct for given type

type cesr_struct
  type (cesr_element_struct) quad(0:n_quad_maxx)
  type (cesr_element_struct) skew_quad(0:n_quad_maxx)
  type (cesr_element_struct) sex(0:n_sex_maxx)
  type (cesr_element_struct) det(0:n_det_maxx)
  type (cesr_element_struct) skew_sex(n_skew_sex_maxx)
  type (cesr_element_struct) oct(n_oct_maxx)
  type (cesr_element_struct) rf(n_rf_maxx)
  type (cesr_element_struct) wig(n_wig_maxx)
  type (cesr_element_struct) sep(n_sep_maxx)
  type (cesr_element_struct) h_steer(0:n_steer_maxx)
  type (cesr_element_struct) v_steer(0:n_steer_maxx)
  type (cesr_element_struct) solenoid       ! solenoid struct
  type (cesr_element_struct) scir_cam_rho(n_scir_cam_maxx)
  type (cesr_element_struct) scir_tilt(n_scir_tilt_maxx)
  type (cesr_element_struct) nir_shuntcur(n_nir_shuntcur_maxx)
  type (cesr_element_struct) scsol_cur(5*n_sc_sol_maxx)

  integer ix_ip_l3                     ! pointer to IP_L3
  integer, pointer :: ix_cesr(:) => null() 
end type

!-----------------------------------------------------------------------------
! For butns.nnnnn files

type detector_struct
  real(rp) x_orb, y_orb
  integer amp(4)
  integer type
  logical ok
end type

type butns_struct
  character(40) lattice
  character(20) date
  character(72) comment(5)
  type (detector_struct) det(0:120)
  integer save_set
  integer turn    ! turn number for injection data
  real(rp) e_cur, p_cur  ! Current in mAmps
end type

!-----------------------------------------------------------------------------
! For cesrv

type a_cesr_freq_struct 
  real(rp) tune
  real(rp) shake
  logical reflection
end type

type cesr_freq_struct 
  type (a_cesr_freq_struct) x, y
  real(rp) rev
end type

type cesr_phase_params_struct 
  integer species                            
  real(rp) current
  character(60) comment
  integer unit, single_unit
  character(60) file_name, raw_file_name
  character(40) lattice
  integer save_set
  logical :: debug
end type        

type cesr_but_struct 
  real(rp) phase
  real(rp) amp
  logical ok
end type

type cesr_det_xy_struct 
  real(rp) amp
  real(rp) phase
end type

type shaking_mode_struct
  type (cesr_det_xy_struct) x, y
  logical ok
end type

type shaking_modes_struct
  type (shaking_mode_struct) a_mode, b_mode
end type

type cesr_det_plane_struct 
  type (cesr_but_struct) but(4)
  type (cesr_det_xy_struct) x, y, z
  real(rp) phase_meas
  real(rp) rms_phase_meas
  real(rp) phase_design
  real(rp) beta_design
  real(rp) cbar11                
  real(rp) cbar12
  real(rp) cbar22
  real(rp) x_amp
  real(rp) x_phase
  real(rp) y_amp
  real(rp) y_phase
  integer n_buts
  logical ok     
  logical shake
  integer system_id   ! 0 = old_system, 1 = new_system
end type

type phase_cbar_data_struct
  real(rp) x_phase, x_cbar22, x_cbar12
  real(rp) y_phase, y_cbar12, y_cbar11
  logical ok_x, ok_y
end type

type cesr_det_dc_position_struct
  real(rp) x, y
  real(rp) signal(4)
end type

type raw_det_struct
  real(rp) phase(4)
  real(rp) amp(4)
  integer system_id
end type

type cesr_xy_data_struct
  real(rp) x, y
  logical good
end type

! New

type cesr_data_params_struct
  character(20) :: data_date = '', data_type = ''
  character(40) :: var_ele_name = '', var_attrib_name = ''
  character(100) :: comment = '', lattice = '', file_name = '', lattice_file_name = ''
  character(40) :: route_name = ''
  integer :: csr_set = 0
  integer :: species = 0
  integer :: ix_data_set = 0       ! Index of the data set. EG: butns.nnnnnn
  real(rp) :: horiz_beta_freq  = 0, vert_beta_freq  = 0  ! Hz
  real(rp) :: dvar = 0
  real(rp) :: chisq = 0
  real(rp) :: ac_z_amp_fit, ac_z_phase_fit
  logical :: horiz_reflection_shake, vert_reflection_shake
end type

type cesr_data1_struct
  real(rp) value
  logical good
end type

type raw_det2_struct
  type (raw_det_struct) a, b
end type

type cesr_all_data_struct
  type (cesr_data1_struct) orbit_x(0:120), phase_a(0:120), eta_x(0:120)
  type (cesr_data1_struct) orbit_y(0:120), phase_b(0:120), eta_y(0:120)
  type (cesr_data1_struct) cbar11_b(0:120), cbar12_a(0:120), cbar12_b(0:120), cbar22_a(0:120) 
  type (cesr_data1_struct) ac_eta_x(0:120), ac_etap_x(0:120), ac_eta_y(0:120), ac_etap_y(0:120)
  type (detector_struct) raw_orbit(0:120)
  type (raw_det2_struct) raw_phase(0:120)
  type (shaking_modes_struct) shake(0:120)
  type (db_struct) db
  type (cesr_data_params_struct) param
end type

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine bmad_to_cesr (lat, cesr)
!
! Subroutine to transfer information from the LAT structure returned from
! BMAD_PARSER to a structure for the CESR lat
!
! WARNING! cesr%ix_cesr is allocated in this subroutine.  It should be
! deallocated in the calling function!
!
! Modules Needed:
!   use cesr_basic_mod
!
! Input:
!   lat      -- lat_struct: Lat to parse.
!   status    -- Common block status structure
!     %type_out    -- If .true. then type error messages if all the elements
!                       are not found
!
! Output:
!   cesr      -- Cesr_struct:
!   status    -- Common block status structure
!     %ok         -- Set .false. if failure to find all the elements.
!
! Notes:
!
! Hardbend steerings are put in CESR.H_STEER(101) through CESR.H_STEER(106)
!-

subroutine bmad_to_cesr (lat, cesr)

  implicit none

  type (lat_struct)  lat
  type (cesr_struct)  cesr
  type (ele_struct)  ele

  character(2) cc2
  character(4) cc4
  character(16) hsteer_name(0:120), vsteer_name(0:99)

  integer i, j, ix

! load names

  if (associated(cesr%ix_cesr)) deallocate (cesr%ix_cesr)
  allocate (cesr%ix_cesr(lat%n_ele_max))

  cesr%quad(:)%ix_lat         = 0
  cesr%skew_quad(:)%ix_lat    = 0
  cesr%scsol_cur(:)%ix_lat    = 0
  cesr%skew_sex(:)%ix_lat     = 0
  cesr%sep(:)%ix_lat          = 0
  cesr%sex(:)%ix_lat          = 0
  cesr%det(:)%ix_lat          = 0
  cesr%oct(:)%ix_lat          = 0
  cesr%wig(:)%ix_lat          = 0
  cesr%rf(:)%ix_lat           = 0
  cesr%h_steer(:)%ix_lat      = 0
  cesr%v_steer(:)%ix_lat      = 0
  cesr%scir_cam_rho(:)%ix_lat = 0
  cesr%scir_tilt(:)%ix_lat    = 0

  cesr%quad(:)%name         = 'DUMMY'            ! assume nothing here
  cesr%skew_quad(:)%name    = 'DUMMY'
  cesr%skew_sex(:)%name     = 'DUMMY'
  cesr%scsol_cur(:)%name    = 'DUMMY'
  cesr%sep(:)%name          = 'DUMMY'
  cesr%sex(:)%name          = 'DUMMY'
  cesr%det(:)%name          = 'DUMMY'
  cesr%oct(:)%name          = 'DUMMY'
  cesr%wig(:)%name          = 'DUMMY'
  cesr%rf(:)%name           = 'DUMMY'
  cesr%h_steer(:)%name      = 'DUMMY'
  cesr%v_steer(:)%name      = 'DUMMY'
  cesr%scir_cam_rho(:)%name = 'DUMMY'            
  cesr%scir_tilt(:)%name    = 'DUMMY'            

!

  do i = 0, 49
    write (cc2, '(i2.2)') i
    cesr%quad(i)%name    = 'Q' // cc2 // 'W'
    cesr%quad(99-i)%name = 'Q' // cc2 // 'E'
    cesr%skew_quad(i)%name    = 'SK_Q' // cc2 // 'W'
    cesr%skew_quad(99-i)%name = 'SK_Q' // cc2 // 'E'
    cesr%sex(i)%name    = 'SEX_' // cc2 // 'W'
    cesr%sex(99-i)%name = 'SEX_' // cc2 // 'E'
    cesr%det(i)%name    = 'DET_' // cc2 // 'W'
    cesr%det(99-i)%name = 'DET_' // cc2 // 'E'
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

  cesr%quad(q49aw$)%name = 'Q49AW'
  cesr%quad(q47aw$)%name = 'Q47AW'
  cesr%quad(q47ae$)%name = 'Q47AE'
  cesr%quad(q49ae$)%name = 'Q49AE'
  cesr%quad(q43aw$)%name = 'Q43AW'
  cesr%quad(q43ae$)%name = 'Q43AE'
  cesr%quad(q08aw$)%name = 'Q08AW'

  cesr%sep(h_sep_08w$)%name = 'H_SEP_08W'
  cesr%sep(h_sep_08e$)%name = 'H_SEP_08E'
  cesr%sep(h_sep_45w$)%name = 'H_SEP_45W'
  cesr%sep(h_sep_45e$)%name = 'H_SEP_45E'
  cesr%sep(v_sep_48w$)%name = 'V_SEP_48W'
  cesr%sep(v_sep_48e$)%name = 'V_SEP_48E'

  cesr%oct(1)%name = 'OCT_45W'
  cesr%oct(2)%name = 'OCT_48W'
  cesr%oct(3)%name = 'OCT_48E'
  cesr%oct(4)%name = 'OCT_45E'

  cesr%rf(rf_w1$)%name = 'RF_W1'
  cesr%rf(rf_w2$)%name = 'RF_W2'
  cesr%rf(rf_e1$)%name = 'RF_E1'
  cesr%rf(rf_e2$)%name = 'RF_E2'

  cesr%skew_sex(1)%name = 'SK_SEX_07W'
  cesr%skew_sex(2)%name = 'SK_SEX_23W'
  cesr%skew_sex(3)%name = 'SK_SEX_23E'
  cesr%skew_sex(4)%name = 'SK_SEX_07E'
  cesr%skew_sex(5)%name = 'SK_SEX_29W'
  cesr%skew_sex(6)%name = 'SK_SEX_29E'
  cesr%skew_sex(7)%name = 'SK_SEX_12W'
  cesr%skew_sex(8)%name = 'SK_SEX_12E'
  cesr%skew_sex(11)%name = 'SK_SEX_02E'

  cesr%wig(wig_w$)%name = 'WIG_W'
  cesr%wig(wig_e$)%name = 'WIG_E'

  cesr%solenoid%name = 'CLEO_SOL'

! phase_iii

  cesr%v_steer(111)%name = 'SC_V01W'
  cesr%v_steer(112)%name = 'SC_V02W'
  cesr%v_steer(113)%name = 'SC_V02E'
  cesr%v_steer(114)%name = 'SC_V01E'

  cesr%skew_quad(111)%name = 'SC_SK_Q01W'
  cesr%skew_quad(112)%name = 'SC_SK_Q02W'
  cesr%skew_quad(113)%name = 'SC_SK_Q02E'
  cesr%skew_quad(114)%name = 'SC_SK_Q01E'

  cesr%quad(111)%name = 'SC_Q01W'
  cesr%quad(112)%name = 'SC_Q02W'
  cesr%quad(113)%name = 'SC_Q02E'
  cesr%quad(114)%name = 'SC_Q01E'

  do i = 1, 5                          
    write (cesr%scir_cam_rho(i)%name,   '(a, i1, a)') 'SC_CAM_', i, 'W'
    write (cesr%scir_cam_rho(i+5)%name, '(a, i1, a)') 'SC_CAM_', i, 'E'
  enddo

  cesr%scir_tilt(scir_tilt_w$)%name    = 'SC_TILT_W'
  cesr%scir_tilt(scir_tilt_e$)%name    = 'SC_TILT_E'
  cesr%scir_tilt(scir_tilt_sk_w$)%name = 'SC_TILT_SK_W'
  cesr%scir_tilt(scir_tilt_sk_e$)%name = 'SC_TILT_SK_E'

  cesr%nir_shuntcur(1)%name = 'NIR_SHUNTCUR___1'
  cesr%nir_shuntcur(2)%name = 'NIR_SHUNTCUR___2'
  cesr%nir_shuntcur(3)%name = 'NIR_SHUNTCUR___3'
  cesr%nir_shuntcur(4)%name = 'NIR_SHUNTCUR___4'

  cesr%scsol_cur(1)%name = 'SCS03W'
  cesr%scsol_cur(6)%name = 'SCS03E'


!-------------------------------------------------------------
! Load elements from LAT to CESR

  ele_loop: do i = 1, lat%n_ele_max

    ele = lat%ele(i)

! quads and skew quads

    if (ele%key == quadrupole$) then

      if (ele%name(1:1) == 'Q' .or. ele%name(1:4) == 'SC_Q') then
        do j = 0, 120
          if (ele%name == cesr%quad(j)%name) then
            cesr%ix_cesr(i) = j
            call insert_info (cesr%quad(j), ele, i)
            cycle ele_loop
          endif
        enddo

      elseif (ele%name(:2) == 'SK' .or. ele%name(1:5) == 'SC_SK') then
        do j = 0, 120
          if (ele%name == cesr%skew_quad(j)%name) then
            cesr%ix_cesr(i) = j
            call insert_info (cesr%skew_quad(j), ele, i)
            cycle ele_loop
          endif
        enddo
      endif

    endif

! sex and skew sex

    if (ele%key == sextupole$) then

      if (ele%name(:3) == 'SEX') then
        do j = 0, 120
          if (ele%name == cesr%sex(j)%name) then
            cesr%ix_cesr(i) = j
            call insert_info (cesr%sex(j), ele, i)
            cycle ele_loop
          endif
        enddo

      elseif (ele%name(:6) == 'SK_SEX') then
        do j = 1, n_skew_sex_maxx
          if (ele%name == cesr%skew_sex(j)%name) then
            cesr%ix_cesr(i) = j
            call insert_info (cesr%skew_sex(j), ele, i)
            cycle ele_loop
          endif
        enddo
      endif

    endif

! octupoles

    if (ele%key == octupole$) then

      do j = 1, n_oct_maxx
        if (ele%name == cesr%oct(j)%name) then
          cesr%oct(j)%ix_lat = i
          cesr%ix_cesr(i) = j
          call insert_info (cesr%oct(j), ele, i)
          cycle ele_loop
        endif
      enddo

    endif

! separators

    if (ele%key == elseparator$) then

      do j = 1, n_sep_maxx
        if (ele%name == cesr%sep(j)%name) then
          cesr%ix_cesr(i) = j
          call insert_info (cesr%sep(j), ele, i)
          cycle ele_loop
        endif
      enddo

    endif

! markers... Detectors or IP_L3
! detector markers

    if (ele%key == marker$) then

      if (ele%name(:3) == 'DET') then
        do j = 0, 120
          if (ele%name == cesr%det(j)%name) then
            cesr%ix_cesr(i) = j
            call insert_info (cesr%det(j), ele, i)
            cycle ele_loop
          endif
        enddo

      elseif (ele%name(:5) == 'IP_L3') then
        cesr%ix_ip_l3 = i
      endif

    endif


! overlays: horz & vert steerings, scir cam & tilts, etc.

    if (ele%key == overlay_lord$) then

      if (ele%name(:1) == 'H') then
        do j = 0, 120
          if (ele%type == hsteer_name(j)) then
            call insert_info (cesr%h_steer(j), ele, i)
            cycle ele_loop
          endif
        enddo

      elseif (ele%name(1:1) == 'V') then
        do j = 0, 99
          if (ele%type == vsteer_name(j)) then
            call insert_info (cesr%v_steer(j), ele, i)
            cycle ele_loop
          endif
        enddo

      elseif (ele%name(1:4) == 'SC_V') then
        do j = 101, ubound(cesr%v_steer, 1)
          if (ele%name == cesr%v_steer(j)%name) then
            call insert_info (cesr%v_steer(j), ele, i)
            cycle ele_loop
          endif
        enddo     

      elseif (ele%name(1:6) == 'SC_CAM') then
        do j = 1, size(cesr%scir_cam_rho)
          if (ele%name == cesr%scir_cam_rho(j)%name) then
            call insert_info (cesr%scir_cam_rho(j), ele, i)
            cycle ele_loop
          endif
        enddo

      elseif (ele%name(1:7) == 'SC_TILT') then
        do j = 1, size(cesr%scir_tilt)
          if (ele%name == cesr%scir_tilt(j)%name) then
            call insert_info (cesr%scir_tilt(j), ele, i)
            cycle ele_loop
          endif
        enddo

      elseif (ele%name(1:12) == 'NIR_SHUNTCUR') then
        do j = 1, size(cesr%nir_shuntcur)
          if (ele%name == cesr%nir_shuntcur(j)%name) then
            call insert_info (cesr%nir_shuntcur(j), ele, i)
            cycle ele_loop
          endif
        enddo     

      endif

    endif

! rf

    if (ele%key == rfcavity$) then
      do j = 1, n_rf_maxx
        if (ele%name == cesr%rf(j)%name) then
          call insert_info (cesr%rf(j), ele, i)
          cycle ele_loop
        endif
      enddo
    endif

! wiggler

    if (ele%key == wiggler$) then
      do j = 1, n_wig_maxx
        if (ele%name == cesr%wig(j)%name) then
          call insert_info (cesr%wig(j), ele, i)
          cycle ele_loop
        endif
      enddo
    endif

! solenoid

    if (ele%name == cesr%solenoid%name) then
      cesr%solenoid%ix_lat = i
      cycle ele_loop
    endif

! SC anti-solenoids
    do j = 1, 5*n_sc_sol_maxx
      if (ele%name == cesr%scsol_cur(j)%name) then
        cesr%scsol_cur(j)%ix_lat = i
        cycle ele_loop
      endif
    enddo

  enddo ele_loop
           

!-------------------------------------------------------------------
! Point to quad overlay instead of quad for nir_shuntcur quads

  do i = lat%n_ele_track+1, lat%n_ele_max
    ele = lat%ele(i)
    if (ele%type(1:12) == 'CSR QUAD CUR') then
      read (ele%type(13:16), *) ix
      call insert_info (cesr%quad(ix), ele, i)
    endif
  enddo

!-------------------------------------------------------------------
! check that we have loaded everything...
! do not check Q01 and Q02's

  if (cesr%quad( 1)%ix_lat == 0) cesr%quad( 1)%name = 'DUMMY'
  if (cesr%quad( 2)%ix_lat == 0) cesr%quad( 2)%name = 'DUMMY'
  if (cesr%quad(97)%ix_lat == 0) cesr%quad(97)%name = 'DUMMY'
  if (cesr%quad(98)%ix_lat == 0) cesr%quad(98)%name = 'DUMMY'

  call bmad_to_cesr_err_type (cesr%quad,        'QUADRUPOLE')
  call bmad_to_cesr_err_type (cesr%sep,         'SEPARATOR')
  call bmad_to_cesr_err_type (cesr%skew_sex,    'SKEW SEXTUPOLE')
  call bmad_to_cesr_err_type (cesr%oct,         'OCTUPOLE')
  call bmad_to_cesr_err_type (cesr%wig,         'WIGGLER')
  call bmad_to_cesr_err_type (cesr%rf,          'RF CAVITY')
  call bmad_to_cesr_err_type (cesr%scir_cam_rho,    'SCIR CAM')
  call bmad_to_cesr_err_type (cesr%scir_tilt,   'SCIR TILT')
  call bmad_to_cesr_err_type (cesr%scsol_cur,   'SC SOL CUR')

!----------------------------------------------------------------
contains

subroutine bmad_to_cesr_err_type (cesr_ele, str)

  type (cesr_element_struct) :: cesr_ele(:)
  integer i
  character(*) str

!

  do i = lbound(cesr_ele, 1), ubound(cesr_ele, 1)
    if (cesr_ele(i)%ix_lat == 0 .and.  &
                                  cesr_ele(i)%name(:5) /= 'DUMMY') then
      bmad_status%ok = .false.
      if (bmad_status%type_out) then
        print *, 'WARNING FROM BMAD_TO_CESR. ELEMENT: ', cesr_ele(i)%name
        print *, '        NOT LOADED INTO CESR STRUCT: ', str
      endif
    endif
  enddo

end subroutine

!------------------------------------------------------------------------

subroutine insert_info (cesr_ele, ele, i_ele)

  implicit none

  type (cesr_element_struct)  cesr_ele
  type (ele_struct)  ele
  integer i_ele, ios

!

  cesr_ele%ix_lat = i_ele
  cesr_ele%db_node_name = ele%type(:12)
  cesr_ele%name = ele%name
  if (ele%type(13:) /= '    ') then
    read (ele%type(13:), *, iostat = ios) cesr_ele%ix_db
    if (ios /= 0) then
      print *, 'ERROR IN INSERT_INFO: READ ERROR FOR NODE INDEX: ', ele%type
    endif
  endif

end subroutine

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine choose_cesr_lattice (lattice, lat_file, current_lat, lat, choice)
!
! Subroutine to let the user choose a lattice. The subroutine will present a
! list to choose from.
!                                                               
! Modules Needed:
!   use cesr_basic_mod
!
! Input:
!   current_lat -- Character(40): Name of current lattice (will be stared in
!                       the list presented to the user).
!                       Use CALL GETLAT (CURRENT_LAT) to get the current name.
!                       Set CURRENT_LAT = ' ' if you do not want to use this
!                       feature.
!                       NOTE: You must be connected to the mpm to use GETLAT.
!   choice      -- Character(*), optional: If present then this will be
!                       used as input instead of querying the user.
!                       ''         -> Query the user.
!                       '0'        -> Use current_lat.
!                       <number>   -> Use N^th lattice in list.
!                       <lat_name> -> Use lattice corresponding to this name.
!
! Output:
!   lattice  -- Character(40): Lattice name choisen. If a file name is given
!                    and lat is not present then lattice = ""
!   lat_file -- Character(*): Name of the lattice file. Typically:
!                    lat_file = '$CESR_MNT/lattice/CESR/bmad/bmad_' // lattice // .lat
!   lat     -- lat_struct, optional: If present then bmad_parser is called
!               to load the lat structure.
!-

subroutine choose_cesr_lattice (lattice, lat_file, current_lat, lat, choice)

  implicit none

  type (lat_struct), optional :: lat

  character(len=*), optional :: choice
  character(*) lat_file, lattice, current_lat
  character(40) lat_list(200), cur_lat
  character(80) lat_choise, lat_dir
   
  integer i, num_lats, i_lat, ix, ios

  logical is_there, ask_for_lat, default_flag

!                   

#ifdef CESR_VMS
  lat_dir = '$CESR_MNT/vms_lattice/cesr/bmad/'
#else
  lat_dir = '$CESR_MNT/lattice/cesr/bmad/'
#endif

  call get_lattice_list (lat_list, num_lats, lat_dir)

  ask_for_lat = .true.

  if (present(choice)) then
    lat_choise = choice
#ifdef CESR_VMS
    if (lat_choise == '*') lat_choise = '0'
#endif
    call downcase_string (lat_choise)
    call string_trim (lat_choise, lat_choise, ix)
    if (ix /= 0) ask_for_lat = .false.
  endif

  cur_lat = current_lat
  call downcase_string (cur_lat)

! loop until we have a valid choice

  do

    if (ask_for_lat) then
      print *
      i_lat = 0
      do i = 1, num_lats
        if (lat_list(i) == cur_lat) then
          print '(1x, a, i3, 3a)', '**', i, ') ', trim(lat_list(i)), &
                                           '   ! Current lattice in Data Base'
          i_lat = i
        else
          print '(i5, 2a)', i, ') ', lat_list(i)
        endif
      enddo
  
      print *, ' [Note: To be in this list a lattice file must have a name of the  ]'
      print *, ' [      form: $CESR_MNT/lattice/cesr/bmad/bmad_<lattice_name>.lat  ]'
      print *, ' [        or: $CESR_MNT/lattice/cesr/bmad/<lattice_name>.lat       ]'

      print *
      print *, 'You can enter a Lattice number or a full file name.'
      if (i_lat == 0) then
        write (*, '(a)', advance = 'no') ' Choice: '
      else
        write (*, '(a, i3, a)', advance = 'no') ' Choice: <CR =', i_lat, '> '
      endif
      read (*, '(a)') lat_choise
    endif

    call string_trim (lat_choise, lat_choise, ix)
    lat_choise = lat_choise(:ix)

    if (ix == 0 .or. (ix == 1 .and. lat_choise == '0')) then
      default_flag = .true.
      do i_lat = 1, num_lats
        if (lat_list(i_lat) == cur_lat) exit
      enddo
    else
      default_flag = .false.
      read (lat_choise, *, iostat = ios) i_lat
    endif

    if (default_flag .or. (ios == 0 .and. index('0123456789', lat_choise(1:1)) /= 0)) then
      if (i_lat < 1 .or. i_lat > num_lats) then
        print *, 'ERROR: WHICH LATTICE? TRY AGAIN...'
        ask_for_lat = .true.
        cycle  ! try again
      endif
      lattice = lat_list(i_lat)
      call lattice_to_bmad_file_name (lattice, lat_file)
    else
      lattice = ""
      lat_file = lat_choise
      ! Look to see if this is a file name.
      inquire (file = lat_file, exist = is_there, name = lat_file) 
      if (.not. is_there) then   ! If not a file name then...
        lattice = lat_choise
        lat_file = trim(lat_dir) // 'bmad_' // lattice
        if (index(lattice, '.') == 0) lat_file = trim(lat_file) // '.lat' 
        ix = index(lattice, '.lat')
        if (ix /= 0) lattice = lattice(:ix-1)
        call FullFileName(lat_file, lat_file)
        inquire (file = lat_file, exist = is_there, name = lat_file)
        if (.not. is_there) then
          print *, 'READ ERROR OR FILE DOES NOT EXIST. TRY AGAIN...'
          ask_for_lat = .true.
          cycle
        endif
      endif
      ix = index(lat_file, ';')
      if (ix /= 0) lat_file = lat_file(:ix-1)
    endif
    exit

  enddo

! load lat if present

  if (present (lat)) then
    call bmad_parser (lat_file, lat)
    if (lattice /= "") then
      lat_choise = lat%lattice
      call downcase_string(lat_choise)
      if (lattice /= lat_choise) print *, &
           'WARNING FROM CHOOSE_CESR_LATTICE: LATTICE NAME IN LAT DOES MATCH FILE NAME!'
    endif
    lattice = lat%lattice
  endif

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine create_vsp_volt_elements (lat, ele_type)
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
!   use cesr_basic_mod
!
! Input:
!   lat     -- lat_struct: Lat to be modified
!   ele_type -- Integer: Type of elements to make
!                   = group$      ! make group controller elements
!                   = overlay$    ! make overlay controller elements
!
! Output:
!   lat -- lat_struct: Modified lat.
!-

subroutine create_vsp_volt_elements (lat, ele_type)

  implicit none

  type (lat_struct)  lat

  type vsp_index_struct
    integer ix(3)
  end type
  type (vsp_index_struct) vsp_index(2)

  integer ele_type
  integer i, j

  character(16) :: vsep_name(2) = (/ 'V_SEP_48W', 'V_SEP_48E' /)

  logical found_vsp(2)

! Setup

  vsp_index(1)%ix = (/ 1, 3, 4 /)
  vsp_index(2)%ix = (/ 2, 5, 6 /)

! find vseps

  found_vsp = .false.

  do i = 1, lat%n_ele_track
    do j = 1, 2
      if (lat%ele(i)%name == vsep_name(j)) then
        found_vsp(j) = .true.
        call do_vsp_eles (lat, i, vsp_index(j)%ix, ele_type)
      endif
    enddo
  enddo

! check that the separators were found.

  do j = 1, 2
    if (.not. found_vsp(j)) then
      print *, 'ERROR IN CREATE_VSP_VOLT_ELEMENTS: CANNOT FIND: ', &
                                                            vsep_name(j)
      if (bmad_status%exit_on_error) call err_exit
      bmad_status%ok = .false.
    endif
  enddo

  bmad_status%ok = .true.

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

subroutine do_vsp_eles (lat, i_vsep, vsp_index, ele_type)

  use bmad_struct
  use bmad_interface

  implicit none

  type (lat_struct)  lat
  type (control_struct)  contl(1)

  integer i_vsep, vsp_index(3), ele_type, i, i_con
  real(rp) vkick
  logical err

!
                 
  if (lat%ele(i_vsep)%control_type /= free$) then
    print *, 'ERROR IN CREATE_VSP_VOLT_ELEMENTS: VSEP NOT FREE!', i_vsep
    return
  endif

  lat%ele(i_vsep)%type = ' '

  contl(1)%ix_attrib = vkick$
  contl(1)%coef = 1.0
  contl(1)%ix_slave = i_vsep
  vkick = lat%ele(i_vsep)%value(vkick$)

  do i = 1, 3

    call new_control (lat, i_con)
    write (lat%ele(i_con)%name, '(a, i1)') 'VSP_VOLT_', vsp_index(i)
    write (lat%ele(i_con)%type, '(a, i4)') 'CSR VSP VOLT', vsp_index(i)

    if (ele_type == group$) then
      call create_group (lat, i_con, contl(1:1), err)
      if (err) call err_exit
    elseif (ele_type == overlay$) then
      call create_overlay (lat, i_con, 'VKICK', contl(1:1), err)
      if (err) call err_exit
      if (i == 2 .or. i == 3) lat%ele(i_con)%value(vkick$) = vkick / 2
    else
      print *, 'ERROR IN CREATE_VSP_VOLT_ELEMENTS: BAD ELE_TYPE: ', ele_type
      call err_exit
    endif

  enddo

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine create_nir_shuntcur_elements (lat, ele_type)
!
! Subroutine to create Bmad lattice elements corresponding to the 
! NIR_SHUNTCUR data base elements.
! For each of the 4 north IR quads two overlays controlling the k1 attribute
! are setup. One for the corresponding CSR_QUAD_CUR node and the other
! for the corresponding NIR_SHUNTCUR node
!
! Use BMAD_TO_DB or BMAD_TO_CESR to find where the elements are located.
!
! Modules Needed:
!   use cesr_basic_mod
!
! Input:
!   lat     -- lat_struct: Lat to be modified
!   ele_type -- Integer: Type of elements to make
!                   = group$      ! make group controller elements
!                   = overlay$    ! make overlay controller elements
!
! Output:
!   lat -- lat_struct: Modified lat.
!-

subroutine create_nir_shuntcur_elements (lat)

  use bmad

  implicit none

  type (lat_struct)  lat

  integer i, j


  character(16) :: nir_quad_name(4) = (/ 'Q48W', 'Q49W', 'Q49E', 'Q48E' /)

  logical found_quad(4) 

! find vseps

  found_quad = .false.

  do i = 1, lat%n_ele_max
    do j = 1, 4
      if (lat%ele(i)%name == nir_quad_name(j)) then
        found_quad(j) = .true.
        call do_setup_nir_shuntcur (lat%ele(i), i, j)
      endif
    enddo
  enddo

  do j = 1, 4
    if (.not. found_quad(j)) then
      print *, 'ERROR IN CREATE_NIR_SHUNTCUR_ELEMENTS: CANNOT FIND: ', &
                                                             nir_quad_name(j)
      if (bmad_status%exit_on_error) call err_exit
      bmad_status%ok = .false.
    endif
  enddo

!-----------------------------------------------------------------------
contains

subroutine do_setup_nir_shuntcur (quad, ix_quad, ix_nir)

  implicit none

  type (ele_struct) quad
  type (control_struct)  contl(1)

  integer ix_quad, ix_nir, ix_lord, i, j, i_con, endj
  real(rp) k1
  logical err

!
                 
  do i = quad%ic1_lord, quad%ic2_lord
    ix_lord = lat%control(lat%ic(i))%ix_lord
    if (lat%ele(ix_lord)%control_type == overlay_lord$ .and. &
        lat%ele(ix_lord)%ix_value == k1$) then
      print *, 'ERROR IN CREATE_NIR_SHUNTCUR_ELEMENTS: VSEP NOT FREE!', ix_quad
      return
    endif
  enddo

  contl(1)%ix_attrib = k1$
  contl(1)%coef = 1.0
  contl(1)%ix_slave = ix_quad
  k1 = lat%ele(ix_quad)%value(k1$)

! create control for CSR_QUAD_CUR

  call new_control (lat, i_con)
  call create_overlay (lat, i_con, 'K1', contl(1:1), err)
  if (err) call err_exit
  lat%ele(i_con)%value(k1$) = k1
  lat%ele(i_con)%type = lat%ele(ix_quad)%type
  lat%ele(i_con)%name = lat%ele(i_con)%type
  endj = len(lat%ele(i_con)%name)
  if (endj > 16) endj = 16
  do j = 1, endj
     if (lat%ele(i_con)%name(j:j) == ' ') lat%ele(i_con)%name(j:j) = '_'
  enddo
  lat%ele(ix_quad)%type = ' '

! Create control for NIR_SHUNTCUR

  call new_control (lat, i_con)
  call create_overlay (lat, i_con, 'K1', contl(1:1), err)
  if (err) call err_exit
  write (lat%ele(i_con)%type, '(a12, i4)') 'NIR SHUNTCUR', ix_nir
  lat%ele(i_con)%name = lat%ele(i_con)%type
  endj = len(lat%ele(i_con)%name)
  if (endj > 16) endj = 16
  do j = 1, endj
     if (lat%ele(i_con)%name(j:j) == ' ') lat%ele(i_con)%name(j:j) = '_'
  enddo

end subroutine

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine db_struct_init_cu_now (db)
! 
! Subroutine to initialize a db_struct variable to zero.
!
! Modules Needed:
!   use cesr_basic_mod
!
! Output:
!   db -- Db_struct: Structure to initialize.
!-

subroutine db_struct_init_cu_now (db)

implicit none

type (db_struct) db
integer i

!

call init_this_db_ele (db%csr_quad_cur)
call init_this_db_ele (db%csr_qadd_cur)
call init_this_db_ele (db%csr_horz_cur)
call init_this_db_ele (db%csr_hbnd_cur)
call init_this_db_ele (db%csr_vert_cur)
call init_this_db_ele (db%csr_hsp_volt)
call init_this_db_ele (db%csr_vsp_volt)
call init_this_db_ele (db%csr_sext_cur)
call init_this_db_ele (db%csr_octu_cur)
call init_this_db_ele (db%csr_sqewquad)
call init_this_db_ele (db%csr_scsolcur)
call init_this_db_ele (db%csr_sqewsext)
call init_this_db_ele (db%scir_quadcur)
call init_this_db_ele (db%scir_skqucur)
call init_this_db_ele (db%scir_vertcur)
call init_this_db_ele (db%nir_shuntcur)
call init_this_db_ele (db%ir_sksxcur)
call init_this_db_ele (db%scir_pos_stp)
call init_this_db_ele (db%scir_enc_cnt)
call init_this_db_ele (db%scir_pos_rd)
call init_this_db_ele (db%quad)
call init_this_db_ele (db%detector)
call init_this_db_ele (db%wiggler)
call init_this_db_ele (db%scir_cam_rho)
call init_this_db_ele (db%scir_tilt)

!-------------------------------------------------------
contains
subroutine init_this_db_ele (db_ele)

type (db_element_struct) db_ele(:)
integer j

!

do i = lbound(db_ele, 1), ubound(db_ele, 1)
  db_ele%cu_now = 0
  db_ele%valid_cu_now = .false.
enddo

end subroutine

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine cesr_all_data_struct_init (data)
! 
! Subroutine to initialize a cesr_all_data_struct variable to zero/False.
!
! Modules Needed:
!   use cesr_basic_mod
!
! Output:
!   data -- Cesr_all_data_struct: Structure to initialize.
!     %...%value     -> 0
!     %...%good      -> False
!     %db%...%cu_now -> 0
!-

subroutine cesr_all_data_struct_init (data)

implicit none

type (cesr_all_data_struct) data

!

data%param%var_attrib_name   = ''
data%param%var_ele_name      = ''
data%param%data_type         = ''
data%param%dvar              = 0
data%param%csr_set           = 0
data%param%lattice           = ''
data%param%species           = 0
data%param%horiz_beta_freq   = 0
data%param%vert_beta_freq    = 0
data%param%data_date         = ''
data%param%comment           = ''
data%param%lattice           = ''
data%param%lattice_file_name = ''

data%orbit_x%value  = 0;  data%orbit_x%good   = .false.
data%orbit_y%value  = 0;  data%orbit_y%good   = .false.
data%phase_a%value  = 0;  data%phase_a%good   = .false.
data%phase_b%value  = 0;  data%phase_b%good   = .false.
data%eta_x%value    = 0;  data%eta_x%good     = .false.
data%eta_y%value    = 0;  data%eta_y%good     = .false.
data%cbar11_b%value = 0;  data%cbar11_b%good  = .false.
data%cbar12_a%value = 0;  data%cbar12_a%good  = .false.
data%cbar12_b%value = 0;  data%cbar12_b%good  = .false.
data%cbar22_a%value = 0;  data%cbar22_a%good  = .false.

call db_struct_init_cu_now (data%db)

end subroutine

end module
