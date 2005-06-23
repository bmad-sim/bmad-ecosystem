module cesr_basic_mod

  use bmad_struct
  use bmad_interface

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
    real beta_mid             ! beta at the midpoint
    real beta_ave             ! beta averaged over the element
    real pos_mid              ! position at midpoint
    real pos_inj              ! position in injection lattice
  end type

  type cesr_element_struct 
    character(16) name              ! bmad name
    type (b_struct)  x, y          ! beta's and positions
    character(12) db_node_name
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

!-----------------------------------------------------------------------------
! For butns.nnnnn files

  type detector_struct
    real x_orb, y_orb
    integer amp(4)
    integer type
    logical ok
  endtype

  type butns_struct
    character(40) lattice
    character(20) date
    character(72) comment(5)
    type (detector_struct) det(0:120)
    integer save_set
    integer file_num
    integer turn    ! turn number for injection data
  end type

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine bmad_to_cesr (ring, cesr)
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
!   ring      -- Ring_struct: Ring to parse.
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

subroutine bmad_to_cesr (ring, cesr)

  implicit none

  type (ring_struct)  ring
  type (cesr_struct)  cesr
  type (ele_struct)  ele

  character(2) cc2
  character(4) cc4
  character(16) hsteer_name(0:120), vsteer_name(0:99)

  integer j, i

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

  cesr%skew_sex_(1)%name = 'SK_SEX_07W'
  cesr%skew_sex_(2)%name = 'SK_SEX_23W'
  cesr%skew_sex_(3)%name = 'SK_SEX_23E'
  cesr%skew_sex_(4)%name = 'SK_SEX_07E'
  cesr%skew_sex_(5)%name = 'SK_SEX_29W'
  cesr%skew_sex_(6)%name = 'SK_SEX_29E'
  cesr%skew_sex_(7)%name = 'SK_SEX_12W'
  cesr%skew_sex_(8)%name = 'SK_SEX_12E'
  cesr%skew_sex_(11)%name = 'SK_SEX_02E'

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
  character(*) str

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

!------------------------------------------------------------------------

subroutine insert_info (cesr_ele, ele, i_ele)

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

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine choose_cesr_lattice (lattice, lat_file, current_lat, ring, choice)
!
! Subroutine to let the user choose a lattice. The subroutine will present a
! list to choose from.
!                                                               
! Modules Needed:
!   use bmad
!
! Input:
!   current_lat -- Character(40): Name of current lattice (will be stared in
!                       the list presented to the user).
!                       Use CALL GETLAT (CURRENT_LAT) to get the current name.
!                       Set CURRENT_LAT = ' ' if you do not want to use this
!                       feature.
!                       NOTE: You must be connected to the mpm to use GETLAT.
!   choice      -- Character(*): [Optional] If present then this will be
!                       used as input instead of querying the user.
!
! Output:
!   lattice  -- Character(40): Lattice name choisen. If a file name is given
!                    and RING is not present then LATTICE = ""
!   lat_file -- Character(*): Name of the lattice file. Typically:
!                    lat_file = 'U:[CESR.BMAD.LAT]BMAD_' // lattice // .LAT
!   ring     -- Ring_struct: OPTIONAL. If present then BMAD_PARSER is called
!               to load the RING structure.
!-

#include "CESR_platform.inc"
               
subroutine choose_cesr_lattice (lattice, lat_file, current_lat, ring, choice)

  implicit none

  type (ring_struct), optional :: ring

  character(len=*), optional :: choice
  character(*) lat_file
  character(40) lattice, current_lat, lat_list(200)
  character(80) line
   
  integer i, num_lats, i_lat, ix, ios

  logical is_there, ask_for_lat, default_flag

!                   

  call get_lattice_list (lat_list, num_lats, 'BMAD_LAT:')

  ask_for_lat = .true.

  if (present(choice)) then
    line = choice
    call string_trim (line, line, ix)
    if (ix /= 0) ask_for_lat = .false.
  endif

! loop until we have a valid choice

  do

    if (ask_for_lat) then
      print *
      i_lat = 0
      do i = 1, num_lats
        if (lat_list(i) == current_lat) then
          print '(1x, a, i3, 3a)', '**', i, ') ', trim(lat_list(i)), &
                                           '   ! Current lattice in Data Base'
          i_lat = i
        else
          print '(i5, 2a)', i, ') ', lat_list(i)
        endif
      enddo
  
      print *, ' [Note: To be in this list a lattice file must have a name   ]'
      print *, ' [      of the form: U:[CESR.BMAD.LAT]bmad_<lattice_name>.lat]'

      print *
      print *, 'You can enter a Lattice number or a full file name.'
      if (i_lat == 0) then
        write (*, '(a)', advance = 'no') ' Choice: '
      else
        write (*, '(a, i3, a)', advance = 'no') ' Choice: <CR =', i_lat, '> '
      endif
      read (*, '(a)') line
    endif

    call string_trim (line, line, ix)
    line = line(:ix)

    if (ix == 0 .or. (ix == 1 .and. line == '*')) then
      default_flag = .true.
      do i_lat = 1, num_lats
        if (lat_list(i_lat) == current_lat) exit
      enddo
    else
      default_flag = .false.
      read (line, *, iostat = ios) i_lat
    endif

    if (default_flag .or. (ios == 0 .and. index('0123456789', line(1:1)) /= 0)) then
      if (i_lat < 1 .or. i_lat > num_lats) then
        print *, 'ERROR: WHICH LATTICE? TRY AGAIN...'
        ask_for_lat = .true.
        cycle  ! try again
      endif
      lattice = lat_list(i_lat)
      call lattice_to_bmad_file_name (lattice, lat_file)
    else
      lattice = ""
      lat_file = line
      inquire (file = lat_file, exist = is_there, name = lat_file)
      if (.not. is_there) then
        lattice = line
        lat_file = 'BMAD_LAT:bmad_' // lattice
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

! load ring if present

  if (present (ring)) then
    call bmad_parser (lat_file, ring)
    if (lattice /= "") then
      if (lattice /= ring%lattice) print *, &
           'WARNING FROM CHOOSE_CESR_LATTICE: LATTICE NAME IN RING DOES MATCH FILE NAME!'
    endif
    lattice = ring%lattice
  endif

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine create_vsp_volt_elements (ring, ele_type)
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
!   use bmad
!
! Input:
!   ring     -- Ring_struct: Ring to be modified
!   ele_type -- Integer: Type of elements to make
!                   = group$      ! make group controller elements
!                   = overlay$    ! make overlay controller elements
!
! Output:
!   ring -- Ring_struct: Modified ring.
!-

#include "CESR_platform.inc"

subroutine create_vsp_volt_elements (ring, ele_type)

  implicit none

  type (ring_struct)  ring

  integer ele_type
  integer :: ix_west(3) = (/ 1, 3, 4 /), ix_east(3) = (/ 2, 5, 6 /)
  integer i

  character(16) :: vsep_west = 'V_SEP_48W', vsep_east = 'V_SEP_48E'

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
        print *, 'ERROR IN CREATE_VSP_VOLT_ELEMENTS: CANNOT FIND WEST VSEP!'
        if (bmad_status%exit_on_error) call err_exit
        bmad_status%ok = .false.
      endif
      return
    endif

  enddo

  print *, 'ERROR IN CREATE_VSP_VOLT_ELEMENTS: CANNOT FIND EAST VSEP!'
  if (bmad_status%exit_on_error) call err_exit
  bmad_status%ok = .false.

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

subroutine do_vsp_eles (ring, i_vsep, ix_, ele_type)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ring_struct)  ring
  type (control_struct)  contl(1)

  integer i_vsep, ix_(3), ele_type, i, i_con

  real(rp) vkick

!
                 
  if (ring%ele_(i_vsep)%control_type /= free$) then
    print *, 'ERROR IN CREATE_VSP_VOLT_ELEMENTS: VSEP NOT FREE!', i_vsep
    return
  endif

  ring%ele_(i_vsep)%type = ' '

  contl(1)%ix_attrib = vkick$
  contl(1)%coef = 1.0
  contl(1)%ix_slave = i_vsep
  vkick = ring%ele_(i_vsep)%value(vkick$)

  do i = 1, 3

    call new_control (ring, i_con)
    write (ring%ele_(i_con)%name, '(a, i1)') 'VSP_VOLT_', ix_(i)
    write (ring%ele_(i_con)%type, '(a, i4)') 'CSR VSP VOLT', ix_(i)

    if (ele_type == group$) then
      call create_group (ring, i_con, contl(1:1))
    elseif (ele_type == overlay$) then
      call create_overlay (ring, i_con, 'VKICK', contl(1:1))
      if (i == 2 .or. i == 3) ring%ele_(i_con)%value(vkick$) = vkick / 2
    else
      print *, 'ERROR IN CREATE_VSP_VOLT_ELEMENTS: BAD ELE_TYPE: ', ele_type
      call err_exit
    endif

  enddo

end subroutine


end module
