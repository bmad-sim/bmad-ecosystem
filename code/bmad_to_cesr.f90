!+
! Subroutine BMAD_TO_CESR (RING, CESR)
!
! Subroutine to transfer information from the RING structure returned from
! BAMD_PARSER to a structure for the CESR ring
!
! Modules Needed:
!   use bmad_struct
!   use bmad_interface
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

!$Id$
!$Log$
!Revision 1.4  2001/10/08 17:18:14  rwh24
!DCS changes to f90 files.
!Bug fixes to c file.
!
!Revision 1.3  2001/10/02 18:49:11  rwh24
!More compatibility updates; also added many explicit variable declarations.
!
!Revision 1.2  2001/09/27 18:31:48  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"



subroutine bmad_to_cesr (ring, cesr)

  use bmad_struct
  implicit none

  type (ring_struct)  ring
  type (cesr_struct)  cesr
  type (ele_struct)  ele

  character cc1*1, cc2*2, cc4*4
  character*16 hsteer_name(0:120), vsteer_name(0:99)

  integer j, i, vnumbr

  logical phase_iii 

! find out if this is a phase_iii lattice

  phase_iii = .false.
  do i = ring%n_ele_ring+1, ring%n_ele_max
    if (ring%ele_(i)%name(1:5) == 'SCIR_') then
      phase_iii = .true.
      exit
    endif
  end do

! load names

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

  if (phase_iii) then

    cesr%v_steer_(101)%name = 'SCIR_V01W'
    cesr%v_steer_(102)%name = 'SCIR_V02W'
    cesr%v_steer_(103)%name = 'SCIR_V02E'
    cesr%v_steer_(104)%name = 'SCIR_V01E'

    do i = 1, 5                          
      write (cesr%scir_cam_rho_(i)%name,   '(a, i1, a)') 'SCIR_CAM_', i, 'W'
      write (cesr%scir_cam_rho_(i+5)%name, '(a, i1, a)') 'SCIR_CAM_', i, 'E'
    enddo

    cesr%scir_tilt_(scir_tilt_w$)%name    = 'SCIR_TILT_W'
    cesr%scir_tilt_(scir_tilt_e$)%name    = 'SCIR_TILT_E'
    cesr%scir_tilt_(scir_tilt_sk_w$)%name = 'SCIR_TILT_SK_W'
    cesr%scir_tilt_(scir_tilt_sk_e$)%name = 'SCIR_TILT_SK_E'

    cesr%skew_quad_( 1)%name = 'DUMMY'
    cesr%skew_quad_( 2)%name = 'SK_Q03W1'
    cesr%skew_quad_( 3)%name = 'SK_Q03W2'
    cesr%skew_quad_(96)%name = 'SK_Q03E2'
    cesr%skew_quad_(97)%name = 'SK_Q03E1'
    cesr%skew_quad_(98)%name = 'DUMMY'

    cesr%skew_quad_(101)%name = 'SK_Q01W'
    cesr%skew_quad_(102)%name = 'SK_Q02W'
    cesr%skew_quad_(103)%name = 'SK_Q02E'
    cesr%skew_quad_(104)%name = 'SK_Q01E'

    cesr%skew_sex_(11)%name = 'SK_SEX_02E'

  endif


!-------------------------------------------------------------
! Load elements from RING to CESR

  ele_loop: do i = 1, ring%n_ele_max

    ele = ring%ele_(i)

! quads and skew quads

    if (ele%key == quadrupole$) then

      if (ele%name(1:1) == 'Q') then
        do j = 0, 120
          if (ele%name == cesr%quad_(j)%name) then
            cesr%ix_cesr(i) = j
            call insert_info (cesr%quad_(j), ele, i)
            cycle ele_loop
          endif
        enddo

      elseif (ele%name(:2) == 'SK') then
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

      elseif (ele%name(1:6) == 'SCIR_V') then
        do j = 101, size(cesr%v_steer_)
          if (ele%name == cesr%v_steer_(j)%name) then
            call insert_info (cesr%v_steer_(j), ele, i)
            cycle ele_loop
          endif
        enddo     

      elseif (ele%name(1:8) == 'SCIR_CAM') then
        do j = 1, size(cesr%scir_cam_rho_)
          if (ele%name == cesr%scir_cam_rho_(j)%name) then
            call insert_info (cesr%scir_cam_rho_(j), ele, i)
            cycle ele_loop
          endif
        enddo

      elseif (ele%name(1:9) == 'SCIR_TILT') then
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

  use bmad_struct

  type (cesr_element_struct) :: cesr_ele(:)
  integer i
  character*(*) str

!

  do i = lbound(cesr_ele, 1), ubound(cesr_ele, 1)
    if (cesr_ele(i)%ix_ring == 0 .and.  &
                                  cesr_ele(i)%name(:5) /= 'DUMMY') then
      bmad_status%ok = .false.
      if (bmad_status%type_out) then
        type *, 'WARNING FROM BMAD_TO_CESR. ELEMENT NOT LOADED INTO'
        type *, '        CESR STRUCT. ', str, ': ', cesr_ele(i)%name
      endif
    endif
  enddo

end subroutine

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------

subroutine insert_info (cesr_ele, ele, i_ele)

  use bmad_struct

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
      type *, 'ERROR IN INSERT_INFO: READ ERROR FOR NODE INDEX: ', ele%type
    endif
  endif

end subroutine
