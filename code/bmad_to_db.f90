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

!$Id$
!$Log$
!Revision 1.12  2003/07/09 01:38:10  dcs
!new bmad with allocatable ring%ele_(:)
!
!Revision 1.11  2003/05/02 15:43:58  dcs
!F90 standard conforming changes.
!
!Revision 1.10  2003/03/04 16:03:28  dcs
!VMS port
!
!Revision 1.9  2003/01/27 14:40:30  dcs
!bmad_version = 56
!
!Revision 1.8  2002/09/14 19:45:24  dcs
!*** empty log message ***
!
!Revision 1.7  2002/06/13 14:54:22  dcs
!Interfaced with FPP/PTC
!
!Revision 1.6  2002/02/23 20:32:11  dcs
!Double/Single Real toggle added
!
!Revision 1.5  2002/01/08 21:44:37  dcs
!Aligned with VMS version  -- DCS
!
!Revision 1.4  2001/10/12 20:53:35  rwh24
!DCS changes and two files added
!
!Revision 1.3  2001/10/02 18:49:11  rwh24
!More compatibility updates; also added many explicit variable declarations.
!
!Revision 1.2  2001/09/27 18:31:48  rwh24
!UNIX compatibility updates
!


subroutine bmad_to_db (ring, db)

  use bmad_struct
  use bmad_interface

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
