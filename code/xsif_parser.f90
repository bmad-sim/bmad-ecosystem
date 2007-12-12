!+
! Subroutine xsif_parser (xsif_file, lat, make_mats6, digested_read_ok, use_line)
!
! Subroutine to parse an XSIF (extended standard input format) lattice file.
! XSIF is used by, for example, LIAR.
! Error messages will be recorded in a local file: 'xsif.err'.
! Standard output messages will be recorded in a local file: 'xsif.out'
!
! Note: The presence of an LCavity element in the input file (even if it is not
! used in the lattice) will make lat%lattice_type = linear_lattice$.
!
! Modules needed:
!   use bmad
!
! Input:
!   xsif_file  -- Character(*): Name of the xsif file.
!   make_mats6 -- Logical, optional: Compute the 6x6 transport matrices for the
!                   Elements? Default is True.
!   use_line   -- Character(*), optional: If present then override the use 
!                   statement in the lattice file and use use_line instead.
!
! Output:
!   lat         -- lat_struct: Structure holding the lattice information.
!   bmad_status  -- Bmad status common block.
!     %ok            -- Set True if parsing is successful. False otherwise.
!-

subroutine xsif_parser (xsif_file, lat, make_mats6, digested_read_ok, use_line)

  use xsif_lat_file_names_mod
  use xsif_inout
  use xsif_interfaces
  use xsif_elements
  use xsif_size_pars
  use xsif_elem_pars
  use bmad_parser_mod, except_dummy => xsif_parser

  implicit none

  type (lat_struct), target :: lat
  type (ele_struct), pointer :: ele

  integer xsif_unit, err_unit, std_out_unit, internal_unit
  integer i, ie, ierr, dat_indx, err_lcl, i_ele, indx, key
  integer ip0, n, it, ix, iep, id, n_names
  integer xsif_io_setup, parchk, ix1, ix2, digested_version

  real(rp) k2, angle, aperture

  character(*) :: xsif_file
  character(*), optional :: use_line
  character(16) name
  character(100) name1, name2, line
  character(200) full_name
  character(200) full_lat_file_name, digested_file
  character(200), allocatable, save :: file_names(:)
  character(16) :: r_name = 'xsif_parser'

  logical, optional :: make_mats6, digested_read_ok
  logical echo_output, err_flag

! See if the digested file is OK

  call form_digested_bmad_file_name (xsif_file, digested_file, full_lat_file_name)
  call read_digested_bmad_file (digested_file, lat, digested_version)

  if (present(use_line)) then
    call str_upcase (name, use_line)
    if (name /= lat%name) bmad_status%ok = .false.
  endif

  if (bmad_status%ok) then
    if (present(digested_read_ok)) digested_read_ok = .true.
    call set_taylor_order (lat%input_taylor_order, .false.)
    call set_ptc (lat%E_TOT, lat%param%particle)
    return
  endif

  if (present(digested_read_ok)) digested_read_ok = .false.

! Init the xsif routines.
! If XSIF_IO_SETUP returns bad status it means that one of the file-open 
!   commands bombed, so we abort execution at this point.

  err_flag = .true.         ! assume the worst
  i = lunget()
  xsif_unit     = i    ! open unit number for reading in the xsif file
  err_unit      = i+1  ! open unit number for error messages
  std_out_unit  = i+2  ! open unit number of standard output messages
  internal_unit = i+3  ! open unit number for a scratch file
  echo_output = .false.     ! lattice will not be echoed to the std out file

  ierr = xsif_io_setup (xsif_file, 'xsif.err', 'xsif.out', xsif_unit, &
               err_unit, std_out_unit, internal_unit,  echo_output, .false.)

  if (ierr /= 0) then
    call xsif_error ('FILE OPEN FAILED FOR OUTPUT FILE.')
    return
  endif


! xsif_cmd_loop parses the file.

  ierr = xsif_cmd_loop ( ) 

  if (present(use_line)) call xuse2 (use_line)

  if (ierr /= 0) then
    call xsif_error ( 'UNABLE TO PARSE LATTICE')
    call xsif_io_close
    return
  endif

  if (.not. line_expanded) then
    call xsif_error ( 'NO "USE" STATEMENT FOUND')
    call xsif_io_close
    return
  endif

! Expand the USEd beamline...
! First perform parameter evaluation

  call parord (err_lcl)

  if (err_lcl /= 0) then
    call xsif_error ('"PARORD" ERROR')
    call xsif_io_close
    return
  endif
  
  call parevl  

  err_lcl = parchk (.true.)

  if (err_lcl /= 0) then
    call xsif_error ('UNDEFINED PARAMETERS')
    call xsif_io_close
    return
  endif

! Init the lattice and look for an LCavity which means this is 
! a linear lattice.

  call init_lat (lat, npos2-npos1+100)
  lat%param%lattice_type = circular_lattice$

  do ie = 1, maxelm
    if (ietyp(ie) == mad_lcav) then
      lat%param%lattice_type = linear_lattice$
      call out_io (s_info$, r_name, &
            'LCavity is present. This implies a Linear Lattice')
      exit
    endif
  enddo

! Allocate elements

  do i = 0, ubound(lat%ele, 1)
    call init_ele (lat%ele(i))
  enddo

  lat%param%particle = positron$
  lat%ele(0)%value(E_TOT$) = 0

  lat%ele(0)%name = 'BEGINNING'     ! Beginning element
  lat%ele(0)%key = init_ele$
  call mat_make_unit (lat%ele(0)%mat6)

! Transfer elements to the lat_struct

  i_ele = 0

  do ie = npos1, npos2-1

    indx = item(ie) 
    if (indx > maxelm) cycle

    aperture = 0
    dat_indx = iedat(indx,1)
    key = ietyp(indx)

    select case (key)

      case (mad_drift)
        call add_ele (drift$)
        ele%value(l$) = pdata(dat_indx)

      case (mad_rbend, mad_sbend)
        call add_ele (sbend$)
        ele%sub_key = sbend$
        ele%value(l$)     = pdata(dat_indx)
        ele%value(angle$) = pdata(dat_indx+1)
        ele%value(k1$)    = pdata(dat_indx+2)
        ele%value(e1$)    = pdata(dat_indx+3)
        ele%value(e2$)    = pdata(dat_indx+4)
        ele%value(tilt$)  = pdata(dat_indx+5)
        k2                = pdata(dat_indx+6)
        ele%value(hgap$)  = pdata(dat_indx+9)
        ele%value(fint$)  = pdata(dat_indx+10)
        ele%value(hgapx$) = pdata(dat_indx+11)
        ele%value(fintx$) = pdata(dat_indx+12)

        if (k2 /= 0) then
          call multipole_init (ele)
          ele%b_pole(2) = k2 / 2
        endif

        if (key == mad_rbend) then  ! transform to an sbend
          angle = ele%value(angle$)
          ele%value(l_chord$) = ele%value(l$)
          ele%value(l$) = ele%value(l_chord$) * angle / (2 * sin(angle/2))
          ele%value(e1$) = ele%value(e1$) + angle / 2
          ele%value(e2$) = ele%value(e2$) + angle / 2
          ele%sub_key = rbend$
        endif

        ele%value(g$) = ele%value(angle$) / ele%value(l$)

      case (mad_quad)
        call add_ele (quadrupole$)
        ele%value(l$)        = pdata(dat_indx)
        ele%value(k1$)       = pdata(dat_indx+1)
        ele%value(tilt$)     = pdata(dat_indx+2)
        aperture = pdata(dat_indx+3)

      case (mad_sext)
        call add_ele (sextupole$)
        ele%value(l$)        = pdata(dat_indx)
        ele%value(k2$)       = pdata(dat_indx+1)
        ele%value(tilt$)     = pdata(dat_indx+2)
        aperture = pdata(dat_indx+3)

      case (mad_octu)
        call add_ele (octupole$)
        ele%value(l$)        = pdata(dat_indx)
        ele%value(k3$)       = pdata(dat_indx+1)
        ele%value(tilt$)     = pdata(dat_indx+2)
        aperture = pdata(dat_indx+3)

      case (mad_multi, mad_dimu)
        call add_ele (multipole$)
        call multipole_init (ele)
        ele%value(l$)        = pdata(dat_indx)
        ele%a_pole(0:20)          = pdata(dat_indx+1:dat_indx+41:2)
        ele%b_pole(0:20)          = pdata(dat_indx+2:dat_indx+42:2)
        aperture = pdata(dat_indx+44)
        ele%value(tilt$)     = pdata(dat_indx+49)

        if (key == mad_dimu) then
          if (ele%value(l$) /= 0) ele%a_pole = ele%a_pole * ele%value(l$)
        endif

        if (abs(pdata(dat_indx+43) - 1) > 1e-6) then
          call xsif_error ('MULTLIPOLE WITH SCALEFAC: ' // ele%name, &
                           'NOT IMPLEMENTED IN BMAD.')
          call err_exit
        endif

        if (pdata(dat_indx+45) /= 0 .or. pdata(dat_indx+46) /= 0) then
          call xsif_error ('MULTLIPOLE WITH SOLENOID: ' // ele%name, &
                           'NOT IMPLEMENTED IN BMAD.')
          call err_exit
        endif

      case (mad_soln)
        call add_ele (solenoid$)
        ele%value(l$)        = pdata(dat_indx)
        ele%value(ks$)       = pdata(dat_indx+1)
        ele%value(k1$)       = pdata(dat_indx+2)
        ele%value(tilt$)     = pdata(dat_indx+3)
        aperture = pdata(dat_indx+4)
        if (ele%value(k1$) /= 0) ele%key = sol_quad$

      case (mad_rfcav)
        call add_ele (rfcavity$)
        ele%value(l$)            = pdata(dat_indx)
        ele%value(voltage$)      = pdata(dat_indx+1) * 1e6
        ele%value(phi0$)         = pdata(dat_indx+2)
        ele%value(harmon$)       = pdata(dat_indx+3)
        ele%value(rf_frequency$) = pdata(dat_indx+10)
        aperture     = pdata(dat_indx+11)

        if (pdata(dat_indx+4) /= 0) then
          call xsif_error ('ENERGY ATTRIBUTE WITH RFCAVITY: ' // ele%name, &
                                'NOT IMPLEMENTED IN BMAD.')
          call err_exit
        endif

        if (pdata(dat_indx+5) /= 0 .or. pdata(dat_indx+6) /= 0) then
          call xsif_error ('WAKEFIELD FILES WITH RFCAVITY: ' // ele%name, &
                                'NOT IMPLEMENTED IN BMAD.')
          call err_exit
        endif

        if (pdata(dat_indx+7) /= 0) then
          call xsif_error ('ELOSS ATTRIBUTE WITH RFCAVITY: ' // ele%name, &
                                'NOT IMPLEMENTED IN BMAD.')
          call err_exit
        endif

      case (mad_sepa)
        call add_ele (elseparator$)
        ele%value(l$)        = pdata(dat_indx)
        ele%value(e_field$)  = pdata(dat_indx+1) * 1e6
        ele%value(tilt$)     = pdata(dat_indx+2)

      case (mad_roll, mad_srot)
        call add_ele (patch$)
        ele%value(tilt$) = pdata(dat_indx)

      case (mad_hkick)
        call add_ele (hkicker$)
        ele%value(l$)     = pdata(dat_indx)
        ele%value(kick$)  = pdata(dat_indx+1)
        ele%value(tilt$)  = pdata(dat_indx+2)
      
      case (mad_vkick)
        call add_ele (vkicker$)
        ele%value(l$)     = pdata(dat_indx)
        ele%value(kick$)  = pdata(dat_indx+1)
        ele%value(tilt$)  = pdata(dat_indx+2)
        
      case (mad_kickmad)
        call add_ele (kicker$)
        ele%value(l$)     = pdata(dat_indx)
        ele%value(hkick$)  = pdata(dat_indx+1)
        ele%value(vkick$)  = pdata(dat_indx+2)
        ele%value(tilt$)  = pdata(dat_indx+3)
  

      case (mad_moni, mad_hmon, mad_vmon)
        call add_ele (monitor$)
        ele%value(l$)     = pdata(dat_indx)

      case (mad_mark)
        call add_ele (marker$)

      case (mad_ecoll)
        call add_ele (ecollimator$)
        ele%value(l$)       = pdata(dat_indx)
        ele%value(x1_limit$) = pdata(dat_indx+1)
        ele%value(x2_limit$) = pdata(dat_indx+1)
        ele%value(y1_limit$) = pdata(dat_indx+2)
        ele%value(y2_limit$) = pdata(dat_indx+2)

      case (mad_rcoll)
        call add_ele (rcollimator$)
        ele%value(l$)       = pdata(dat_indx)
        ele%value(x1_limit$) = pdata(dat_indx+1)
        ele%value(x2_limit$) = pdata(dat_indx+1)
        ele%value(y1_limit$) = pdata(dat_indx+2)
        ele%value(y2_limit$) = pdata(dat_indx+2)


      case (mad_quse)  ! Quad/Sextupole
        call add_ele (quadrupole$)
        ele%value(l$)        = pdata(dat_indx)
        ele%value(k1$)       = pdata(dat_indx+1)
        k2                   = pdata(dat_indx+2)
        ele%value(tilt$)     = pdata(dat_indx+3)
        aperture = pdata(dat_indx+4)

        if (k2 /= 0) then
          call multipole_init (ele)
          ele%b_pole(2) = k2 / 2
        endif

      case (mad_gkick)
        call add_ele (patch$)
        ele%value(x_offset$) = pdata(dat_indx+1)
        ele%value(x_pitch$)  = pdata(dat_indx+2)
        ele%value(y_offset$) = pdata(dat_indx+3)
        ele%value(y_pitch$)  = pdata(dat_indx+4)
        ele%value(l$)        = pdata(dat_indx+5)
        ele%value(pz_offset$)= pdata(dat_indx+6)
        ele%value(tilt$)     = pdata(dat_indx+7)
        ele%value(z_offset$) = pdata(dat_indx+8)

      case (mad_arbit)
        call add_ele (custom$)
        ele%value(l$)       = pdata(dat_indx)
        ele%value(val1$:val12$) = pdata(dat_indx+1:dat_indx+12)
        if (any (pdata(dat_indx+13:dat_indx+20) /= 0)) then
          call xsif_error ('NON-ZERO Pn WITH n > 12 FOR: ' // ele%name)
          call err_exit
        endif

      case (mad_mtwis)
        call add_ele (taylor$)

        allocate (ele%taylor(5)%term(1))
        ele%taylor(5)%term(1) = &
                        taylor_term_struct(1.0_rp, (/ 0, 0, 0, 0, 1, 0 /))

        allocate (ele%taylor(6)%term(1))
        ele%taylor(6)%term(1) = &
                        taylor_term_struct(1.0_rp, (/ 0, 0, 0, 0, 0, 1 /))

        call twiss_to_taylor (dat_indx+0, 0)
        call twiss_to_taylor (dat_indx+3, 2)

      case (mad_matr)
        call add_ele (taylor$)
        do i = 1, 6
          ip0 = dat_indx + 27 * (i-1)
          n = count(pdata(ip0+1:ip0+27) /= 0)
          if (count(pdata(ip0+1:ip0+6) /= 0) == 0) n = n + 1
          allocate (ele%taylor(i)%term(n))
          it = 0
          call add_t_term (ele%taylor(i)%term, it, 0, ip0+1)  ! Rij terms
          if (it == 0) then
            ele%taylor(i)%term(1)%coef = 1
            ele%taylor(i)%term(1)%exp = 0
            ele%taylor(i)%term(1)%exp(i) = 1
            it = 1
          endif
          call add_t_term (ele%taylor(i)%term, it, 1, ip0+7)  ! Ti1k terms
          call add_t_term (ele%taylor(i)%term, it, 2, ip0+13) ! Ti2k terms
          call add_t_term (ele%taylor(i)%term, it, 3, ip0+18) ! Ti3k terms
          call add_t_term (ele%taylor(i)%term, it, 4, ip0+22) ! Ti4k terms
          call add_t_term (ele%taylor(i)%term, it, 5, ip0+25) ! Ti5k terms
          call add_t_term (ele%taylor(i)%term, it, 6, ip0+27) ! Ti6k terms
        enddo

      case (mad_lcav)
        call add_ele (lcavity$)
        ele%value(l$)            = pdata(dat_indx)
        ele%value(gradient$)     = pdata(dat_indx+2) * 1e6 / ele%value(l$)
        ele%value(phi0$)         = pdata(dat_indx+3)
        ele%value(rf_frequency$) = pdata(dat_indx+4) * 1e6
        ele%value(e_loss$)       = pdata(dat_indx+9) 
        aperture     = pdata(dat_indx+12)

        ix1 = nint(pdata(dat_indx+7))
        ix2 = nint(pdata(dat_indx+8))

        if ((ix1 /= 0) .neqv. (ix2 /= 0)) then
          call xsif_error ( &
            'LCAVITY DOES NOT HAVE BOTH L AND T WAKE FILES: ' // ele%name)
          call err_exit
        endif

        if ((ix1 /= 0) .and. (ix2 /= 0)) then
          allocate (ele%wake)
          allocate (ele%wake%sr_mode_long(0), ele%wake%sr_mode_trans(0), ele%wake%lr(0))
          name1 = arr_to_str(lwake_file(ix1)%fnam_ptr)
          name2 = arr_to_str(twake_file(ix2)%fnam_ptr)
          ele%wake%sr_file = 'xsif:: ' // trim(name1) // ' | ' // trim(name2)
          call read_xsif_wake (ele%wake%sr_table, name1, 'LONG')
          call read_xsif_wake (ele%wake%sr_table, name2, 'TRANS')
          ele%wake%z_sr_mode_max = 0
        endif

        lat%param%lattice_type = linear_lattice$

      case (mad_inst, mad_blmo, mad_prof, mad_wire, mad_slmo, mad_imon)
        call add_ele (instrument$)
        ele%value(l$) = pdata(dat_indx)

      case (mad_zrot, mad_yrot)
        call add_ele (patch$)
        ele%value(x_pitch$) = -pdata(dat_indx)

      case default
        write (line, '(a, i8)') 'UNKNOWN ELEMENT TYPE:', key
        call xsif_error (line)
        call err_exit

    end select

    iep = errptr(ie)
    if (iep /= 0) then
      id = iedat(iep, 1)
      ele%value(y_pitch$)  = -pdata(id+3)   ! phi
      ele%value(x_pitch$)  =  pdata(id+4)   ! theta
      ele%value(tilt$)     =  pdata(id+5)   ! psi
      ele%value(x_offset$) = pdata(id)   + ele%value(l$) * ele%value(x_pitch$) / 2
      ele%value(y_offset$) = pdata(id+1) + ele%value(l$) * ele%value(y_pitch$) / 2
      ele%value(s_offset$) = pdata(id+2)
    endif

    ele%value(x1_limit$) = aperture
    ele%value(x2_limit$) = aperture
    ele%value(y1_limit$) = aperture
    ele%value(y2_limit$) = aperture

  enddo

! beam and beta0

  if (ibeta0_ptr /= 0) then
    dat_indx = iedat(ibeta0_ptr, 1)
    ele => lat%ele(0)
    ele%a%beta  = pdata(dat_indx)
    ele%a%alpha = pdata(dat_indx+1)
    ele%a%phi   = pdata(dat_indx+2)
    ele%b%beta  = pdata(dat_indx+3)
    ele%b%alpha = pdata(dat_indx+4)
    ele%b%phi   = pdata(dat_indx+5)
    ele%a%eta   = pdata(dat_indx+6)
    ele%a%etap  = pdata(dat_indx+7)
    ele%b%eta   = pdata(dat_indx+8)
    ele%b%etap  = pdata(dat_indx+9)
    ele%x%eta   = ele%a%eta  
    ele%x%etap  = ele%a%etap 
    ele%y%eta   = ele%b%eta  
    ele%y%etap  = ele%b%etap 
    lat%beam_start%vec(1) =  pdata(dat_indx+10)
    lat%beam_start%vec(2) =  pdata(dat_indx+11)
    lat%beam_start%vec(3) =  pdata(dat_indx+12)
    lat%beam_start%vec(4) =  pdata(dat_indx+13)
    lat%beam_start%vec(5) =  pdata(dat_indx+14)
    lat%beam_start%vec(6) =  pdata(dat_indx+15)
    ele%value(e_tot$) = pdata(dat_indx+26) * 1e9

    if (ele%a%beta /= 0) ele%a%gamma = (1 + ele%a%alpha**2) / ele%a%beta
    if (ele%b%beta /= 0) ele%b%gamma = (1 + ele%b%alpha**2) / ele%b%beta
  endif


  if (ibeam_ptr /= 0) then
    dat_indx = iedat(ibeam_ptr, 1)
    ele => lat%ele(0)
    select case (nint(pdata(dat_indx)))  ! particle type
    case (0) ! default
      lat%param%particle = positron$
    case (1)
      lat%param%particle = positron$
    case (2)
      lat%param%particle = electron$
    case (3)
      lat%param%particle = proton$
    case (4)
      lat%param%particle = antiproton$
    case default
      write (line, '(a, i8)') 'UNKNOWN PARTICLE TYPE:', nint(pdata(dat_indx))
      call xsif_error (line)
      call err_exit
    end select

    lat%param%n_part = pdata(dat_indx+14)
    lat%a%emit = pdata(dat_indx+6)
    lat%b%emit = pdata(dat_indx+8)

    if (pdata(dat_indx+3) /= 0) ele%value(E_TOT$) = &
                                                  pdata(dat_indx+3) * 1e9
  endif

! Global stuff

  inquire (file = xsif_file, name = full_name) 
  lat%input_file_name = full_name    

  lat%name = ' '
  lat%lattice = ' '

  lat%n_ele_track  = i_ele
  lat%n_ele_track = i_ele
  lat%n_ele_max  = i_ele

  lat%version            = bmad_inc_version$
  lat%param%aperture_limit_on  = .true.
  lat%n_ic_max           = 0                     
  lat%n_control_max      = 0    
  lat%name = ktext  ! ktext is global xsif variable

  call set_taylor_order (lat%input_taylor_order, .false.)
  call set_ptc (lat%E_TOT, lat%param%particle)

! Element cleanup

  call compute_reference_energy (lat)
  do i = 1, lat%n_ele_max
    if (ele%key == elseparator$) then
      if (ele%value(E_TOT$) == 0) cycle
      ele%value(vkick$) = ele%value(e_field$) * ele%value(l$) / &
                                                      ele%value(E_TOT$)
    endif
  enddo

! last

  call s_calc (lat)
  call lat_geometry (lat)
  if (logic_option (.true., make_mats6)) call lat_make_mat6 (lat, -1)
  err_flag = .false.
  call xsif_io_close

! Make a digested file

  call file_name_list_show (file_names, n_names)
  do i = 1, n_names
    call fullfilename (file_names(i), file_names(i))
    inquire (file = file_names(i), name = file_names(i))
  enddo
  call write_digested_bmad_file (digested_file, lat, n_names, file_names) 
  deallocate (file_names)

!------------------------------------------------------------------------
contains

subroutine add_ele (key)

  integer key

  i_ele = i_ele + 1
  ele => lat%ele(i_ele)
  ele%key = key
  ele%name = kelem(indx)
  ele%type = ketyp(indx)

  ele%value(ds_step$) = bmad_com%default_ds_step

end subroutine

!------------------------------------------------------------------------
!+
! Subroutine add_t_term (term, it, i1, ip1)
!
! Subroutine to add a taylor term.
! This subroutine makes specific use of the ordering of pdata.
!-

subroutine add_t_term (term, it, i1, ip1)

  type (taylor_term_struct) term(:)
  integer j, it, i1, ip1, ip2

!

  ip2 = ip1 + (6-i1)
  if (i1 == 0) ip2 = ip2 - 1
  
  do j = ip1, ip2
    if (pdata(j) == 0) cycle
    it = it + 1
    term(it)%coef = pdata(j)
    term(it)%exp = 0
    if (i1 /= 0) term(it)%exp(i1) = term(it)%exp(i1) + 1
    ix = 6 - ip2 + j
    term(it)%exp(ix) = term(it)%exp(ix) + 1
  enddo

  

end subroutine

!------------------------------------------------------------------------

subroutine read_xsif_wake (wake, file_name, this)

  type (sr_table_wake_struct), pointer :: wake(:)

  real(rp) s_pos(1000), field(1000), ds

  integer iu, k, n, ios, ix

  character(*) file_name, this
  character(80) line

!

  iu = lunget()
  open (iu, file = file_name, status = 'old', iostat = ios)
  if (ios /= 0) then
    call xsif_error ('CANNOT OPEN WAKE FILE: ' // file_name)
    call err_exit
  endif

  n = 0
  do
    read (iu, '(a)', iostat = ios) line
    if (ios /= 0) exit
    if (line(1:1) == '(') cycle
    call string_trim(line, line, ix)
    if (ix == 0) cycle
    n = n + 1
    read (line, *, iostat = ios) k, s_pos(n), field(n)
    if (ios /= 0) then
      call xsif_error ('CANNOT READ WAKEFILE: ' // file_name, &
                                      'CANNOT READ LINE: ' // line)
      call err_exit
    endif
  enddo
  close (iu)

  if (n < 1) then
    call xsif_error ('WAKE FILE SEEMS TO BE EMPTY: ' // file_name)
    call err_exit
  endif

! Signs: In Bmad +z is pointing forward. In LIAR this is the opposite.

  s_pos = -s_pos

! Check that the wake is ordered correctly.

  ds = s_pos(n) / (n - 1)
  do i = 1, n
    if (abs(s_pos(i) - ds * (i - 1)) > -1e-4 * ds) then
      write (line, '(a, i5)') 'NOT UNIFORMLY ASSENDING. PROBLEM IS INDEX:', i
      call xsif_error ('"S" VALUES IN WAKE FILE: ' // file_name, line)
      call err_exit
    endif
  enddo

! If we have already read in one wake file check for consistancey between
! the current file and the previous one.

  if (associated(wake)) then
    if (size(wake) /= n) then
      call xsif_error ('WAKE FILES HAVE UNEQUAL LENGTHS!', &
                                    'FOR ELEMENT: ' // ele%name)
      call err_exit
    endif
    if (abs(ds - wake(n-1)%z / (n-1)) > -1e-4 * ds) then
      call xsif_error ('WAKE FILES HAVE DIFFERENT dZ BETWEEN POINTS.', &
                                                  'FOR ELEMENT: ' // ele%name)
      call err_exit
    endif
  else
    allocate (wake(0:n-1))
  endif

! Load wake structure.
! Signs: Bmad has a different sign convention for the transverse wake.

  wake%z = s_pos(1:n)  
  if (this == 'LONG') then
    wake%long = field(1:n)
  elseif (this == 'TRANS') then
    wake%trans = -field(1:n)
  else
    call xsif_error ('INTERNAL ERROR!')
    call err_exit
  endif

end subroutine

!------------------------------------------------------------------------

subroutine twiss_to_taylor (d_ix, i0)

  type (twiss_struct) twiss
  real(rp) phi, mat2(2,2)
  integer d_ix, i0

!

  phi         = pdata(d_ix+1)
  twiss%beta  = pdata(d_ix+2)
  twiss%alpha = pdata(d_ix+3)
  

  if (twiss%beta == 0) &
            call xsif_error ('BETA = 0 FOR MTWISS ELEMENT: ' // ele%name)

  call twiss_to_1_turn_mat (twiss, phi, mat2)

  allocate (ele%taylor(i0+1)%term(2))

  ele%taylor(i0+1)%term(1)%exp = 0
  ele%taylor(i0+1)%term(1)%exp(i0+1) = 1
  ele%taylor(i0+1)%term(1)%coef = mat2(1,1)

  ele%taylor(i0+1)%term(2)%exp = 0
  ele%taylor(i0+1)%term(2)%exp(i0+2) = 1
  ele%taylor(i0+1)%term(2)%coef = mat2(1,2)

  allocate (ele%taylor(i0+2)%term(2))

  ele%taylor(i0+2)%term(1)%exp = 0
  ele%taylor(i0+2)%term(1)%exp(i0+1) = 1
  ele%taylor(i0+2)%term(1)%coef = mat2(2,1)

  ele%taylor(i0+2)%term(2)%exp = 0
  ele%taylor(i0+2)%term(2)%exp(i0+2) = 1
  ele%taylor(i0+2)%term(2)%coef = mat2(2,2)

end subroutine

!------------------------------------------------------------------------

subroutine xsif_error (line1, line2)

  character(*) line1
  character(*), optional :: line2

!

  print *, 'ERROR IN XSIF_PARSER: ', trim(line1)
  if (present(line2)) print *, '      ', trim(line2)
  print *, '      FOR XSIF FILE: ', trim(xsif_file)

end subroutine

end subroutine
