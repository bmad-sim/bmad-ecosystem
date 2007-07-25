!+
! Subroutine read_digested_bmad_file (digested_name, lat, version)
!
! Subroutine to read in a digested file. The subroutine will check that
! the version of the digested file is up to date and that the digested file
! is current with respect to the original BMAD files that were used. [See
! write_digested_bmad_file.]
!
! Note: This subroutine also reads in the common structures for BMAD_PARSER2
!
! Modules Needed:
!   use bmad
!
! Input:
!   digested_name -- Character(*): Name of the digested file.
!
! Output:
!   lat         -- lat_struct: Output lattice structure
!   version     -- Integer: bmad_inc_version number stored in the lattice file.
!                   If the file is current this number should be the same 
!                   as the global parameter bmad_inc_version$.
!   bmad_status -- Bmad_status_struct: Common block status structure.
!     %ok           -- Set .false. if read failure.
!-

#include "CESR_platform.inc"

subroutine read_digested_bmad_file (digested_name, lat, version)

  use bmad_struct
  use bmad_interface, except_dummy => read_digested_bmad_file
  use multipole_mod

  implicit none


  type old_twiss_struct
    real(rp) beta, alpha, gamma, phi, eta, etap
    real(rp) eta_lab, etap_lab   ! dispersion along the x or y axis
    real(rp) sigma
  end type  

  type (lat_struct), target, intent(inout) :: lat
  type (ele_struct), pointer :: ele
  type (old_twiss_struct) old_a, old_b, old_z

  integer d_unit, n_files, version, i, j, k, ix
  integer ix_wig, ix_const, ix_r(4), ix_d, ix_m, ix_t(6), ios
  integer ix_sr_table, ix_sr_mode_long, ix_sr_mode_trans, ix_lr, ierr, stat
  integer stat_b(12), idate_old

  character(*) digested_name
  character(200) fname1, fname2, fname3, input_file_name, full_digested_name
  character(200), allocatable :: file_names(:)
  character(25) :: r_name = 'read_digested_bmad_file'

  logical found_it, v80, v81, v82, v83, v_old, v_now

! init all elements in lat

  call init_ele (lat%ele_init)  ! init pointers
  call deallocate_lat_pointers (lat)

! Read the digested file.
! Some old versions can be read even though they are not the current version.

  d_unit = lunget()
  bmad_status%ok = .true.
  lat%n_ele_track = 0

  call fullfilename (digested_name, full_digested_name)
  inquire (file = full_digested_name, name = full_digested_name)
  open (unit = d_unit, file = full_digested_name, status = 'old',  &
                     form = 'unformatted', action = 'READ', err = 9000)

  read (d_unit, err = 9100) n_files, version

  v80 = (version == 80)
  v81 = (version == 81)
  v82 = (version == 82)
  v83 = (version == 83)
  v_old = v80 .or. v81 .or. v82 .or. v83
  v_now = (version == bmad_inc_version$)  ! v84

  if (version < bmad_inc_version$) then
    if (bmad_status%type_out) call out_io (s_warn$, r_name, &
           (/ 'DIGESTED FILE VERSION OUT OF DATE \i4\ > \i4\ ' /),  &
            i_array = (/ bmad_inc_version$, version /) )
    if (v_old) then 
      allocate (file_names(n_files))
      bmad_status%ok = .false.
    else
      close (d_unit)
      bmad_status%ok = .false.
      return
    endif
  endif

  if (version > bmad_inc_version$) then
    if (bmad_status%type_out) call out_io (s_warn$, r_name, &
       'DIGESTED FILE HAS VERSION: \i4\ ', &
       'GREATER THAN VERSION OF THIS PROGRAM: \i4\ ', &
       'WILL NOT USE THE DIGESTED FILE. YOU SHOULD RECOMPILE THIS PROGRAM.', &
       i_array = (/ version, bmad_inc_version$ /) )
    close (d_unit)
    bmad_status%ok = .false.
    return
  endif

! if the digested file is out of date then we still read in the file since
! we can possibly reuse the taylor series.

  do i = 1, n_files
    read (d_unit, err = 9100) fname1, idate_old
      
    if (fname1(1:10) == '!DIGESTED:') then
      fname1 = fname1(11:)
      if (fname1 /= full_digested_name) then
        if (bmad_status%type_out .and. bmad_status%ok) call out_io(s_warn$, &
                                          r_name, ' NOTE: MOVED DIGESTED FILE.')
        bmad_status%ok = .false.
      endif
      cycle
   endif

    if (fname1 == '!RAN FUNCTION WAS CALLED') then
      if (bmad_status%type_out) call out_io(s_warn$, r_name, &
                'NOTE: THE RANDOM NUMBER FUNCTION WAS USED IN THE LATTICE FILE SO THIS', &
                '      LATTICE WILL DIFFER FROM OTHER LATTICES GENERATED FROM THE SAME FILE.')
      bmad_status%ok = .false.
      cycle
    endif

    call simplify_path (fname1, fname1)
    if (v_old) file_names(i) = fname1  ! fake out
    ix = index(fname1, ';')
    stat_b = 0
    if (ix > 0) then    ! has VMS version number
      fname2 = fname1(:ix-1)
    else
      fname2 = fname1
#ifndef CESR_VMS 
      ierr = stat(fname2, stat_b)
#endif
     endif
     inquire (file = fname2, exist = found_it, name = fname3)
     call simplify_path (fname3, fname3)
     if (.not. found_it .or. fname1 /= fname3 .or. &
                                             stat_b(10) /= idate_old) then
       if (bmad_status%type_out .and. bmad_status%ok) call out_io(s_warn$, &
                                      r_name, 'NOTE: DIGESTED FILE OUT OF DATE.')
       bmad_status%ok = .false.
     endif

   enddo

! we read (and write) the lat in pieces since it is
! too big to write in one piece

  read (d_unit, err = 9100)  &   
          lat%name, lat%lattice, lat%input_file_name, lat%title, &
          lat%a, lat%b, lat%z, lat%param, lat%version, lat%n_ele_track, &
          lat%n_ele_track, lat%n_ele_max, &
          lat%n_control_max, lat%n_ic_max, lat%input_taylor_order

  call allocate_lat_ele(lat, lat%n_ele_max+100)
  allocate (lat%control(lat%n_control_max+100))
  allocate (lat%ic(lat%n_ic_max+100))

!

  do i = 0, lat%n_ele_max

    ele => lat%ele(i)

    if (v80) then
      read (d_unit, err = 9100) ix_wig, ix_const, ix_r, ix_d, ix_m, ix_t, &
            ix_sr_table, ix_sr_mode_long, ix_sr_mode_trans, ix_lr, &
            ele%name(1:16), ele%type(1:16), ele%alias(1:16), ele%attribute_name(1:16), &
            old_a, old_b, old_z, ele%value(1:55), ele%gen0, ele%vec0, ele%mat6, &
            ele%c_mat, ele%gamma_c, ele%s, ele%key, ele%floor, &
            ele%is_on, ele%sub_key, ele%control_type, ele%ix_value, &
            ele%n_slave, ele%ix1_slave, ele%ix2_slave, ele%n_lord, &
            ele%ic1_lord, ele%ic2_lord, ele%ix_pointer, ele%ixx, &
            ele%ix_ele, ele%mat6_calc_method, ele%tracking_method, &
            ele%num_steps, ele%integrator_order, ele%ptc_kind, &
            ele%taylor_order, ele%symplectify, ele%mode_flip, &
            ele%multipoles_on, ele%map_with_offsets, ele%Field_master, &
            ele%logic, ele%old_is_on, ele%field_calc, ele%aperture_at, &
            ele%coupler_at, ele%on_an_girder
    elseif (v81) then
      read (d_unit, err = 9100) ix_wig, ix_const, ix_r, ix_d, ix_m, ix_t, &
            ix_sr_table, ix_sr_mode_long, ix_sr_mode_trans, ix_lr, &
            ele%name, ele%type, ele%alias, ele%attribute_name, old_a, &
            old_b, old_z, ele%value(1:55), ele%gen0, ele%vec0, ele%mat6, &
            ele%c_mat, ele%gamma_c, ele%s, ele%key, ele%floor, &
            ele%is_on, ele%sub_key, ele%control_type, ele%ix_value, &
            ele%n_slave, ele%ix1_slave, ele%ix2_slave, ele%n_lord, &
            ele%ic1_lord, ele%ic2_lord, ele%ix_pointer, ele%ixx, &
            ele%ix_ele, ele%mat6_calc_method, ele%tracking_method, &
            ele%num_steps, ele%integrator_order, ele%ptc_kind, &
            ele%taylor_order, ele%symplectify, ele%mode_flip, &
            ele%multipoles_on, ele%map_with_offsets, ele%Field_master, &
            ele%logic, ele%old_is_on, ele%field_calc, ele%aperture_at, &
            ele%coupler_at, ele%on_an_girder
    elseif (v82) then
      read (d_unit, err = 9100) ix_wig, ix_const, ix_r, ix_d, ix_m, ix_t, &
            ix_sr_table, ix_sr_mode_long, ix_sr_mode_trans, ix_lr, &
            ele%name, ele%type, ele%alias, ele%attribute_name, old_a, &
            old_b, old_z, ele%value(1:55), ele%gen0, ele%vec0, ele%mat6, &
            ele%c_mat, ele%gamma_c, ele%s, ele%key, ele%floor, &
            ele%is_on, ele%sub_key, ele%control_type, ele%ix_value, &
            ele%n_slave, ele%ix1_slave, ele%ix2_slave, ele%n_lord, &
            ele%ic1_lord, ele%ic2_lord, ele%ix_pointer, ele%ixx, &
            ele%ix_ele, ele%mat6_calc_method, ele%tracking_method, &
            ele%num_steps, ele%integrator_order, ele%ptc_kind, &
            ele%taylor_order, ele%symplectify, ele%mode_flip, &
            ele%multipoles_on, ele%map_with_offsets, ele%Field_master, &
            ele%logic, ele%old_is_on, ele%field_calc, ele%aperture_at, &
            ele%coupler_at, ele%on_an_girder, ele%csr_calc_on
    elseif (v83) then
      read (d_unit, err = 9100) ix_wig, ix_const, ix_r, ix_d, ix_m, ix_t, &
            ix_sr_table, ix_sr_mode_long, ix_sr_mode_trans, ix_lr, &
            ele%name, ele%type, ele%alias, ele%attribute_name, ele%x, ele%y, &
            ele%a, ele%b, ele%z, ele%value(1:55), ele%gen0, ele%vec0, ele%mat6, &
            ele%c_mat, ele%gamma_c, ele%s, ele%key, ele%floor, &
            ele%is_on, ele%sub_key, ele%control_type, ele%ix_value, &
            ele%n_slave, ele%ix1_slave, ele%ix2_slave, ele%n_lord, &
            ele%ic1_lord, ele%ic2_lord, ele%ix_pointer, ele%ixx, &
            ele%ix_ele, ele%mat6_calc_method, ele%tracking_method, &
            ele%num_steps, ele%integrator_order, ele%ptc_kind, &
            ele%taylor_order, ele%symplectify, ele%mode_flip, &
            ele%multipoles_on, ele%map_with_offsets, ele%Field_master, &
            ele%logic, ele%old_is_on, ele%field_calc, ele%aperture_at, &
            ele%coupler_at, ele%on_an_girder, ele%csr_calc_on
    elseif (v_now) then
      read (d_unit, err = 9100) ix_wig, ix_const, ix_r, ix_d, ix_m, ix_t, &
            ix_sr_table, ix_sr_mode_long, ix_sr_mode_trans, ix_lr, &
            ele%name, ele%type, ele%alias, ele%attribute_name, ele%x, ele%y, &
            ele%a, ele%b, ele%z, ele%value, ele%gen0, ele%vec0, ele%mat6, &
            ele%c_mat, ele%gamma_c, ele%s, ele%key, ele%floor, &
            ele%is_on, ele%sub_key, ele%control_type, ele%ix_value, &
            ele%n_slave, ele%ix1_slave, ele%ix2_slave, ele%n_lord, &
            ele%ic1_lord, ele%ic2_lord, ele%ix_pointer, ele%ixx, &
            ele%ix_ele, ele%mat6_calc_method, ele%tracking_method, &
            ele%num_steps, ele%integrator_order, ele%ptc_kind, &
            ele%taylor_order, ele%symplectify, ele%mode_flip, &
            ele%multipoles_on, ele%map_with_offsets, ele%Field_master, &
            ele%logic, ele%old_is_on, ele%field_calc, ele%aperture_at, &
            ele%coupler_at, ele%on_an_girder, ele%csr_calc_on, &
            ele%ref_orb_in, ele%ref_orb_out

    endif

    ! Transfer twiss

    if (v80 .or. v81 .or. v82) then
      ele%x%eta   = old_a%eta_lab
      ele%x%etap  = old_a%etap_lab

      ele%a%eta   = old_a%eta
      ele%a%etap  = old_a%etap
      ele%a%beta  = old_a%beta
      ele%a%alpha = old_a%alpha
      ele%a%phi   = old_a%phi
      ele%a%gamma = old_a%gamma
      ele%a%sigma = old_a%sigma

      ele%y%eta   = old_b%eta_lab
      ele%y%etap  = old_b%etap_lab

      ele%b%eta   = old_b%eta
      ele%b%etap  = old_b%etap
      ele%b%beta  = old_b%beta
      ele%b%alpha = old_b%alpha
      ele%b%phi   = old_b%phi
      ele%b%gamma = old_b%gamma
      ele%b%sigma = old_b%sigma

      ele%z%eta   = old_z%eta
      ele%z%etap  = old_z%etap
      ele%z%beta  = old_z%beta
      ele%z%alpha = old_z%alpha
      ele%z%phi   = old_z%phi
      ele%z%gamma = old_z%gamma
      ele%z%sigma = old_z%sigma
    endif

    ! l_pole$ attribute did not exist before now.
    if (v80) then
      if (ele%key == wiggler$ .and. ele%sub_key == periodic_type$) then
        if (ele%value(n_pole$) /= 0) ele%value(l_pole$) = ele%value(l$) / ele%value(n_pole$)
      endif
    endif

    if (ix_wig /= 0) then
      allocate (ele%wig_term(ix_wig))
      do j = 1, ix_wig
        read (d_unit, err = 9200) ele%wig_term(j)
      enddo
    endif

    ! This is to cover up the change where periodic_type wigglers now have a single wig_term
    ! where before they did not have any.
    if (ele%key == wiggler$ .and. ele%sub_key == periodic_type$ .and. &
                                                       .not. associated(ele%wig_term)) then
      allocate (ele%wig_term(1))
    endif

    if (ix_const /= 0) then
      allocate (ele%const(ix_const))
      read (d_unit, err = 9300) ele%const
    endif

    if (any (ix_r /= 0)) then
      allocate (ele%r(ix_r(1):ix_r(3), ix_r(2):ix_r(4)))
      read (d_unit, err = 9400) ele%r
    endif

    if (ix_d /= 0) then
      allocate (ele%descrip)
      read (d_unit, err = 9500) ele%descrip
    endif

    if (ix_m /= 0) then
      call multipole_init (ele)
      read (d_unit, err = 9600) ele%a_pole, ele%b_pole
    endif
    
    do j = 1, 6
      if (ix_t(j) == 0) cycle
      read (d_unit) ele%taylor(j)%ref
      allocate (ele%taylor(j)%term(ix_t(j)))
      do k = 1, ix_t(j)
        read (d_unit, err = 9700) ele%taylor(j)%term(k)
      enddo
    enddo

    ! If ix_lr is negative then it is a pointer to a previously read wake. 
    ! See write_digested_bmad_file.

    if (ix_sr_table /= 0 .or. ix_sr_mode_long /= 0 .or. ix_sr_mode_trans /= 0 .or. ix_lr /= 0) then
      if (ix_lr < 0) then
        call transfer_wake (lat%ele(abs(ix_lr))%wake, ele%wake)

      else
        call init_wake (ele%wake, ix_sr_table, ix_sr_mode_long, ix_sr_mode_trans, ix_lr)
        read (d_unit, err = 9800) ele%wake%sr_file
        read (d_unit, err = 9810) ele%wake%sr_table
        read (d_unit, err = 9840) ele%wake%sr_mode_long
        read (d_unit, err = 9850) ele%wake%sr_mode_trans
        read (d_unit, err = 9820) ele%wake%lr_file
        read (d_unit, err = 9830) ele%wake%lr
        read (d_unit, err = 9860) ele%wake%z_sr_mode_max
      endif
    endif

  enddo

!

  do i = 1, lat%n_control_max
    read (d_unit, err = 9100) lat%control(i)
  enddo

  do i = 1, lat%n_ic_max
    read (d_unit, err = 9100) lat%ic(i)
  enddo

  read (d_unit, iostat = ios) lat%beam_start

  close (d_unit)

  lat%param%stable = .true.  ! Assume this 

  return

!------------------

9000  continue
  if (bmad_status%type_out) then
     call out_io (s_warn$, r_name, 'DIGESTED FILE DOES NOT EXIST.')
  endif
  close (d_unit)
  bmad_status%ok = .false.
  version = -1
  return

!--------------------------------------------------------------

9100  continue
  if (bmad_status%type_out) then
     call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.')
  endif
  close (d_unit)
  bmad_status%ok = .false.
  return

9200  continue
  if (bmad_status%type_out) then
     call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING WIGGLER TERM FOR ELEMENT: ' // ele%name)
  endif
  close (d_unit)
  bmad_status%ok = .false.
  return

9300  continue
  if (bmad_status%type_out) then
     call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING IX_CONST TERM FOR ELEMENT: ' // ele%name)
  endif
  close (d_unit)
  bmad_status%ok = .false.
  return

9400  continue
  if (bmad_status%type_out) then
     call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING R TERM FOR ELEMENT: ' // ele%name)
  endif
  close (d_unit)
  bmad_status%ok = .false.
  return

9500  continue
  if (bmad_status%type_out) then
     call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING DESCRIP TERM FOR ELEMENT: ' // ele%name)
  endif
  close (d_unit)
  bmad_status%ok = .false.
  return

9600  continue
  if (bmad_status%type_out) then
     call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING AN,BN TERMS FOR ELEMENT: ' // ele%name)
  endif
  close (d_unit)
  bmad_status%ok = .false.
  return

9700  continue
  if (bmad_status%type_out) then
     call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING TAYLOR TERMS FOR ELEMENT: ' // ele%name)
  endif
  close (d_unit)
  bmad_status%ok = .false.
  return

9800  continue
  if (bmad_status%type_out) then
     call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING WAKE%SR_FILE FOR ELEMENT: ' // ele%name)
  endif
  close (d_unit)
  bmad_status%ok = .false.
  return

9810  continue
  if (bmad_status%type_out) then
     call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING WAKE%sr_table FOR ELEMENT: ' // ele%name)
  endif
  close (d_unit)
  bmad_status%ok = .false.
  return

9820  continue
  if (bmad_status%type_out) then
     call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING WAKE%LR_FILE FOR ELEMENT: ' // ele%name)
  endif
  close (d_unit)
  bmad_status%ok = .false.
  return

9830  continue
  if (bmad_status%type_out) then
     call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING WAKE%LR FOR ELEMENT: ' // ele%name)
  endif
  close (d_unit)
  bmad_status%ok = .false.
  return

9840  continue
  if (bmad_status%type_out) then
     call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING WAKE%sr_mode_long FOR ELEMENT: ' // ele%name)
  endif
  close (d_unit)
  bmad_status%ok = .false.
  return

9850  continue
  if (bmad_status%type_out) then
     call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING WAKE%sr_mode_trans FOR ELEMENT: ' // ele%name)
  endif
  close (d_unit)
  bmad_status%ok = .false.
  return

9860  continue
  if (bmad_status%type_out) then
     call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING WAKE%Z_CUT_SR FOR ELEMENT: ' // ele%name)
  endif
  close (d_unit)
  bmad_status%ok = .false.
  return

!--------------------------------------------------------------------
!--------------------------------------------------------------------
contains

subroutine simplify_path (name_in, name_out)

  implicit none

  character(*) name_in, name_out
  integer i, ix

! 

  name_out = name_in
  out_loop: do 
    ix = index(name_out, '/..')
    if (ix == 0) return
    do i = ix-1, 1, -1
      if (name_out(i:i) == '/') then
        name_out = name_out(:i-1) // name_out(ix+3:)
        cycle out_loop
      endif
    enddo
    name_out = name_out(ix+3:)
  enddo out_loop

end subroutine

end subroutine
