!+
! Subroutine READ_DIGESTED_BMAD_FILE (IN_FILE_NAME, RING, VERSION)
!
! Subroutine to read in a digested file. The subroutine will check that
! the version of the digested file is up to date and that the digested file
! is current with respect to the original BMAD files that were used. [See
! WRITE_DIGESTED_BMAD_FILE.F77]
!
! Note: This subroutine also reads in the common structures for BMAD_PARSER2
!
! Modules Needed:
!   use bmad
!
! Input:
!     IN_FILE_NAME -- Character*(*): Name of the digested file
!
! Output:
!     RING      -- Ring_struct: Output structure
!     VERSION   -- Integer: Version number of RING.
!     STATUS    -- Common block status structure
!       .OK       -- Set .false. if read failure.
!-

#include "CESR_platform.inc"

subroutine read_digested_bmad_file (in_file_name, ring, version)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ring_struct), target, intent(out) :: ring
  type (ele_struct), pointer :: ele
  
  integer d_unit, lunget, n_files, version, i, j, k, ix
  integer ix_w, ix_d, ix_m, ix_t(6)
  integer stat_b(12), stat, ierr, idate_old

  character*(*) in_file_name
  character*200 fname(3)

  logical found_it

! init all elements in ring

  call init_ele (ring%ele_init)  ! init pointers

  do i = 0, n_ele_maxx
    call init_ele (ring%ele_(i))
  enddo

! read the digested file

  d_unit = lunget()
  bmad_status%ok = .true.
  ring%n_ele_ring = 0

  open (unit = d_unit, file = in_file_name, status = 'old',  &
                     form = 'unformatted', action = 'READ', err = 9000)

  read (d_unit, err = 9100) n_files, version

  if (version < bmad_inc_version$) then
    if (bmad_status%type_out) print '(1x, a, i4, a, i4)',  &
           'READ_DIGESTED_BMAD_FILE: DIGESTED FILE VERSION OUT OF DATE',  &
            version, ' <', bmad_inc_version$
    close (d_unit)
    bmad_status%ok = .false.
    return
  endif

  if (version > bmad_inc_version$) then
    if (bmad_status%type_out) then
      print *, 'READ_DIGESTED_BMAD_FILE: DIGESTED FILE HAS VERSION:',  &
                                                              version
      print *, '     GREATER THAN VERSION OF THIS PROGRAM:',  &
                                                  bmad_inc_version$
      print *, '     WILL NOT USE THE DIGESTED FILE.'
      print *, '     YOU SHOULD RECOMPILE THIS PROGRAM.'
    endif
    close (d_unit)
    bmad_status%ok = .false.
    return
  endif

! if the digested file is out of date then we still read in the file since
! we can possibly reuse the taylor series.

  do i = 1, n_files
    read (d_unit, err = 9100) fname(1), idate_old
    ix = index(fname(1), ';')
    stat_b = 0
    if (ix > 0) then    ! has VMS version number
      fname(2) = fname(1)(:ix-1)
    else
      fname(2) = fname(1)
#ifdef CESR_UNIX
      ierr = stat(fname(2), stat_b)
#endif
    endif
    inquire (file = fname(2), exist = found_it, name = fname(3))
    if (.not. found_it .or. fname(1) /= fname(3) .or. &
                                             stat_b(10) /= idate_old) then
      if (bmad_status%type_out .and. bmad_status%ok) print *, &
              'READ_DIGESTED_BMAD_FILE: WARNING: DIGESTED FILE OUT OF DATE.'
      bmad_status%ok = .false.
    endif
    if (i == 1 .and. fname(2) /= ring%input_file_name) then
      if (bmad_status%type_out .and. bmad_status%ok) print *, &
                    'READ_DIGESTED_BMAD_FILE: WARNING: MOVED DIGESTED FILE.'
      bmad_status%ok = .false.
    endif
   enddo

! we read (and write) the ring in pieces since it is
! too big to write in one piece

  read (d_unit, err = 9100)  &   
          ring%name, ring%lattice, ring%input_file_name, ring%title, &
          ring%x, ring%y, ring%z, ring%param, ring%version, ring%n_ele_ring, &
          ring%n_ele_symm, ring%n_ele_use, ring%n_ele_max, &
          ring%n_control_array, ring%n_ic_array, ring%input_taylor_order

  if (ring%n_ele_max > n_ele_maxx) then
    print *, 'ERROR IN READ_DIGESTED_BMAD_FILE: NUMBER OF ELEMENTS:',  &
                                                            ring%n_ele_max
    print *, '      IS GREATER THAN THE ELEMENT ARRAY SIZE:', n_ele_maxx
    print *, '      YOU NEED TO RECOMPILE'
    call err_exit
  endif

!

  do i = 0, ring%n_ele_max

    ele => ring%ele_(i)
    read (d_unit, err = 9100) ix_w, ix_d, ix_m, ix_t, &
            ele%name, ele%type, ele%alias, ele%attribute_name, ele%x, &
            ele%y, ele%z, ele%value, ele%gen0, ele%vec0, ele%mat6, &
            ele%c_mat, ele%gamma_c, ele%s, ele%x_position, ele%y_position, &
            ele%z_position, ele%theta_position, ele%phi_position, ele%key, &
            ele%is_on, ele%sub_key, ele%control_type, ele%ix_value, &
            ele%n_slave, ele%ix1_slave, ele%ix2_slave, ele%n_lord, &
            ele%ic1_lord, ele%ic2_lord, ele%ix_pointer, ele%ixx, &
            ele%iyy, ele%mat6_calc_method, ele%tracking_method, &
            ele%num_steps, ele%integration_order, ele%ptc_kind, &
            ele%taylor_order, ele%symplectify, ele%mode_flip, &
            ele%multipoles_on, ele%exact_rad_int_calc, ele%B_field_master

    ele%pointer_init = 0  ! signal that pointers have garbage
    call deallocate_ele_pointers (ele)  ! and deallocate 

    if (ix_w /= 0) then
      allocate (ele%wig_term(ix_w))
      do j = 1, ix_w
        read (d_unit) ele%wig_term(j)
      enddo
    endif

    if (ix_d /= 0) then
      allocate (ele%descrip)
      read (d_unit) ele%descrip
    endif

    if (ix_m /= 0) then
      allocate (ele%a(0:n_pole_maxx), ele%b(0:n_pole_maxx))
      read (d_unit) ele%a, ele%b
    endif
    
    do j = 1, 6
      if (ix_t(j) /= 0) then
        allocate (ele%taylor(j)%term(ix_t(j)))
        do k = 1, ix_t(j)
          read (d_unit) ele%taylor(j)%term(k)
        enddo
      endif
    enddo
    
  enddo

!

  if (ring%n_control_array > n_control_maxx) then
    print *, 'ERROR IN READ_DIGESTED_BMAD_FILE: NUMBER OF ELEMENTS:',  &
                                                      ring%n_control_array
    print *, '      IS GREATER THAN THE CONTROL ARRAY SIZE:', n_control_maxx
    print *, '      YOU NEED TO RECOMPILE'
    call err_exit
  endif

  do i = 1, ring%n_control_array
    read (d_unit, err = 9100) ring%control_(i)
  enddo

  do i = 1, ring%n_ic_array
    read (d_unit, err = 9100) ring%ic_(i)
  enddo

  close (d_unit)
  return

!------------------

9000  continue
  if (bmad_status%type_out) then
    print *, 'READ_DIGESTED_BMAD_FILE: DIGESTED FILE DOES NOT EXIST.'
  endif
  close (d_unit)
  bmad_status%ok = .false.
  version = -1
  return

9100  continue
  if (bmad_status%type_out) then
    print *, 'READ_DIGESTED_BMAD_FILE: ERROR READING DIGESTED FILE.'
  endif
  close (d_unit)
  bmad_status%ok = .false.
  return

end subroutine
