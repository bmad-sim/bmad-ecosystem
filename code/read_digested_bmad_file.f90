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

!$Id$
!$Log$
!Revision 1.7  2002/10/21 20:29:53  dcs
!Init all elements
!
!Revision 1.6  2002/09/14 19:45:24  dcs
!*** empty log message ***
!
!Revision 1.5  2002/07/31 14:32:41  dcs
!Modified so moved digested file handled correctly.
!
!Revision 1.4  2002/06/13 14:54:28  dcs
!Interfaced with FPP/PTC
!
!Revision 1.3  2002/02/23 20:32:23  dcs
!Double/Single Real toggle added
!
!Revision 1.2  2001/09/27 18:31:56  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine read_digested_bmad_file (in_file_name, ring, version)

  use bmad

  implicit none

  type (ring_struct), target, intent(out) :: ring
  type (ele_struct), pointer :: ele
  type (ele_digested_struct) :: u_ele
  
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
                     form = 'unformatted', readonly, shared, err = 9000)

  read (d_unit, err = 9100) n_files, version

  if (version < bmad_inc_version$) then
    if (bmad_status%type_out) type '(1x, a, i4, a, i4)',  &
            'READ_DIGESTED_BMAD_FILE: DIGESTED FILE VERSION OUT OF DATE',  &
            version, ' <', bmad_inc_version$
    close (d_unit)
    bmad_status%ok = .false.
    return
  endif

  if (version > bmad_inc_version$) then
    if (bmad_status%type_out) then
      type *, 'READ_DIGESTED_BMAD_FILE: DIGESTED FILE HAS VERSION:',  &
                                                              version
      type *, '     GREATER THAN VERSION OF THIS PROGRAM:',  &
                                                  bmad_inc_version$
      type *, '     WILL NOT USE THE DIGESTED FILE.'
      type *, '     YOU SHOULD RECOMPILE THIS PROGRAM.'
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
      if (bmad_status%type_out) type *, &
                        'READ_DIGESTED_BMAD_FILE: DIGESTED FILE OUT OF DATE.'
      bmad_status%ok = .false.
    endif
    if (i == 1 .and. fname(2) /= ring%input_file_name) then
      if (bmad_status%type_out) type *, &
                        'READ_DIGESTED_BMAD_FILE: MOVED DIGESTED FILE.'
      bmad_status%ok = .false.
    endif
   enddo

! we read (and write) the ring in pieces since it is
! too big to write in one piece

  read (d_unit, err = 9100) ring%parameters

  if (ring%n_ele_max > n_ele_maxx) then
    type *, 'ERROR IN READ_DIGESTED_BMAD_FILE: NUMBER OF ELEMENTS:',  &
                                                            ring%n_ele_max
    type *, '      IS GREATER THAN THE ELEMENT ARRAY SIZE:', n_ele_maxx
    type *, '      YOU NEED TO RECOMPILE'
    call err_exit
  endif

!

  do i = 0, ring%n_ele_max

    read (d_unit, err = 9100) u_ele%digested, ix_w, ix_d, ix_m, ix_t
    u_ele%ele%pointer_init = 0              ! signal that pointers have garbage
    call deallocate_ele_pointers (u_ele%ele)  ! and deallocate 

    ele => ring%ele_(i)
    ele = u_ele%ele

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
    type *, 'ERROR IN READ_DIGESTED_BMAD_FILE: NUMBER OF ELEMENTS:',  &
                                                      ring%n_control_array
    type *, '      IS GREATER THAN THE CONTROL ARRAY SIZE:', n_control_maxx
    type *, '      YOU NEED TO RECOMPILE'
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
    type *, 'READ_DIGESTED_BMAD_FILE: DIGESTED FILE DOES NOT EXIST.'
  endif
  close (d_unit)
  bmad_status%ok = .false.
  return

9100  continue
  if (bmad_status%type_out) then
    type *, 'READ_DIGESTED_BMAD_FILE: ERROR READING DIGESTED FILE.'
  endif
  close (d_unit)
  bmad_status%ok = .false.
  return

end subroutine
