!+
! Subroutine write_digested_bmad_file (digested_name, ring, n_files, file_names)
!
! Subroutine to write a digested file. The names of the original files used
! to create the RING structure are also put in the digested file and are used
! by other routines to check if the digested file is out of date.
!
! Modules Needed:
!   use bmad
!
! Input:
!     digested_name -- Character(*): Name for the digested file.
!     ring          -- Ring_struct: Input ring structure.
!     n_files       -- Number of original files
!     file_names(*) -- Character(*), optional: Names of the original 
!                       files used to create the ring structure.
!-

#include "CESR_platform.inc"

subroutine write_digested_bmad_file (digested_name, ring,  &
                                                  n_files, file_names)

  use bmad_struct
  use bmad_interface

  implicit none

  type ele_digested_struct
    union
      map
        integer(rdef) dummy(1000)
      end map
      map
        type (ele_struct) ele
      end map
    end union
  end type
  
  type (ring_struct), target, intent(in) :: ring
  type (ele_digested_struct) :: u_ele
  type (ele_struct), pointer :: ele
  type (taylor_struct), pointer :: tt(:)
  
  integer d_unit, lunget, n_files, i, j, k, ix_w, ix_d, ix_m, ix_t(6)
  integer stat_b(12), stat, ierr, i_write

  character(*) digested_name
  character(*), optional :: file_names(:)
  character(200) fname

  external stat

! Find out minimum size to write

  u_ele%dummy = 0
  u_ele%ele%is_on = .true.
  do i = 1, size(u_ele%dummy)
    if (u_ele%dummy(i) /= 0) exit
    if (i == size(u_ele%dummy)) then
      print *, 'ERROR IN WRITE_DIGESTED_BMAD_FILE: END OF ELE_STRUCT NOT FOUND!'
      call err_exit
    endif
  enddo

  i_write = i + 4
!!  print *, 'WRITE_DIGESTED_BMAD_FILE: Ele_struct size is:', i_write

! write input file names to the digested file

  d_unit = lunget()
  open (unit = d_unit, file = digested_name, form = 'unformatted', err = 9000)

  write (d_unit, err = 9010) n_files, bmad_inc_version$
  write (d_unit, err = 9010) i_write

  do j = 1, n_files
    fname = file_names(j)
    stat_b = 0
#ifdef CESR_UNIX
    ierr = stat(fname, stat_b)
#endif
    write (d_unit) fname, stat_b(10)  ! stat_b(10) = Modification date
  enddo

! write the ring structure to the digested file. We do this in pieces
! since the whole structure is too big to write in 1 statement.
! also: because of the pointers we need to disguise that we are writting ele's
!  by using the digested_ele_struct

  write (d_unit) ring%parameters
  
  do i = 0, ring%n_ele_max
  
    ele => ring%ele_(i)
    tt => ele%taylor
    u_ele%ele = ele
    
    ix_w = 0; ix_d = 0; ix_m = 0; ix_t = 0

    if (ele%pointer_init == has_been_inited$) then
      if (associated(ele%wig_term)) ix_w = size(ele%wig_term)
      if (associated(ele%descrip))  ix_d = 1
      if (associated(ele%a))        ix_m = 1
      if (associated(tt(1)%term))   ix_t = (/ (size(tt(j)%term), j = 1, 6) /)
    endif

    write (d_unit) ix_w, ix_d, ix_m, ix_t, u_ele%dummy(1:i_write)
    do j = 1, ix_w
      write (d_unit) ele%wig_term(j)
    enddo

    if (associated(ele%descrip))  write (d_unit) ele%descrip
    if (associated(ele%a))        write (d_unit) ele%a, ele%b
    
    do j = 1, 6
      do k = 1, ix_t(j)
        write (d_unit) tt(j)%term(k)
      enddo
    enddo

  enddo

! write the control info

  do i = 1, ring%n_control_array
    write (d_unit) ring%control_(i)
  enddo

  do i = 1, ring%n_ic_array
    write (d_unit) ring%ic_(i)
  enddo

  close (d_unit)

  return

! Errors

9000  type *
  type *, 'WRITE_DIGESTED_BMAD_FILE: NOTE: CANNOT OPEN FILE FOR OUTPUT:'
  type *, '    ', trim(digested_name)
  type *, '     [This does not affect program operation]'
  return

9010  type *
  type *, 'WRITE_DIGESTED_BMAD_FILE: NOTE: CANNOT WRITE TO FILE FOR OUTPUT:'
  type *, '    ', trim(digested_name)
  type *, '     [This does not affect program operation]'
  close (d_unit)
  return

end subroutine
