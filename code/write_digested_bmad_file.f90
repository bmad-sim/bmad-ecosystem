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

  type (ring_struct), target, intent(in) :: ring
  type (ele_struct), pointer :: ele
  type (taylor_struct), pointer :: tt(:)
  
  integer d_unit, lunget, n_files, i, j, k
  integer ix_wig, ix_const, ix_d, ix_m, ix_t(6)
  integer stat_b(12), stat, ierr

  character(*) digested_name
  character(*), optional :: file_names(:)
  character(200) fname

  external stat

! write input file names to the digested file

  d_unit = lunget()
  open (unit = d_unit, file = digested_name, form = 'unformatted', err = 9000)

  write (d_unit, err = 9010) n_files, bmad_inc_version$

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

  write (d_unit) &
          ring%name, ring%lattice, ring%input_file_name, ring%title, &
          ring%x, ring%y, ring%z, ring%param, ring%version, ring%n_ele_ring, &
          ring%n_ele_symm, ring%n_ele_use, ring%n_ele_max, &
          ring%n_control_array, ring%n_ic_array, ring%input_taylor_order
  
  do i = 0, ring%n_ele_max
  
    ele => ring%ele_(i)
    tt => ele%taylor
    
    ix_wig = 0; ix_d = 0; ix_m = 0; ix_t = 0; ix_const = 0

    if (ele%pointer_init == has_been_inited$) then
      if (associated(ele%wig_term)) ix_wig = size(ele%wig_term)
      if (associated(ele%const))    ix_const = size(ele%const)
      if (associated(ele%descrip))  ix_d = 1
      if (associated(ele%a))        ix_m = 1
      if (associated(tt(1)%term))   ix_t = (/ (size(tt(j)%term), j = 1, 6) /)
    endif

    write (d_unit) ix_wig, ix_const, ix_d, ix_m, ix_t, &
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

    do j = 1, ix_wig
      write (d_unit) ele%wig_term(j)
    enddo

    if (associated(ele%const))    write (d_unit) ele%const
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

9000  print *
  print *, 'WRITE_DIGESTED_BMAD_FILE: NOTE: CANNOT OPEN FILE FOR OUTPUT:'
  print *, '    ', trim(digested_name)
  print *, '     [This does not affect program operation]'
  return

9010  print *
  print *, 'WRITE_DIGESTED_BMAD_FILE: NOTE: CANNOT WRITE TO FILE FOR OUTPUT:'
  print *, '    ', trim(digested_name)
  print *, '     [This does not affect program operation]'
  close (d_unit)
  return

end subroutine
