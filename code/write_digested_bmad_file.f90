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
!     n_files       -- Integer, optional: Number of original files
!     file_names(*) -- Character(*), optional: Names of the original 
!                       files used to create the ring structure.
!-

#include "CESR_platform.inc"

subroutine write_digested_bmad_file (digested_name, ring,  &
                                                  n_files, file_names)

  use bmad_struct
  use bmad_interface, except => write_digested_bmad_file
  use equality_mod, only: operator(==)

  implicit none

  type (ring_struct), target, intent(in) :: ring
  type (ele_struct), pointer :: ele
  type (taylor_struct), pointer :: tt(:)
  type (wake_struct), pointer :: wake
  
  integer, intent(in), optional :: n_files
  integer d_unit, i, j, k, n_file
  integer ix_wig, ix_const, ix_r(4), ix_d, ix_m, ix_t(6)
  integer stat_b(12), stat, n_wake
  integer ix_sr1, ix_sr2_long, ix_sr2_trans, ix_lr, ierr
  integer, automatic :: ix_wake(ring%n_ele_max)

  character(*) digested_name
  character(*), optional :: file_names(:)
  character(200) fname
  character(32) :: r_name = 'write_digested_bmad_file'

  logical write_wake

  external stat

! write input file names to the digested file

  n_file = 0
  if (present(n_files)) n_file = n_files

  d_unit = lunget()

  call fullfilename (digested_name, fname)
  open (unit = d_unit, file = fname, form = 'unformatted', err = 9000)

  write (d_unit, err = 9010) n_file, bmad_inc_version$

  do j = 1, n_file
    fname = file_names(j)
    stat_b = 0
#ifndef CESR_VMS 
    ierr = stat(fname, stat_b)
#endif
    write (d_unit) fname, stat_b(10)  ! stat_b(10) = Modification date
  enddo

! write the ring structure to the digested file. We do this in pieces
! since the whole structure is too big to write in 1 statement.

  write (d_unit) &
          ring%name, ring%lattice, ring%input_file_name, ring%title, &
          ring%x, ring%y, ring%z, ring%param, ring%version, ring%n_ele_use, &
          ring%n_ele_use, ring%n_ele_max, &
          ring%n_control_max, ring%n_ic_max, ring%input_taylor_order
  
  n_wake = 0  ! number of wakes written to the digested file.

  do i = 0, ring%n_ele_max
  
    ele => ring%ele_(i)
    tt => ele%taylor
    
    ix_wig = 0; ix_d = 0; ix_m = 0; ix_t = 0; ix_const = 0; ix_r = 0
    ix_sr1 = 0; ix_sr2_long = 0; ix_sr2_trans = 0; ix_lr = 0

    if (associated(ele%wig_term)) ix_wig = size(ele%wig_term)
    if (associated(ele%const))    ix_const = size(ele%const)
    if (associated(ele%r))        ix_r = (/ lbound(ele%r), ubound(ele%r) /)
    if (associated(ele%descrip))  ix_d = 1
    if (associated(ele%a))        ix_m = 1
    if (associated(tt(1)%term))   ix_t = (/ (size(tt(j)%term), j = 1, 6) /)

    ! Since some large lattices with a large number of wakes can take a lot of time writing 
    ! the wake info we only write a wake when needed.
    ! The idea is that ix_lr serves as a pointer to a previously written wake.

    write_wake = .true.
    if (associated(ele%wake)) then
      do j = 1, n_wake
        if (.not. ring%ele_(ix_wake(j))%wake == ele%wake) cycle
        write_wake = .false.
        ix_lr = -ix_wake(j)        
      enddo

      if (write_wake) then
        if (associated(ele%wake%sr1))       ix_sr1       = size(ele%wake%sr1)
        if (associated(ele%wake%sr2_long))  ix_sr2_long  = size(ele%wake%sr2_long)
        if (associated(ele%wake%sr2_trans)) ix_sr2_trans = size(ele%wake%sr2_trans)
        if (associated(ele%wake%lr))        ix_lr        = size(ele%wake%lr)
        n_wake = n_wake + 1
        ix_wake(n_wake) = i
      endif
    endif

    ! Now write the element info

    write (d_unit) ix_wig, ix_const, ix_r, ix_d, ix_m, ix_t, &
            ix_sr1, ix_sr2_long, ix_sr2_trans, ix_lr, &
            ele%name, ele%type, ele%alias, ele%attribute_name, ele%x, &
            ele%y, ele%z, ele%value, ele%gen0, ele%vec0, ele%mat6, &
            ele%c_mat, ele%gamma_c, ele%s, ele%key, ele%floor, &
            ele%is_on, ele%sub_key, ele%control_type, ele%ix_value, &
            ele%n_slave, ele%ix1_slave, ele%ix2_slave, ele%n_lord, &
            ele%ic1_lord, ele%ic2_lord, ele%ix_pointer, ele%ixx, &
            ele%ix_ele, ele%mat6_calc_method, ele%tracking_method, &
            ele%num_steps, ele%integrator_order, ele%ptc_kind, &
            ele%taylor_order, ele%symplectify, ele%mode_flip, &
            ele%multipoles_on, ele%exact_rad_int_calc, ele%Field_master, &
            ele%logic, ele%internal_logic, ele%field_calc, ele%aperture_at, &
            ele%on_an_i_beam

    do j = 1, ix_wig
      write (d_unit) ele%wig_term(j)
    enddo

    if (associated(ele%const))    write (d_unit) ele%const
    if (associated(ele%r))        write (d_unit) ele%r
    if (associated(ele%descrip))  write (d_unit) ele%descrip
    if (associated(ele%a))        write (d_unit) ele%a, ele%b
    
    do j = 1, 6
      if (ix_t(j) == 0) cycle
      write (d_unit) tt(j)%ref
      do k = 1, ix_t(j)
        write (d_unit) tt(j)%term(k)
      enddo
    enddo

    if (associated(ele%wake) .and. write_wake) then
      write (d_unit) ele%wake%sr_file
      write (d_unit) ele%wake%sr1
      write (d_unit) ele%wake%sr2_long
      write (d_unit) ele%wake%sr2_trans
      write (d_unit) ele%wake%lr_file
      write (d_unit) ele%wake%lr
      write (d_unit) ele%wake%z_sr2_max
    endif

  enddo

! write the control info

  do i = 1, ring%n_control_max
    write (d_unit) ring%control_(i)
  enddo

  do i = 1, ring%n_ic_max
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
