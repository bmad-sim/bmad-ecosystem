!+
! Subroutine read_digested_bmad_file (digested_name, ring, version)
!
! Subroutine to read in a digested file. The subroutine will check that
! the version of the digested file is up to date and that the digested file
! is current with respect to the original BMAD files that were used. [See
! write_digested_bmad_file]
!
! Note: This subroutine also reads in the common structures for BMAD_PARSER2
!
! Modules Needed:
!   use bmad
!
! Input:
!   digested_name -- Character(*): Name of the digested file
!
! Output:
!   ring        -- Ring_struct: Output structure
!   version     -- Integer: Version number of RING.
!   bmad_status -- Bmad_status_struct: Common block status structure.
!     %ok           -- Set .false. if read failure.
!-

#include "CESR_platform.inc"

subroutine read_digested_bmad_file (digested_name, ring, version)

  use bmad_struct
  use bmad_interface, except => read_digested_bmad_file

  implicit none

  type (ring_struct), target, intent(inout) :: ring
  type (ele_struct), pointer :: ele
  
  integer d_unit, n_files, version, i, j, k, ix, dum1, dum2
  integer ix_wig, ix_const, ix_r(4), ix_d, ix_m, ix_t(6)
  integer ix_sr1, ix_sr2_long, ix_sr2_trans, ix_lr, ierr, stat
  integer stat_b(12), idate_old

  character(*) digested_name
  character(200) fname(3), input_file_name
  character(200), allocatable :: file_names(:)
  character(25) :: r_name = 'read_digested_bmad_file'

  logical found_it, v75, v76, v_old, v_now

  type old_param_struct
    real(rp) garbage2           ! Saved for future use.
    real(rp) n_part             ! Particles/bunch (for BeamBeam elements).
    real(rp) garbage            ! Saved for future use.
    real(rp) total_length       ! total_length of ring
    real(rp) growth_rate        ! growth rate/turn if not stable
    real(rp) t1_with_RF(6,6)    ! Full 1-turn matrix with RF on.
    real(rp) t1_no_RF(6,6)      ! Full 1-turn matrix with RF off.
    integer particle            ! positron$, electron$, etc.
    integer ix_lost             ! Index of element particle was lost at.
    integer end_lost_at         ! entrance_end$ or exit_end$
    integer lattice_type        ! linear_lattice$, etc...
    integer ixx                 ! Integer for general use
    integer ran_seed            ! Random number seed used.
    logical stable              ! is closed ring stable?
    logical aperture_limit_on   ! use apertures in tracking?
    logical lost                ! for use in tracking
  end type

  type (old_param_struct) old_param

! init all elements in ring

  call init_ele (ring%ele_init)  ! init pointers
  call deallocate_ring_pointers (ring)

! read the digested file.
! some old versions can be read even though they are not the current version

  d_unit = lunget()
  bmad_status%ok = .true.
  ring%n_ele_use = 0

  call fullfilename (digested_name, input_file_name)
  open (unit = d_unit, file = input_file_name, status = 'old',  &
                     form = 'unformatted', action = 'READ', err = 9000)

  read (d_unit, err = 9100) n_files, version

  v75 = (version == 75)
  v76 = (version == 76)
  v_old = v75 .or. v76
  v_now = (version == bmad_inc_version$)  ! v77

  if (version < bmad_inc_version$) then
    if (bmad_status%type_out) call out_io (s_warn$, r_name, &
           (/ 'DIGESTED FILE VERSION OUT OF DATE \i4\ < \i4\ ' /),  &
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

  call simplify_path(ring%input_file_name, input_file_name)

  do i = 1, n_files
    read (d_unit, err = 9100) fname(1), idate_old
    call simplify_path (fname(1), fname(1))
    if (v_old) file_names(i) = fname(1)  ! fake out
    ix = index(fname(1), ';')
    stat_b = 0
    if (ix > 0) then    ! has VMS version number
      fname(2) = fname(1)(:ix-1)
    else
      fname(2) = fname(1)
#ifndef CESR_VMS 
      ierr = stat(fname(2), stat_b)
#endif
    endif
    inquire (file = fname(2), exist = found_it, name = fname(3))
    call simplify_path (fname(3), fname(3))
    if (.not. found_it .or. fname(1) /= fname(3) .or. &
                                             stat_b(10) /= idate_old) then
      if (bmad_status%type_out .and. bmad_status%ok) call out_io(s_warn$ ,r_name, &
              'NOTE: DIGESTED FILE OUT OF DATE.')

      bmad_status%ok = .false.
    endif
    if (i == 1 .and. fname(2) /= input_file_name) then
      if (bmad_status%type_out .and. bmad_status%ok) call out_io(s_warn$, r_name, &
                    ' NOTE: MOVED DIGESTED FILE.')

      bmad_status%ok = .false.
    endif
   enddo

! we read (and write) the ring in pieces since it is
! too big to write in one piece

  if (v75) then
    read (d_unit, err = 9100)  &   
          ring%name, ring%lattice, ring%input_file_name, ring%title, &
          ring%x, ring%y, ring%z, old_param, ring%version, ring%n_ele_use, &
          ring%n_ele_ring, ring%n_ele_max, &
          ring%n_control_max, ring%n_ic_max, ring%input_taylor_order
    ring%param%n_part            = old_param%n_part
    ring%param%total_length      = old_param%total_length
    ring%param%growth_rate       = old_param%growth_rate
    ring%param%t1_with_RF(6,6)   = old_param%t1_with_RF(6,6)
    ring%param%t1_no_RF(6,6)     = old_param%t1_no_RF(6,6)
    ring%param%particle          = old_param%particle
    ring%param%ix_lost           = old_param%ix_lost
    ring%param%end_lost_at       = old_param%end_lost_at
    ring%param%lattice_type      = old_param%lattice_type
    ring%param%ixx               = old_param%ixx
    ring%param%stable            = old_param%stable
    ring%param%aperture_limit_on = old_param%aperture_limit_on
    ring%param%lost              = old_param%lost
    
  else
    read (d_unit, err = 9100)  &   
          ring%name, ring%lattice, ring%input_file_name, ring%title, &
          ring%x, ring%y, ring%z, ring%param, ring%version, ring%n_ele_use, &
          ring%n_ele_ring, ring%n_ele_max, &
          ring%n_control_max, ring%n_ic_max, ring%input_taylor_order
  endif

  call allocate_ring_ele_(ring, ring%n_ele_max+100)
  allocate (ring%control_(ring%n_control_max+100))
  allocate (ring%ic_(ring%n_ic_max+100))

!

  do i = 0, ring%n_ele_max

    ele => ring%ele_(i)
    if (v75) then
      read (d_unit, err = 9100) ix_wig, ix_const, ix_r, ix_d, ix_m, ix_t, &
                              dum1, ix_sr1, dum2, ix_lr, &
            ele%name, ele%type, ele%alias, ele%attribute_name, ele%x, &
            ele%y, ele%z, ele%value(1:41), ele%gen0, ele%vec0, ele%mat6, &
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
      ix_sr2_long = 0
      ix_sr2_trans = 0

    elseif (v76) then
      read (d_unit, err = 9100) ix_wig, ix_const, ix_r, ix_d, ix_m, ix_t, &
            ix_sr1, ix_sr2_long, ix_sr2_trans, ix_lr, &
            ele%name, ele%type, ele%alias, ele%attribute_name, ele%x, &
            ele%y, ele%z, ele%value(1:41), ele%gen0, ele%vec0, ele%mat6, &
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
      ele%value(38:44) = ele%value(35:41)
      ele%value(36:37) = ele%value(34)
      ele%value(25:35) = ele%value(23:33)

    elseif (v_now) then
      read (d_unit, err = 9100) ix_wig, ix_const, ix_r, ix_d, ix_m, ix_t, &
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
    endif

    if (ix_wig /= 0) then
      allocate (ele%wig_term(ix_wig))
      do j = 1, ix_wig
        read (d_unit, err = 9200) ele%wig_term(j)
      enddo
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
      allocate (ele%a(0:n_pole_maxx), ele%b(0:n_pole_maxx))
      read (d_unit, err = 9600) ele%a, ele%b
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

    if (ix_sr1 /= 0 .or. ix_sr2_long /= 0 .or. ix_sr2_trans /= 0 .or. ix_lr /= 0) then
      if (ix_lr < 0) then
        call transfer_wake (ring%ele_(abs(ix_lr))%wake, ele%wake)

      else
        call init_wake (ele%wake, ix_sr1, ix_sr2_long, ix_sr2_trans, ix_lr)
        if (v75) then
          if (ix_sr1 /= 0) then
            allocate (ele%wake%sr1(0:ix_sr1-1))
            read (d_unit, err = 9800) ele%wake%sr_file
            read (d_unit, err = 9810) ele%wake%sr1
          endif
          if (ix_lr /= 0) then
            allocate (ele%wake%lr(0:ix_lr-1))
            read (d_unit, err = 9820) ele%wake%lr_file
            read (d_unit, err = 9830) ele%wake%lr
          endif
        elseif (v76 .or. v_now) then
          read (d_unit, err = 9800) ele%wake%sr_file
          read (d_unit, err = 9810) ele%wake%sr1
          read (d_unit, err = 9840) ele%wake%sr2_long
          read (d_unit, err = 9850) ele%wake%sr2_trans
          read (d_unit, err = 9820) ele%wake%lr_file
          read (d_unit, err = 9830) ele%wake%lr
          read (d_unit, err = 9860) ele%wake%z_sr2_max
        endif
      endif
    endif

  enddo

!

  do i = 1, ring%n_control_max
    read (d_unit, err = 9100) ring%control_(i)
  enddo

  do i = 1, ring%n_ic_max
    read (d_unit, err = 9100) ring%ic_(i)
  enddo

  close (d_unit)

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
          'ERROR READING WAKE%SR1 FOR ELEMENT: ' // ele%name)
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
          'ERROR READING WAKE%SR2_LONG FOR ELEMENT: ' // ele%name)
  endif
  close (d_unit)
  bmad_status%ok = .false.
  return

9850  continue
  if (bmad_status%type_out) then
     call out_io(s_error$, r_name, 'ERROR READING DIGESTED FILE.', &
          'ERROR READING WAKE%SR2_TRANS FOR ELEMENT: ' // ele%name)
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
