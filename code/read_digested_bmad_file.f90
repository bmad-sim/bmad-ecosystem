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
!   digested_name -- Character*(*): Name of the digested file
!
! Output:
!   ring      -- Ring_struct: Output structure
!   version   -- Integer: Version number of RING.
!   status    -- Common block status structure
!     %ok       -- Set .false. if read failure.
!-

#include "CESR_platform.inc"

subroutine read_digested_bmad_file (digested_name, ring, version)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ring_struct), target, intent(inout) :: ring
  type (ele_struct), pointer :: ele
  
  type v69_param_struct
    real(rp) beam_energy        ! beam energy in eV
    real(rp) n_part             ! Particles/bunch (for BeamBeam elements).
    real(rp) charge             ! Charge of a bunch (used by LCavities).
    real(rp) total_length       ! total_length of ring
    real(rp) growth_rate        ! growth rate/turn if not stable
    real(rp) t1_with_RF(6,6)    ! Full 1-turn matrix with RF on.
    real(rp) t1_no_RF(4,4)      ! Full 1-turn matrix with RF off.
    integer particle            ! +1 = positrons, -1 = electrons
    integer iy                  ! Not currently used.
    integer ix_lost             ! If lost at what element?
    integer lattice_type        ! linac_lattice$, etc...
    integer ixx                 ! Integer for general use
    logical stable              ! is closed ring stable?
    logical aperture_limit_on   ! use apertures in tracking?
    logical lost                ! for use in tracking
  end type
  type (v69_param_struct) param69

  integer d_unit, n_files, version, i, j, k, ix
  integer ix_wig, ix_const, ix_r, ix_d, ix_m, ix_t(6)
  integer ix_srf, ix_sr, ix_lrf, ix_lr, i_garbage
  integer stat_b(12), stat, ierr, idate_old

  character(*) digested_name
  character(200) fname(3), input_file_name
  character(200), allocatable :: file_names(:)

  logical found_it, v69, v_now

! init all elements in ring

  call init_ele (ring%ele_init)  ! init pointers
  call deallocate_ring_pointers (ring)

! read the digested file
! versions 66, 67, and 68 can be read even though it is not the current version

  d_unit = lunget()
  bmad_status%ok = .true.
  ring%n_ele_ring = 0

  open (unit = d_unit, file = digested_name, status = 'old',  &
                     form = 'unformatted', action = 'READ', err = 9000)

  read (d_unit, err = 9100) n_files, version

  v69 = (version == 69)
  v_now = (version == bmad_inc_version$)

  if (version < bmad_inc_version$) then
!    if (bmad_status%type_out) print '(1x, a, i4, a, i4)',  &
    if (bmad_status%type_out) print *,  &
           'READ_DIGESTED_BMAD_FILE: DIGESTED FILE VERSION OUT OF DATE',  &
            version, ' <', bmad_inc_version$
    if (v69) then 
      allocate (file_names(n_files))
      bmad_status%ok = .false.
    else
      close (d_unit)
      bmad_status%ok = .false.
      return
    endif
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

  call simplify_path(ring%input_file_name, input_file_name)

  do i = 1, n_files
    read (d_unit, err = 9100) fname(1), idate_old
    call simplify_path (fname(1), fname(1))
    if (v69) file_names(i) = fname(1)  ! fake out
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
    call simplify_path (fname(3), fname(3))
    if (.not. found_it .or. fname(1) /= fname(3) .or. &
                                             stat_b(10) /= idate_old) then
      if (bmad_status%type_out .and. bmad_status%ok) print *, &
              'READ_DIGESTED_BMAD_FILE: NOTE: DIGESTED FILE OUT OF DATE.'
      bmad_status%ok = .false.
    endif
    if (i == 1 .and. fname(2) /= input_file_name) then
      if (bmad_status%type_out .and. bmad_status%ok) print *, &
                    'READ_DIGESTED_BMAD_FILE: NOTE: MOVED DIGESTED FILE.'
      bmad_status%ok = .false.
    endif
   enddo

! we read (and write) the ring in pieces since it is
! too big to write in one piece

  if (v_now) then
    read (d_unit, err = 9100)  &   
          ring%name, ring%lattice, ring%input_file_name, ring%title, &
          ring%x, ring%y, ring%z, ring%param, ring%version, ring%n_ele_ring, &
          ring%n_ele_use, ring%n_ele_max, &
          ring%n_control_max, ring%n_ic_max, ring%input_taylor_order
  elseif (v69) then
    read (d_unit, err = 9100)  &   
          ring%name, ring%lattice, ring%input_file_name, ring%title, &
          ring%x, ring%y, ring%z, param69, ring%version, ring%n_ele_ring, &
          i_garbage, ring%n_ele_use, ring%n_ele_max, &
          ring%n_control_max, ring%n_ic_max, ring%input_taylor_order
    ring%param%beam_energy  = param69%beam_energy
    ring%param%particle     = param69%particle
    ring%param%lattice_type = param69%lattice_type
  else
    print *, 'ERROR IN READ_DIGESTED_BMAD_FILE: INTERNAL ERROR: RING.'
    print *, '      PLEASE GET EXPERT HELP!'
    call err_exit
  endif

  call allocate_ring_ele_(ring, ring%n_ele_max+100)
  allocate (ring%control_(ring%n_control_max+100))
  allocate (ring%ic_(ring%n_ic_max+100))

!

  do i = 0, ring%n_ele_max

    ele => ring%ele_(i)
    if (v69 .or. v_now) then
      read (d_unit, err = 9100) ix_wig, ix_const, ix_r, ix_d, ix_m, ix_t, &
                              ix_srf, ix_sr, ix_lrf, ix_lr, &
            ele%name, ele%type, ele%alias, ele%attribute_name, ele%x, &
            ele%y, ele%z, ele%value, ele%gen0, ele%vec0, ele%mat6, &
            ele%c_mat, ele%gamma_c, ele%s, ele%key, ele%position, &
            ele%is_on, ele%sub_key, ele%control_type, ele%ix_value, &
            ele%n_slave, ele%ix1_slave, ele%ix2_slave, ele%n_lord, &
            ele%ic1_lord, ele%ic2_lord, ele%ix_pointer, ele%ixx, &
            ele%iyy, ele%mat6_calc_method, ele%tracking_method, &
            ele%num_steps, ele%integration_ord, ele%ptc_kind, &
            ele%taylor_order, ele%symplectify, ele%mode_flip, &
            ele%multipoles_on, ele%exact_rad_int_calc, ele%Field_master, &
            ele%logic, ele%internal_logic, ele%field_calc
    else
      print *, 'ERROR IN READ_DIGESTED_BMAD_FILE: INTERNAL ERROR: ELE.'
      print *, '      PLEASE GET EXPERT HELP!'
      call err_exit
    endif

    if (ix_wig /= 0) then
      allocate (ele%wig_term(ix_wig))
      do j = 1, ix_wig
        read (d_unit) ele%wig_term(j)
      enddo
    endif

    if (ix_const /= 0) then
      allocate (ele%const(ix_const))
      read (d_unit) ele%const
    endif

!    if (ix_r /= 0) then
!      allocate (ele%r(ix_r))
!      read (d_unit) ele%r
!    endif

    if (ix_d /= 0) then
      allocate (ele%descrip)
      read (d_unit) ele%descrip
    endif

    if (ix_m /= 0) then
      allocate (ele%a(0:n_pole_maxx), ele%b(0:n_pole_maxx))
      read (d_unit) ele%a, ele%b
    endif
    
    do j = 1, 6
      if (ix_t(j) == 0) cycle
      read (d_unit) ele%taylor(j)%ref
      allocate (ele%taylor(j)%term(ix_t(j)))
      do k = 1, ix_t(j)
        read (d_unit) ele%taylor(j)%term(k)
      enddo
    enddo

    if (ix_srf /= 0) then
      allocate (ele%wake%sr_file)
      read (d_unit) ele%wake%sr_file
    endif

    if (ix_sr /= 0) then
      allocate (ele%wake%sr(0:ix_sr-1))
      read (d_unit) ele%wake%sr
    endif

    if (ix_lrf /= 0) then
      allocate (ele%wake%lr_file)
      read (d_unit) ele%wake%lr_file
    endif

    if (ix_lr /= 0) then
      allocate (ele%wake%lr(0:ix_lr-1))
      read (d_unit) ele%wake%lr
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
