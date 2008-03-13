!+
! Subroutine tao_init_lattice (lat_file_name)
!
! Subroutine to initialize the design lattices.
!
! Input:
!   lat_file_name  -- character(*): file name containing lattice file
!                                             namestructs
!
! Output:
!    %u(:)%design -- Initialized design lattices.
!-

subroutine tao_init_lattice (lat_file_name)

use tao_mod
use tao_input_struct
use ptc_interface_mod

implicit none

type (tao_design_lat_input) design_lattice(0:200)
type (tao_universe_struct), pointer :: u

character(*) lat_file_name
character(200) complete_file_name, file_name
character(40) unique_name_suffix, suffix
character(20) :: r_name = 'tao_init_lattice'
character(16) aperture_limit_on

integer i, j, n, iu, ios, version, taylor_order, ix, key, n_universes

logical custom_init, override, combine_consecutive_elements_of_like_name
logical unified_lattices
logical err, is_set

namelist / tao_design_lattice / design_lattice, taylor_order, &
       combine_consecutive_elements_of_like_name, unique_name_suffix, &
       aperture_limit_on, unified_lattices, n_universes

! Read lattice info

call tao_hook_init_read_lattice_info (lat_file_name, is_set)

if (.not. is_set) then
  call tao_open_file ('TAO_INIT_DIR', lat_file_name, iu, complete_file_name)
  call out_io (s_blank$, r_name, '*Init: Opening File: ' // complete_file_name)
  if (iu == 0) then
    call out_io (s_fatal$, r_name, 'ERROR OPENING PLOTTING FILE. WILL EXIT HERE...')
    call err_exit
  endif

  ! Defaults

  design_lattice%file   = ''
  design_lattice%file2  = ''
  design_lattice%parser = ''
  unified_lattices = .false.
  n_universes = tao_com%n_universes
  taylor_order = 0
  combine_consecutive_elements_of_like_name = .false.
  unique_name_suffix = ''
  aperture_limit_on = ''

  read (iu, nml = tao_design_lattice, iostat = ios)
  if (ios /= 0) then
    call out_io (s_abort$, r_name, 'TAO_DESIGN_LATTICE NAMELIST READ ERROR.')
    rewind (iu)
    do
      read (iu, nml = tao_design_lattice)  ! force printing of error message
    enddo
  endif

  if (taylor_order /= 0) call set_taylor_order (taylor_order)
  tao_com%combine_consecutive_elements_of_like_name = combine_consecutive_elements_of_like_name
  tao_com%unified_lattices = unified_lattices
  tao_com%n_universes = n_universes
  tao_com%aperture_limit_on = aperture_limit_on
  tao_com%unique_name_suffix = unique_name_suffix
endif

!

if (tao_com%unified_lattices) then
  allocate (s%u(tao_com%n_universes))
  tao_com%u_common = 1
  tao_com%u_working = 2
else
  allocate (s%u(tao_com%n_universes))
endif

! Are we using a custom initialization?

custom_init = .false.
call tao_hook_init_lattice (design_lattice, custom_init)
if (custom_init) then
  close (iu)
  return
endif

! See what type of lattice file we have and issue a warning if old style syntax

override = .false.
if (tao_com%init_lat_file(1) /= '') override = .true.

do i = lbound(s%u, 1), ubound(s%u, 1)

  u => s%u(i)
  u%is_on = .true.          ! turn universe on

  if (design_lattice(i)%parser /= '') then
    call out_io (s_error$, r_name, (/ &
        '************************************************************', &
        '***** OLD STYLE "DESIGN_LATTICE()%PARSER" SYNTAX USED! *****', &
        '*****         PLEASE CONVERT TO THE NEW STYLE!         *****', &
        '************************************************************' /) )
  endif

  ! If unified then point back to the common lattice (universe #1) and the working lattice (#2)

  if (tao_com%unified_lattices) then
    if (i == tao_com%u_common) then
      u%common_uni = .true.
    else
      u%common => s%u(tao_com%u_common)
      u%ele => s%u(tao_com%u_common)%ele
      if (i == tao_com%u_working) then   ! Working lattices
        allocate (u%design, u%base, u%model)
        u%design%lat = u%common%design%lat
        u%base%lat   = u%common%base%lat
        u%model%lat  = u%common%model%lat
      else
        u%design => s%u(tao_com%u_working)%design
        u%model  => s%u(tao_com%u_working)%model
        u%base   => s%u(tao_com%u_working)%base
      endif
      cycle
    endif
  endif

  ! Else not unified...
  ! We need to read in the lattice

  if (override) then
    file_name = tao_com%init_lat_file(i)
  else
    file_name = design_lattice(i)%file
    if (design_lattice(i)%parser /= '') file_name = &
              trim(design_lattice(i)%parser) // '::' // trim(file_name)
  endif

  ix = index(file_name(1:10), '::')
  if (ix == 0) then
    design_lattice(i)%file = file_name
    design_lattice(i)%parser = 'bmad'
  else
    design_lattice(i)%file = file_name(ix+2:)
    design_lattice(i)%parser = file_name(1:ix-1)
  endif

  if (design_lattice(i)%file == ' ') design_lattice(i) = design_lattice(1)
  allocate (u%design, u%base, u%model)
  select case (design_lattice(i)%parser)
  case ('bmad')
    call bmad_parser (design_lattice(i)%file, u%design%lat)
  case ('xsif')
    call xsif_parser (design_lattice(i)%file, u%design%lat)
  case ('aml')
    call aml_parser (design_lattice(i)%file, u%design%lat)
  case ('digested')
    call out_io (s_blank$, r_name, &
                "Reading digested BMAD file " // trim(design_lattice(i)%file))
    call read_digested_bmad_file (design_lattice(i)%file, u%design%lat, version)
  case default
    call out_io (s_abort$, r_name, 'PARSER NOT RECOGNIZED: ' // &
                                              design_lattice(i)%parser)
    call err_exit
  end select

  u%design%modes%a%emittance = u%design%lat%a%emit
  u%design%modes%b%emittance = u%design%lat%b%emit

  if (tao_com%combine_consecutive_elements_of_like_name) &
                              call combine_consecutive_elements(u%design%lat)

  if (tao_com%unique_name_suffix /= '') then
    call tao_string_to_element_id (tao_com%unique_name_suffix, key, suffix, err, .true.)
    if (err) call err_exit
    call create_unique_ele_names (u%design%lat, key, suffix)
  endif

  ! Call bmad_parser2 if wanted

  if (design_lattice(i)%file2 /= '') then
    call bmad_parser2 (design_lattice(i)%file2, u%design%lat)
  endif

  ! Aperture limit

  if (tao_com%aperture_limit_on /= '') then
   read (tao_com%aperture_limit_on, *) u%design%lat%param%aperture_limit_on
  endif

  if (u%design%lat%param%lattice_type .eq. circular_lattice$) then
    call out_io (s_blank$, r_name, "RFCavities will be turned off in lattices")
    call set_on_off (rfcavity$, u%design%lat, off$)
  endif

  ! Init model, base, and u%ele

  n = u%design%lat%n_ele_max
  allocate (u%model%orb(0:n), u%model%bunch_params(0:n))
  allocate (u%design%orb(0:n), u%design%bunch_params(0:n))
  allocate (u%base%orb(0:n), u%base%bunch_params(0:n))
  call reallocate_coord (u%design%orb, n) 
  call reallocate_coord (u%model%orb, n)
  call reallocate_coord (u%base%orb, n)
  u%model = u%design
  u%base  = u%design
  allocate (u%ele(0:u%design%lat%n_ele_max))

enddo

close (iu)

end subroutine tao_init_lattice
