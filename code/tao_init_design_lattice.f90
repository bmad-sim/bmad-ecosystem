
!+
! Subroutine tao_init_design_lattice (tao_design_lattice_file)
!
! Subroutine to initialize the design lattices.
!
! Input:
!   tao_design_lattice_file  -- character(*): file name containing lattice file
!                                             namestructs
!
! Output:
!    %u(:)%design -- Initialized design lattices.
!-

subroutine tao_init_design_lattice (tao_design_lattice_file)

  use tao_mod
  use tao_input_struct
  use ptc_interface_mod

  implicit none

  type (tao_design_lat_input)  design_lattice(100)
  type (tao_universe_struct), pointer :: u

  character(*) tao_design_lattice_file
  character(200) complete_file_name, file_name
  character(40) :: r_name = 'tao_init_design_lattice'
  integer i, j, iu, ios, version, taylor_order, ix

  logical custom_init, override

  namelist / tao_design_lattice / design_lattice, taylor_order

!

  call tao_open_file ('TAO_INIT_DIR', tao_design_lattice_file, iu, complete_file_name)
  call out_io (s_blank$, r_name, '*Init: Opening File: ' // complete_file_name)

  design_lattice%file = ''
  design_lattice%parser = ''
  taylor_order = 0

  read (iu, nml = tao_design_lattice, iostat = ios)
  if (ios /= 0) then
    call out_io (s_abort$, r_name, 'TAO_DESIGN_LATTICE NAMELIST READ ERROR.')
    rewind (iu)
    do
      read (iu, nml = tao_design_lattice)  ! force printing of error message
    enddo
  endif
  close (iu)

  if (taylor_order /= 0) call set_taylor_order (taylor_order)

  ! are we using a custom initialization?

  custom_init = .false.
  call tao_hook_init_design_lattice (design_lattice, custom_init)
  if (custom_init) return

  ! TAO does its own bookkeeping

  bmad_com%auto_bookkeeper = .false.

  ! See what type of lattice file we have and issue a warning if old style syntax

  override = .false.
  if (tao_com%init_lat_file(1) /= '') override = .true.

  do i = 1, size(s%u)
    if (design_lattice(i)%parser /= '') then
      call out_io (s_error$, r_name, &
          '***************************************************************************************', &
          '***** OLD STYLE "DESIGN_LATTICE%PARSER" SYNTAX USED. PLEASE CONVERT TO NEW STYLE! *****', &
          '***************************************************************************************')
    endif

    if (override) then
      file_name = tao_com%init_lat_file(i)
    else
      file_name = design_lattice(i)%file
      if (design_lattice(i)%parser /= '') file_name = trim(design_lattice(i)%parser) // '::' // trim(file_name)
    endif

    ix = index(file_name(1:10), '::')
    if (ix == 0) then
      design_lattice(i)%file = file_name
      design_lattice(i)%parser = 'bmad'
    else
      design_lattice(i)%file = file_name(ix+2:)
      design_lattice(i)%parser = file_name(1:ix-1)
    endif

    u => s%u(i)

    if (design_lattice(i)%file == ' ') design_lattice(i) = design_lattice(1)
    select case (design_lattice(i)%parser)
    case ('bmad')
      call bmad_parser (design_lattice(i)%file, u%design%lat)
      u%design%modes%a%emittance = u%design%lat%a%emit
      u%design%modes%b%emittance = u%design%lat%b%emit
    case ('xsif')
      call xsif_parser (design_lattice(i)%file, u%design%lat)
    case ('digested')
      call out_io (s_blank$, r_name, &
                  "Reading digested BMAD file " // trim(design_lattice(i)%file))
      call read_digested_bmad_file (design_lattice(i)%file, u%design%lat, version)
    case default
      call out_io (s_abort$, r_name, 'PARSER NOT RECOGNIZED: ' // &
                                                design_lattice(i)%parser)
      call err_exit
    end select

    ! Initialize calibration array
    ! This must be performed or tao_read_bpm and tao_do_wire_scan will crash.
    ! r(1,:) is for bpm and steering callibration
    ! r(2,:) is for saving ele parameters

    do j = 1, u%design%lat%n_ele_max
        allocate(u%design%lat%ele(j)%r(2,5))
        u%design%lat%ele(j)%r = 0.0
    enddo

    if (u%design%lat%param%lattice_type .eq. circular_lattice$) then
      call out_io (s_blank$, r_name, &
                  "RFCavities will be turned off in lattices")
      call set_on_off (rfcavity$, u%design%lat, off$)
    endif

    ! Init beam save

    allocate (u%beam_at_element(0:u%design%lat%n_ele_track))

  enddo

end subroutine tao_init_design_lattice
