
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

  type (tao_design_lat_input)  design_lattice_file(20)

  character(*) tao_design_lattice_file
  character(200) complete_file_name
  character(40) :: r_name = 'tao_init_design_lattice'
  integer i, iu, ios, version, taylor_order

  logical custom_init

  namelist / tao_design_lattice / design_lattice_file, taylor_order

!

  call tao_open_file ('TAO_INIT_DIR', tao_design_lattice_file, iu, complete_file_name)

  design_lattice_file%file = ' '
  design_lattice_file%parser = 'bmad'
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
  call tao_hook_init_design_lattice (design_lattice_file, custom_init)
  if (custom_init) return

  
  do i = 1, size(s%u)
    select case (design_lattice_file(i)%parser)
    case ('bmad')
      call bmad_parser (design_lattice_file(i)%file, s%u(i)%design)
    case ('xsif')
      call xsif_parser (design_lattice_file(i)%file, s%u(i)%design)
    case ('digested')
      call out_io (s_blank$, r_name, &
                  "Reading digested BMAD file " // trim(design_lattice_file(i)%file))
      call read_digested_bmad_file (design_lattice_file(i)%file, s%u(i)%design, version)
    case default
      call out_io (s_abort$, r_name, 'PARSER NOT RECOGNIZED: ' // &
                                                design_lattice_file(i)%parser)
      call err_exit
    end select
  enddo

  ! TAO does its own bookkeeping
  bmad_com%auto_bookkeeper = .false.

  ! turn off rfcavities in rings
  do i = 1, size(s%u)
    if (s%u(i)%design%param%lattice_type .eq. circular_lattice$) then
      call out_io (s_blank$, r_name, &
                  "RFCavities will be turned off in rings")
      call set_on_off (rfcavity$, s%u(i)%design, off$)
    endif
  enddo

end subroutine tao_init_design_lattice
