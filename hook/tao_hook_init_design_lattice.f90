!+
! Subroutine tao_hook_init_design_lattice (design_lattice_file, custom_init)
!
! Subroutine to initialize custom design lattices. If this routine is used to do
! a custom lattice initialization then set custom_init to TRUE or else, the
! standard intialization will overwrite whatever you do here!
!
! For reference:
!  type tao_design_lat_input
!    character(200) file
!    character(16) :: parser = 'bmad'
!  end type
!
! Input:
!   design_lattice_file(:)  -- tao_design_lat_input: lattice file namestructs
!   custom_init             -- Logical: is this hook used to init lattice?
!
! Output:
!    %u(:)%design -- Initialized design lattices.
!-

subroutine tao_hook_init_design_lattice (design_lattice_file, custom_init)

  use tao_mod
  use tao_input_struct

  implicit none

  type (tao_design_lat_input)  design_lattice_file(:)

  character(40) :: r_name = 'tao_hook_init_design_lattice'
  integer i, iu

  logical custom_init
  
!

  custom_init = .false.
  
  return
  
!--------------------------------------------------------------------  
! For reference, this is what's called in tao_init_design_lattice:

! do i = 1, size(s%u)
!   select case (design_lattice_file(i)%parser)
!   case ('bmad')
!     call bmad_parser (design_lattice_file(i)%file, s%u(i)%design)
!   case ('xsif')
!     call xsif_parser (design_lattice_file(i)%file, s%u(i)%design)
!   case default
!     call out_io (s_abort$, r_name, 'PARSER NOT RECOGNIZED: ' // &
!                                               design_lattice_file(i)%parser)
!     call err_exit
!   end select
! enddo

end subroutine 
