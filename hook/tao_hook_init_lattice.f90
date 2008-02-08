!+
! Subroutine tao_hook_init_lattice (lattice_file, custom_init)
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
!   lattice_file(:)  -- tao_design_lat_input: lattice file namestructs
!   custom_init      -- Logical: is this hook used to init lattice?
!
! Output:
!    %u(:)%design -- Initialized design lattices.
!-

subroutine tao_hook_init_lattice (lattice_file, custom_init)

  use tao_mod
  use tao_input_struct

  implicit none

  type (tao_design_lat_input)  lattice_file(:)

  character(40) :: r_name = 'tao_hook_init_lattice'
  integer i, iu

  logical custom_init
  
!

  custom_init = .false.
  
end subroutine 
