!+
! Subroutine tao_init (init_file)
!
! Subroutine to initialize the tao structures.
!
! Input:
!   init_file -- Character(*): Initialization file.
! Output:
!-

subroutine tao_init (init_file)

  use tao_mod

  implicit none

  character(*) init_file
  character(200) lattice_file, plot_file, data_and_var_file, file_name
  character(200) single_mode_file
  character(20) :: r_name = 'tao_init'
  integer i, n_universes, iu

  namelist / tao_start / lattice_file, &
               data_and_var_file, plot_file, n_universes

! Find namelist files

  lattice_file = ' '      ! set default
  plot_file = ' '         ! set default
  data_and_var_file = ' ' ! set default
  single_mode_file = ' '
  n_universes = 1         ! set default
  call tao_open_file ('TAO_INIT_DIR', init_file, iu, file_name)
  read (iu, nml = tao_start)
  close (iu)

  allocate (s%u(n_universes))

  if (lattice_file == ' ')      lattice_file      = init_file
  if (plot_file == ' ')         plot_file         = init_file
  if (data_and_var_file == ' ') data_and_var_file = init_file
  if (single_mode_file == ' ')  single_mode_file  = init_file

! Init

  call tao_hook_init_design_lattice (lattice_file) 
  call tao_init_global_and_universes (data_and_var_file)
  call tao_init_single_mode (single_mode_file)
  call tao_hook_init ()

  call tao_lattice_calc ()
  do i = 1, size(s%u)
    if (associated (s%u(i)%data)) then
      s%u(i)%data%design_value = s%u(i)%data%model_value
      s%u(i)%data%base_value = s%u(i)%data%model_value
    endif
  enddo

  call tao_set_data_useit_opt ()
  call tao_set_var_useit_opt ()

  call tao_lattice_calc()      ! calculate Twiss parameters, closed orbit
  if (s%global%plot_on) call tao_init_plotting (plot_file)
  call tao_plot_data_setup ()  ! transfer data to the plotting structures
  call tao_plot_out ()         ! Update the plotting window

end subroutine tao_init
