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
  use tao_lattice_calc_mod

  implicit none

  character(*) init_file
  character(200) lattice_file, plot_file, data_and_var_file, file_name
  character(200) single_mode_file, startup_file
  character(n_universe_maxx) :: r_name = 'tao_init'
  integer i, n_universes, iu

  namelist / tao_start / lattice_file, startup_file, &
               data_and_var_file, plot_file, n_universes

! Find namelist files

  lattice_file = ' '      ! set default
  plot_file = ' '         ! set default
  data_and_var_file = ' ' ! set default
  single_mode_file = ' '
  startup_file = 'tao.startup'
  n_universes = 1         ! set default
  call tao_open_file ('TAO_INIT_DIR', init_file, iu, file_name)
  read (iu, nml = tao_start)
  close (iu)

  if (n_universes .gt. n_universe_maxx) then
    call out_io (s_abort$, r_name, "Too many universe! Maximum number of &
                    & universes: \i\ ", n_universe_maxx)
    call err_exit
  endif
  
  if (associated(s%u)) call deallocate_everything ()
  allocate (s%u(n_universes))

  if (lattice_file == ' ')      lattice_file      = init_file
  if (plot_file == ' ')         plot_file         = init_file
  if (data_and_var_file == ' ') data_and_var_file = init_file
  if (single_mode_file == ' ')  single_mode_file  = init_file

! Init

  call tao_init_design_lattice (lattice_file) 
  call tao_init_global_and_universes (data_and_var_file)
  call tao_init_single_mode (single_mode_file)
  call tao_hook_init ()

  ! set up design and base lattices
  do i = 1, size(s%u)
    s%u(i)%model = s%u(i)%design; s%u(i)%model_orb = s%u(i)%design_orb
  enddo
  s%global%lattice_recalc = .true.
  call tao_lattice_calc (.true.) ! .true. => init design lattice
  do i = 1, size(s%u)
    s%u(i)%design = s%u(i)%model; s%u(i)%design_orb = s%u(i)%model_orb
    s%u(i)%base  = s%u(i)%design; s%u(i)%base_orb  = s%u(i)%design_orb
    s%u(i)%data%design_value = s%u(i)%data%model_value
    s%u(i)%data%base_value = s%u(i)%data%model_value
  enddo

  call tao_init_plotting (plot_file)
  call tao_plot_data_setup ()  ! transfer data to the plotting structures
  call tao_plot_out ()         ! Update the plotting window

! Look for a startup file

  call tao_open_file ('TAO_INIT_DIR', startup_file, iu, file_name)
  if (iu /= 0) then
    call out_io (s_blank$, r_name, 'Using startup file: ' // file_name)
    call tao_call_cmd (file_name)
  endif

contains
!------------------------------------------------------------------------------
! every pointer and allocatable needs to be deallocated now before the universe
! is reallocated.

subroutine deallocate_everything ()

implicit none

integer i, j, k, istat

! Variables  
  if (associated (s%v1_var)) then
    do i = 1, size(s%v1_var)
      nullify(s%v1_var(i)%v)
    enddo
    deallocate(s%v1_var, stat=istat)
  endif
  
  if (associated (s%var)) then
    do i = lbound(s%var,1), ubound(s%var,1)
      do j = 1, size(s%var(i)%this)
        nullify(s%var(i)%this(j)%model_ptr)
        nullify(s%var(i)%this(j)%base_ptr)
      enddo
      nullify(s%var(i)%v1)
      deallocate(s%var(i)%this, stat=istat)
    enddo
    deallocate(s%var, stat=istat)
  endif
 
! Keytable 
  if (associated (s%key)) deallocate(s%key, stat=istat)

! plotting  
  do i = 1, size(s%template_plot)
    if (.not. associated (s%template_plot(i)%graph)) cycle
    do j = 1, size(s%template_plot(i)%graph)
      do k = 1, size(s%template_plot(i)%graph(j)%curve)
        deallocate(s%template_plot(i)%graph(j)%curve(k)%x_line, stat=istat)
        deallocate(s%template_plot(i)%graph(j)%curve(k)%y_line, stat=istat)
        deallocate(s%template_plot(i)%graph(j)%curve(k)%x_symb, stat=istat)
        deallocate(s%template_plot(i)%graph(j)%curve(k)%y_symb, stat=istat)
        deallocate(s%template_plot(i)%graph(j)%curve(k)%ix_symb, stat=istat)
      enddo
      deallocate(s%template_plot(i)%graph(j)%curve, stat=istat)
    enddo
    deallocate(s%template_plot(i)%graph, stat=istat)
  enddo

  nullify(s%plot_page%plot)

! Universes 
  if (associated (s%u)) then
    do i = 1, size(s%u)
      ! Orbits
      deallocate(s%u(i)%model_orb, stat=istat)
      deallocate(s%u(i)%design_orb, stat=istat)
      deallocate(s%u(i)%base_orb, stat=istat)
      
      ! Beams
      deallocate(s%u(i)%macro_beam%ix_lost, stat=istat)
  
      ! d2_data
      do j = 1, size(s%u(i)%d2_data)
        do k = 1, size(s%u(i)%d2_data(j)%d1)
          nullify(s%u(i)%d2_data(j)%d1(k)%d2)
          nullify(s%u(i)%d2_data(j)%d1(k)%d)
        enddo
        deallocate(s%u(i)%d2_data(j)%d1, stat=istat)
      enddo
      deallocate(s%u(i)%d2_data, stat=istat)
 
 
      ! Data
      do j = lbound(s%u(i)%data,1), ubound(s%u(i)%data,1)
        nullify(s%u(i)%data(j)%d1)
      enddo
      deallocate(s%u(i)%data, stat=istat)
 
      ! dModel_dVar
      deallocate(s%u(i)%dmodel_dvar, stat=istat)
 
      ! Lattices
      call deallocate_lattice_internals(s%u(i)%model)
      call deallocate_lattice_internals(s%u(i)%design)
      call deallocate_lattice_internals(s%u(i)%base)
    enddo
  endif
    
end subroutine deallocate_everything
    

!------------------------------------------------------------------------------

subroutine deallocate_lattice_internals (lat)

implicit none

type (ring_struct) :: lat

integer j, istat

  ! Lattice elements
  do j = 0, size(lat%ele_)
    if (associated(lat%ele_(j)%r)) deallocate(lat%ele_(j)%r, stat=istat)
    if (associated(lat%ele_(j)%a)) deallocate(lat%ele_(j)%a, stat=istat)
    if (associated(lat%ele_(j)%b)) deallocate(lat%ele_(j)%b, stat=istat)
    if (associated(lat%ele_(j)%const)) deallocate(lat%ele_(j)%const, stat=istat)
    if (associated(lat%ele_(j)%descrip)) deallocate(lat%ele_(j)%descrip, stat=istat)
    if (associated(lat%ele_(j)%gen_field)) deallocate(lat%ele_(j)%gen_field, stat=istat)
    if (associated(lat%ele_(j)%wake)) then
      if (associated(lat%ele_(j)%wake%sr1)) deallocate(lat%ele_(j)%wake%sr1, stat=istat)  
      if (associated(lat%ele_(j)%wake%sr2_long)) deallocate(lat%ele_(j)%wake%sr2_long, stat=istat)  
      if (associated(lat%ele_(j)%wake%sr2_trans)) deallocate(lat%ele_(j)%wake%sr2_trans, stat=istat)  
      if (associated(lat%ele_(j)%wake%lr)) deallocate(lat%ele_(j)%wake%lr, stat=istat)  
      deallocate(lat%ele_(j)%wake, stat=istat)
    endif
    if (associated(lat%ele_(j)%wig_term)) deallocate(lat%ele_(j)%wig_term, stat=istat)
    nullify (lat%ele_(j)%r)
    nullify (lat%ele_(j)%a)
    nullify (lat%ele_(j)%b)
    nullify (lat%ele_(j)%const)
    nullify (lat%ele_(j)%descrip)
    nullify (lat%ele_(j)%gen_field)
    nullify (lat%ele_(j)%wake)
    nullify (lat%ele_(j)%wig_term)
  enddo
  deallocate(lat%ele_, stat=istat)
  nullify(lat%ele_)
  
  ! other stuff in lattice
  if (associated(lat%control_)) deallocate(lat%control_, stat=istat)
  if (associated(lat%ic_)) deallocate(lat%ic_, stat=istat)
  nullify(lat%control_)
  nullify(lat%ic_)
  nullify(lat%beam_energy)

end subroutine deallocate_lattice_internals

end subroutine tao_init
