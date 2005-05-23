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
  use tao_plot_mod

  implicit none

  type (tao_universe_struct), pointer :: u
  type (tao_var_struct), pointer :: var_ptr
  type (tao_this_var_struct), pointer :: this

  character(*) init_file
  character(200) lattice_file, plot_file, data_file, var_file, file_name
  character(200) single_mode_file, startup_file
  character(16) :: r_name = 'tao_init'
  integer i, j, n_universes, iu, ix

  namelist / tao_start / lattice_file, startup_file, &
               data_file, var_file, plot_file, single_mode_file, n_universes

! Find namelist files

  lattice_file     = init_file      ! set default
  plot_file        = init_file      ! set default
  data_file        = init_file      ! set default
  var_file         = init_file      ! set default
  single_mode_file = init_file      ! set default
  startup_file     = 'tao.startup'  ! set default
  n_universes = 1                   ! set default
  call tao_open_file ('TAO_INIT_DIR', init_file, iu, file_name)
  read (iu, nml = tao_start)
  close (iu)

  if (associated(s%u)) call deallocate_everything ()
  allocate (s%u(n_universes))

! Init

  call tao_init_design_lattice (lattice_file) 
  call tao_init_global_and_universes (init_file, data_file, var_file)
  call tao_init_single_mode (single_mode_file)
  call tao_hook_init (init_file)

! check variables
! check if vars are good

  call tao_set_var_useit_opt

  do i = 1, size(s%var)
    var_ptr => s%var(i)
    if (.not. var_ptr%exists) cycle
    do j = 1, size(var_ptr%this)
      this => var_ptr%this(j)
      u => s%u(this%ix_uni)
      ix = attribute_index (u%model%ele_(this%ix_ele), var_ptr%attrib_name)
      if (ix < 1) then
        call out_io (s_abort$, r_name, &
                'Error: Attribute not recognized: ' // var_ptr%attrib_name, &
                '       For element: ' // u%model%ele_(this%ix_ele)%name, &
                '       Which is a: ' // key_name(u%model%ele_(this%ix_ele)%key))
        call err_exit
      endif
      if (.not. attribute_free (u%model%ele_(this%ix_ele), ix, u%model)) then
        call out_io (s_abort$, r_name, &
                'Error: Variable trying to control an attribute that is not free to vary.', &
                '       Variable:  ' // var_ptr%name, &
                '       Element:   ' // var_ptr%ele_name, &
                '       Attribute: ' // var_ptr%attrib_name)
        call err_exit
      endif
    enddo
  enddo

! set up model and base lattices

  ! must first transfer to model lattice for tao_lattice_calc to run
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
  if (s%global%lattice_recalc) call tao_lattice_calc 

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

  type (tao_plot_struct), pointer :: plot
  type (tao_curve_struct), pointer :: curve

  integer i, j, k, istat

! Variables  

  if (associated (s%v1_var)) then
    deallocate(s%v1_var, stat=istat)
  endif
  
  if (associated (s%var)) then
    do i = lbound(s%var,1), ubound(s%var,1)
      deallocate(s%var(i)%this, stat=istat)
    enddo
    deallocate(s%var, stat=istat)
  endif

! Keytable 

  if (associated (s%key)) deallocate(s%key, stat=istat)

! plotting  

  nullify(s%plot_page%region)

  do i = 1, size(s%template_plot)
    plot => s%template_plot(i)
    if (.not. associated (plot%graph)) cycle
    do j = 1, size(plot%graph)
      if (.not. associated (plot%graph(j)%curve)) cycle
      do k = 1, size(plot%graph(j)%curve)
        curve => plot%graph(j)%curve(k)
        if (associated(curve%x_line)) deallocate(curve%x_line, stat=istat)
        if (associated(curve%y_line)) deallocate(curve%y_line, stat=istat)
        if (associated(curve%x_symb)) deallocate(curve%x_symb, stat=istat)
        if (associated(curve%y_symb)) deallocate(curve%y_symb, stat=istat)
        if (associated(curve%ix_symb)) deallocate(curve%ix_symb, stat=istat)
      enddo
      deallocate(plot%graph(j)%curve, stat=istat)
    enddo
    deallocate(plot%graph, stat=istat)
  enddo

! Universes 

  if (associated (s%u)) then
    do i = 1, size(s%u)
      ! Orbits
      deallocate(s%u(i)%model_orb, stat=istat)
      deallocate(s%u(i)%design_orb, stat=istat)
      deallocate(s%u(i)%base_orb, stat=istat)
      
      ! Beams
      deallocate (s%u(i)%macro_beam%ix_lost, stat=istat)
      call reallocate_macro_beam (s%u(i)%macro_beam%beam, 0, 0, 0)
      call reallocate_beam (s%u(i)%beam%beam, 0, 0)
  
      ! Coupling
      call deallocate_ele_pointers (s%u(i)%coupling%coupling_ele)
      call reallocate_macro_beam (s%u(i)%coupling%injecting_macro_beam, 0, 0, 0)
      call reallocate_beam (s%u(i)%coupling%injecting_beam, 0, 0)
      
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
 
      ! ix_data
      do j = 0, ubound(s%u(i)%ix_data,1)
        if (associated(s%u(i)%ix_data(j)%ix_datum)) deallocate(s%u(i)%ix_data(j)%ix_datum)
      enddo
      deallocate(s%u(i)%ix_data)
      
      ! dModel_dVar
      deallocate(s%u(i)%dmodel_dvar, stat=istat)
 
      ! Lattices
      call deallocate_ring_pointers (s%u(i)%model)
      call deallocate_ring_pointers (s%u(i)%design)
      call deallocate_ring_pointers (s%u(i)%base)
    enddo
  endif
    
end subroutine deallocate_everything
    
end subroutine tao_init


