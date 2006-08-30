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
  use tao_command_mod
  use tao_plot_mod

  implicit none

  type (tao_universe_struct), pointer :: u
  type (tao_var_struct), pointer :: var_ptr
  type (tao_this_var_struct), pointer :: this
  type (tao_plot_struct), pointer :: p
  real(rp), pointer :: ptr_attrib

  character(*) init_file
  character(200) lattice_file, plot_file, data_file, var_file, file_name
  character(200) single_mode_file, startup_file
  character(40) name1, name2
  character(16) :: r_name = 'tao_init'
  character(16) init_name
  integer i, j, i2, j2, n_universes, iu, ix, ix_attrib

  logical err, calc_ok

  namelist / tao_start / lattice_file, startup_file, &
               data_file, var_file, plot_file, single_mode_file, &
               n_universes, init_name

! Find namelist files

  lattice_file     = init_file      ! set default
  plot_file        = init_file      ! set default
  data_file        = init_file      ! set default
  var_file         = init_file      ! set default
  single_mode_file = init_file      ! set default
  startup_file     = 'tao.startup'  ! set default
  n_universes      = 1              ! set default
  init_name        = "Tao"          ! set default
  call tao_open_file ('TAO_INIT_DIR', init_file, iu, file_name)
  read (iu, nml = tao_start)
  close (iu)
  tao_com%init_name = init_name

  if (associated(s%u)) call deallocate_everything ()
  allocate (s%u(n_universes))

! Init

  call tao_init_design_lattice (lattice_file) 
  call tao_init_global_and_universes (init_file, data_file, var_file)
  call tao_init_single_mode (single_mode_file)
  call tao_hook_init (init_file)

  bmad_status%exit_on_error = .false.

! check variables
! check if vars are good

  call tao_set_var_useit_opt

  do i = 1, size(s%var)
    var_ptr => s%var(i)
    if (.not. var_ptr%exists) cycle
    do j = 1, size(var_ptr%this)
      this => var_ptr%this(j)
      u => s%u(this%ix_uni)
      call pointer_to_attribute (u%model%lat%ele_(this%ix_ele), var_ptr%attrib_name, &
                                 .false., ptr_attrib, err, .false., ix_attrib)
      if (err) then
        call out_io (s_abort$, r_name, &
                'Error: Attribute not recognized: ' // var_ptr%attrib_name, &
                '       For element: ' // u%model%lat%ele_(this%ix_ele)%name, &
                '       Which is a: ' // key_name(u%model%lat%ele_(this%ix_ele)%key))
        call err_exit
      endif
      if (.not. attribute_free (this%ix_ele, ix_attrib, u%model%lat)) then
        call out_io (s_abort$, r_name, &
                'Error: Variable trying to control an attribute that is not free to vary.', &
                '       Variable:  ' // var_ptr%name, &
                '       Element:   ' // var_ptr%ele_name, &
                '       Attribute: ' // var_ptr%attrib_name)
        call err_exit
      endif
    enddo
  enddo

! make sure two variables do not vary the same attribute

  do i = 1, size(s%var)
    do j = 1, size(s%var(i)%this)
      do i2 = i, size(s%var)
        do j2 = 1, size(s%var(i2)%this)
          if (i == i2 .and. j == j2) cycle
          if (associated (s%var(i)%this(j)%model_ptr, &
                            s%var(i2)%this(j2)%model_ptr)) then
            write (name1, '(2a, i0, a)') trim(s%var(i)%v1%name), '[', s%var(i)%ix_v1, ']'  
            write (name2, '(2a, i0, a)') trim(s%var(i2)%v1%name), '[', s%var(i2)%ix_v1, ']'  
            call out_io (s_error$, r_name, &
               'ERROR: VARIABLE:     ' // name1, &
               '       AND VARIABLE: ' // name2, &
               '       CONTROL THE SAME EXACT THING!', &
               '       YOU HAVE BEEN WARNED!!!')            
          endif
        enddo
      enddo
    enddo
  enddo

! plotting

  call tao_init_plotting (plot_file)
  
! set up model and base lattices

  ! must first transfer to model lattice for tao_lattice_calc to run
  do i = 1, size(s%u)
    s%u(i)%model = s%u(i)%design; s%u(i)%model%orb = s%u(i)%design%orb
  enddo
  s%global%lattice_recalc = .true.
  call tao_lattice_calc (calc_ok, .true.) ! .true. => init design lattice
  do i = 1, size(s%u)
    s%u(i)%design = s%u(i)%model; s%u(i)%design%orb = s%u(i)%model%orb
    s%u(i)%base  = s%u(i)%design; s%u(i)%base%orb  = s%u(i)%design%orb
    s%u(i)%data%design_value = s%u(i)%data%model_value
    s%u(i)%data%base_value = s%u(i)%data%model_value
  enddo

  call tao_plot_data_setup ()  ! transfer data to the plotting structures
  call tao_plot_out ()         ! Update the plotting window

! Look for a startup file

  call tao_open_file ('TAO_INIT_DIR', startup_file, iu, file_name)
  if (iu /= 0) then
    call out_io (s_blank$, r_name, 'Using startup file: ' // file_name)
    tao_com%cmd_from_cmd_file = .false.
    call tao_cmd_history_record ('call ' // startup_file)
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
  type (tao_universe_struct), pointer :: u

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

! Plotting  

  nullify(s%plot_page%region)

  do i = 1, size(s%template_plot)
    plot => s%template_plot(i)
    if (.not. allocated (plot%graph)) cycle
    do j = 1, size(plot%graph)
      if (.not. allocated (plot%graph(j)%curve)) cycle
      do k = 1, size(plot%graph(j)%curve)
        curve => plot%graph(j)%curve(k)
        if (allocated(curve%x_line)) deallocate(curve%x_line, stat=istat)
        if (allocated(curve%y_line)) deallocate(curve%y_line, stat=istat)
        if (allocated(curve%x_symb)) deallocate(curve%x_symb, stat=istat)
        if (allocated(curve%y_symb)) deallocate(curve%y_symb, stat=istat)
        if (allocated(curve%ix_symb)) deallocate(curve%ix_symb, stat=istat)
      enddo
      deallocate(plot%graph(j)%curve, stat=istat)
    enddo
    deallocate(plot%graph, stat=istat)
  enddo

! Universes 

  if (associated (s%u)) then
    do i = 1, size(s%u)

      u => s%u(i)
      ! radiation integrals cache
      if (u%ix_rad_int_cache /= 0) call release_rad_int_cache(u%ix_rad_int_cache)


      ! Orbits
      deallocate(u%model%orb, stat=istat)
      deallocate(u%design%orb, stat=istat)
      deallocate(u%base%orb, stat=istat)
      
      ! Beams
      deallocate (u%macro_beam%ix_lost, stat=istat)
      call reallocate_macro_beam (u%macro_beam%beam, 0, 0, 0)
      call reallocate_beam (u%beam%beam, 0, 0)
  
      ! Coupling
      call deallocate_ele_pointers (u%coupling%coupling_ele)
      call reallocate_macro_beam (u%coupling%injecting_macro_beam, 0, 0, 0)
      call reallocate_beam (u%coupling%injecting_beam, 0, 0)
      
      ! d2_data
      do j = 1, size(u%d2_data)
        do k = 1, size(u%d2_data(j)%d1)
          nullify(u%d2_data(j)%d1(k)%d2)
          nullify(u%d2_data(j)%d1(k)%d)
        enddo
        deallocate(u%d2_data(j)%d1, stat=istat)
      enddo
      deallocate(u%d2_data, stat=istat)
 
      ! Data
      do j = lbound(u%data,1), ubound(u%data,1)
        nullify(u%data(j)%d1)
      enddo
      deallocate(u%data, stat=istat)
 
      ! ix_data
      do j = 0, ubound(u%ix_data,1)
        if (associated(u%ix_data(j)%ix_datum)) deallocate(u%ix_data(j)%ix_datum)
      enddo
      deallocate(u%ix_data)
      
      ! dModel_dVar
      deallocate(u%dmodel_dvar, stat=istat)
 
      ! Lattices
      call deallocate_ring_pointers (u%model%lat)
      call deallocate_ring_pointers (u%design%lat)
      call deallocate_ring_pointers (u%base%lat)
    enddo
  endif
    
end subroutine deallocate_everything
    
end subroutine tao_init


