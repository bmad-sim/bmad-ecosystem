!+
! Subroutine pointers_to_attribute (lat, ele_name, attrib_name, do_allocation,
!                     ptr_array, err_flag, err_print_flag, eles, ix_attrib, do_unlink)
!
! Returns an array of pointers to an attribute with name attrib_name within elements with name ele_name.
!
! Note: If, for example, ele_name = 'BMAD_COM', there is not corresponding lattice element and
!   therefore eles will have size 0 on output.
! Note: ele_name = 'PARTICLE_START' corresponds to the lat%particle_start substructure. 
! Note: ele_name can be a list of element indices. For example:
!           ele_name = "3:5"
!  This sets elements 3, 4, and 5 in the lat%ele(:) array.
! Note: ele_name can be in key:name format and include the wild card characters "*" and "%".
!       For example: "quad:q*"
! Note: Use attribute_free to see if the attribute may be varied independently.
! Note: When using wild cards, it is *not* an error if some of the matched elements do 
!       not have the the required attribute as long as at least one does.
! Note: Alternatively consider the routines:
!     set_ele_attribute
!     value_of_attribute
!
! Input:
!   lat             -- lat_struct: Lattice.
!   ele_name        -- Character(*): Element name. Must be uppercase
!   attrib_name     -- Character(*): Attribute name. Must be uppercase.
!                       For example: "HKICK".
!   do_allocation   -- Logical: If True then do an allocation if needed.
!                       EG: The multipole An and Bn arrays need to be allocated
!                       before their use.
!   err_print_flag  -- Logical, optional: If present and False then suppress
!                       printing of an error message on error.
!
! Output:
!   ptr_array(:) -- all_pointer_struct, allocatable: Pointer to the attribute.
!                     Size of ptr_array will be set to 0 if there is a problem.
!   err_flag     -- Logical: Set True if attribtute not found.
!   eles(:)      -- Ele_pointer_struct, optional, allocatable: Array of element pointers.
!                     size(eles) = size(ptr_array). If there are no associated lattice 
!                     elements (EG if ele_name = 'PARTICLE_START'), eles(i)%ele will be null.
!                     
!   ix_attrib    -- Integer, optional: If applicable then this is the index to the 
!                     attribute in the ele%value(:) array.
!-

Subroutine pointers_to_attribute (lat, ele_name, attrib_name, do_allocation, &
                        ptr_array, err_flag, err_print_flag, eles, ix_attrib, do_unlink)

use bmad_interface, except_dummy => pointers_to_attribute

implicit none

type (lat_struct), target :: lat
type (ele_struct), target :: particle_start
type (all_pointer_struct), allocatable :: ptr_array(:)
type (all_pointer_struct), allocatable :: ptrs(:)
type (ele_pointer_struct), optional, allocatable :: eles(:)
type (ele_pointer_struct), allocatable :: eles2(:)
type (all_pointer_struct) a_ptr

integer, optional :: ix_attrib
integer n, i, ix, key, ix_a, n_loc

character(*) ele_name, attrib_name
character(20) str
character(*), parameter :: r_name = 'pointers_to_attribute'

logical err_flag, do_allocation, do_print
logical, optional :: err_print_flag, do_unlink

! init

err_flag = .false.
do_print = logic_option (.true., err_print_flag)
call re_allocate(ptr_array, 0)
if (present(eles)) call re_allocate_eles (eles, 0)

!---

select case (ele_name)

! Selected parameters in bmad_com

case ('BMAD_COM')
  select case (attrib_name)
  case ('D_ORB')
    call re_allocate (ptr_array, 6)
    do i = 1, 6
      ptr_array(i)%r => bmad_com%d_orb(i)
    enddo
    return
  end select

  if (attrib_name(1:5) == 'D_ORB') then
    str = attrib_name(6:)
    if (str(1:1) /= '(') err_flag = .true.
    call string_trim(str(2:), str, ix)
    n = index('123456', str(1:1))
    if (n == 0) err_flag = .true.
    call string_trim(str(2:), str, ix)
    if (str /= ')') err_flag = .true.
    if (.not. err_flag) then
      call re_allocate (ptr_array, 1)
      ptr_array(1)%r => bmad_com%d_orb(n)
      return
    endif
  endif    

  call re_allocate (ptr_array, 1)

  select case(attrib_name)
  case ('MAX_APERTURE_LIMIT');              ptr_array(1)%r => bmad_com%max_aperture_limit
  case ('DEFAULT_DS_STEP');                 ptr_array(1)%r => bmad_com%default_ds_step
  case ('SIGNIFICANT_LENGTH');              ptr_array(1)%r => bmad_com%significant_length
  case ('REL_TOL_TRACKING');                ptr_array(1)%r => bmad_com%rel_tol_tracking
  case ('ABS_TOL_TRACKING');                ptr_array(1)%r => bmad_com%abs_tol_tracking
  case ('REL_TOL_ADAPTIVE_TRACKING');       ptr_array(1)%r => bmad_com%rel_tol_adaptive_tracking
  case ('ABS_TOL_ADAPTIVE_TRACKING');       ptr_array(1)%r => bmad_com%abs_tol_adaptive_tracking
  case ('INIT_DS_ADAPTIVE_TRACKING');       ptr_array(1)%r => bmad_com%init_ds_adaptive_tracking
  case ('MIN_DS_ADAPTIVE_TRACKING');        ptr_array(1)%r => bmad_com%min_ds_adaptive_tracking
  case ('FATAL_DS_ADAPTIVE_TRACKING');      ptr_array(1)%r => bmad_com%fatal_ds_adaptive_tracking
  case ('AUTOSCALE_AMP_ABS_TOL');           ptr_array(1)%r => bmad_com%autoscale_amp_abs_tol
  case ('AUTOSCALE_AMP_REL_TOL');           ptr_array(1)%r => bmad_com%autoscale_amp_rel_tol
  case ('AUTOSCALE_PHASE_TOL');             ptr_array(1)%r => bmad_com%autoscale_phase_tol
  case ('ELECTRIC_DIPOLE_MOMENT');          ptr_array(1)%r => bmad_com%electric_dipole_moment
  case ('SYNCH_RAD_SCALE');                 ptr_array(1)%r => bmad_com%synch_rad_scale
  case ('SAD_EPS_SCALE');                   ptr_array(1)%r => bmad_com%sad_eps_scale
  case ('SAD_AMP_MAX');                     ptr_array(1)%r => bmad_com%sad_amp_max

  case ('SAD_N_DIV_MAX');                   ptr_array(1)%i => bmad_com%sad_n_div_max
  case ('TAYLOR_ORDER');                    ptr_array(1)%i => bmad_com%taylor_order
  case ('DEFAULT_INTEG_ORDER');             ptr_array(1)%i => bmad_com%default_integ_order
  case ('RUNGE_KUTTA_ORDER');               ptr_array(1)%i => bmad_com%runge_kutta_order
  case ('MAX_NUM_RUNGE_KUTTA_STEP');        ptr_array(1)%i => bmad_com%max_num_runge_kutta_step

  case ('RF_PHASE_BELOW_TRANSITION_REF');   ptr_array(1)%l => bmad_com%rf_phase_below_transition_ref
  case ('SR_WAKES_ON');                     ptr_array(1)%l => bmad_com%sr_wakes_on
  case ('LR_WAKES_ON');                     ptr_array(1)%l => bmad_com%lr_wakes_on
  case ('AUTO_BOOKKEEPER');                 ptr_array(1)%l => bmad_com%auto_bookkeeper
  case ('HIGH_ENERGY_SPACE_CHARGE_ON');     ptr_array(1)%l => bmad_com%high_energy_space_charge_on
  case ('CSR_AND_SPACE_CHARGE_ON');         ptr_array(1)%l => bmad_com%csr_and_space_charge_on
  case ('SPIN_TRACKING_ON');                ptr_array(1)%l => bmad_com%spin_tracking_on
  case ('SPIN_SOKOLOV_TERNOV_FLIPPING_ON'); ptr_array(1)%l => bmad_com%spin_sokolov_ternov_flipping_on
  case ('RADIATION_DAMPING_ON');            ptr_array(1)%l => bmad_com%radiation_damping_on
  case ('RADIATION_ZERO_AVERAGE');          ptr_array(1)%l => bmad_com%radiation_zero_average
  case ('RADIATION_FLUCTUATIONS_ON');       ptr_array(1)%l => bmad_com%radiation_fluctuations_on
  case ('CONSERVE_TAYLOR_MAPS');            ptr_array(1)%l => bmad_com%conserve_taylor_maps
  case ('ABSOLUTE_TIME_TRACKING');          ptr_array(1)%l => bmad_com%absolute_time_tracking
  case ('ABSOLUTE_TIME_REF_SHIFT');         ptr_array(1)%l => bmad_com%absolute_time_ref_shift
  case ('CONVERT_TO_KINETIC_MOMENTUM');     ptr_array(1)%l => bmad_com%convert_to_kinetic_momentum
  case ('APERTURE_LIMIT_ON');               ptr_array(1)%l => bmad_com%aperture_limit_on
  case ('DEBUG');                           ptr_array(1)%l => bmad_com%debug

  case default
    if (do_print) call out_io (s_error$, r_name, &
             'INVALID ATTRIBUTE: ' // attrib_name, 'FOR ELEMENT: ' // ele_name)
    call re_allocate(ptr_array, 0)
    err_flag = .true.
  end select

  return

! space_charge_com

case ('SPACE_CHARGE_COM')
  call re_allocate (ptr_array, 1)

  select case (attrib_name)

  case ('SPACE_CHARGE_MESH_SIZE')
    call re_allocate (ptr_array, 3)
    do i = 1, 3
      ptr_array(i)%i => space_charge_com%space_charge_mesh_size(i)
    enddo

  case ('CSR3D_MESH_SIZE')
    call re_allocate (ptr_array, 3)
    do i = 1, 3
      ptr_array(i)%i => space_charge_com%csr3d_mesh_size(i)
    enddo

  case ('DS_TRACK_STEP');                   ptr_array(1)%r => space_charge_com%ds_track_step
  case ('DT_TRACK_STEP');                   ptr_array(1)%r => space_charge_com%dt_track_step
  case ('CATHODE_STRENGTH_CUTOFF');         ptr_array(1)%r => space_charge_com%cathode_strength_cutoff
  case ('REL_TOL_TRACKING');                ptr_array(1)%r => space_charge_com%rel_tol_tracking
  case ('ABS_TOL_TRACKING');                ptr_array(1)%r => space_charge_com%abs_tol_tracking
  case ('BEAM_CHAMBER_HEIGHT');             ptr_array(1)%r => space_charge_com%beam_chamber_height
  case ('LSC_SIGMA_CUTOFF');                ptr_array(1)%r => space_charge_com%lsc_sigma_cutoff
  case ('PARTICLE_SIGMA_CUTOFF');           ptr_array(1)%r => space_charge_com%particle_sigma_cutoff
  case ('N_BIN');                           ptr_array(1)%i => space_charge_com%n_bin
  case ('PARTICLE_BIN_SPAN');               ptr_array(1)%i => space_charge_com%particle_bin_span
  case ('N_SHIELD_IMAGES');                 ptr_array(1)%i => space_charge_com%n_shield_images
  case ('SC_MIN_IN_BIN');                   ptr_array(1)%i => space_charge_com%sc_min_in_bin
  case ('LSC_KICK_TRANSVERSE_DEPENDENCE');  ptr_array(1)%l => space_charge_com%lsc_kick_transverse_dependence
  case ('DEBUG');                           ptr_array(1)%l => space_charge_com%debug

  case default
    if (do_print) call out_io (s_error$, r_name, &
             'INVALID ATTRIBUTE: ' // attrib_name, 'FOR ELEMENT: ' // ele_name)
    call re_allocate(ptr_array, 0)
    err_flag = .true.
  end select

  return

! particle_start

case ('PARTICLE_START')

  call re_allocate (ptr_array, 1)
  if (present(eles)) call re_allocate_eles (eles, 1)

  ix = attribute_index (particle_start, attrib_name)
  if (ix < 1) then
    if (do_print) call out_io (s_error$, r_name, &
           'INVALID ATTRIBUTE: ' // attrib_name, 'FOR ELEMENT: ' // ele_name)
    call re_allocate(ptr_array, 0)
    if (present(eles)) call re_allocate_eles (eles, 0)
    err_flag = .true.
    return
  endif
  if (present(ix_attrib)) ix_attrib = ix
  select case (ix)
  case (x$);                             ptr_array(1)%r => lat%particle_start%vec(1)
  case (px$);                            ptr_array(1)%r => lat%particle_start%vec(2)
  case (y$);                             ptr_array(1)%r => lat%particle_start%vec(3)
  case (py$);                            ptr_array(1)%r => lat%particle_start%vec(4)
  case (z$);                             ptr_array(1)%r => lat%particle_start%vec(5)
  case (pz$);                            ptr_array(1)%r => lat%particle_start%vec(6)
  case (field_x$);                       ptr_array(1)%r => lat%particle_start%field(1)
  case (field_y$);                       ptr_array(1)%r => lat%particle_start%field(2)
  case (phase_x$);                       ptr_array(1)%r => lat%particle_start%phase(1)
  case (phase_y$);                       ptr_array(1)%r => lat%particle_start%phase(2)
  case (t$);                             ptr_array(1)%q => lat%particle_start%t    ! Quad precision
  case (e_photon$);                      ptr_array(1)%r => lat%particle_start%p0c
  case (spin_x$);                        ptr_array(1)%r => lat%particle_start%spin(1)
  case (spin_y$);                        ptr_array(1)%r => lat%particle_start%spin(2)
  case (spin_z$);                        ptr_array(1)%r => lat%particle_start%spin(3)
  case (emittance_a$);                   ptr_array(1)%r => lat%a%emit
  case (emittance_b$);                   ptr_array(1)%r => lat%b%emit
  case (emittance_z$);                   ptr_array(1)%r => lat%z%emit
  case (sig_x$);                         ptr_array(1)%r => lat%a%sigma
  case (sig_y$);                         ptr_array(1)%r => lat%b%sigma
  case (sig_z$);                         ptr_array(1)%r => lat%z%sigma
  case (sig_e$);                         ptr_array(1)%r => lat%z%sigmap

  case default
    if (do_print) call out_io (s_error$, r_name, &
             'INVALID ATTRIBUTE: ' // attrib_name, 'FOR ELEMENT: ' // ele_name)
    call re_allocate(ptr_array, 0)
    if (present(eles)) call re_allocate_eles (eles, 0)
    err_flag = .true.
  end select

  return

! This section is for things like "parameter[n_part]" whose value is stored
! in a non-standard location. Things like "parameter[e_tot]" are handled in 
! the usual way in the code section after this one.

case ('PARAMETER')

  select case(attrib_name)
  case ('N_PART')
    call re_allocate (ptr_array, 1)
    if (present(eles)) call re_allocate_eles (eles, 1)
    ptr_array(1)%r => lat%param%n_part
    return
  end select

end select

! Locate elements

call lat_ele_locator (ele_name, lat, eles2, n_loc, err_flag)
if (n_loc == 0) then
  if (do_print) call out_io (s_error$, r_name, 'ELEMENT NOT FOUND: ' // ele_name)
  call re_allocate(ptr_array, 0)
  if (present(eles)) call re_allocate_eles (eles, 0)
  err_flag = .true.
  return  
endif

! Locate attributes

call re_allocate (ptrs, n_loc)
n = 0
do i = 1, n_loc
  call pointer_to_attribute (eles2(i)%ele, attrib_name, do_allocation, a_ptr, err_flag, .false., ix_a, do_unlink)
  if (err_flag) cycle
  n = n + 1
  ptrs(n) = a_ptr
  eles2(n)%ele => eles2(i)%ele
  if (present(ix_attrib)) ix_attrib = ix_a
enddo

if (n == 0) then
  if (do_print) call out_io (s_error$, r_name, 'ATTRIBUTE: ' // attrib_name, &
                                               'NOT FOUND FOR: ' // ele_name)
  call re_allocate(ptr_array, 0)
  if (present(eles)) call re_allocate_eles (eles, 0)
  err_flag = .true.
  return  
endif

! Transfer pointers to ptr_array

if (present(eles)) then
  call re_allocate_eles (eles, n)
  eles = eles2(1:n)
endif

call re_allocate (ptr_array, n)
do i = 1, n
  ptr_array(i) = ptrs(i)
enddo

err_flag = .false.

end subroutine

