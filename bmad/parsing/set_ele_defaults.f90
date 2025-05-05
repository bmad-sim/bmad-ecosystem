!+
! Subroutine set_ele_defaults (ele, do_allocate)
!
! Subroutine set the defaults for an element of a given type.
! For example, the default aperture type for an ecollimator$
!   element is ele%aperture_type = elliptical$.
!
! Input:
!   ele           -- ele_struct: Element to init.
!     %key          -- Type of element.
!   do_allocate   -- logical, optional: Do default allocation of element components? Default is True.
!
! Output:
!   ele   -- ele_struct: Initialized element.
!-

subroutine set_ele_defaults (ele, do_allocate)

use bmad_interface, dummy => set_ele_defaults

implicit none

type (ele_struct) ele
integer i, j
logical, optional :: do_allocate
real(rp) g

!

ele%component_name = ''
ele%taylor_map_includes_offsets = .true.

! Some overall defaults.

if (has_attribute(ele, 'FRINGE_AT'))            ele%value(fringe_at$) = both_ends$
if (has_attribute(ele, 'FRINGE_TYPE'))          ele%value(fringe_type$) = none$
if (has_attribute(ele, 'SPIN_FRINGE_ON'))       ele%value(spin_fringe_on$) = true$
if (has_attribute(ele, 'PTC_CANONICAL_COORDS')) ele%value(ptc_canonical_coords$) = true$
if (has_attribute(ele, 'MAT6_CALC_METHOD'))     ele%mat6_calc_method = auto$
call init_multipole_cache(ele)

! Other inits.

select case (ele%key)

case (ac_kicker$)
  if (logic_option(.true., do_allocate)) then
    if (associated(ele%ac_kick)) deallocate(ele%ac_kick)
    allocate (ele%ac_kick)
  endif
  ele%value(interpolation$) = cubic$

case (beambeam$)
  ele%value(charge$) = -1
  ele%value(n_slice$) = 1
  ele%value(species_strong$) = real_garbage$
  ele%value(E_tot_strong$) = -1
  ele%value(pc_strong$) = -1

case (beginning_ele$)
  ele%value(e_tot$) = -1
  ele%value(p0c$) = -1
  ele%value(inherit_from_fork$) = real_garbage$
  ele%value(deta_ds_master$) = false$
  ele%z%etap = 1
  ele%z%deta_ds = 1
  call mat_make_unit (ele%mat6)

case (capillary$)
  ele%offset_moves_aperture = .true.

case (converter$)
  if (logic_option(.true., do_allocate)) then
    if (associated(ele%converter)) deallocate(ele%converter)
    allocate (ele%converter)
  endif

case (crab_cavity$)
  ele%value(field_autoscale$) = 1
  ele%value(num_steps$) = 10

case (crystal$)
  ele%value(is_mosaic$) = false$
  ele%value(mosaic_angle_rms_out_plane$) = -1
  ele%value(ref_orbit_follows$) = bragg_diffracted$
  ele%aperture_at = surface$
  ele%offset_moves_aperture = .true.
  ele%value(use_reflectivity_table$) = false$
  if (logic_option(.true., do_allocate)) then
    ! Avoid "ele%photon = photon_element_struct()" to get around Ifort bug. 4/10/2019
    if (associated(ele%photon)) deallocate(ele%photon)
    allocate(ele%photon)
  endif

case (custom$)  
  ele%tracking_method  = custom$
  ele%field_calc       = custom$

case (def_particle_start$)

case (def_line$)
  g = real_garbage$
  ele%value = g
  ele%s     = g
  ele%ref_time = g
  ele%a = twiss_struct(g, g, g, g, g, g, g, g, g, g, g, g, g)
  ele%b = twiss_struct(g, g, g, g, g, g, g, g, g, g, g, g, g)
  ele%z = twiss_struct(g, g, g, g, g, g, g, g, g, g, g, g, g)
  ele%x = xy_disp_struct(g, g, g, g)
  ele%y = xy_disp_struct(g, g, g, g)
  ele%floor = floor_position_struct([g,g,g], mat3_unit$, g,g,g)
  ele%value(inherit_from_fork$) = g

case (def_mad_beam$)
  ele%ref_species = not_set$

case (def_parameter$)
  ele%value(geometry$) = real_garbage$
  ele%value(live_branch$) = real_garbage$
  ele%value(high_energy_space_charge_on$) = real_garbage$
  ele%ref_species = positron$
  ele%value(default_tracking_species$) = real_garbage$
  ele%value(ix_branch$) = -1

case (detector$)
  ele%aperture_at = surface$
  ele%aperture_type = auto_aperture$
  if (logic_option(.true., do_allocate)) then
    ! Avoid "ele%photon = photon_element_struct()" to get around ifort bug. 4/10/2019
    if (associated(ele%photon)) deallocate(ele%photon)
    allocate(ele%photon)
  endif

case (diffraction_plate$)
  ele%aperture_at = surface$
  ele%aperture_type = auto_aperture$
  ele%offset_moves_aperture = .true.
  ele%value(mode$) = transmission$
  if (logic_option(.true., do_allocate)) then
    ! Avoid "ele%photon = photon_element_struct()" to get around ifort bug. 4/10/2019
    if (associated(ele%photon)) deallocate(ele%photon)
    allocate(ele%photon)
  endif

case (e_gun$)
  ele%tracking_method = time_runge_kutta$
  ele%value(field_autoscale$) = 1
  ele%value(fringe_at$) = exit_end$
  ele%value(fringe_type$) = full$
  ele%value(autoscale_amplitude$) = true$
  ele%value(autoscale_phase$) = true$

case (ecollimator$)
  ele%aperture_type = elliptical$
  ele%offset_moves_aperture = .true.

case (em_field$)
  ele%tracking_method = runge_kutta$
  ele%value(fringe_type$) = full$
  ele%value(field_autoscale$) = 1
  ele%value(constant_ref_energy$) = true$
  ele%value(polarity$) = 1.0

case (feedback$)

case (fiducial$)
  ele%value(origin_ele_ref_pt$) = center_pt$

case (floor_shift$)
  ele%value(origin_ele_ref_pt$) = exit_end$
  ele%value(upstream_coord_dir$) = 1
  ele%value(downstream_coord_dir$) = 1

case (foil$)
  ele%value(num_steps$) = 10
  ele%value(final_charge$) = real_garbage$
  ele%value(x1_edge$) = -99.0_rp
  ele%value(y1_edge$) = -99.0_rp
  ele%value(x2_edge$) =  99.0_rp
  ele%value(y2_edge$) =  99.0_rp
  ele%value(f_factor$)       = 0.98_rp
  ele%value(scatter_method$) = highland$
  if (associated(ele%foil)) deallocate(ele%foil)
  allocate(ele%foil)

case (fork$, photon_fork$)
  ele%value(direction$) = 1
  ele%ref_species = not_set$   ! Used when element is after a converter
  ele%value(new_branch$) = true$

case (girder$)
  ele%value(origin_ele_ref_pt$) = center_pt$

case (hybrid$) 
  ele%tracking_method = taylor$
  ele%mat6_calc_method = taylor$

case (lcavity$)
  ele%value(coupler_at$) = exit_end$
  ele%value(field_autoscale$) = 1
  ele%value(n_cell$) = 1
  ele%value(cavity_type$) = standing_wave$
  ele%value(fringe_type$) = full$
  ele%value(autoscale_amplitude$) = true$
  ele%value(autoscale_phase$) = true$
  ele%value(longitudinal_mode$) = 0
  ! So to not affect the changeover when the step loop was finally implemented in 2/2024
  ele%value(num_steps$) = 1       

case (marker$)
  ele%ref_species = not_set$    ! Used when element is after a converter

case (mask$)
  ele%aperture_at = surface$
  ele%aperture_type = auto_aperture$
  ele%offset_moves_aperture = .true.
  ele%value(mode$) = transmission$
  if (logic_option(.true., do_allocate)) then
    ! Avoid "ele%photon = photon_element_struct()" to get around ifort bug. 4/10/2019
    if (associated(ele%photon)) deallocate(ele%photon)
    allocate(ele%photon)
  endif

case (match$)
  ele%value(matrix$) = standard$
  ele%value(kick0$) = standard$
  ele%value(recalc$) = true$
  ele%value(deta_ds_master$) = false$
  ele%spin_tracking_method = off$

case (mirror$)
  ele%aperture_at = surface$
  ele%offset_moves_aperture = .true.
  ele%value(use_reflectivity_table$) = false$
  if (logic_option(.true., do_allocate)) then
    ! Avoid "ele%photon = photon_element_struct()" to get around ifort bug. 4/10/2019
    if (associated(ele%photon)) deallocate(ele%photon)
    allocate(ele%photon)
  endif

case (multilayer_mirror$)
  ele%aperture_at = surface$
  ele%offset_moves_aperture = .true.
  if (logic_option(.true., do_allocate)) then
    ! Avoid "ele%photon = photon_element_struct()" to get around Ifort bug. 4/10/2019
    if (associated(ele%photon)) deallocate(ele%photon)
    allocate(ele%photon)
  endif

case (ab_multipole$)
  if (logic_option(.true., do_allocate)) then
    call multipole_init (ele, magnetic$, .true.)
  endif
  ele%scale_multipoles = .false.

case (multipole$)
  if (logic_option(.true., do_allocate)) then
    call multipole_init (ele, magnetic$, .true.)
    call multipole_init (ele, electric$, .true.)
  endif
  ele%scale_multipoles = .false.

case (patch$)
  ele%value(flexible$) = false$ 
  ele%value(ref_coords$)= exit_end$
  ele%value(upstream_coord_dir$) = 1
  ele%value(downstream_coord_dir$) = 1

case (photon_init$)
  ele%value(ds_slice$) = 0.01
  ele%value(velocity_distribution$) = gaussian$
  ele%value(energy_distribution$) = gaussian$
  ele%value(spatial_distribution$) = gaussian$
  ele%value(transverse_sigma_cut$) = 3
  ele%value(E_center_relative_to_ref$) = true$
  if (logic_option(.true., do_allocate)) then
    ! Avoid "ele%photon = photon_element_struct()" to get around Ifort bug. 4/10/2019
    if (associated(ele%photon)) deallocate(ele%photon)
    allocate(ele%photon)
  endif

case (rf_bend$)
  ele%tracking_method = runge_kutta$
  ele%field_calc = fieldmap$
  ele%value(num_steps$) = 10
  ele%value(fiducial_pt$) = none_pt$
  ele%value(init_needed$) = true$

case (rbend$, sbend$)
  ele%value(fintx$) = real_garbage$
  ele%value(hgapx$) = real_garbage$
  ele%value(fringe_type$) = basic_bend$
  ele%value(ptc_fringe_geometry$) = x_invariant$
  ele%value(exact_multipoles$) = off$
  ele%value(ptc_field_geometry$) = sector$
  ele%value(fiducial_pt$) = none_pt$
  ele%value(init_needed$) = true$

case (rcollimator$)
  ele%offset_moves_aperture = .true.

case (rfcavity$)
  ele%value(coupler_at$) = exit_end$
  ele%value(field_autoscale$) = 1
  ele%value(n_cell$) = 1
  ele%value(cavity_type$) = standing_wave$
  ele%value(fringe_type$) = full$
  ele%value(autoscale_amplitude$) = true$
  ele%value(autoscale_phase$) = true$
  ele%value(longitudinal_mode$) = 0
  ele%value(num_steps$) = 10

case (sad_mult$)
  ele%value(eps_step_scale$) = 1
  ele%scale_multipoles = .false.
  if (logic_option(.true., do_allocate)) then
    call multipole_init (ele, magnetic$, .true.)
  endif

case (sample$)
  ele%aperture_at = surface$
  ele%value(mode$) = reflection$
  if (logic_option(.true., do_allocate)) then
    ! Avoid "ele%photon = photon_element_struct()" to get around Ifort bug. 4/10/2019
    if (associated(ele%photon)) deallocate(ele%photon)
    allocate(ele%photon)
  endif

case (taylor$)   ! start with unit matrix
  ele%tracking_method = taylor$  
  ele%mat6_calc_method = taylor$ 
  ele%taylor_map_includes_offsets = .false.
  if (logic_option(.true., do_allocate)) then
    call taylor_make_unit (ele%taylor)
  
    do i = 0, 3
      ele%spin_taylor(i)%ref = 0
      call init_taylor_series (ele%spin_taylor(i), 0)
    enddo
  endif

case (wiggler$, undulator$) 
  ele%field_calc = int_garbage$
  ele%value(polarity$) = 1.0
  ele%value(delta_ref_time_user_set$) = false$

end select

! %bookkeeping_state inits
! Note: Groups, for example, do not have a reference energy, etc. so set the bookkeeping
! state to OK$ for these categories.

call set_status_flags (ele%bookkeeping_state, stale$)

select case (ele%key)

case (group$)
  ele%bookkeeping_state%attributes     = ok$
  ele%bookkeeping_state%s_position     = ok$
  ele%bookkeeping_state%floor_position = ok$
  ele%bookkeeping_state%ref_energy     = ok$
  ele%bookkeeping_state%mat6           = ok$
  ele%bookkeeping_state%rad_int        = ok$
  ele%bookkeeping_state%ptc            = ok$
  ele%field_calc            = no_field$
  ele%value(gang$)          = true$
  ele%value(interpolation$) = cubic$

case (overlay$)
  ele%bookkeeping_state%attributes     = ok$
  ele%bookkeeping_state%mat6           = ok$
  ele%bookkeeping_state%rad_int        = ok$
  ele%bookkeeping_state%ptc            = ok$
  ele%field_calc            = no_field$
  ele%value(gang$)          = true$
  ele%value(interpolation$) = cubic$

case (ramper$)
  ele%bookkeeping_state%attributes     = ok$
  ele%bookkeeping_state%mat6           = ok$
  ele%bookkeeping_state%rad_int        = ok$
  ele%bookkeeping_state%ptc            = ok$
  ele%field_calc            = no_field$
  ele%value(interpolation$) = cubic$

case (girder$)
  ele%bookkeeping_state%attributes     = ok$
  ele%bookkeeping_state%ref_energy     = ok$
  ele%bookkeeping_state%mat6           = ok$
  ele%bookkeeping_state%rad_int        = ok$
  ele%bookkeeping_state%ptc            = ok$
  ele%field_calc = no_field$

case (beginning_ele$)
  ele%bookkeeping_state%rad_int        = ok$
  ele%bookkeeping_state%control        = ok$
  ele%field_calc = no_field$

case default
  ele%bookkeeping_state%control        = ok$

end select

end subroutine set_ele_defaults

