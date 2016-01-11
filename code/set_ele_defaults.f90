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

! Default fringe set for non bend elements

if (attribute_index(ele, 'FRINGE_AT') /= 0)        ele%value(fringe_at$) = both_ends$
if (attribute_index(ele, 'FRINGE_TYPE') /= 0)      ele%value(fringe_type$) = none$
if (attribute_index(ele, 'SPIN_FRINGE_ON') /= 0)   ele%value(spin_fringe_on$) = true$

! %value() inits

select case (ele%key)

case (beambeam$)
  ele%value(charge$) = -1
  ele%value(n_slice$) = 1

case (beginning_ele$)
  ele%value(e_tot$) = -1
  ele%value(p0c$) = -1
  ele%value(spinor_polarization$) = 1
  ele%value(spin_z$) = 1

case (bend_sol_quad$) 
  ele%mat6_calc_method = symp_lie_bmad$
  ele%tracking_method  = symp_lie_bmad$

case (fork$, photon_fork$)
  ele%value(direction$) = 1
  ele%value(particle$) = real_garbage$
  ele%value(geometry$) = open$

case (capillary$)
  ele%offset_moves_aperture = .true.
  
case (crystal$)
  ele%value(ref_orbit_follows$) = bragg_diffracted$
  ele%aperture_at = surface$
  ele%offset_moves_aperture = .true.
  if (logic_option(.true., do_allocate)) then
    if (.not. associated(ele%photon)) allocate(ele%photon)
    !!! Due to ifort bug cannot do:  ele%photon = photon_element_struct()
    call init_photon_element_struct(ele%photon)
  endif

case (custom$)  
  ele%mat6_calc_method = custom$
  ele%tracking_method  = custom$
  ele%field_calc       = custom$

case (def_beam_start$)
  ele%value(spinor_polarization$) = 1
  ele%value(spin_z$) = 1

case (def_mad_beam$)
  ele%value(particle$) = positron$

case (def_parameter$)
  ele%value(geometry$) = -1
  ele%value(particle$) = positron$
  ele%value(default_tracking_species$) = real_garbage$

case (detector$)
  ele%aperture_type = auto_aperture$
  if (logic_option(.true., do_allocate)) then
    if (.not. associated(ele%photon)) allocate(ele%photon)
    !!! Due to ifort bug cannot do:  ele%photon = photon_element_struct()
    call init_photon_element_struct(ele%photon)
  endif

case (diffraction_plate$)
  ele%aperture_at = surface$
  ele%aperture_type = auto_aperture$
  ele%offset_moves_aperture = .true.
  ele%value(geometry$) = transmission$
  if (logic_option(.true., do_allocate)) then
    if (.not. associated(ele%photon)) allocate(ele%photon)
    !!! Due to ifort bug cannot do: ele%photon = photon_element_struct()
    call init_photon_element_struct(ele%photon)
  endif

case (e_gun$)
  ele%tracking_method = time_runge_kutta$
  ele%mat6_calc_method = tracking$
  ele%value(field_factor$) = 1
  ele%value(fringe_at$) = exit_end$
  ele%value(fringe_type$) = full$
  ele%value(autoscale_amplitude$) = true$
  ele%value(autoscale_phase$) = true$

case (ecollimator$)
  ele%aperture_type = elliptical$
  ele%offset_moves_aperture = .true.

case (em_field$)
  ele%tracking_method = runge_kutta$
  ele%mat6_calc_method = tracking$
  ele%value(fringe_type$) = full$
  ele%value(field_factor$) = 1

case (fiducial$)
  ele%value(origin_ele_ref_pt$) = center_pt$

case (floor_shift$)
  ele%value(origin_ele_ref_pt$) = exit_end$
  ele%value(upstream_ele_dir$) = 1
  ele%value(downstream_ele_dir$) = 1

case (girder$)
  ele%value(origin_ele_ref_pt$) = center_pt$

case (hybrid$)   ! start with unit matrix
  if (logic_option(.true., do_allocate)) then
    do i = 1, 3
      ele%taylor(i)%ref = 0
      call init_taylor_series (ele%taylor(i), 0)
    enddo
    do i = 1, 3; do j = 1, 3
      ele%spin_taylor(i,j)%ref = 0
      call init_taylor_series (ele%spin_taylor(i,j), 0)
    enddo; enddo
  endif

case (lcavity$)
  ele%value(coupler_at$) = exit_end$
  ele%value(field_factor$) = 1
  ele%value(n_cell$) = 1
  ele%value(cavity_type$) = standing_wave$
  ele%value(fringe_type$) = full$
  ele%value(autoscale_amplitude$) = true$
  ele%value(autoscale_phase$) = true$

case (line_ele$)
  ele%value(particle$) = real_garbage$
  ele%value(geometry$) = real_garbage$
  ele%value(default_tracking_species$) = real_garbage$
  ele%value(e_tot$) = -1
  ele%value(p0c$) = -1

case (mask$)
  ele%aperture_at = surface$
  ele%aperture_type = auto_aperture$
  ele%offset_moves_aperture = .true.
  ele%value(geometry$) = transmission$
  if (logic_option(.true., do_allocate)) then
    if (.not. associated(ele%photon)) allocate(ele%photon)
    !!! Due to ifort bug cannot do: ele%photon = photon_element_struct()
    call init_photon_element_struct(ele%photon)
  endif

case (mirror$)
  ele%aperture_at = surface$
  ele%offset_moves_aperture = .true.
  if (logic_option(.true., do_allocate)) then
    if (.not. associated(ele%photon)) allocate(ele%photon)
    !!! Due to ifort bug cannot do:  ele%photon = photon_element_struct()
    call init_photon_element_struct(ele%photon)
  endif

case (multilayer_mirror$)
  ele%aperture_at = surface$
  ele%offset_moves_aperture = .true.
  if (logic_option(.true., do_allocate)) then
    if (.not. associated(ele%photon)) allocate(ele%photon)
    !!! Due to ifort bug cannot do:  ele%photon = photon_element_struct()
    call init_photon_element_struct(ele%photon)
  endif

case (multipole$, ab_multipole$)
  if (logic_option(.true., do_allocate)) then
    call multipole_init (ele, .true.)
  endif
  ele%scale_multipoles = .false.

case (patch$)
  ele%value(flexible$) = false$ 
  ele%value(new_branch$) = true$
  ele%value(ref_coordinates$)= exit_end$
  ele%value(upstream_ele_dir$) = 1
  ele%value(downstream_ele_dir$) = 1

case (photon_init$)
  ele%value(ds_slice$) = 0.01
  ele%value(velocity_distribution$) = gaussian$
  ele%value(energy_distribution$) = gaussian$
  ele%value(spatial_distribution$) = gaussian$
  ele%value(transverse_sigma_cut$) = 3
  ele%value(E_center_relative_to_ref$) = true$
  if (logic_option(.true., do_allocate)) then
    if (.not. associated(ele%photon)) allocate(ele%photon)
    !!! Due to ifort bug cannot do:  ele%photon = photon_element_struct()
    call init_photon_element_struct(ele%photon)
  endif

case (rbend$, sbend$)
  ele%value(fintx$) = real_garbage$
  ele%value(hgapx$) = real_garbage$
  ele%value(fringe_type$) = basic_bend$
  ele%value(higher_order_fringe_type$) = none$
  ele%value(ptc_field_geometry$) = sector$
  ele%value(ptc_fringe_geometry$) = x_invariant$

case (rcollimator$)
  ele%offset_moves_aperture = .true.

case (rfcavity$)
  ele%value(coupler_at$) = exit_end$
  ele%value(field_factor$) = 1
  ele%value(n_cell$) = 1
  ele%value(cavity_type$) = standing_wave$
  ele%value(fringe_type$) = full$
  ele%value(autoscale_amplitude$) = true$
  ele%value(autoscale_phase$) = true$

case (sad_mult$)
  ele%value(eps_step_scale$) = 1
  ele%scale_multipoles = .false.
  if (logic_option(.true., do_allocate)) then
    call multipole_init (ele, .true.)
  endif

case (sample$)
  ele%aperture_at = surface$
  ele%value(geometry$) = reflection$
  if (logic_option(.true., do_allocate)) then
    if (.not. associated(ele%photon)) allocate(ele%photon)
    !!! Due to ifort bug cannot do:  ele%photon = photon_element_struct()
    call init_photon_element_struct(ele%photon)
  endif

case (taylor$)   ! start with unit matrix
  ele%tracking_method = taylor$  
  ele%mat6_calc_method = taylor$ 
  ele%taylor_map_includes_offsets = .false.
  if (logic_option(.true., do_allocate)) then
    call taylor_make_unit (ele%taylor)
  
    do i = 1, 3; do j = 1, 3
      ele%spin_taylor(i,j)%ref = 0
      call init_taylor_series (ele%spin_taylor(i,j), 0)
    enddo; enddo
  endif

case (wiggler$, undulator$) 
  ele%sub_key = periodic_type$   
  ele%value(polarity$) = 1.0     

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
  ele%field_calc = no_field$

case (overlay$)
  ele%bookkeeping_state%attributes     = ok$
  ele%bookkeeping_state%mat6           = ok$
  ele%bookkeeping_state%rad_int        = ok$
  ele%bookkeeping_state%ptc            = ok$
  ele%field_calc = no_field$

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

!------------------------------
contains

subroutine init_photon_element_struct(photon_element)
type (photon_element_struct) photon_element

photon_element%surface = photon_surface_struct()
photon_element%target  = photon_target_struct() 

end subroutine init_photon_element_struct

end subroutine set_ele_defaults

