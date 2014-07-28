!+
! Subroutine set_ele_defaults (ele)
!
! Subroutine set the defaults for an element of a given type.
! For example, the default aperture type for an ecollimator$
!   element is ele%aperture_type = elliptical$.
!
! Input:
!   ele   -- ele_struct: Element to init.
!     %key -- Type of element.
!
! Output:
!   ele   -- ele_struct: Initialized element.
!-

subroutine set_ele_defaults (ele)

use bmad_interface, dummy => set_ele_defaults

implicit none

type (ele_struct) ele

! Default fringe set for non bend elements

if (attribute_index(ele, 'FRINGE_AT') /= 0)   ele%value(fringe_at$) = both_ends$
if (attribute_index(ele, 'FRINGE_TYPE') /= 0) ele%value(fringe_type$) = none$

! %value() inits

select case (ele%key)

case (beambeam$)
  ele%value(charge$) = -1

case (beginning_ele$)
  ele%value(e_tot$) = -1
  ele%value(p0c$) = -1

case (bend_sol_quad$) 
  ele%mat6_calc_method = symp_lie_bmad$
  ele%tracking_method  = symp_lie_bmad$

case (fork$, photon_fork$)
  ele%value(direction$) = 1
  ele%value(particle$) = real_garbage$
  ele%value(geometry$) = open$

case (crystal$)
  ele%value(ref_orbit_follows$) = bragg_diffracted$
  ele%aperture_at = surface$
  ele%offset_moves_aperture = .true.
  if (.not. associated(ele%photon)) allocate(ele%photon)
!!! Due to ifort bug:  ele%photon = photon_element_struct()
  call init_photon_element_struct(ele%photon)

case (custom$)  
  ele%mat6_calc_method = custom$
  ele%tracking_method  = custom$
  ele%field_calc       = custom$

case (def_beam$)
  ele%value(particle$) = positron$

case (def_parameter$)
  ele%value(geometry$) = -1
  ele%value(particle$) = positron$  

case (detector$)
  if (.not. associated(ele%photon)) allocate(ele%photon)
!!! Due to ifort bug:  ele%photon = photon_element_struct()
  call init_photon_element_struct(ele%photon)

case (diffraction_plate$)
  ele%aperture_at = surface$
  ele%aperture_type = auto_aperture$
  ele%value(geometry$) = transmission$
  if (.not. associated(ele%photon)) allocate(ele%photon)
!!! Due to ifort bug:  ele%photon = photon_element_struct()
  call init_photon_element_struct(ele%photon)

case (e_gun$)
  ele%tracking_method = time_runge_kutta$
  ele%mat6_calc_method = tracking$
  ele%value(field_factor$) = 1
  ele%value(fringe_at$) = exit_end$
  ele%value(fringe_type$) = full_straight$

case (ecollimator$)
  ele%aperture_type = elliptical$
  ele%offset_moves_aperture = .true.

case (em_field$)
  ele%value(fringe_type$) = full_straight$

case (fiducial$)
  ele%value(origin_ele_ref_pt$) = center_pt$

case (floor_shift$)
  ele%value(origin_ele_ref_pt$) = exit_end$
  ele%value(upstream_ele_dir$) = 1
  ele%value(downstream_ele_dir$) = 1

case (girder$)
  ele%value(origin_ele_ref_pt$) = center_pt$

case (lcavity$)
  ele%value(coupler_at$) = exit_end$
  ele%value(field_factor$) = 1
  ele%value(n_cell$) = 1
  ele%value(traveling_wave$) = false$
  ele%value(fringe_type$) = full_straight$

case (line_ele$)
  ele%value(particle$) = real_garbage$
  ele%value(geometry$) = real_garbage$
  ele%value(rel_tracking_charge$) = real_garbage$
  ele%value(e_tot$) = -1
  ele%value(p0c$) = -1

case (mirror$)
  ele%aperture_at = surface$
  ele%offset_moves_aperture = .true.
  if (.not. associated(ele%photon)) allocate(ele%photon)
!!! Due to ifort bug:  ele%photon = photon_element_struct()
  call init_photon_element_struct(ele%photon)

case (multilayer_mirror$)
  ele%aperture_at = surface$
  ele%offset_moves_aperture = .true.
  if (.not. associated(ele%photon)) allocate(ele%photon)
!!! Due to ifort bug:  ele%photon = photon_element_struct()
  call init_photon_element_struct(ele%photon)

case (multipole$, ab_multipole$)
  call multipole_init (ele, .true.)

case (patch$)
  ele%value(flexible$) = false$ 
  ele%value(new_branch$) = true$
  ele%value(ref_coordinates$)= exit_end$
  ele%value(upstream_ele_dir$) = 1
  ele%value(downstream_ele_dir$) = 1

case (rbend$, sbend$)
  ele%value(fintx$) = real_garbage$
  ele%value(hgapx$) = real_garbage$
  ele%value(fringe_type$) = basic_bend$
  ele%value(ptc_field_geometry$) = sector$

case (rcollimator$)
  ele%offset_moves_aperture = .true.

case (rfcavity$)
  ele%value(coupler_at$) = exit_end$
  ele%value(field_factor$) = 1
  ele%value(n_cell$) = 1
  ele%value(traveling_wave$) = false$
  ele%value(fringe_type$) = full_straight$

case (sad_mult$)
  ele%value(eps_step_scale$) = 1
  call multipole_init (ele, .true.)

case (sample$)
  ele%aperture_at = surface$
  if (.not. associated(ele%photon)) allocate(ele%photon)
!!! Due to ifort bug:  ele%photon = photon_element_struct()
  call init_photon_element_struct(ele%photon)
  ele%value(geometry$) = reflection$

case (taylor$)   ! start with unit matrix
  ele%tracking_method = taylor$  
  ele%mat6_calc_method = taylor$ 
  ele%taylor_map_includes_offsets = .false.
  call taylor_make_unit (ele%taylor)

case (wiggler$, undulator$) 
  ele%sub_key = periodic_type$   
  ele%value(polarity$) = 1.0     

case (x_ray_init$)
  if (.not. associated(ele%photon)) allocate(ele%photon)
!!! Due to ifort bug:  ele%photon = photon_element_struct()
  call init_photon_element_struct(ele%photon)

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

case (overlay$)
  ele%bookkeeping_state%attributes     = ok$
  ele%bookkeeping_state%mat6           = ok$
  ele%bookkeeping_state%rad_int        = ok$
  ele%bookkeeping_state%ptc            = ok$

case (girder$)
  ele%bookkeeping_state%attributes     = ok$
  ele%bookkeeping_state%ref_energy     = ok$
  ele%bookkeeping_state%mat6           = ok$
  ele%bookkeeping_state%rad_int        = ok$
  ele%bookkeeping_state%ptc            = ok$

case (beginning_ele$)
  ele%bookkeeping_state%rad_int        = ok$
  ele%bookkeeping_state%control        = ok$

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

