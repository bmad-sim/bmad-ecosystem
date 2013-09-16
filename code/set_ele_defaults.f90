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

! %value() inits

select case (ele%key)

case (beambeam$)
  ele%value(charge$) = -1

case (bend_sol_quad$) 
  ele%mat6_calc_method = symp_lie_bmad$
  ele%tracking_method  = symp_lie_bmad$

case (branch$, photon_branch$)
  ele%value(direction$) = 1
  ele%value(particle$) = real_garbage$
  ele%value(geometry$) = open$

case (crystal$)
  ele%value(ref_orbit_follows$) = bragg_diffracted$
  ele%value(ref_polarization$) = sigma_polarization$ 
  ele%aperture_at = surface$
  ele%offset_moves_aperture = .true.
  if (.not. associated(ele%surface)) allocate(ele%surface)
  ele%surface = photon_surface_struct()

case (custom$)  
  ele%mat6_calc_method = custom$
  ele%tracking_method  = custom$
  ele%field_calc       = custom$

case (ecollimator$)
  ele%aperture_type = elliptical$
  ele%offset_moves_aperture = .true.

case (fiducial$)
  ele%value(origin_ele_ref_pt$) = center_pt$

case (floor_shift$)
  ele%value(origin_ele_ref_pt$) = exit_end$

case (girder$)
  ele%value(origin_ele_ref_pt$) = center_pt$

case (lcavity$)
  ele%value(coupler_at$) = exit_end$
  ele%value(field_scale$) = 1
  ele%value(n_cell$) = 1
  ele%value(traveling_wave$) = 0

case (mirror$)
  ele%aperture_at = surface$
  ele%offset_moves_aperture = .true.
  if (.not. associated(ele%surface)) allocate(ele%surface)
  ele%surface = photon_surface_struct()

case (multilayer_mirror$)
  ele%value(ref_polarization$) = sigma_polarization$  
  ele%aperture_at = surface$
  ele%offset_moves_aperture = .true.
  if (.not. associated(ele%surface)) allocate(ele%surface)
  ele%surface = photon_surface_struct()

case (multipole$, ab_multipole$)
  call multipole_init (ele, .true.)

case (patch$)
  ele%value(flexible$) = 0 ! False
  ele%value(new_branch$) = 1    ! True
  ele%value(ref_coordinates$)= exit_end$

case (rbend$, sbend$)
  ele%value(fintx$) = real_garbage$
  ele%value(hgapx$) = real_garbage$
  ele%value(kill_fringe$) = no_end$
  ele%value(fringe_type$) = basic_bend$
  ele%value(ptc_field_geometry$) = sector$

case (rcollimator$)
  ele%offset_moves_aperture = .true.

case (rfcavity$)
  ele%value(coupler_at$) = exit_end$
  ele%value(field_scale$) = 1
  ele%value(n_cell$) = 1
  ele%value(traveling_wave$) = 0

case (taylor$)   ! start with unit matrix
  ele%tracking_method = taylor$  
  ele%mat6_calc_method = taylor$ 
  ele%map_with_offsets = .false.
  call taylor_make_unit (ele%taylor)

case (wiggler$, undulator$) 
  ele%sub_key = periodic_type$   
  ele%value(polarity$) = 1.0     

case (e_gun$)
  ele%tracking_method = time_runge_kutta$
  ele%mat6_calc_method = tracking$
  ele%value(field_scale$) = 1

end select

! Fringe set for non bend elements

if (ele%key /= sbend$ .and. ele%key /= rbend$ .and. attribute_index(ele, 'KILL_FRINGE') /= 0) THEN
  ele%value(kill_fringe$) = no_end$
  ele%value(fringe_type$) = none$
endif

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

case (init_ele$)
  ele%bookkeeping_state%rad_int        = ok$
  ele%bookkeeping_state%control        = ok$

case default
  ele%bookkeeping_state%control        = ok$

end select

end subroutine set_ele_defaults

