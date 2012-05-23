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

use bmad_struct

implicit none

type (ele_struct) ele

!

select case (ele%key)

case (bend_sol_quad$) 
  ele%mat6_calc_method = symp_lie_bmad$
  ele%tracking_method  = symp_lie_bmad$

case (branch$, photon_branch$)
  ele%value(direction$) = 1
  ele%value(particle$) = real_garbage$
  ele%value(lattice_type$) = linear_lattice$

case (crystal$)
  ele%value(follow_diffracted_beam$) = 1  ! True
  ele%value(ref_polarization$) = sigma_polarization$ 
  ele%aperture_at = surface$
  ele%offset_moves_aperture = .true.

case (custom$)  
  ele%mat6_calc_method = custom$
  ele%tracking_method  = custom$
  ele%field_calc       = custom$

case (ecollimator$)
  ele%aperture_type = elliptical$
  ele%offset_moves_aperture = .true.

case (lcavity$)
  ele%value(coupler_at$) = exit_end$
  ele%value(field_scale$) = 1
  ele%value(n_cell$) = 1

case (mirror$)
  ele%aperture_at = surface$
  ele%offset_moves_aperture = .true.

case (multilayer_mirror$)
  ele%value(ref_polarization$) = sigma_polarization$  
  ele%aperture_at = surface$
  ele%offset_moves_aperture = .true.

case (rbend$, sbend$)
  ele%value(fintx$) = real_garbage$
  ele%value(hgapx$) = real_garbage$

case (rcollimator$)
  ele%offset_moves_aperture = .true.

case (rfcavity$)
  ele%value(coupler_at$) = exit_end$
  ele%value(field_scale$) = 1
  ele%value(n_cell$) = 1

case (taylor$)   ! start with unit matrix
  ele%tracking_method = taylor$  
  ele%mat6_calc_method = taylor$ 
  ele%map_with_offsets = .false.
  call taylor_make_unit (ele%taylor)

case (wiggler$) 
  ele%sub_key = periodic_type$   
  ele%value(polarity$) = 1.0     

case (e_gun$)
  ele%tracking_method = time_runge_kutta$
  ele%mat6_calc_method = tracking$
  ele%value(field_scale$) = 1

end select

end subroutine set_ele_defaults

