!+
! Subroutine init_lrbbi(lat, oppos_lat, lrbbi_ele, ix_lrbbi, ix_oppos)
!
! Subroutine to calculate the basic parameters of a beambeam element. 
! Initializes the element and establishes the following values:
!     key, value(n_slice), value(charge), value(current -- provided 
!     lat%param%n_part has been specified), x_offset, y_offset,  
!     sigma_x, sigma_y, sigma_z. 
! Names the element 'LRBBI'.
! DOES NOT establish any other values, including %s (s_positions)!!
!                                                                     
!
! Modules needed:
!   use bmad
!
! Input:
!   lat            -- lat_struct: lat in which LRBBI's are inserted
!   oppos_lat      -- Lat_sturct: lat for the oppositely circulating beam
!   LRBBI_ele       -- Ele_struct: ele_struct for the lrbbi element
!   ix_LRBBI        -- Integer: index of LRBBI element in lat (probably the 
!                        same as ix_split from split_lat since the
!                        LRBBI element has zero length)
!   ix_oppos
!
! Output:
!   LRBBI_ele       -- Ele_struct: updated lrbbi structure.
!-  

#include "CESR_platform.inc"

subroutine init_LRBBI(lat, oppos_lat, LRBBI_ele, ix_LRBBI, ix_oppos)

  use bmad_struct
  use bmad_interface, except_dummy => init_LRBBI
  
  implicit none

  type (lat_struct)  lat
  type (lat_struct) oppos_lat
  type (ele_struct)     LRBBI_ele
  integer, intent(in) :: ix_LRBBI, ix_oppos 

  type (normal_modes_struct)  mode
  type (coord_struct), allocatable, save ::  orbit(:)
  type (coord_struct), allocatable, save ::  oppos_orbit(:)
  
  real(rp) :: beta_a, beta_b, eta_a, eta_b, sigma_z
  real(rp) :: sigma_x(lat%n_ele_max), sigma_y(lat%n_ele_max)

  real(rp) :: e_spread, emit_x, emit_y, current
  integer :: i

  character(40) :: call_it

! init

  call init_ele(LRBBI_ele)

  call_it = 'lrbbi'
  LRBBI_ele%name = call_it
  LRBBI_ele%key = beambeam$
  LRBBI_ele%value(charge$) = -1
  LRBBI_ele%value(n_slice$) = 1

  current = lat%param%n_part*e_charge*c_light/(lat%param%total_length)

  call reallocate_coord (orbit, lat%n_ele_max)
  call reallocate_coord (oppos_orbit, lat%n_ele_max)

!

  call twiss_at_start(lat)
  call twiss_propagate_all(lat)
  call closed_orbit_calc(lat, orbit, 4)

! Get additional needed parameter values.

  call radiation_integrals(lat, orbit, mode)

  e_spread = mode%sigE_E
  sigma_z = mode%sig_z
  emit_x = mode%a%emittance !a=horizontal-like. Use %b% for vertical-like.
  emit_y = emit_x * 0.01

! Calculate sigmas for the beambeam element.

  do i=1, lat%n_ele_track
    beta_a = lat%ele(i)%a%beta
    beta_b = lat%ele(i)%b%beta
    eta_a = lat%ele(i)%a%eta
    eta_b = lat%ele(i)%b%eta
    sigma_x(i) = sqrt(emit_x * beta_a + (eta_a**2) * ((e_spread)**2))
    sigma_y(i) = sqrt(emit_y * beta_b + (eta_b**2) *((e_spread)**2))
  enddo

  LRBBI_ele%value(sig_x$) = sigma_x(ix_LRBBI)
  LRBBI_ele%value(sig_y$) = sigma_y(ix_LRBBI)
  LRBBI_ele%value(sig_z$) = sigma_z

  call twiss_at_start(oppos_lat)
  call twiss_propagate_all(oppos_lat)
  call closed_orbit_calc(oppos_lat, oppos_orbit, 4)

  LRBBI_ele%value(x_offset$) = oppos_orbit(ix_oppos)%vec(1)
  LRBBI_ele%value(y_offset$) = oppos_orbit(ix_oppos)%vec(3)

end subroutine
