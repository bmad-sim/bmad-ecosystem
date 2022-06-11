!+
! Subroutine convert_particle_coordinates_s_to_t  (particle, s_body, orientation)
!
! Subroutine to convert particle coordinates from s-based to t-based system. 
!
! Note: t coordinates are:            
!     vec(1) = x                              [m]
!     vec(2) = c*p_x = m c^2 gamma beta_x     [eV]
!     vec(3) = y                              [m]
!     vec(4) = c*p_y = m c^2 gamma beta_y     [eV]
!     vec(5) = z                              [m]
!     vec(6) = c*p_s = m c^2 gamma beta_s     [eV]
!
! Input:
!   particle    -- coord_struct: Particle with %vec(:) in s-coords.
!   s_body      -- real(rp): s-position in element body coords.
!   orientation -- integer: ele%orientation for vec(6).
!
! Output:
!   particle    -- coord_struct: Particle with %vec(:) in t-coords.
!-

subroutine convert_particle_coordinates_s_to_t (particle, s_body, orientation)

use bmad_struct

implicit none

type (coord_struct), intent(inout), target :: particle
real(rp) s_body
real(rp), pointer :: vec(:)
integer :: orientation

! Convert s to t. vec(1) & vec(3) are unchanged

vec => particle%vec

if (particle%beta == 0) then
  particle%dt_ref = 0
else
  particle%dt_ref = -vec(5) / (c_light * particle%beta)
endif

vec(6) = particle%direction * orientation * particle%p0c * sqrt(((1+vec(6)))**2 - vec(2)**2 -vec(4)**2)
vec(2) = vec(2) * particle%p0c
vec(4) = vec(4) * particle%p0c
vec(5) = s_body

end subroutine convert_particle_coordinates_s_to_t

