!+
! Subroutine convert_particle_coordinates_t_to_s (particle, ele, s_body)
!
! Subroutine to convert particle coordinates from t-based to s-based system. 
!
! Input:
!   particle    -- coord_struct: Particle with %vec(:) in t-coords.
!   ele         -- ele_sturct, optional: Element particle is going through.
!
! Output:
!   particle    -- coord_struct: Particle with %vec(:) in s-coords.
!   s_body      -- real(rp), optional: s-position in element body coords.
!-

subroutine convert_particle_coordinates_t_to_s (particle, ele, s_body)

use bmad_struct

implicit none

type (coord_struct), intent(inout), target :: particle
type (ele_struct), optional :: ele
real(rp), optional :: s_body
real(rp) :: p0c, pctot
real(rp), pointer :: vec(:)
integer edir

!

vec => particle%vec
p0c = particle%p0c
pctot = sqrt (vec(2)**2 + vec(4)**2 + vec(6)**2)

if (present(s_body)) s_body = vec(5)

! If vec(6) = 0 then leave %direction as is.

edir = 1
if (present(ele)) edir = ele%orientation

if (vec(6)*edir > 0) then
  particle%direction = 1
elseif (vec(6)*edir < 0) then
  particle%direction = -1
endif

! Convert t to s. vec(1) and vec(3) are unchanged.

vec(2) = vec(2)/p0c
vec(4) = vec(4)/p0c
vec(5) = -c_light * particle%beta * particle%dt_ref
vec(6) = pctot/p0c - 1.0_rp

end subroutine convert_particle_coordinates_t_to_s

