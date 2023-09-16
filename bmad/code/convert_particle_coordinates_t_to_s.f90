!+
! Subroutine convert_particle_coordinates_t_to_s (particle, ele, s_body, use_downstream_p0c)
!
! Subroutine to convert particle coordinates from t-based to s-based system. 
!
! Input:
!   particle            -- coord_struct: Particle with %vec(:) in t-coords.
!   ele                 -- ele_struct: Element particle is going through.
!   use_downstream_p0c  -- logical, optional: If True (the default), use ele%value(p0c$) as the reference momentum.
!                           If False, use ele%value(p0c_start$) as the reference.
!
! Output:
!   particle    -- coord_struct: Particle with %vec(:) in s-coords.
!   s_body      -- real(rp), optional: s-position in element body coords.
!-

subroutine convert_particle_coordinates_t_to_s (particle, ele, s_body, use_downstream_p0c)

use bmad_struct

implicit none

type (coord_struct), target :: particle
type (ele_struct) :: ele
real(rp), optional :: s_body
real(rp) :: p0c, pctot
real(rp), pointer :: vec(:)
logical, optional :: use_downstream_p0c

!

vec => particle%vec
pctot = sqrt (vec(2)**2 + vec(4)**2 + vec(6)**2)

if (present(s_body)) s_body = vec(5)

! If vec(6) = 0 then leave %direction as is.

if (vec(6)*ele%orientation > 0) then
  particle%direction = 1
elseif (vec(6)*ele%orientation < 0) then
  particle%direction = -1
endif

! Convert t to s. vec(1) and vec(3) are unchanged.

if (logic_option(.true., use_downstream_p0c)) then
  p0c = ele%value(p0c$)
else
  p0c = ele%value(p0c_start$)
endif

vec(2) = vec(2)/p0c
vec(4) = vec(4)/p0c
vec(5) = -c_light * particle%beta * particle%dt_ref
vec(6) = pctot/p0c - 1.0_rp
particle%p0c = p0c

end subroutine convert_particle_coordinates_t_to_s

