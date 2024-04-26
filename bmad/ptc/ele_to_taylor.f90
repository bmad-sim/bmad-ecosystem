!+
! Subroutine ele_to_taylor (ele, param, orb0, taylor_map_includes_offsets, include_damping, orbital_taylor, spin_taylor)
!
! Subroutine to make orbital and spin (if spin tracking is on) taylor maps for an element. 
! The order of the map is set by set_ptc
!
! Input:
!   ele               -- Element_struct: Element to construct map for.
!   orb0              -- Coord_struct, optional: Starting coords around which the Taylor map is evaluated.
!                         Default is the zero orbit.
!   param             -- lat_param_struct: 
!   taylor_map_includes_offsets 
!                     -- Logical, optional: If present then value overrides ele%taylor_map_includes_offsets.
!   include_damping   -- logical, optional: Sets if radiation damping is included. Default is what is set in ptc_private%base_state.
!
! Output:
!   orbital_taylor(6) -- taylor_struct, optional: Orbital taylor map.
!                         If not present then the map is put in ele%taylor.
!   spin_taylor(0:3)  -- taylor_struct, optional: Spin taylor map. 
!                         If not present then the map is put in ele%spin_taylor.
!-

subroutine ele_to_taylor (ele, param, orb0, taylor_map_includes_offsets, include_damping, orbital_taylor, spin_taylor)

use ptc_interface_mod, dummy => ele_to_taylor, dummy2 => dp
use s_tracking
use mad_like, only: real_8, fibre, ring_l, survey, CONVERSION_XPRIME_IN_ABELL, internal_state
use ptc_spin, only: track_probe_x, track_probe
use ptc_multiparticle, only: survey
use madx_ptc_module, only: bmadl

implicit none

type (ele_struct), target :: ele
type (lat_param_struct) :: param
type (coord_struct), optional, intent(in) :: orb0
type (coord_struct) c0
type (taylor_struct), optional, target :: orbital_taylor(6), spin_taylor(0:3)
type (taylor_struct), pointer :: orb_tylr(:), spin_tylr(:)
type (probe) ptc_probe
type (probe_8) ptc_probe8
type (fibre), pointer :: ptc_fibre
type (real_8) y2(6)
type (c_damap) ptc_cdamap
type (internal_state) ptc_state

real(dp) x(6), beta
integer i, print12

logical, optional :: taylor_map_includes_offsets, include_damping
logical use_offsets, err_flag

character(16) :: r_name = 'ele_to_taylor'

!

CONVERSION_XPRIME_IN_ABELL = (.not. bmad_com%convert_to_kinetic_momentum) ! Only affects cylindrical map eles

if (present(orbital_taylor)) then
  orb_tylr => orbital_taylor
else
  orb_tylr => ele%taylor
endif

if (present(spin_taylor)) then
  spin_tylr => spin_taylor
else
  spin_tylr => ele%spin_taylor
endif

ptc_state = ptc_private%base_state
if (present(include_damping)) then
  select case (include_damping)
  case (.true.);  ptc_state = ptc_state + radiation0
  case default;   ptc_state = ptc_state - radiation0
  end select
endif

! Match elements and helical wiggler/undulators without a map are not implemented in PTC so just use the matrix.

if (ele%key == match$ .or. &
            (.not. associated(ele%cylindrical_map) .and. .not. associated(ele%cartesian_map) .and. &
            .not. associated(ele%gen_grad_map) .and. &
            (ele%key == wiggler$ .or. ele%key == undulator$) .and. ele%field_calc == helical_model$)) then
  call mat6_to_taylor (ele%vec0, ele%mat6, orb_tylr)
  call taylor_make_quaternion_unit (spin_tylr)
  if (.not. present(spin_taylor)) ele%spin_taylor_ref_orb_in = 0
  return
endif

! Init. 

call ptc_set_taylor_order_if_needed()

use_offsets = logic_option(ele%taylor_map_includes_offsets, taylor_map_includes_offsets)

call ele_to_fibre (ele, ptc_fibre, param, use_offsets, err_flag, ref_in = orb0)
if (err_flag) return

call alloc(ptc_cdamap)
call alloc(ptc_probe8)

! Initial map

if (present(orb0)) then
  orb_tylr(:)%ref = orb0%vec
  x = orb0%vec
else
  orb_tylr(:)%ref = 0
  x = 0
endif

ptc_probe = 0
ptc_probe = x

ptc_cdamap = 1
ptc_probe8 = ptc_cdamap + ptc_probe ! = IdentityMap + const

! It does not make sense to do spin tracking if spin_taylor is not present but orbital_taylor is.

if (bmad_com%spin_tracking_on .and. (present(spin_taylor) .or. (.not. present(spin_taylor) .and. .not. present(orbital_taylor)))) then
  call track_probe (ptc_probe8, ptc_state+SPIN0, fibre1 = bmadl%start)

  do i = 0, 3
    spin_tylr(i) = ptc_probe8%q%x(i)%t
  enddo

  if (.not. present(spin_taylor)) ele%spin_taylor_ref_orb_in = x

else
  call track_probe (ptc_probe8, ptc_state-SPIN0, fibre1 = bmadl%start)
endif

! take out the offset

!if (any(x /= 0)) then
!  call alloc(y2)
!  y2 = -x  ! y2 = IdentityMap - x
!  call concat_real_8 (y2, y0, y0)
!  call kill(y2)
!endif

! convert to orb_tylr_struct

orb_tylr = ptc_probe8%x

call kill (ptc_probe8)
call kill (ptc_cdamap)

call set_ele_status_stale (ele, mat6_group$)

CONVERSION_XPRIME_IN_ABELL = .true. ! Reset to normal.

end subroutine ele_to_taylor
