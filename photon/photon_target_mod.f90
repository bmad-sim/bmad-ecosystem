module photon_target_mod

use geometry_mod
use photon_utils_mod

contains

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine photon_target_setup (ele)
!
! Routine to calculate and store the parmeters needed for photon targeting.
!
! Input:
!   ele       -- ele_struct: Source element to setup.
!
! Output:
!   ele       -- ele_struct: Source element with target parameters setup.
!-

subroutine photon_target_setup (ele)

implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: ap_ele
type (photon_target_struct), pointer :: target

real(rp), pointer :: val(:)
real(rp) z
logical :: is_bending_element, follow_fork
character(*), parameter :: r_name = 'photon_target_setup '

! Init

if (ele%lord_status == super_lord$) then
  ap_ele => pointer_to_slave (ele, 1)
else
  ap_ele => ele
endif

is_bending_element = .false.

if (ap_ele%branch%param%particle == photon$) then
  follow_fork = .false.
else
  follow_fork = .true.
endif

! Find next element with an aperture

do 

  ap_ele => pointer_to_next_ele (ap_ele, skip_beginning = .true., follow_fork = follow_fork)

  if (ap_ele%value(x1_limit$) /= 0) exit

  select case (ap_ele%key)
  case (diffraction_plate$, crystal$, capillary$, mirror$, multilayer_mirror$, sample$, patch$)
    is_bending_element = .true.
  end select

  if (is_bending_element .or. ap_ele%ix_ele == ap_ele%branch%n_ele_track) then
    call out_io (s_fatal$, r_name, 'NO ELEMENT WITH APERTURE FOUND DOWNSTREAM FROM: ' // ele%name)
    if (global_com%exit_on_error) call err_exit
    return
  endif

enddo

! get aperture corners 

if (.not. associated(ele%photon)) allocate(ele%photon)
target => ele%photon%target
val => ap_ele%value


z = 0
if (stream_ele_end (ap_ele%aperture_at, ap_ele%orientation) == upstream_end$) z = -ap_ele%value(l$)

call photon_target_corner_calc (ap_ele,  0.0_rp,          0.0_rp,         z, ele, target%center)

call photon_target_corner_calc (ap_ele, -val(x1_limit$), -val(y1_limit$), z, ele, target%corner(1))
call photon_target_corner_calc (ap_ele, -val(x1_limit$),  val(y2_limit$), z, ele, target%corner(2))
call photon_target_corner_calc (ap_ele,  val(x2_limit$), -val(y1_limit$), z, ele, target%corner(3))
call photon_target_corner_calc (ap_ele,  val(x2_limit$),  val(y1_limit$), z, ele, target%corner(4))
target%n_corner = 4

! If there is surface curvature then the aperture rectangle becomes a 3D aperture box.

if (associated(ap_ele%photon)) then
  if (ap_ele%photon%surface%has_curvature .and. ap_ele%aperture_at == surface$) then
    z = z_at_surface(ele, -val(x1_limit$), -val(y1_limit$))
    call photon_target_corner_calc (ap_ele, -val(x1_limit$), -val(y1_limit$), z, ele, target%corner(5))

    z = z_at_surface(ele, -val(x1_limit$), val(y1_limit$))
    call photon_target_corner_calc (ap_ele, -val(x1_limit$),  val(y2_limit$), z, ele, target%corner(6))

    z = z_at_surface(ele, val(x1_limit$), -val(y1_limit$))
    call photon_target_corner_calc (ap_ele,  val(x2_limit$), -val(y1_limit$), z, ele, target%corner(7))

    z = z_at_surface(ele, val(x1_limit$), val(y1_limit$))
    call photon_target_corner_calc (ap_ele,  val(x2_limit$),  val(y1_limit$), z, ele, target%corner(8))
    target%n_corner = 8
  endif
endif

target%enabled = .true.

end subroutine photon_target_setup 

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine photon_target_corner_calc (aperture_ele, x_lim, y_lim, z_lim, source_ele, corner)
!
! Routine to calculate the corner coords in the source_ele ref frame.
!
! Input:
!   aperture_ele  -- ele_struct: Element containing the aperture
!   x_lim, y_lim  -- real(rp): Transverse corner points in aperture_ele coord frame.
!   source_ele    -- ele_struct: Photon source element.
!
! Output:
!   corner        -- target_point_struct: Corner coords in source_ele ref frame.
!-

subroutine photon_target_corner_calc (aperture_ele, x_lim, y_lim, z_lim, source_ele, corner)

implicit none

type (ele_struct), target :: aperture_ele, source_ele
type (ele_struct), pointer :: ele0
type (target_point_struct) corner
type (floor_position_struct) floor
type (coord_struct) orb

real(rp) x_lim, y_lim, z_lim

! Corner in aperture_ele coords

corner%r = [x_lim, y_lim, z_lim]

select case (stream_ele_end (aperture_ele%aperture_at, aperture_ele%orientation))
case (upstream_end$, downstream_end$)

case (surface$) 
  orb%vec = 0
  orb%vec(1:5:2) = corner%r
  call offset_photon (aperture_ele, orb, unset$, offset_position_only = .true.)
  corner%r = orb%vec(1:5:2)

case default
  call err_exit
end select

! Convert to floor coords and then to source_ele coords

floor = coords_relative_to_floor (aperture_ele%floor, corner%r)
ele0 => pointer_to_next_ele (source_ele, -1)
floor = coords_floor_to_relative (ele0%floor, floor, .false.)

orb%vec = 0
orb%vec(1:5:2) = floor%r
call offset_photon (source_ele, orb, set$, offset_position_only = .true.)

corner%r = orb%vec(1:5:2)

end subroutine photon_target_corner_calc

end module
