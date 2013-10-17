module photon_target_mod

use lat_geometry_mod

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
!   ele       -- ele_struct: Element to setup.
!
! Output:
!   ele       -- ele_struct: Element with target parameters setup.
!-

subroutine photon_target_setup (ele)

implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: ap_ele
type (branch_struct), pointer :: branch
type (photon_target_struct), pointer :: target

real(rp), pointer :: val(:)
logical :: good
character(*), parameter :: r_name = 'photon_target_setup '

! Find next element with an aperture

branch => ele%branch

ap_ele => pointer_to_next_ele (ele)
good = .true.

do 
  if (ap_ele%value(x1_limit$) /= 0) exit

  select case (ap_ele%key)
  case (diffraction_plate$, crystal$, capillary$, mirror$, multilayer_mirror$, sample$)
    good = .false.
  end select

  if (ap_ele%ix_ele == branch%n_ele_track) good = .false.

  if (.not. good) then
    call out_io (s_fatal$, r_name, 'NO APERTURE ELEMENT FOUND DOWNSTREAM FROM: ' // ele%name)
    if (global_com%exit_on_error) call err_exit
    return
  endif

  ap_ele => pointer_to_next_ele (ap_ele)

enddo

! get aperture corners 

if (.not. associated(ele%photon)) allocate(ele%photon)
target => ele%photon%target
val => ap_ele%value

call photon_target_corner_calc (ap_ele, -val(x1_limit$), -val(y1_limit$), ele, target%corner(1))
call photon_target_corner_calc (ap_ele, -val(x1_limit$),  val(y2_limit$), ele, target%corner(2))
call photon_target_corner_calc (ap_ele,  val(x2_limit$), -val(y1_limit$), ele, target%corner(3))
call photon_target_corner_calc (ap_ele,  val(x2_limit$),  val(y1_limit$), ele, target%corner(4))
call photon_target_corner_calc (ap_ele,  0.0_rp,          0.0_rp,         ele, target%center)

target%enabled = .true.

end subroutine photon_target_setup 

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine photon_target_corner_calc (aperture_ele, x_lim, y_lim, source_ele, corner)
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

subroutine photon_target_corner_calc (aperture_ele, x_lim, y_lim, source_ele, corner)

implicit none

type (ele_struct), target :: aperture_ele, source_ele
type (ele_struct), pointer :: ele0
type (target_point_struct) corner
type (floor_position_struct) floor
type (coord_struct) orb

real(rp) x_lim, y_lim

! Corner in aperture_ele coords

corner%r = [x_lim, y_lim, 0.0_rp]

select case (stream_ele_end (aperture_ele%aperture_at, aperture_ele%orientation))
case (upstream_end$)
  corner%r(3) = -aperture_ele%value(l$)

case (downstream_end$)

case (surface$) 
  orb%vec = 0
  orb%vec(1:5:2) = corner%r
  call offset_photon (aperture_ele, orb, unset$, offset_position_only = .true.)
  corner%r = orb%vec(1:5:2)

case default
  call err_exit
end select

! Convert to floor coords and then to source_ele coords

floor = local_to_floor (aperture_ele%floor, corner%r)
ele0 => pointer_to_next_ele (source_ele, -1)
floor = floor_to_local (ele0%floor, floor, .false.)

orb%vec = 0
orb%vec(1:5:2) = floor%r
call offset_photon (source_ele, orb, set$, offset_position_only = .true.)

corner%r = orb%vec(1:5:2)

end subroutine photon_target_corner_calc

end module
