module photon_target_mod

use photon_utils_mod

implicit none

contains

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine photon_target_setup (ele)
!
! Routine to calculate and store the parmeters needed for photon targeting.
! This routine is called by Bmad parsing routines and is not meant for general use.
!
! Photon initialization with targeting is done by the routine init_a_photon_from_a_photon_init_ele
! Which is called by init_coord. 
!
! Input:
!   ele       -- ele_struct: Source element to setup. 
!                 Element will be of type: sample, diffraction_plate or photon_init.
!
! Output:
!   ele       -- ele_struct: Source element with target parameters setup in ele%photon%target.
!-

subroutine photon_target_setup (ele)

type (ele_struct), target :: ele
type (ele_struct), pointer :: ap_ele
type (photon_target_struct), pointer :: target
type (surface_grid_struct), pointer :: gr
type (branch_struct), pointer :: branch

real(rp), pointer :: val(:)
real(rp) z
logical :: is_bending_element, follow_fork, grid_defined, err_flag
character(*), parameter :: r_name = 'photon_target_setup'

! Init

if (.not. associated(ele%photon)) allocate(ele%photon)
target => ele%photon%target

if (ele%lord_status == super_lord$) then
  ap_ele => pointer_to_slave (ele, 1)
else
  ap_ele => ele
endif

branch => pointer_to_branch(ap_ele)

is_bending_element = .false.

if (branch%param%particle == photon$) then
  follow_fork = .false.
else
  follow_fork = .true.
endif

! Find next element with an aperture

do 
  ap_ele => pointer_to_next_ele (ap_ele, skip_beginning = .true., follow_fork = follow_fork)

  if (ap_ele%value(x1_limit$) /= 0 .or. ap_ele%key == detector$) exit

  select case (ap_ele%key)
  case (diffraction_plate$, mask$, crystal$, capillary$, mirror$, multilayer_mirror$, sample$)
    is_bending_element = .true.
  end select

  if (is_bending_element .or. ap_ele%ix_ele == branch%n_ele_track) then
    target%type = off$
    return
  endif
enddo

! Get aperture corners...
! Target info is stored in ele%photon%target so allocate ele%photon if needed.

grid_defined = .false.
if (associated(ap_ele%photon)) then
  gr => ap_ele%photon%surface%grid
  grid_defined = (gr%dr(1) /= 0 .or. gr%dr(2) /= 0)
endif

! If a grid has been defined use that as the target

if (grid_defined) then 
  target%type = grided$
  target%ele_loc = lat_ele_loc_struct(ap_ele%ix_ele, ap_ele%ix_branch)

  z = 0
  
  call photon_target_corner_calc (ap_ele, gr%r0(1),          gr%r0(2),          z, ele, target%center)
  call photon_target_corner_calc (ap_ele, gr%r0(1)+gr%dr(1), gr%r0(2),          z, ele, target%corner(1))
  call photon_target_corner_calc (ap_ele, gr%r0(1),          gr%r0(2)+gr%dr(2), z, ele, target%corner(2))

! If no grid defined use the element limits.

else

  target%type = rectangular$

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
      z = z_at_surface(ele, -val(x1_limit$), -val(y1_limit$), err_flag)
      call photon_target_corner_calc (ap_ele, -val(x1_limit$), -val(y1_limit$), z, ele, target%corner(5))

      z = z_at_surface(ele, -val(x1_limit$), val(y1_limit$), err_flag)
      call photon_target_corner_calc (ap_ele, -val(x1_limit$),  val(y2_limit$), z, ele, target%corner(6))

      z = z_at_surface(ele, val(x1_limit$), -val(y1_limit$), err_flag)
      call photon_target_corner_calc (ap_ele,  val(x2_limit$), -val(y1_limit$), z, ele, target%corner(7))

      z = z_at_surface(ele, val(x1_limit$), val(y1_limit$), err_flag)
      call photon_target_corner_calc (ap_ele,  val(x2_limit$),  val(y1_limit$), z, ele, target%corner(8))
      target%n_corner = 8
    endif
  endif

endif

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

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine photon_add_to_detector_statistics (orbit0, orbit, ele, ix_pt, iy_pt, grid_pt)
!
! Routine to add photon statistics to the appropriate pixel of a "detector" grid.
!
! Input:
!   orbit0  -- coord_struct: Photon coords at beginning of lattice
!   orbit   -- coord_struct: Photon coords at the detector.
!   ele     -- ele_struct: Element with grid.
!   grid_pt -- surface_grid_pt_struct, optional: If present then use this grid point
!                 instead of the grid point determined by the (x, y) coords of the photon
!   
!
! Output:
!   ele           -- ele_struct: Element with updatted grid.
!   ix_pt, iy_pt  -- integer, optional: Index of upgraded ele%photon%surface%grid%pt(:,:) point.
!                       These arguments are not set if the grid_pt argument is present.
!-

subroutine photon_add_to_detector_statistics (orbit0, orbit, ele, ix_pt, iy_pt, grid_pt)

type (coord_struct) orbit0, orbit, orb
type (ele_struct), target :: ele
type (surface_grid_pt_struct), optional, target :: grid_pt
type (surface_grid_struct), pointer :: grid
type (surface_grid_pt_struct), pointer :: pix
type (branch_struct), pointer :: branch

real(rp) phase, intens_x, intens_y, intens, dE
integer, optional :: ix_pt, iy_pt
integer nx, ny

!  Convert to detector and then to angle coords.

orb = to_photon_angle_coords (orbit, ele, .true.)

! Find grid pt to update.
! Note: dr(i) can be zero for 1-dim grid

if (present(grid_pt)) then
  pix => grid_pt

else
  grid => ele%photon%surface%grid

  nx = 0; ny = 0
  if (grid%dr(1) /= 0) nx = nint((orb%vec(1) - grid%r0(1)) / grid%dr(1))
  if (grid%dr(2) /= 0) ny = nint((orb%vec(3) - grid%r0(2)) / grid%dr(2))

  if (present(ix_pt)) ix_pt = nx
  if (present(iy_pt)) iy_pt = ny

  ! If outside of detector region then do nothing.

  if (nx < lbound(grid%pt, 1) .or. nx > ubound(grid%pt, 1) .or. &
      ny < lbound(grid%pt, 2) .or. ny > ubound(grid%pt, 2)) return

  pix => grid%pt(nx,ny)
endif

! Add to det stat

branch => pointer_to_branch(ele)
pix%n_photon  = pix%n_photon + 1
if (branch%lat%photon_type == coherent$) then
  phase = orb%phase(1) 
  pix%E_x = pix%E_x + orb%field(1) * cmplx(cos(phase), sin(phase), rp)
  phase = orb%phase(2) 
  pix%E_y = pix%E_y + orb%field(2) * cmplx(cos(phase), sin(phase), rp)
else
  intens_x = orbit%field(1)**2
  intens_y = orbit%field(2)**2
  intens = intens_x + intens_y
  pix%intensity_x     = pix%intensity_x    + intens_x
  pix%intensity_y     = pix%intensity_y    + intens_y
  pix%intensity       = pix%intensity      + intens
  pix%orbit           = pix%orbit          + intens * orb%vec
  pix%orbit_rms       = pix%orbit_rms      + intens * orb%vec**2

  orb = to_photon_angle_coords (orbit0, ele, .false.)
  pix%init_orbit      = pix%init_orbit     + intens * orb%vec
  pix%init_orbit_rms  = pix%init_orbit_rms + intens * orb%vec**2
endif

end subroutine photon_add_to_detector_statistics

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Function to_photon_angle_coords (orb_in, ele, to_ele_coords) result (orb_out)
!
! Routine to convert from standard photon coords to "angle" coords defined as:
!       x, angle_x, y, angle_y, z, E-E_ref
!
! Input:
!   orb_in        -- coord_struct: orbit in standard photon coords.
!   ele           -- ele_struct: Reference element (generally the detector element.)
!   to_ele_coords -- logical: Transform from lab to element coordinates?
!
! Output:
!   orb_out       -- coord_struct: Transformed coordinates.
!-

function to_photon_angle_coords (orb_in, ele, to_ele_coords) result (orb_out)

type (coord_struct) orb_in, orb_out
type (ele_struct) ele
logical to_ele_coords

! To element coords?

orb_out = orb_in
if (to_ele_coords) call offset_photon (ele, orb_out, set$)  ! Go to coordinates of the detector

! To angle coords

orb_out%vec(2) = atan2(orb_out%vec(2), orb_out%vec(6))
orb_out%vec(4) = atan2(orb_out%vec(4), orb_out%vec(6))
orb_out%vec(6) = orb_out%p0c - ele%value(E_tot$)

end function to_photon_angle_coords

end module
