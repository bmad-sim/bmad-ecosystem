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
type (photon_element_struct), pointer :: ph
type (branch_struct), pointer :: branch

real(rp), pointer :: val(:)
real(rp) z, x_lim(2), y_lim(2)
logical :: is_bending_element, follow_fork, err_flag
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

val => ap_ele%value
x_lim = [-val(x1_limit$), val(x2_limit$)]
y_lim = [-val(y1_limit$), val(y2_limit$)]

if (associated(ap_ele%photon)) then
  ph => ap_ele%photon
  if (allocated(ph%displacement%pt)) then
    x_lim(1) = min(x_lim(1), ph%displacement%r0(1) + ph%displacement%dr(1) * lbound(ph%displacement%pt,1)-0.5_rp)
    x_lim(2) = max(x_lim(2), ph%displacement%r0(1) + ph%displacement%dr(1) * ubound(ph%displacement%pt,1)+0.5_rp)
    y_lim(1) = min(y_lim(1), ph%displacement%r0(2) + ph%displacement%dr(2) * lbound(ph%displacement%pt,2)-0.5_rp)
    y_lim(2) = max(y_lim(2), ph%displacement%r0(2) + ph%displacement%dr(2) * ubound(ph%displacement%pt,2)+0.5_rp)
  endif
endif

if (associated(ap_ele%photon)) then
  ph => ap_ele%photon
  if (allocated(ph%h_misalign%pt)) then
    x_lim(1) = min(x_lim(1), ph%h_misalign%r0(1) + ph%h_misalign%dr(1) * lbound(ph%h_misalign%pt,1)-0.5_rp)
    x_lim(2) = max(x_lim(2), ph%h_misalign%r0(1) + ph%h_misalign%dr(1) * ubound(ph%h_misalign%pt,1)+0.5_rp)
    y_lim(1) = min(y_lim(1), ph%h_misalign%r0(2) + ph%h_misalign%dr(2) * lbound(ph%h_misalign%pt,2)-0.5_rp)
    y_lim(2) = max(y_lim(2), ph%h_misalign%r0(2) + ph%h_misalign%dr(2) * ubound(ph%h_misalign%pt,2)+0.5_rp)
  endif
endif

if (associated(ap_ele%photon)) then
  ph => ap_ele%photon
  if (allocated(ph%segmented%pt)) then
    x_lim(1) = min(x_lim(1), ph%segmented%r0(1) + ph%segmented%dr(1) * lbound(ph%segmented%pt,1)-0.5_rp)
    x_lim(2) = max(x_lim(2), ph%segmented%r0(1) + ph%segmented%dr(1) * ubound(ph%segmented%pt,1)+0.5_rp)
    y_lim(1) = min(y_lim(1), ph%segmented%r0(2) + ph%segmented%dr(2) * lbound(ph%segmented%pt,2)-0.5_rp)
    y_lim(2) = max(y_lim(2), ph%segmented%r0(2) + ph%segmented%dr(2) * ubound(ph%segmented%pt,2)+0.5_rp)
  endif
endif

! 

target%type = rectangular$
target%ele_loc = lat_ele_loc_struct(ap_ele%ix_ele, ap_ele%ix_branch)

z = 0
if (stream_ele_end (ap_ele%aperture_at, ap_ele%orientation) == upstream_end$) z = -ap_ele%value(l$)

call photon_target_corner_calc (ap_ele, (x_lim(1)+x_lim(2))/2, (y_lim(1)+y_lim(2))/2, z, ele, target%center)

call photon_target_corner_calc (ap_ele, x_lim(1), y_lim(1), z, ele, target%corner(1))
call photon_target_corner_calc (ap_ele, x_lim(1), y_lim(2), z, ele, target%corner(2))
call photon_target_corner_calc (ap_ele, x_lim(2), y_lim(1), z, ele, target%corner(3))
call photon_target_corner_calc (ap_ele, x_lim(2), y_lim(2), z, ele, target%corner(4))

target%n_corner = 4

! If there is surface curvature then the aperture rectangle becomes a 3D aperture box.

if (associated(ap_ele%photon)) then
  if (has_curvature(ap_ele%photon) .and. ap_ele%aperture_at == surface$) then
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
! Subroutine photon_add_to_detector_statistics (orbit0, orbit, ele, ix_pt, iy_pt, pixel_pt)
!
! Routine to add photon statistics to the appropriate pixel of a "detector" grid.
!
! It is assumed that track_to_surface has been called so that the photon is at the
! detector surface and that orbit%vec(1) and %vec(3) are in surface coords (needed for curved detectors).
!
! Input:
!   orbit0    -- coord_struct: Photon coords at beginning of lattice
!   orbit     -- coord_struct: Photon coords at the detector.
!   ele       -- ele_struct: Element with grid.
!   pixel_pt  -- pixel_pt_struct, optional: If present then use this grid point
!                 instead of the grid point determined by the (x, y) coords of the photon
!   
!
! Output:
!   ele           -- ele_struct: Element with updatted grid.
!   ix_pt, iy_pt  -- integer, optional: Index of upgraded ele%photon%surface%grid%pt(:,:) point.
!                       These arguments are not set if the pixel_pt argument is present.
!-

subroutine photon_add_to_detector_statistics (orbit0, orbit, ele, ix_pt, iy_pt, pixel_pt)

type (coord_struct) orbit0, orbit, orb
type (ele_struct), target :: ele
type (pixel_pt_struct), optional, target :: pixel_pt
type (pixel_pt_struct), pointer :: pix
type (pixel_detec_struct), pointer :: pixel
type (branch_struct), pointer :: branch

real(rp) phase, intens_x, intens_y, intens, dE

integer, optional :: ix_pt, iy_pt
integer ix_pix(2)

! Find grid pt to update.
! Note: dr(i) can be zero for 1-dim grid

if (present(pixel_pt)) then
  pix => pixel_pt

else
  pixel => ele%photon%pixel
  ix_pix = detector_pixel_pt(orbit, ele)

  if (present(ix_pt)) ix_pt = ix_pix(1)
  if (present(iy_pt)) iy_pt = ix_pix(2)

  ! If outside of detector region then do nothing.

  pixel%n_hit_detec = pixel%n_hit_detec + 1
  if (ix_pix(1) < lbound(pixel%pt, 1) .or. ix_pix(1) > ubound(pixel%pt, 1) .or. &
      ix_pix(2) < lbound(pixel%pt, 2) .or. ix_pix(2) > ubound(pixel%pt, 2)) return

  pix => pixel%pt(ix_pix(1), ix_pix(2))
  pixel%n_hit_pixel = pixel%n_hit_pixel + 1
endif

! Add to det stat

branch => pointer_to_branch(ele)
pix%n_photon  = pix%n_photon + 1
if (branch%lat%photon_type == coherent$) then
  phase = orbit%phase(1) 
  pix%E_x = pix%E_x + orbit%field(1) * cmplx(cos(phase), sin(phase), rp)
  phase = orbit%phase(2) 
  pix%E_y = pix%E_y + orbit%field(2) * cmplx(cos(phase), sin(phase), rp)
else
  intens_x = orbit%field(1)**2
  intens_y = orbit%field(2)**2
  intens = intens_x + intens_y
  pix%intensity_x     = pix%intensity_x    + intens_x
  pix%intensity_y     = pix%intensity_y    + intens_y
  pix%intensity       = pix%intensity      + intens

  orb = to_photon_angle_coords (orbit, ele)
  pix%orbit           = pix%orbit          + intens * orb%vec
  pix%orbit_rms       = pix%orbit_rms      + intens * orb%vec**2

  orb = to_photon_angle_coords (orbit0, ele)
  pix%init_orbit      = pix%init_orbit     + intens * orb%vec
  pix%init_orbit_rms  = pix%init_orbit_rms + intens * orb%vec**2
endif

end subroutine photon_add_to_detector_statistics

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Function detector_pixel_pt (orbit, ele) result (ix_pix)
!
! Routine to return the pixel a particle is hitting.
!
! Input:
!   orbit     -- coord_struct: Orbit at surface.
!   ele       -- ele_struct: Detector element.
!
! Output:
!   ix_pix(2) -- integer: index of ele%photon%pixel%pt(:,:) the particle is in.
!- 

function detector_pixel_pt (orbit, ele) result (ix_pix)

type (coord_struct) orbit
type (ele_struct), target :: ele
type (pixel_detec_struct), pointer :: pixel

integer ix_pix(2)

!

ix_pix = 0
pixel => ele%photon%pixel
if (pixel%dr(1) /= 0) ix_pix(1) = nint((orbit%vec(1) - pixel%r0(1)) / pixel%dr(1))
if (pixel%dr(2) /= 0) ix_pix(2) = nint((orbit%vec(3) - pixel%r0(2)) / pixel%dr(2))


end function detector_pixel_pt

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Function to_photon_angle_coords (orb_in, ele) result (orb_out)
!
! Routine to convert from standard photon coords to "angle" coords defined as:
!       x, angle_x, y, angle_y, z, E-E_ref
!
! Input:
!   orb_in        -- coord_struct: orbit in standard photon coords.
!   ele           -- ele_struct: Reference element (generally the detector element.)
!
! Output:
!   orb_out       -- coord_struct: Transformed coordinates.
!-

function to_photon_angle_coords (orb_in, ele) result (orb_out)

type (coord_struct) orb_in, orb_out
type (ele_struct) ele

! To angle coords

orb_out = orb_in

orb_out%vec(2) = atan2(orb_out%vec(2), orb_out%vec(6))
orb_out%vec(4) = atan2(orb_out%vec(4), orb_out%vec(6))
orb_out%vec(6) = orb_out%p0c - ele%value(E_tot$)

end function to_photon_angle_coords

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine to_surface_coords (lab_orbit, ele, surface_orbit)
!
! Routine to convert lab_orbit laboratory coordinates to surface body coordinates.
!
! Input:
!   lab_orbit     -- coord_struct: Photon position in laboratory coords.
!   ele           -- ele_struct: Detector element.
!
! Output:
!   surface_orbit -- coord_struct: Photon position in element body coordinates.
!     %state      -- Set to lost$ if orbit outside of surface (can happen with sperical surface).
!-

subroutine to_surface_coords (lab_orbit, ele, surface_orbit)

use super_recipes_mod, only: super_qromb

type (coord_struct) lab_orbit, surface_orbit
type (ele_struct) ele

real(rp) w_surf(3,3), z
integer idim
logical err_flag

!

surface_orbit = lab_orbit
call offset_photon (ele, surface_orbit, set$)  ! Go to coordinates of the detector
call track_to_surface (ele, surface_orbit, ele%branch%param, w_surf)

if (.not. has_curvature(ele%photon)) return

z = z_at_surface(ele, surface_orbit%vec(1), surface_orbit%vec(3), err_flag)
if (err_flag) then
  surface_orbit%state = lost$
  return
endif

!

idim = 1
surface_orbit%vec(1) = super_qromb(qfunc, 0.0_rp, surface_orbit%vec(1), 1e-6_rp, 1e-8_rp, 4, err_flag)

idim = 2
surface_orbit%vec(3) = super_qromb(qfunc, 0.0_rp, surface_orbit%vec(3), 1e-6_rp, 1e-8_rp, 4, err_flag)

surface_orbit%vec(5) = z_at_surface(ele, surface_orbit%vec(1), surface_orbit%vec(3), err_flag, .true.)

!------------------------------------------------------
contains

function qfunc (x) result (value)

real(rp), intent(in) :: x(:)
real(rp) :: value(size(x))
real(rp) z, xy(2), dz_dxy(2)
logical err_flag

integer i

!

xy = 0

do i = 1, size(x)
  xy(idim) = x(i)
  z = z_at_surface(ele, xy(1), xy(2), err_flag, .true., dz_dxy)
  value(i) = sqrt(1.0_rp + dz_dxy(idim)**2)
enddo

end function qfunc

end subroutine to_surface_coords

end module
