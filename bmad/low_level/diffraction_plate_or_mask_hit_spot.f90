!+
! Function diffraction_plate_or_mask_hit_spot (ele, orbit) result (ix_section)
!
! Routine to determine where a particle hits on a diffraction_plate or mask element.
!
! Note: It is assumed that orbit is in the frame of reference of the element.
! That is, offset_photon/offset_particle needs to be called before this routine.
!
! Input:
!   ele     -- ele_struct: diffraction_plate or mask element.
!   orbit   -- coord_struct: particle position.
!
! Output:
!   ix_section -- integer, Set to index of clear section hit. Set to zero if
!                   photon is outside all clear areas.
!-

function diffraction_plate_or_mask_hit_spot (ele, orbit) result (ix_section)

use wall3d_mod, dummy => diffraction_plate_or_mask_hit_spot

implicit none

type (ele_struct), target :: ele
type (coord_struct) orbit
type (wall3d_struct), pointer :: wall3d

integer :: ix_section
integer i, ix_sec

! Logic: A particle is in a clear section if it is inside the section and outside
! all subsiquent opaque sections up to the next clear section.

wall3d => ele%wall3d(1)
ix_section = 0
ix_sec = 0

section_loop: do 

  ! Skip any opaque sections

  do
    ix_sec = ix_sec + 1
    if (ix_sec > size(wall3d%section)) exit section_loop
    if (wall3d%section(ix_sec)%type == clear$) exit
  enddo

  ! Section must be clear. Check if photon is within the section

  if (.not. in_section(wall3d%section(ix_sec))) cycle
  ix_section = ix_sec

  ! Now check if photon is within an opaque section

  do
    ix_sec = ix_sec + 1
    if (ix_sec > size(wall3d%section)) return           ! In this clear area
    if (wall3d%section(ix_sec)%type == clear$) return   ! In this clear area
    if (in_section(wall3d%section(ix_sec))) cycle section_loop  ! Is opaque
  enddo

enddo section_loop

! Not in a clear area...

ix_section = 0

!------------------------------------------------------------
contains

function in_section(sec) result (is_in)

type (wall3d_section_struct) sec
logical is_in

real(rp) x, y, norm, r_wall, dr_dtheta

!

x = orbit%vec(1) - sec%r0(1);  y = orbit%vec(3) - sec%r0(2)

if (x == 0 .and. y == 0) then
  is_in = .true.
  return
endif

norm = norm2([x, y])
call calc_wall_radius (sec%v, x/norm, y/norm, r_wall, dr_dtheta)
is_in = (norm <= r_wall)

end function in_section

end function diffraction_plate_or_mask_hit_spot

