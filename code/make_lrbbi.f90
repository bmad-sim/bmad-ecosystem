!+
! subroutine MAKE_LRBBI(master_lat_oppos, lat, ix_LRBBI, master_ix_LRBBI_oppos)
!
! This subroutine takes lattices with markers for LRBBI's (like those created
!   in MARK_LRBBI) and turns the marked locations into LRBBI elements. The
!   master_lat and master_lat_oppos supply the data for calculating sigmas
!   and offsets.
!
! Modules Needed:
!   use bmad
!
! Input:
!   master_lat_oppos    -- Lat struct: Lat for oppositely circulating
!                           particles, with markers at parasitic crossings.
!   lat(:)              -- Lat struct: Each bunch has its own lat(i)
!                           with markers at all crossings that bunch sees.
!   ix_LRBBI(:,:)        -- Real(rp): First index (i) is the index of the lat
!                           (i.e., bunch), second index (j) is the index of a
!                           beam-beam element's position in lat(i).
!   master_ix_LRBBI_oppos(:,:) -- Real(rp): First index (i) is the index of the
!                           lat, second index (j) is the index of a beam-beam
!                           element (seen by the ith bunch) in the master_lat_oppos.
!                           This index is used to calculate sigmas and offsets.
!
! Output:
!   lat(:)              -- Lat struct: Lat for each bunch with
!                           markers updated to beam-beam elements.
!
! NOTE: The following attributes of the LRBBI elements are changed:
!         %name = "LRBBI"
!         %key = beambeam
!         %value(charge) = -1
!         %value(n_slice) = 1
!         %value(sig_x), (sig_y), (sig_z), (x_offset), (y_offset)
!           (These are assigned values from master_lat and master_lat_oppos.)
!    All other values (i.e., %s_pos) must be set elsewhere!! Also, this does
!     NOT call lat_make_mat6.
!-

#include "CESR_platform.inc"

subroutine MAKE_LRBBI(master_lat_oppos, lat, ix_LRBBI, master_ix_LRBBI_oppos)

  use bmad_struct
  use bmad_interface, except => MAKE_LRBBI

  implicit none

  type (lat_struct), dimension(:) :: lat
  type (lat_struct) :: master_lat_oppos
  type (coord_struct), allocatable, save :: orbit_oppos(:)
  type (normal_modes_struct) :: mode_oppos

  real(rp) :: beta_a, beta_b, eta_a, eta_b
  real(rp) :: sigma_z
  real(rp) :: e_spread, emit_x, emit_y
  real(rp) :: sigma_x(master_lat_oppos%n_ele_track), sigma_y(master_lat_oppos%n_ele_track)
  real(rp) :: n_part

  integer :: i,j
  integer, dimension(:,:) :: ix_LRBBI
  integer, dimension(:,:) :: master_ix_LRBBI_oppos
  
  character*12 type

! init

  call reallocate_coord (orbit_oppos, master_lat_oppos%n_ele_max)

! calc master_lat_opps parameters

  n_part = master_lat_oppos%param%n_part
!  master_lat_oppos%param%n_part = 0.
  call twiss_at_start(master_lat_oppos)
  call twiss_propagate_all(master_lat_oppos)
  call closed_orbit_at_start(master_lat_oppos, orbit_oppos(0), 4, .true.)
  if (.not. bmad_status%ok) return
  call track_all(master_lat_oppos, orbit_oppos)
  call twiss_at_start(master_lat_oppos)
  if (.not. bmad_status%ok) return
  call twiss_propagate_all(master_lat_oppos)
  call radiation_integrals(master_lat_oppos, orbit_oppos, mode_oppos)

!

  do i = 1, size(lat)
    if(lat(i)%param%particle == positron$) type = 'POSITRON'
    if(lat(i)%param%particle == electron$) type = 'ELECTRON'
    do j = 1, size(ix_LRBBI,2)
      if (ix_LRBBI(i,j) == 0) cycle
      lat(i)%ele(ix_LRBBI(i,j))%name = "LRBBI"
      lat(i)%ele(ix_LRBBI(i,j))%key = beambeam$
      lat(i)%ele(ix_LRBBI(i,j))%value(charge$) = -1
      lat(i)%ele(ix_LRBBI(i,j))%value(n_slice$) = 1
      lat(i)%ele(ix_LRBBI(i,j))%type = type
    enddo
  enddo

  e_spread = mode_oppos%sigE_E
  sigma_z = mode_oppos%sig_z

  emit_x = mode_oppos%a%emittance
  emit_y = emit_x * 0.01
  do j = 1, master_lat_oppos%n_ele_track
    beta_a = master_lat_oppos%ele(j)%a%beta
    beta_b = master_lat_oppos%ele(j)%b%beta
    eta_a = master_lat_oppos%ele(j)%a%eta
    eta_b = master_lat_oppos%ele(j)%b%eta
    sigma_x(j) = sqrt(emit_x * beta_a + (eta_a**2) * ((e_spread)**2))
    sigma_y(j) = sqrt(emit_y * beta_b + (eta_b**2) *((e_spread)**2))
  enddo

  do i = 1, size(lat)

    do j = 1, size(ix_LRBBI, 2)
      if (ix_lrbbi(i,j) == 0) cycle
      lat(i)%ele(ix_lrbbi(i,j))%value(sig_x$) = sigma_x(master_ix_LRBBI_oppos(i,j))
      lat(i)%ele(ix_lrbbi(i,j))%value(sig_y$) = sigma_y(master_ix_LRBBI_oppos(i,j))
      lat(i)%ele(ix_lrbbi(i,j))%value(sig_z$) = sigma_z
    enddo
     
    do j = 1, size(ix_LRBBI, 2)
      if (ix_lrbbi(i,j) == 0) cycle
      lat(i)%ele(ix_lrbbi(i,j))%value(x_offset$) = orbit_oppos(master_ix_lrbbi_oppos(i,j))%vec(1)
      lat(i)%ele(ix_lrbbi(i,j))%value(y_offset$) = orbit_oppos(master_ix_lrbbi_oppos(i,j))%vec(3)
    enddo
        
  enddo


end subroutine
