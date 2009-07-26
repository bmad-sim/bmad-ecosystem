!+
! Program synrad3d
!
! Program to calculate photoelectron distributions in a lattice
!-

program synrad3d

use synrad_3d

implicit none

type (ele_struct) ele_here
type (ele_struct), pointer :: ele
type (lat_struct), target :: lat
type (coord_struct), allocatable :: orb(:)
type (coord_struct) orbit_here
type (rad_int_common_struct) rad_int_ele
type (photon_track_struct), allocatable :: photons(:)

real(rp) ds_step_min, d_i2, i2_tot, ds, gx, gy

integer ix_ele
integer ix_ele_track_start, ix_ele_track_end
integer photon_direction, num_photons, n_pt

character(100) lattice_file

namelist / parameters / ix_ele_track_start, ix_ele_track_end, &
            photon_direction, num_photons, lattice_file, ds_step_min

! Get parameters.
! Radiation is produced from the end of ix_ele_track_start to the end of ix_ele_track_end.

ix_ele_track_start = 0    
ix_ele_track_end = -1
ds_step_min = 0.01

open (1, file = 'photoelectron.init', status = 'old')
read (1, nml = parameters)
close (1)

! Get lattice

call bmad_parser (lattice_file, lat)
call twiss_and_track (lat, orb, ok)
if (.not. ok) stop

if (ix_ele_track_end == -1) ix_ele_track_end = lat%n_ele_track

! Find out much radiation is produced

call radiation_integrals (lat, orb, mode, rad_int_by_ele = rad_int_ele)
i2_tot = sum(rad_int_ele(ix_ele_track_start+1:ix_ele_track_end)%i2)
d_i2 = i2_tot / num_photons

! Track through the elements and generate photons.

n_photon_tot = 0
allocate (photons(1.1*num_photons))   ! Allow for some slop

ix_ele = ix_ele_track_start
do 

  ix_ele = ix_ele + 1
  if (ix_ele > lat%n_ele_track) ix_ele = 1

  ele => lat%ele(ix_ele)

  n_pt = nint(rad_int_ele(i)%i2 / d_i2)
  if (n_pt == 0) cycle
  ds = max(ds_step_min, ele%value(l$) / n_pt)

  ! Loop over all photon generating points

  s_offset = ds / 2
  i2_ele = 0         ! Integrated i2 for this element
  n_photon_ele = 0   

  do
    call get_emission_pt_params (lat, orb, ix_ele, s_offset, ele_here, orb_here, gx, gy)
    
    ! Generate photons, track to the wall 

    n_photon_here = nint((gx**2 + gy**2) * ds / d_i2)
    do j = 1, n_photon_here
      photon => photons(n_photon_tot + j)
      call emit_photon (ele_here, orb_here, gx, gy, &
                             gx, gy, emit_a, emit_b, photon_direction, photon%init)
      photon%ix_source = i
      call track_photon (photon, lat)
    enddo
    n_photon_tot = n_photon_tot + n_photon_here

    s_offset = s_offset + ds
    if (s_offset > ele%value(l$)) exit

  enddo

  if (ix_ele == ix_ele_track_end) exit

enddo

! Write results to a file

end program
