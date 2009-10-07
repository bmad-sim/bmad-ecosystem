!+
! Program synrad3d
!
! Program to calculate photoelectron distributions in a lattice
!-

program synrad3d

use synrad3d_track_mod
use synrad3d_utils

implicit none

type (ele_struct) ele_here
type (ele_struct), pointer :: ele
type (lat_struct), target :: lat
type (coord_struct), allocatable :: orb(:)
type (coord_struct) orbit_here
type (rad_int_common_struct) rad_int_ele
type (normal_modes_struct) modes
type (photon3d_track_struct), allocatable, target :: photons(:)
type (photon3d_track_struct), pointer :: photon
type (wall3d_struct) wall
type (wall3d_pt_struct) wall_pt(0:100)

real(rp) ds_step_min, d_i0, i0_tot, ds, gx, gy, s_offset
real(rp) emit_a, emit_b, sig_e, g, gamma

integer i, j, n_wall_pt_max, random_seed
integer ix_ele, n_photon_tot, i0_ele, n_photon_ele, n_photon_here
integer ix_ele_track_start, ix_ele_track_end
integer photon_direction, num_photons, n_pt

character(100) lattice_file, dat_file, wall_file
character(16) :: r_name = 'synrad3d'

logical ok

namelist / synrad3d_parameters / ix_ele_track_start, ix_ele_track_end, &
            photon_direction, num_photons, lattice_file, ds_step_min, &
            emit_a, emit_b, sig_e, synrad3d_params, wall_file, dat_file, random_seed

namelist / synrad3d_wall / wall_pt, n_wall_pt_max

! Get parameters.
! Radiation is produced from the end of ix_ele_track_start to the end of ix_ele_track_end.

random_seed = 0
ix_ele_track_start = 0    ! defaults
ix_ele_track_end = -1
ds_step_min = 0.01
emit_a = -1
emit_b = -1
sig_e  = -1
dat_file = 'synrad3d.dat'
wall_file = 'synrad3d.wall'
photon_direction = 1

open (1, file = 'synrad3d.init', status = 'old')
read (1, nml = synrad3d_parameters)
close (1)

open (1, file = wall_file, status = 'old')
read (1, nml = synrad3d_wall)
close (1)

! Get lattice

call bmad_parser (lattice_file, lat)
call twiss_and_track (lat, orb, ok)
if (.not. ok) stop
  
if (ix_ele_track_end == -1) ix_ele_track_end = lat%n_ele_track

allocate (wall%pt(0:n_wall_pt_max))
wall%pt = wall_pt(0:n_wall_pt_max)
wall%n_pt_max = n_wall_pt_max
wall%pt(n_wall_pt_max)%s = lat%ele(lat%n_ele_track)%s

call ran_seed_put (random_seed)

! Find out much radiation is produced

call radiation_integrals (lat, orb, modes, rad_int_by_ele = rad_int_ele)
i0_tot = sum(rad_int_ele%i0(ix_ele_track_start+1:ix_ele_track_end))
if (i0_tot == 0) then
  call out_io (s_fatal$, r_name, 'No bends in region of interest')
  call err_exit
endif

d_i0 = i0_tot / num_photons

if (emit_a < 0) then
  emit_a = modes%a%emittance
  print *, 'Using emit_a =', emit_a
endif

if (emit_b < 0) then
  emit_b = modes%b%emittance
  print *, 'Using emit_b =', emit_b
endif

if (sig_e < 0) then
  sig_e  = modes%sige_e
  print *, 'Using sig_e =', sig_e
endif

! Track through the elements and generate photons.

open (1, file = dat_file)
print *, 'Data file is: ', trim(dat_file)

n_photon_tot = 0
allocate (photons(nint(1.1*num_photons)))   ! Allow for some slop

ix_ele = ix_ele_track_start
do 

  if (ix_ele == ix_ele_track_end) exit
  ix_ele = ix_ele + 1
  if (ix_ele > lat%n_ele_track) ix_ele = 1

  ele => lat%ele(ix_ele)

  n_pt = nint(rad_int_ele%i0(ix_ele) / d_i0)
  if (n_pt == 0) cycle

  ds = max(ds_step_min, ele%value(l$) / n_pt)

  ! Loop over all photon generating points

  s_offset = ds / 2
  i0_ele = 0         ! Integrated i0 for this element
  n_photon_ele = 0   

  do
    call get_emission_pt_params (lat, orb, ix_ele, s_offset, ele_here, orbit_here, gx, gy)
    g = sqrt(gx**2 + gy**2) 
    call convert_total_energy_to (ele%value(e_tot$),  lat%param%particle, gamma)
    ! Generate photons, track to the wall 

    n_photon_here = nint(g * gamma * ds / d_i0)
    do j = 1, n_photon_here
      n_photon_tot = n_photon_tot + 1
      photon => photons(n_photon_tot)
      call emit_photon (ele_here, orbit_here, gx, gy, &
                             emit_a, emit_b, sig_e, photon_direction, photon%start)
      photon%n_reflect = 0
      photon%start%ix_ele = ix_ele
      photon%start%intensity = 5 * sqrt(3.0) * r_e * mass_of(lat%param%particle) * i0_tot / &
                                                      (4 * h_bar_planck * c_light * num_photons)
      call track_photon (photon, lat, wall)

    enddo

    s_offset = s_offset + ds
    if (s_offset > ele%value(l$)) exit

  enddo

enddo

! Write results

do i = 1, n_photon_tot      
  photon => photons(i)
  write (1, '(2i8, f10.2, a)') i, photon%n_reflect, photon%start%energy, '  ! index, n_reflect, eV'
  write (1, '(10x, 4f12.6, f12.3, f12.6, a)') photon%start%vec, '  ! Starting position vector'
  write (1, '(10x, 4f12.6, f12.3, f12.6, a)') photon%now%vec,   '  ! End position vector'
enddo

end program
