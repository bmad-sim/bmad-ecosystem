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
type (sr3d_params_struct) synrad3d_params

real(rp) ds_step_min, d_i0, i0_tot, ds, gx, gy, s_offset
real(rp) emit_a, emit_b, sig_e, g, gamma, radius

integer i, j, iu, n_wall_pt_max, random_seed
integer ix_ele, n_photon_tot, i0_ele, n_photon_ele, n_photon_here
integer ix_ele_track_start, ix_ele_track_end
integer photon_direction, num_photons, n_phot

character(100) lattice_file, dat_file, dat2_file, wall_file, param_file
character(16) :: r_name = 'synrad3d'

logical ok

namelist / synrad3d_parameters / ix_ele_track_start, ix_ele_track_end, &
            photon_direction, num_photons, lattice_file, ds_step_min, &
            emit_a, emit_b, sig_e, synrad3d_params, wall_file, dat_file, random_seed

namelist / synrad3d_wall / wall_pt, n_wall_pt_max

! Get parameter file name

if (cesr_iargc() > 1) then
  print *, 'TOO MANY ARGUMENTS ON THE COMMAND LINE!'
  stop
endif

param_file = 'synrad3d.init'
if (cesr_iargc() == 1) call cesr_getarg(1, param_file)

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

print *, 'Input parameter file: ', trim(param_file)
open (1, file = param_file, status = 'old')
read (1, nml = synrad3d_parameters)
close (1)

sr3d_params = synrad3d_params

n_wall_pt_max = -1
wall_pt%type = ''
wall_pt%antechamber_height2_plus = -1
wall_pt%antechamber_height2_minus = -1
wall_pt%width2_plus = -1
wall_pt%width2_minus = -1

open (1, file = wall_file, status = 'old')
read (1, nml = synrad3d_wall)
close (1)

if (n_wall_pt_max > 0) then
  print *, 'NOTE: YOU DO NOT NEED TO SPECIFY N_WALL_PT_MAX IN YOUR WALL FILE!'
  print *, '      THIS SET WILL BE IGNORED.'
endif

do i = 1, ubound(wall_pt, 1)
  if (wall_pt(i)%type == '') then
    n_wall_pt_max = i - 1
    exit
  endif
enddo

print *, 'n_wall_pt_max:', n_wall_pt_max

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

print *, 'I0 Radiation Integral of entire lattice:', modes%synch_int(0)
print *, 'I0 Radiation Integral over emission region:', i0_tot

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

if (sr3d_params%stop_if_hit_antechamber .and. &
    (any(wall_pt%antechamber_height2_plus > 0) .or. &
     any(wall_pt%antechamber_height2_minus > 0))) then
  dat2_file = trim(dat_file) // '.antechamber'
  open (2, file = dat2_file)
  print *, 'Data file for photons hitting the antechamber: ', trim(dat_file)
endif

n_photon_tot = 0
allocate (photons(nint(1.1*num_photons)))   ! Allow for some slop

ix_ele = ix_ele_track_start
do 

  if (ix_ele == ix_ele_track_end) exit
  ix_ele = ix_ele + 1
  if (ix_ele > lat%n_ele_track) ix_ele = 1

  ele => lat%ele(ix_ele)

  n_phot = nint(rad_int_ele%i0(ix_ele) / d_i0)
  if (n_phot == 0) cycle

  ds = ele%value(l$) / n_phot
  if (ds < ds_step_min) ds = (1+int(ds_step_min/ds)) * ds

  ! Loop over all photon generating points

  s_offset = ds / 2
  i0_ele = 0         ! Integrated i0 for this element
  n_photon_ele = 0   

  do
    call sr3d_get_emission_pt_params (lat, orb, ix_ele, s_offset, ele_here, orbit_here, gx, gy)
    g = sqrt(gx**2 + gy**2) 
    call convert_total_energy_to (ele%value(e_tot$),  lat%param%particle, gamma)
    ! Generate photons, track to the wall 

    n_photon_here = nint(g * gamma * ds / d_i0)
    do j = 1, n_photon_here
      n_photon_tot = n_photon_tot + 1
      photon => photons(n_photon_tot)
      photon%ix_photon = n_photon_tot
      call sr3d_emit_photon (ele_here, orbit_here, gx, gy, &
                             emit_a, emit_b, sig_e, photon_direction, photon%start)
      photon%n_reflect = 0
      photon%start%ix_ele = ix_ele

      call sr3d_photon_radius (photon%start, wall, radius)
      if (radius > 1) then
        print *,              'ERROR: INITIALIZED PHOTON IS OUTSIDE THE WALL!', n_photon_tot
        print '(a, 6f10.4)', '        INITIALIZATION PT: ', photon%start%vec      
        cycle
      endif

      photon%intensity = 5 * sqrt(3.0) * r_e * mass_of(lat%param%particle) * i0_tot / &
                                                      (6 * h_bar_planck * c_light * num_photons)
      call sr3d_track_photon (photon, lat, wall)

    enddo

    s_offset = s_offset + ds
    if (s_offset > ele%value(l$)) exit

  enddo

enddo

! Write results

write (1, *) 'ix_ele_track_start =', ix_ele_track_start
write (1, *) 'ix_ele_track_end   =', ix_ele_track_end
write (1, *) 'photon_direction   =', photon_direction
write (1, *) 'num_photons        =', num_photons
write (1, *) 'lattice_file       =', lattice_file
write (1, *) 'ds_step_min        =', ds_step_min
write (1, *) 'emit_a             =', emit_a
write (1, *) 'emit_b             =', emit_b
write (1, *) 'sig_e              =', sig_e
write (1, *) 'wall_file          =', wall_file
write (1, *) 'dat_file           =', dat_file
write (1, *) 'random_seed        =', random_seed
write (1, *) 'synrad3d_params%allow_reflections =', synrad3d_params%allow_reflections
write (1, *)

do i = 1, n_photon_tot      
  photon => photons(i)
  iu = 1
  if (sr3d_params%stop_if_hit_antechamber .and. photon%hit_antechamber) iu = 2
  write (iu, '(2i8, f10.2, es11.3, 2x, a)') i, photon%n_reflect, photon%start%energy, photon%intensity, &
                                             '! index, n_reflect, eV, intensity'
  write (iu, '(4f12.6, f12.3, f12.6, a)') photon%start%vec, '  ! Start position'
  write (iu, '(4f12.6, f12.3, f12.6, a)') photon%now%vec,   '  ! End position'
  j = photon%now%ix_ele
  write (iu, '(i8, 3x, 2a)') j, key_name(lat%ele(j)%key), '  ! Lat ele index and class'
  do j = 0, photon%n_reflect + 1
    ! photon%reflect(j)%vec 
    ! photon%reflect(j)%track_len
  enddo 
enddo

close (1)

end program
