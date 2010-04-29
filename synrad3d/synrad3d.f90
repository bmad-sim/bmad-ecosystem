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
type (wall3d_struct), target :: wall
type (wall3d_pt_struct) wall_pt(0:100)
type (photon3d_coord_struct) p
type (random_state_struct) ran_state
type (photon3d_wall_hit_struct), allocatable :: wall_hit(:)
type (polygon_vertex_struct) v(100)
type (wall3d_polygon_struct), pointer :: poly

real(rp) ds_step_min, d_i0, i0_tot, ds, gx, gy, s_offset
real(rp) emit_a, emit_b, sig_e, g, gamma, radius, r
real(rp) e_filter_min, e_filter_max, s_filter_min, s_filter_max

integer i, j, n, nn, iu, n_wall_pt_max, random_seed, iu_start
integer ix_ele, n_photon_generated, n_photon_array, i0_ele, n_photon_ele, n_photon_here
integer ix_ele_track_start, ix_ele_track_end, iu_hit_file, iu_lat_file
integer photon_direction, num_photons, num_photons_per_pass, n_phot, ios, ix_polygon
integer n_photons_per_pass, ix_generated_warn

character(200) lattice_file, wall_hit_file, reflect_file, lat_ele_file
character(200) photon_start_input_file, photon_start_output_file

character(100) dat_file, dat2_file, wall_file, param_file
character(16) :: r_name = 'synrad3d'

logical ok, filter_on, s_wrap_on, filter_this
logical is_outside

namelist / synrad3d_parameters / ix_ele_track_start, ix_ele_track_end, &
            photon_direction, num_photons, lattice_file, ds_step_min, num_photons_per_pass, &
            emit_a, emit_b, sig_e, sr3d_params, wall_file, dat_file, random_seed, &
            e_filter_min, e_filter_max, s_filter_min, s_filter_max, wall_hit_file, &
            photon_start_input_file, photon_start_output_file, reflect_file, lat_ele_file, &
            ix_generated_warn

namelist / synrad3d_wall / wall_pt
namelist / polygon_def / ix_polygon, v

namelist / start / p, ran_state

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
e_filter_min = -1
e_filter_max = -1
s_filter_min = -1
s_filter_max = -1
wall_hit_file = ''
reflect_file = ''
lat_ele_file = ''
iu_lat_file = 0
iu_hit_file = 0
sr3d_params%debug_on = .false.
photon_start_input_file = '' 
photon_start_output_file = ''
num_photons = -1
num_photons_per_pass = -1
ix_generated_warn = -1

print *, 'Input parameter file: ', trim(param_file)
open (1, file = param_file, status = 'old')
read (1, nml = synrad3d_parameters)
close (1)

if (reflect_file /= '') wall_hit_file = reflect_file  ! Accept old syntax.

! When a filter parameter is set, only photons that satisfy the filter criteria are kept

filter_on = (e_filter_min > 0) .or. (e_filter_max > 0) .or. (s_filter_min >= 0) .or. (s_filter_max >= 0)
s_wrap_on = (s_filter_min >= 0) .and. (s_filter_max >= 0) .and. (s_filter_min > s_filter_max)

! Get wall info

n_wall_pt_max = -1
wall_pt%basic_shape = ''
wall_pt%ante_height2_plus = -1
wall_pt%ante_height2_minus = -1
wall_pt%width2_plus = -1
wall_pt%width2_minus = -1
wall_pt%ix_polygon = -1

open (1, file = wall_file, status = 'old')
read (1, nml = synrad3d_wall)

n = 0
do i = 0, ubound(wall_pt, 1)
  if (wall_pt(i)%basic_shape == 'polygon') then
    wall_pt(i)%ix_polygon = nint(wall_pt(i)%width2)
    n = max (n, wall_pt(i)%ix_polygon)
  endif
  if (wall_pt(i)%basic_shape == '') then
    n_wall_pt_max = i - 1
    exit
  endif
enddo

if (n > 0) then
  allocate (wall%polygon(n))
  do
    v = polygon_vertex_struct(0.0_rp, 0.0_rp, 0.0_rp)
    read (1, nml = polygon_def, iostat = ios)
    if (ios > 0) then ! If error
      print *, 'ERROR READING POLYGON_DEF NAMELIST.'
      rewind (1)
      do
        read (1, nml = polygon_def) ! Generate error message
      enddo
    endif
    if (ios < 0) exit  ! End of file
    if (ix_polygon > n .or. ix_polygon < 1) then
      print *, 'BAD IX_POLYGON VALUE IN WALL FILE: ', ix_polygon
      call err_exit
    endif

    ! Count number of vertices and calc angles.

    do n = 1, size(v)
      v(n)%angle = atan2(v(n)%y, v(n)%x)
      if (n > 1) then
        if (v(n)%angle <= v(n-1)%angle) v(n)%angle = v(n)%angle + twopi
        if (v(n)%angle >= v(n-1)%angle + pi .or. v(n)%angle <= v(n-1)%angle) then
          print *, 'POLYGON SHAPE IS BAD.'
          print *, '  FOR IX_POLYGON =', ix_polygon
          call err_exit
        endif
      endif
      if (v(n+1)%x == 0 .and. v(n+1)%y == 0) exit
    enddo

    if (v(1)%angle < 0 .or. v(n)%angle > twopi) then
      print *, 'FIRST VERTEX CANNOT HAVE ANGLE < 0 AND LAST VERTEX CANNOT HAVE ANGLE > 0.'
      print *, '  FOR IX_POLYGON =', ix_polygon
      call err_exit
    endif

    poly => wall%polygon(ix_polygon)
    nn = n  ! Total number of vertices

    ! If all y >= 0, only half the polygon has been specified.
    ! In this case, assume up/down symmetry.

    if (all(v(1:n)%y >= 0)) then
      nn = 2 * n ! Total number of vetices
      if (v(n)%y == 0) then  ! Do not duplicate v(n) vertex
        nn = nn - 1
        v(n+1:nn) = v(n-1:1:-1)
      else
        v(n+1:nn) = v(n:1:-1)
      endif
      v(n+1:nn)%y     = -v(n+1:nn)%y
      v(n+1:nn)%angle = twopi - v(n+1:nn)%angle
      if (v(1)%y == 0) nn = nn - 1  ! Do not duplicate v(1) vertex
    endif

    ! Transfer the information to the poly%v array.
    ! The last vertex is the first vertex and closes the polygon

    allocate(poly%v(nn+1))
    poly%v(1:nn) = v(1:nn)
    poly%v(nn+1) = v(1)
    poly%v(nn+1)%angle = v(1)%angle + twopi

  enddo
endif

close (1)

! Get lattice

if (lattice_file(1:6) == 'xsif::') then
  call xsif_parser(lattice_file(7:), lat)
else
  call bmad_parser (lattice_file, lat)
endif

call twiss_and_track (lat, orb, ok)
if (.not. ok) stop
  
if (ix_ele_track_end < 0) ix_ele_track_end = lat%n_ele_track

call ran_seed_put (random_seed)

! Transfer info from wall_pt to wall%pt

print *, 'n_wall_pt_max:', n_wall_pt_max

allocate (wall%pt(0:n_wall_pt_max))
wall%pt = wall_pt(0:n_wall_pt_max)
wall%n_pt_max = n_wall_pt_max
wall%pt(n_wall_pt_max)%s = lat%ele(lat%n_ele_track)%s

call sr3d_check_wall (wall)

! Find out much radiation is produced

call radiation_integrals (lat, orb, modes, rad_int_by_ele = rad_int_ele)

if (ix_ele_track_end > ix_ele_track_start) then
  i0_tot = sum(rad_int_ele%i0(ix_ele_track_start+1:ix_ele_track_end))
else
  i0_tot = sum(rad_int_ele%i0(ix_ele_track_start+1:lat%n_ele_track)) + &
           sum(rad_int_ele%i0(1:ix_ele_track_end))
endif

if (i0_tot == 0) then
  call out_io (s_fatal$, r_name, 'No bends in region of interest')
  call err_exit
endif

print *, 'I0 Radiation Integral of entire lattice:', modes%synch_int(0)
print *, 'I0 Radiation Integral over emission region:', i0_tot

! d_i0 determines the number of photons to generatie per unit i0 integral.

n_photons_per_pass = num_photons_per_pass
if (n_photons_per_pass < 1) n_photons_per_pass = num_photons
d_i0 = i0_tot / n_photons_per_pass

! Determine the emittance

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

! Open files

if (wall_hit_file /= '') then
  iu_hit_file = lunget()
  open (iu_hit_file, file = wall_hit_file)
  print *, 'Creating photon hit point output file: ', trim(wall_hit_file)
endif

if (lat_ele_file /= '') then
  iu_lat_file = lunget()
  open (iu_lat_file, file = lat_ele_file, recl = 120)
  print *, 'Creating lattice element output file: ', trim(lat_ele_file)
    write (iu_lat_file, *) 'I0 Radiation Integral of entire lattice:', modes%synch_int(0)
    write (iu_lat_file, *) 'I0 Radiation Integral over emission region:', i0_tot
    write (iu_lat_file, *) ''
    write (iu_lat_file, *) &
        'Index  Name                Type                  S       L          I0    N_phot  ds_step'
endif

! Track through the elements and generate photons.

bmad_com%auto_bookkeeper = .false.  ! Since we are not changing any element params.

n_photon_generated = 0
n_photon_array = 0

allocate (wall_hit(10))

!--------------------------------------------------------------------------
! If the photon_start input file exists then use that

if (photon_start_input_file /= '') then

  ! Open photon start input file and count the number of photons

  print *, 'Opening photon starting position input file: ', trim(photon_start_input_file)
  open (1, file = photon_start_input_file, status = 'old')

  allocate (photons(1000))

  ! Now read the file and initialize photons and the random number generator.
  ! Only set the random number generator if ran_state is set in the file.

  do
    ran_state%iy = -1  ! To see if ran_state is set by the read.
    read (1, nml = start, iostat = ios)
    if (ios < 0) exit
    if (ios > 0) then
      print *, 'Error reading photon starting position at photon index:', n_photon_generated
      cycle
    endif
    n_photon_generated = n_photon_generated + 1
    n_photon_array = n_photon_array + 1
    if (n_photon_array > size(photons)) call reallocate_photon_array (photons, 2*size(photons))
    photon => photons(n_photon_array)
    photon%start = p
    photon%n_wall_hit = 0
    if (ran_state%iy > 0) call ran_seed_put (state = ran_state)
    call check_if_photon_init_coords_outside_wall (is_outside)
    if (is_outside) cycle
    call sr3d_track_photon (photon, lat, wall, wall_hit)
    call check_filter_restrictions(ok)
    if (ok) call print_hit_points (iu_hit_file, photon, wall_hit)
  enddo

  close (1)

! Regular photon generation

else

  ! Open photon start output file

  if (photon_start_output_file /= '') then
    iu_start = lunget()
    open (iu_start, file = photon_start_output_file, recl = 140)
    print *, 'Creating photon start output file: ', trim(photon_start_output_file)
  endif

  ! Allocate photons array

  n = max(num_photons, n_photons_per_pass)
  if (filter_on) then
    allocate (photons(nint(2.1*n)))   
  else
    allocate (photons(nint(1.1*n)))   ! Allow for some slop
  endif

  ix_ele = ix_ele_track_start
  do 

    if (ix_ele == ix_ele_track_end) then
      if (.not. filter_on .or. n_photon_array > 0.9 * num_photons) exit
      ix_ele = ix_ele_track_start
      if (iu_lat_file > 0) close (iu_lat_file)
      iu_lat_file = 0 ! To stop further output
    endif

    ix_ele = ix_ele + 1
    if (ix_ele > lat%n_ele_track) ix_ele = 0

    ele => lat%ele(ix_ele)

    n_phot = nint(rad_int_ele%i0(ix_ele) / d_i0)
    if (n_phot == 0) cycle

    ds = ele%value(l$) / n_phot
    if (ds < ds_step_min) ds = (1+int(ds_step_min/ds)) * ds

    ! Write info to lat_ele_file

    if (iu_lat_file > 0) then
      write (iu_lat_file, '(i6, 2x, a20, a16, f10.3, f8.3, f10.1, i9, f8.3)') ix_ele, ele%name, &
              key_name(ele%key), ele%s, ele%value(l$), rad_int_ele%i0(ix_ele), n_phot, ds
    endif

    ! Loop over all photon generating points.
    ! First point is random in the range [0, ds] to avoid correlations between passes when
    ! there are multiple passes.

    call ran_uniform (r)
    s_offset = r * ds 
    i0_ele = 0         ! Integrated i0 for this element
    n_photon_ele = 0   

    do
      call sr3d_get_emission_pt_params (lat, orb, ix_ele, s_offset, ele_here, orbit_here, gx, gy)
      g = sqrt(gx**2 + gy**2) 
      call convert_total_energy_to (ele%value(e_tot$),  lat%param%particle, gamma)
      ! Generate photons, track to the wall 

      n_photon_here = nint(g * gamma * ds / d_i0)
      do j = 1, n_photon_here
        n_photon_generated = n_photon_generated + 1
        n_photon_array = n_photon_array + 1
        if (n_photon_array > size(photons)) then
          print *, 'INTERNAL ERROR: NUMBER OF PHOTONS GENERATED TOO LARGE!'
          call err_exit
        endif
        photon => photons(n_photon_array)
        photon%ix_photon = n_photon_array
        photon%ix_photon_generated = n_photon_generated

        if (n_photon_generated == ix_generated_warn) then
          print *, 'Note: At ix_generated_warn:', ix_generated_warn ! For debug purposes.
        endif

        call sr3d_emit_photon (ele_here, orbit_here, gx, gy, &
                               emit_a, emit_b, sig_e, photon_direction, photon%start)
        photon%n_wall_hit = 0
        photon%start%ix_ele = ix_ele

        if (photon_start_output_file /= '') then
          call ran_seed_get (state = ran_state)
          write (iu_start, '(a)')           '&start'
          write (iu_start, '(a, 6es20.12)') '  p%vec = ', photon%start%vec
          write (iu_start, '(a, es20.12)')  '  p%energy =', photon%start%energy
          write (iu_start, *)               '  ran_state = ', ran_state
          write (iu_start, '(a)')           '/'
        endif

        call check_if_photon_init_coords_outside_wall (is_outside)
        if (is_outside) cycle
        call sr3d_track_photon (photon, lat, wall, wall_hit)
        call check_filter_restrictions (ok)
        if (ok) call print_hit_points (iu_hit_file, photon, wall_hit)

      enddo

      s_offset = s_offset + ds
      if (s_offset > ele%value(l$)) exit

    enddo

  enddo

endif

!

photons(:)%intensity = 5 * sqrt(3.0) * r_e * mass_of(lat%param%particle) * i0_tot / &
                                             (6 * h_bar_planck * c_light * n_photon_generated)

! Write results

open (1, file = dat_file)
print *, 'Data file is: ', trim(dat_file)

if (sr3d_params%stop_if_hit_antechamber .and. &
    (any(wall_pt%ante_height2_plus > 0) .or. &
     any(wall_pt%ante_height2_minus > 0))) then
  dat2_file = trim(dat_file) // '.antechamber'
  open (2, file = dat2_file)
  print *, 'Data file for photons hitting the antechamber: ', trim(dat_file)
endif

write (1, *) 'ix_ele_track_start   =', ix_ele_track_start
write (1, *) 'ix_ele_track_end     =', ix_ele_track_end
write (1, *) 'photon_direction     =', photon_direction
write (1, *) 'num_photons          =', num_photons
write (1, *) 'num_photons_per_pass =', num_photons_per_pass
write (1, *) 'lattice_file         =', trim(lattice_file)
write (1, *) 'ds_step_min          =', ds_step_min
write (1, *) 'emit_a               =', emit_a
write (1, *) 'emit_b               =', emit_b
write (1, *) 'sig_e                =', sig_e
write (1, *) 'wall_file            =', trim(wall_file)
write (1, *) 'dat_file             =', trim(dat_file)
write (1, *) 'random_seed          =', random_seed
write (1, *) 'sr3d_params%allow_reflections =', sr3d_params%allow_reflections
write (1, *) 'e_filter_min     =', e_filter_min
write (1, *) 'e_filter_max     =', e_filter_max
write (1, *) 's_filter_min     =', s_filter_min
write (1, *) 's_filter_max     =', s_filter_max
write (1, *)

do i = 1, n_photon_array   
  photon => photons(i)
  iu = 1
  if (sr3d_params%stop_if_hit_antechamber .and. photon%hit_antechamber) iu = 2
  write (iu, '(2i8, f12.4, es11.3, 2x, a)') i, photon%n_wall_hit, photon%start%energy, photon%intensity, &
                                             '! index, n_wall_hit, eV, intensity'
  write (iu, '(4f12.6, f12.3, f12.6, a)') photon%start%vec, '  ! Start position'
  write (iu, '(4f12.6, f12.3, f12.6, a)') photon%now%vec,   '  ! End position'
  j = photon%now%ix_ele
  write (iu, '(i8, 3x, 2a)') j, key_name(lat%ele(j)%key), '  ! Lat ele index and class'
enddo

close (1)

!--------------------------------------------------------------------------------------------
contains

!+
! Subroutine check_filter_restrictions (ok)
!
! Routine to check if a photon has passed the filter requirements.
!
! Output:
!   ok -- logical: Set True if passed. False otherwise.
!-

subroutine check_filter_restrictions (ok)

logical ok

! Check filter restrictions

ok = .true.

if (filter_on) then
  filter_this = .false.
  if (e_filter_min > 0 .and. photon%now%energy < e_filter_min) filter_this = .true.
  if (e_filter_max > 0 .and. photon%now%energy > e_filter_max) filter_this = .true.
  if (s_wrap_on) then
    if (photon%now%vec(5) > s_filter_max .and. photon%now%vec(5) < s_filter_min) filter_this = .true.
  else
    if (s_filter_min > 0 .and. photon%now%vec(5) < s_filter_min) filter_this = .true.
    if (s_filter_max > 0 .and. photon%now%vec(5) > s_filter_max) filter_this = .true.
  endif

  if (filter_this) then
    n_photon_array = n_photon_array - 1  ! Delete photon from the array.
    ok = .false.
  endif

endif

end subroutine

!--------------------------------------------------------------------------------------------
! contains

subroutine check_if_photon_init_coords_outside_wall (is_outside)

logical is_outside

call sr3d_photon_radius (photon%start, wall, radius)
if (radius > 1) then
  print *,              'ERROR: INITIALIZED PHOTON IS OUTSIDE THE WALL!', n_photon_generated
  print '(a, 6f10.4)', '        INITIALIZATION PT: ', photon%start%vec      
  n_photon_array = n_photon_array - 1
endif

end subroutine

!--------------------------------------------------------------------------------------------
! contains

!+
! Subroutine reallocate_photon_array (photon_array, n_size)
! 
! Routine to enlarge an array of photons while saving the original information.
! The final size will be at least n_size.
!
! Input:
!   photon_array(:) -- photon3d_track_struct, allocatable: Array of photons
!   n_size          -- Integer: Minimum size
!
! Output:
!   photon_array(:) -- photon3d_track_struct, allocatable: Array with size at least n_size.
!-

subroutine reallocate_photon_array (photon_array, n_size)

implicit none

type (photon3d_track_struct), allocatable :: photon_array(:), temp(:)
integer n, n_size, n_old, j

!

if (.not. allocated(photon_array)) then
  allocate (photon_array(n_size))
  return
endif

n_old = size(photon_array)
if (n_old >= n_size) return

allocate(temp(n_old))
temp = photon_array

deallocate(photon_array)
allocate(photon_array(n_size))

photon_array(1:n_old) = temp

deallocate (temp)

end subroutine

end program
