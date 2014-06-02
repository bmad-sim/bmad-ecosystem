program lux

use lux_plot_mod

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: detec_ele
type (ele_pointer_struct), allocatable :: eles(:)
type (photon_init_struct) photon_init
type (lux_params_struct) lux_param
type (lux_photon_struct), target :: photon
type (surface_grid_struct), pointer :: detector
type (coord_struct), pointer :: end_orb, this_orb
type (surface_grid_pt_struct), pointer :: pix
type (surface_grid_pt_struct) :: pixel
type (branch_struct), pointer :: s_branch, d_branch

real(rp) intensity, intens_max, intensity_tot, dtime
real(rp) cut, normalization, intensity_normalization_coef, area
real(rp) total_dead_intens, pix_in_file_intensity, x, phase
real(rp) x_sum, y_sum, x2_sum, y2_sum, x_ave, y_ave, e_rms, e_ave, e_ref

integer i, j, n, ie, nt, n_track_at_energy, n_track_tot, n_live, track_state
integer nx, ny, nx_min, nx_max, ny_min, ny_max, ie_max, n_loc
integer nx_active_min, nx_active_max, ny_active_min, ny_active_max
integer random_seed, n_photon1_file, n_throw, ix_ele_photon1_file

character(3) num_str
character(16) random_engine
character(40) arg, plotting, number_file
character(100) param_file, lattice_file, photon1_out_file, det_pix_out_file

logical ok, is_there, reject_dead_at_det_photon1, accept, err

namelist / params / lattice_file, random_seed, photon1_out_file, det_pix_out_file, lux_param, &
    intensity_normalization_coef, random_engine, ix_ele_photon1_file, photon_init, &
    reject_dead_at_det_photon1

! Get inputs

param_file = 'lux.init'
plotting = ''

i = 0
do while (i < cesr_iargc())
  i = i + 1
  call cesr_getarg(i, arg)
  select case (arg)
  case ('-plot')
    i = i + 1
    call cesr_getarg (i, plotting)
  case default
    if (arg(1:1) == '-') then
      print *, 'I DO NOT UNDERSTAND: ', trim(arg)
      print *
      ok = .false.
    endif
    param_file = arg
  end select
enddo

photon1_out_file = ''
ix_ele_photon1_file = -1  ! Use detector
reject_dead_at_det_photon1 = .false.
det_pix_out_file = ''
random_seed = 0
random_engine = 'pseudo'
intensity_normalization_coef = 1e6

open (1, file = param_file, status = 'old')
read (1, nml = params)
close (1) 

if (lux_param%debug) then
  print *, 'Note: lux_param%debug = True'
endif

call ran_seed_put (random_seed)

if (random_engine == 'quasi') then
  if (random_seed == 0) then
    call ran_uniform(x)
    n_throw = 10000 * x
  else
    n_throw = modulo(random_seed-1, 10000)
  endif
  call ran_engine(set = random_engine)
  do i = 1, n_throw
    call ran_uniform(x)
  enddo
endif

! Add number to file name

if (index(photon1_out_file, '#') /= 0 .or. index(det_pix_out_file, '#') /= 0) then

  number_file = 'lux_out_file.number'
  inquire(file = number_file, exist = is_there)
  if (.not. is_there) then
    open (1, file = number_file)
    write (1, *) 0
    close (1)
  endif

  call increment_file_number (number_file, 3, n, num_str)
  call sub_in (photon1_out_file)
  call sub_in (det_pix_out_file)

endif

! Read in lattice

call bmad_parser (lattice_file, lat)

select case (lux_param%simulation_type)
case ('ROCKING_CURVE')
  s_branch => lat%branch(0)  ! Source branch
  d_branch => s_branch       ! detector branch
  lux_com%source_ele => lat%ele(1)
  if (photon_init%e_field_x == 0 .and. photon_init%e_field_y == 0) then
    print *, 'photon_init%e_field_x and/or photon_init%e_field_y must be non-zero for a rocking curve.'
    stop
  endif

case ('NORMAL')
  call lat_ele_locator (lux_param%source_element, lat, eles, n_loc, err)
  if (n_loc == 0) then
    print *, 'NO SOURCE ELEMENTS FOUND MATCHING NAME: "', trim(lux_param%source_element), '"'
    stop
  elseif (n_loc > 1) then
    print *, 'MULTIPLE SOURCE ELEMENTS FOUND MATCHING NAME: "', trim(lux_param%source_element), '"'
    stop
  endif

  lux_com%source_ele => eles(1)%ele
  s_branch => lux_com%source_ele%branch

  select case (lux_com%source_ele%key)
  case (x_ray_init$)
    d_branch => s_branch       ! detector branch
    if (photon_init%e_field_x == 0 .and. photon_init%e_field_y == 0) then
      print *, 'photon_init%e_field_x and/or photon_init%e_field_y must be non-zero for a rocking curve.'
      stop
    endif

  case (sbend$, wiggler$, undulator$)
    nullify (lux_com%fork_ele)
    do i = 1, s_branch%n_ele_track
      if (s_branch%ele(i)%key /= photon_fork$) cycle
      lux_com%fork_ele => s_branch%ele(i)
      exit
    enddo
    if (.not. associated(lux_com%fork_ele)) then
      print *, 'NO APPROPRIATE FORK/PHOTON_FORK ELEMENT FOUND IN LATTICE!'
      stop
    endif
    d_branch => lat%branch(nint(lux_com%fork_ele%value(ix_to_branch$)))       ! detector branch

  case default
    print *, 'CANNOT SIMULATE PHOTONS GENERATED IN ELEMENT OF TYPE: ', trim(key_name(lux_com%source_ele%key))
    stop
  end select

case default
  print *, 'UNKNOWN LUX_PARAM%SIMULATION_TYPE: ' // lux_param%simulation_type
  stop
end select

lux_com%s_branch => s_branch
lux_com%d_branch => d_branch

!

call reallocate_coord (photon%orb, d_branch%n_ele_track)
nt = d_branch%n_ele_track
detec_ele => d_branch%ele(nt)
if (detec_ele%name == 'END' .and. detec_ele%key == marker$) then
  nt = nt - 1
  detec_ele => d_branch%ele(nt)
endif
lux_com%det_ele => detec_ele
detector => detec_ele%photon%surface%grid

if (detec_ele%value(x1_limit$) == 0 .or. detec_ele%value(x2_limit$) == 0 .or. &
    detec_ele%value(y1_limit$) == 0 .or. detec_ele%value(y2_limit$) == 0) then
  print *, 'LIMITS NOT SET AT DETECTOR!'
endif

! Some init

call run_timer('START')

call lux_setup (photon, lat, photon_init, lux_param)

if (photon1_out_file /= '') then
  open (1, file = photon1_out_file, recl = 240)
  write (1, '(a)') '#      |                               Start (x 1e3)                    |                                    End (x 1e3)                     |                 End'
  write (1, '(a)') '#   Ix |          x           vx            y           vy            z |             x           vx            y           vy            z  |     Energy     Intens_x     Intens_y'
endif

call reallocate_coord (photon%orb, lat, d_branch%ix_branch)

nx_min = lbound(detector%pt, 1)
nx_max = ubound(detector%pt, 1)
ny_min = lbound(detector%pt, 2)
ny_max = ubound(detector%pt, 2)

n_photon1_file = 0

x_sum = 0; x2_sum = 0
y_sum = 0; y2_sum = 0
e_ave = 0; e_rms = 0

!------------------------------------------

ie_max = 1
if (d_branch%param%photon_type == coherent$) then
  ie_max = lux_param%n_energy_pts
  if (photon_init%e_field_x == 0 .and. photon_init%e_field_y == 0) then
    print *, 'WARNING: INPUT E_FIELD IS ZERO SO RANDOM FILED WILL BE GENERATED WITH COHERENT PHOTONS!'
  endif
endif

n_live = 0
n_track_tot = 0
intensity_tot = 0

energy_loop: do ie = 1, ie_max

  n_track_at_energy = 0
  intensity_tot = 0

  do 
    if (lux_param%simulation_type /= 'ROCKING_CURVE') then
      if (lux_param%stop_total_intensity > 0 .and. intensity_tot >= lux_param%stop_total_intensity) exit
    endif
    if (lux_param%stop_num_photons > 0 .and. n_track_at_energy >= lux_param%stop_num_photons) exit
    n_track_at_energy = n_track_at_energy + 1
    n_track_tot = n_track_tot + 1
    photon%n_photon_generated = n_track_tot

    call lux_generate_photon (photon, lat, photon_init, ie, lux_param)
    if (lux_param%debug) then
      call init_coord (photon%orb(0), lat%beam_start, d_branch%ele(0), .true., photon$, &
                                          1, d_branch%ele(0)%value(E_tot$) * (1 + lat%beam_start%vec(6)))
    endif

    call track_all (lat, photon%orb, d_branch%ix_branch, track_state)

    end_orb => photon%orb(nt)
    call offset_photon (detec_ele, end_orb, set$)

    intensity = end_orb%field(1)**2 + end_orb%field(2)**2

    ! Write results

    if (photon1_out_file /= '' .and. intensity >= lux_param%intensity_min_photon1_cutoff) then
      if (ix_ele_photon1_file < 1) then
        this_orb => end_orb
      else
        this_orb => photon%orb(ix_ele_photon1_file)
      endif

      accept = .true.
      if (reject_dead_at_det_photon1 .and. track_state /= moving_forward$) accept = .false.
      if (this_orb%state /= alive$) accept = .false.
      if (accept) then
        write (1, '(i6, 5f13.6, 3x, 5f13.6, 3x, f11.3, 2es13.4)') n_track_tot, 1d3*photon%orb(1)%vec(1:5), &
                        1d3*this_orb%vec(1:5), this_orb%p0c, end_orb%field(1)**2, end_orb%field(2)**2
        n_photon1_file = n_photon1_file + 1
      endif
    endif

    !

    if (track_state /= moving_forward$) cycle

    ! Go to coordinates of the detector

    n_live = n_live + 1
    intensity_tot = intensity_tot + intensity

    E_ref = detec_ele%value(e_tot$)
    x_sum  = x_sum  + intensity * end_orb%vec(1)
    x2_sum = x2_sum + intensity * end_orb%vec(1)**2
    y_sum  = y_sum  + intensity * end_orb%vec(3)
    y2_sum = y2_sum + intensity * end_orb%vec(3)**2
    e_ave  = e_ave  + intensity * (end_orb%p0c - e_ref)
    e_rms = e_rms + intensity * (end_orb%p0c - e_ref)**2

    nx = nint((end_orb%vec(1) - detector%r0(1)) / detector%dr(1))
    ny = nint((end_orb%vec(3) - detector%r0(2)) / detector%dr(2))

    if (nx_min <= nx .and. nx <= nx_max .and. ny_min <= ny .and. ny <= ny_max) then
      pix => detector%pt(nx,ny)
      pix%n_photon  = pix%n_photon + 1
      if (d_branch%param%photon_type == coherent$) then
        phase = end_orb%phase(1) 
        pix%E_x = pix%E_x + end_orb%field(1) * [cos(phase), sin(phase)]
        phase = end_orb%phase(2) 
        pix%E_y = pix%E_y + end_orb%field(2) * [cos(phase), sin(phase)]
      else
        pix%intensity = pix%intensity + intensity
        pix%energy_ave  = pix%energy_ave  + intensity * (end_orb%p0c - e_ref)
        pix%energy_rms = pix%energy_rms + intensity * (end_orb%p0c - e_ref)**2
      endif
    endif

  enddo

  if (d_branch%param%photon_type == coherent$) then
    do nx= nx_min, nx_max; do ny= ny_min, ny_max
      pix => detector%pt(nx,ny)
      intensity = pix%E_x(1)**2 + pix%E_x(2)**2 + pix%E_y(1)**2 + pix%E_y(2)**2
      pix%intensity = pix%intensity + intensity
      pix%energy_ave  = pix%energy_ave  + intensity * (end_orb%p0c - e_ref)
      pix%energy_rms = pix%energy_rms + intensity * (end_orb%p0c - e_ref)**2
    enddo; enddo
  endif

enddo energy_loop

area = fourpi
normalization = intensity_normalization_coef * area / (n_track_tot * fourpi)

close(1)

!------------------------------------------
! det_pix_out_file

pix_in_file_intensity = 0
intensity_tot = 0

if (det_pix_out_file /= '') then
  open (3, file = det_pix_out_file, recl = 160)
  intens_max = maxval(detector%pt%intensity)
  cut = intens_max * lux_param%intensity_min_det_pixel_cutoff

  nx_active_min = nx_max;    nx_active_max = nx_min
  ny_active_min = ny_max;    ny_active_max = ny_min

  do i = nx_min, nx_max
  do j = ny_min, ny_max
    if (detector%pt(i,j)%intensity <= cut) cycle
    nx_active_min = min(nx_active_min, i);  nx_active_max = max(nx_active_max, i)
    ny_active_min = min(ny_active_min, j);  ny_active_max = max(ny_active_max, j)
  enddo
  enddo

  write (3, '(3a)')        'master_input_file = "', trim(param_file), '"'
  write (3, '(3a)')        'lattice_file      = "', trim(lattice_file), '"'
  write (3, '(a, es14.6)') 'normalization     =', normalization
  write (3, '(a, f10.6)')  'dx_pixel          =', detector%dr(1)
  write (3, '(a, f10.6)')  'dy_pixel          =', detector%dr(2)
  write (3, '(a, i8)')     'nx_active_min     =', nx_active_min
  write (3, '(a, i8)')     'nx_active_max     =', nx_active_max
  write (3, '(a, i8)')     'ny_active_min     =', ny_active_min
  write (3, '(a, i8)')     'ny_active_max     =', ny_active_max
  write (3, '(a)')         '#-----------------------------------------------------'
  write (3, '(a)')         '#     ix      iy        x_pix        y_pix      Intensity  N_photn     E_ave     E_rms'

  do i = nx_min, nx_max
  do j = ny_min, ny_max
    pix => detector%pt(i,j)
    intensity_tot = intensity_tot + pix%intensity 
    if (pix%intensity <= cut .or. pix%n_photon == 0) cycle
    pix%energy_ave = pix%energy_ave / pix%intensity
    pix%energy_rms = sqrt(max(0.0_rp, pix%energy_rms / pix%intensity - pix%energy_ave**2)) 
    write (3, '(2i8, 2f13.8, es16.5, i8, 2f10.3)') i, j, [i,j]*detector%dr+detector%r0, &
           pix%intensity * normalization, pix%n_photon, pix%energy_ave, pix%energy_rms
    pix_in_file_intensity = pix_in_file_intensity + pix%intensity 
  enddo
  enddo

  close(3)

  open (3, file = trim(det_pix_out_file) // '.x')
  write (3, '(a)')        '#-----------------------------------------------------'
  write (3, '(a)')        '#     ix        x_pix      Intensity  N_photn     E_ave     E_rms'
  do i = nx_min, nx_max
    pixel = surface_grid_pt_struct()
    pixel%intensity = sum(detector%pt(i,:)%intensity)
    pixel%n_photon  = sum(detector%pt(i,:)%n_photon)
    if (pixel%intensity <= cut .or. pixel%n_photon == 0) cycle
    pixel%energy_ave = sum(detector%pt(i,:)%energy_ave) / pixel%intensity
    pixel%energy_rms = sqrt(max(0.0_rp, sum(detector%pt(i,:)%energy_rms) / pixel%intensity - pixel%energy_ave**2))
    write (3, '(i8, f13.8, es16.5, i8, 2f10.3)') i, i*detector%dr(1)+detector%r0(1), &
                       pixel%intensity * normalization, pixel%n_photon, pixel%energy_ave, pixel%energy_rms
  enddo
  close(3)

  open (3, file = trim(det_pix_out_file) // '.y')
  write (3, '(a)')        '#-----------------------------------------------------'
  write (3, '(a)')        '#     iy        y_pix      Intensity  N_photn     E_ave     E_rms'
  do j = ny_min, ny_max
    pixel = surface_grid_pt_struct()
    pixel%intensity = sum(detector%pt(:,j)%intensity)
    pixel%n_photon  = sum(detector%pt(:,j)%n_photon)
    if (pixel%intensity <= cut .or. pixel%n_photon == 0) cycle
    pixel%energy_ave = sum(detector%pt(:,j)%energy_ave) / pixel%intensity
    pixel%energy_rms = sqrt(max(0.0_rp, sum(detector%pt(:,j)%energy_rms) / pixel%intensity - pixel%energy_ave**2)) 
    write (3, '(i8, f13.8, es16.5, i8, 2f10.3)') j, j*detector%dr(2)+detector%r0(2), &
                       pixel%intensity * normalization, pixel%n_photon, pixel%energy_ave, pixel%energy_rms
  enddo
  close(3)

endif

!------------------------------------------

call run_timer ('READ', dtime)

print *, 'Photons Tracked:                       ', n_track_tot
print *, 'Photons at detector:                   ', n_live
print *, 'Normalization factor:                  ', normalization
print *, 'Total intensity (unnormalized):        ', intensity_tot
print *, 'Total intensity (normalized):          ', intensity_tot * normalization
print *
print *, 'Photons at detector making photon1 cut:', n_photon1_file
if (intensity_tot /= 0) then
  x_ave = x_sum / intensity_tot; y_ave = y_sum / intensity_tot
  e_ave = e_ave / intensity_tot
  print '(a, f10.3)', &
          ' %Intensity at pixels not in det_pix file:', 100 * (intensity_tot - pix_in_file_intensity) / intensity_tot
  print '(a, 2f12.5)', 'Average position at det (x, y) (mm):     ', 1000 * x_ave, 1000 * y_ave
  print '(a, 2f12.5)', 'RMS at det (x, y) (mm):                  ', 1000 * sqrt(x2_sum / intensity_tot - x_ave**2), &
                                                                    1000 * sqrt(y2_sum / intensity_tot - y_ave**2)
  print '(a, 2f12.5)', 'Average energy deviation at det (eV):    ', e_ave
  print '(a, 2f12.5)', 'RMS energy at det (eV):                  ', sqrt(max(0.0_rp, e_rms / intensity_tot - e_ave**2))
endif
print '(a, f10.2)', &
        ' Simulation time (min):               ', dtime/60
print *
print *, 'Photon1 data file:        ', trim(photon1_out_file)
print *, 'Detector pixel data file: ', trim(det_pix_out_file)

! End plotting

if (plotting /= '') then
  if (index('detector', trim(plotting)) == 1) then
    call lux_plot_detector (lat, lux_param, detector)
  endif
endif

!--------------------------------------------------------------------
contains

subroutine sub_in (file_name)

character(*) file_name
integer ix

ix = index(file_name, '#')
if (ix == 0) return
file_name = file_name(:ix-1) // num_str // trim(file_name(ix+1:))

end subroutine

end program
