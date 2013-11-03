program lux

use lux_plot_mod

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: detec_ele
type (lux_source_struct), target :: source
type (lux_params_struct) lux_param
type (lux_photon_struct), target :: photon
type (lux_direction_tile_struct), pointer :: tt
type (surface_grid_struct), pointer :: detector
type (coord_struct), pointer :: end_orb, this_orb
type (surface_grid_point_struct), pointer :: pix
type (surface_grid_point_struct) :: pixel
type (branch_struct), pointer :: branch

real(rp) intensity, intens_max, total_intensity, dtime
real(rp) cut, normalization, intensity_normalization_coef, area
real(rp) total_dead_intens, pix_in_file_intensity, x
real(rp) x_sum, y_sum, x2_sum, y2_sum, x_ave, y_ave, e_sum, e2_sum, e_ave, e_ref

integer i, j, n, ie, nt, n_track, n_live, track_state, ix_tile
integer nx, ny, nx_min, nx_max, ny_min, ny_max, ix_tracking_mode
integer nx_active_min, nx_active_max, ny_active_min, ny_active_max
integer random_seed, n_photon1_file, n_throw, ix_ele_photon1, 

character(3) num_str
character(16) random_engine, tracking_mode
character(40) arg, plotting, number_file
character(100) param_file, lattice_file, photon1_out_file, det_pix_out_file
character(100) emission_tile_out_file

logical ok, is_there, reject_dead_at_det_photon1, accept

namelist / params / lattice_file, random_seed, photon1_out_file, det_pix_out_file, lux_param, &
    emission_tile_out_file, intensity_normalization_coef, random_engine, ix_ele_photon1, &
    reject_dead_at_det_photon1, tracking_mode

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

tracking_mode = 'incoherent'
photon1_out_file = ''
ix_ele_photon1 = -1  ! Use detector
reject_dead_at_det_photon1 = .false.
det_pix_out_file = ''
emission_tile_out_file = ''
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

call match_word (tracking_mode, tracking_mode_name, ix_tracking_mode)
if (ix_tracking_mode < 1) then
  print *, 'Unknown tracking mode: ' // tracking_mode
  stop
endif

! Add number to file name

if (index(photon1_out_file, '#') /= 0 .or. index(det_pix_out_file, '#') /= 0 .or. &
    index(emission_tile_out_file, '#') /= 0) then

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
  call sub_in (emission_tile_out_file)

endif

! Read in lattice

call bmad_parser (lattice_file, lat)

select case (lux_param%source_type)
case ('SPHERICAL', 'PLANAR')
  branch => lat%branch(0)
  source%source_ele => lat%ele(1)

case ('BEND')
  branch => lat%branch(1)
  source%source_ele => lat%ele(1)
  if (source%source_ele%slave_status == super_slave$) source%source_ele => pointer_to_lord(source%source_ele, 1)
  nullify (source%branch_ele)
  do i = 1, lat%n_ele_track
    if (lat%ele(i)%key /= photon_branch$) cycle
    source%branch_ele => lat%ele(i)
    exit
  enddo
  if (.not. associated(source%branch_ele)) then
    print *, 'NO BRANCH_ELE ELEMENT FOUND IN LATTICE!'
    stop
  endif    

case default
  print *, 'UNKNOWN LUX_PARAM%SOURCE_TYPE: ' // lux_param%source_type
  stop
end select

call reallocate_coord (photon%orb, branch%n_ele_track)
nt = branch%n_ele_track
detec_ele => branch%ele(nt)
source%det_ele => detec_ele
detector => detec_ele%surface%grid

if (detec_ele%value(x1_limit$) == 0 .or. detec_ele%value(x2_limit$) == 0 .or. &
    detec_ele%value(y1_limit$) == 0 .or. detec_ele%value(y2_limit$) == 0) then
  print *, 'LIMITS NOT SET AT DETECTOR!'
endif

! Divide the sphere of initial photon propagation directions into 
! a grid of roughly rectangular tiles 

call run_timer('START')

call lux_setup (photon, lat, lux_param, source)

! Plot tiles

if (plotting /= '') then
  if (index('tiles', trim(plotting)) == 1) then
    call lux_plot_tiles (lux_param, source)
  endif
endif

! Some init

if (photon1_out_file /= '') then
  open (1, file = photon1_out_file, recl = 200)
endif

n_live = 0
n_track = 0
call reallocate_coord (photon%orb, lat, branch%ix_branch)

nx_min = lbound(detector%pt, 1)
nx_max = ubound(detector%pt, 1)
ny_min = lbound(detector%pt, 2)
ny_max = ubound(detector%pt, 2)

total_intensity = 0
n_photon1_file = 0

x_sum = 0; x2_sum = 0
y_sum = 0; y2_sum = 0
e_sum = 0; e2_sum = 0

!------------------------------------------
! Track photons

do 
  if (lux_param%stop_num_photons > 0 .and. n_track >= lux_param%stop_num_photons) exit
  n_track = n_track + 1

  call generate_photon (photon, lat, lux_param, source)
  if (lux_param%debug) then
    call init_coord (photon%orb(0), lat%beam_start, branch%ele(0), .true., photon$, &
                                        1, branch%ele(0)%value(E_tot$) * (1 + lat%beam_start%vec(6)))
  endif

  call track_all (lat, photon%orb, branch%ix_branch, track_state)

  end_orb => photon%orb(nt)
  call offset_photon (detec_ele, end_orb, set$)

  intensity = end_orb%field(1)**2 + end_orb%field(2)**2

  ! Write results

  if (photon1_out_file /= '' .and. intensity >= lux_param%intensity_min_photon1_cutoff) then
    if (ix_ele_photon1 < 1) then
      this_orb => end_orb
    else
      this_orb => photon%orb(ix_ele_photon1)
    endif

    accept = .true.
    if (reject_dead_at_det_photon1 .and. track_state /= moving_forward$) accept = .false.
    if (this_orb%state /= alive$) accept = .false.
    if (accept) then
      write (1, '(i9, 6f15.10, a, 6f10.5, f11.3, f10.6)') n_live, photon%orb(1)%vec, ' : ', &
                      this_orb%vec, this_orb%p0c, intensity
      n_photon1_file = n_photon1_file + 1
    endif
  endif

  !

  if (track_state /= moving_forward$) cycle

  ! Go to coordinates of the detector

  n_live = n_live + 1
  total_intensity = total_intensity + intensity

  E_ref = detec_ele%value(e_tot$)
  x_sum  = x_sum  + intensity * end_orb%vec(1)
  x2_sum = x2_sum + intensity * end_orb%vec(1)**2
  y_sum  = y_sum  + intensity * end_orb%vec(3)
  y2_sum = y2_sum + intensity * end_orb%vec(3)**2
  e_sum  = e_sum  + intensity * (end_orb%p0c - e_ref)
  e2_sum = e2_sum + intensity * (end_orb%p0c - e_ref)**2

  source%tile(photon%ix_tile)%n_photon_live = source%tile(photon%ix_tile)%n_photon_live + 1
  source%tile(photon%ix_tile)%det_intensity = source%tile(photon%ix_tile)%det_intensity + intensity

  nx = nint((end_orb%vec(1) - detector%r0(1)) / detector%dr(1))
  ny = nint((end_orb%vec(3) - detector%r0(2)) / detector%dr(2))

  if (nx_min <= nx .and. nx <= nx_max .and. ny_min <= ny .and. ny <= ny_max) then
    pix => detector%pt(nx,ny)
    pix%intensity = pix%intensity + intensity
    pix%n_photon  = pix%n_photon + 1
    pix%e_sum  = pix%e_sum  + intensity * (end_orb%p0c - e_ref)
    pix%e2_sum = pix%e2_sum + intensity * (end_orb%p0c - e_ref)**2
  endif

  if (lux_param%stop_total_intensity > 0 .and. total_intensity >= lux_param%stop_total_intensity) exit
enddo

if (lux_param%source_type == 'SPHERICAL') then
  area = size(source%tile) * lux_param%del_phi * lux_param%del_y
else
  area = fourpi
endif

normalization = intensity_normalization_coef * area / (n_track * fourpi)

close(1)

!------------------------------------------
! det_pix_out_file

pix_in_file_intensity = 0

if (det_pix_out_file /= '') then
  open (3, file = det_pix_out_file)
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
  write (3, '(a, i8)')     'n_tile_tot        =', source%n_tile_tot
  write (3, '(a, i8)')     'n_used_tiles      =', size(source%tile)
  write (3, '(a, es14.6)') 'area_used_tiles   =', area
  write (3, '(a, es14.6)') 'normalization     =', normalization
  write (3, '(a, f10.6)')  'dx_pixel          =', detector%dr(1)
  write (3, '(a, f10.6)')  'dy_pixel          =', detector%dr(2)
  write (3, '(a, i8)')     'nx_active_min     =', nx_active_min
  write (3, '(a, i8)')     'nx_active_max     =', nx_active_max
  write (3, '(a, i8)')     'ny_active_min     =', ny_active_min
  write (3, '(a, i8)')     'ny_active_max     =', ny_active_max
  write (3, '(a)')         '#-----------------------------------------------------'
  write (3, '(a)')         '#     ix      iy      x_pix      y_pix      Intensity  N_photn     E_ave     E_rms'

  do i = nx_min, nx_max
  do j = ny_min, ny_max
    pix => detector%pt(i,j)
    if (pix%intensity <= cut .or. pix%n_photon == 0) cycle
    pix%e_ave = pix%e_sum / pix%intensity
    pix%e_rms = sqrt(max(0.0_rp, pix%e2_sum / pix%intensity - pix%e_ave**2)) 
    write (3, '(2i8, 2f11.6, es16.5, i8, 2f10.3)') i, j, [i,j]*detector%dr+detector%r0, &
           pix%intensity * normalization, pix%n_photon, pix%e_ave, pix%e_rms
    pix_in_file_intensity = pix_in_file_intensity + pix%intensity 
  enddo
  enddo

  close(3)

  open (3, file = trim(det_pix_out_file) // '.x')
  write (3, '(a)')        '#-----------------------------------------------------'
  write (3, '(a)')        '#     ix      x_pix      Intensity  N_photn     E_ave     E_rms'
  do i = nx_min, nx_max
    pixel = lux_detector_pixel_struct()
    pixel%intensity = sum(detector%pt(i,:)%intensity)
    pixel%n_photon  = sum(detector%pt(i,:)%n_photon)
    if (pixel%intensity <= cut .or. pixel%n_photon == 0) cycle
    pixel%e_ave = sum(detector%pt(i,:)%e_sum) / pixel%intensity
    pixel%e_rms = sqrt(max(0.0_rp, sum(detector%pt(i,:)%e2_sum) / pixel%intensity - pixel%e_ave**2))
    write (3, '(i8, f11.6, es16.5, i8, 2f10.3)') i, i*detector%dr(1)+detector%r0(1), &
                       pixel%intensity * normalization, pixel%n_photon, pixel%e_ave, pixel%e_rms
  enddo
  close(3)

  open (3, file = trim(det_pix_out_file) // '.y')
  write (3, '(a)')        '#-----------------------------------------------------'
  write (3, '(a)')        '#     iy      y_pix      Intensity  N_photn     E_ave     E_rms'
  do j = ny_min, ny_max
    pixel = lux_detector_pixel_struct()
    pixel%intensity = sum(detector%pt(:,j)%intensity)
    pixel%n_photon  = sum(detector%pt(:,j)%n_photon)
    if (pixel%intensity <= cut .or. pixel%n_photon == 0) cycle
    pixel%e_ave = sum(detector%pt(:,j)%e_sum) / pixel%intensity
    pixel%e_rms = sqrt(max(0.0_rp, sum(detector%pt(:,j)%e2_sum) / pixel%intensity - pixel%e_ave**2)) 
    write (3, '(i8, f11.6, es16.5, i8, 2f10.3)') j, j*detector%dr(2)+detector%r0(2), &
                       pixel%intensity * normalization, pixel%n_photon, pixel%e_ave, pixel%e_rms
  enddo
  close(3)

endif

!------------------------------------------
! emission_tile_out_file

if (emission_tile_out_file /= '') then
  open (10, file = emission_tile_out_file)
  do i = 1, size(source%tile)
    tt => source%tile(i)
    write (10, '(i6, 2x, 2i5, 2x, 2f8.3, 2x, 2i5, es11.3, l4)') &
                      i, tt%iphi, tt%iy, tt%iphi * lux_param%del_phi, tt%iy * lux_param%del_y, &
                      tt%n_photon, tt%n_photon_live, tt%det_intensity * normalization, tt%alive
    if (tt%alive .or. tt%n_photon_live == 0) cycle
    if (tt%det_intensity > 1d-5 * total_intensity) then
      print '(a, 2i5, f10.5)', '*** NOTE! PHOTONS ARE REACHING DETECTOR FROM DEAD TILE I_PHI, I_Y: ', &
              tt%iphi, tt%iy, tt%det_intensity / total_intensity
    endif
  enddo
  close(10)
endif

!------------------------------------------

total_dead_intens = 0
do i = 1, size(source%tile)
  tt => source%tile(i)
  if (tt%alive) cycle
  total_dead_intens = total_dead_intens + tt%det_intensity
enddo
if (total_dead_intens > 1e-4 * total_intensity) then
  print '(a, f10.5, a)', 'TOTAL INTENSITY FROM DEAD TILES:', 100 * total_dead_intens / total_intensity, '%'
elseif (total_intensity == 0) then
  print '(a, f10.5, a)', 'Total intensity from dead tiles:', 0.0_rp, '%'
else
  print '(a, f10.5, a)', 'Total intensity from dead tiles:', 100 * total_dead_intens / total_intensity, '%'
endif

!

call run_timer ('READ', dtime)

print *, 'Photons Tracked:                       ', n_track
print *, 'Photons at detector:                   ', n_live
print *, 'Normalization factor:                  ', normalization
print *, 'Total intensity (unnormalized):        ', total_intensity
print *, 'Total intensity (normalized):          ', total_intensity * normalization
print *
print *, 'Photons at detector making photon1 cut:', n_photon1_file
if (total_intensity /= 0) then
  x_ave = x_sum / total_intensity; y_ave = y_sum / total_intensity
  e_ave = e_sum / total_intensity
  print '(a, f10.3)', &
          ' %Intensity at pixels not in det_pix file:', 100 * (total_intensity - pix_in_file_intensity) / total_intensity
  print '(a, 2f12.5)', 'Average position at det (x, y) (mm):     ', 1000 * x_ave, 1000 * y_ave
  print '(a, 2f12.5)', 'RMS at det (x, y) (mm):                  ', 1000 * sqrt(x2_sum / total_intensity - x_ave**2), &
                                                                    1000 * sqrt(y2_sum / total_intensity - y_ave**2)
  print '(a, 2f12.5)', 'Average energy deviation at det (eV):    ', e_ave
  print '(a, 2f12.5)', 'RMS energy at det (eV):                  ', sqrt(max(0.0_rp, e2_sum / total_intensity - e_ave**2))
endif
print '(a, f10.2)', &
        ' Simulation time (min):               ', dtime/60
print *
print *, 'Photon1 data file:        ', trim(photon1_out_file)
print *, 'Detector pixel data file: ', trim(det_pix_out_file)
print *, 'Emission tile data file:  ', trim(emission_tile_out_file)

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
