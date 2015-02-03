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

real(rp) intens, intens_x, intens_y, intens_max, intens_tot, dtime
real(rp) cut, normalization, area, intens_tot_x, intens_tot_y
real(rp) total_dead_intens, pix_in_file_intens, x, phase, phase_x, phase_y
real(rp) x_sum, y_sum, x2_sum, y2_sum, x_ave, y_ave, e_rms, e_ave, e_ref

integer i, j, n, ix, ie, nt, n_track_at_energy, n_track_tot, n_live, track_state
integer nx, ny, nx_min, nx_max, ny_min, ny_max, n_loc
integer nx_active_min, nx_active_max, ny_active_min, ny_active_max
integer random_seed, n_photon1_file, n_throw

character(3) num_str
character(16) random_engine
character(40) arg, plotting, number_file
character(100) param_file, lattice_file, photon1_out_file, det_pix_out_file

logical ok, is_there, accept, err

namelist / params / lattice_file, random_seed, photon1_out_file, det_pix_out_file, lux_param, &
    random_engine, photon_init

! Get inputs

param_file = 'lux.init'
plotting = ''
lux_com%verbose = .true.

i = 0
do while (i < cesr_iargc())
  i = i + 1
  call cesr_getarg(i, arg)
  select case (arg)
  case ('-silent')
    lux_com%verbose = .false.
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
det_pix_out_file = ''
random_seed = 0
random_engine = 'pseudo'
lux_param%intensity_normalization_coef = 1e6

open (1, file = param_file, status = 'old')
read (1, nml = params)
close (1) 

if (lux_param%debug) then
  print *, 'Note: lux_param%debug = True'
endif

if (.not. lux_com%verbose) call output_direct (0, .false., max_level = s_success$) ! Suppress bmad_parser info output.

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
lux_com%lat => lat

! locate source element

call lat_ele_locator (lux_param%source_element, lat, eles, n_loc, err)
if (n_loc == 0) then
  print *, 'NO SOURCE ELEMENT FOUND MATCHING NAME: "', trim(lux_param%source_element), '"'
  stop
elseif (n_loc > 1) then
  print *, 'MULTIPLE SOURCE ELEMENTS FOUND MATCHING NAME: "', trim(lux_param%source_element), '"'
  stop
endif

lux_com%source_ele => eles(1)%ele
s_branch => lux_com%source_ele%branch

! Locate detector element

call lat_ele_locator (lux_param%detector_element, lat, eles, n_loc, err)
if (n_loc == 0) then
  print *, 'NO DETECTOR ELEMENT FOUND MATCHING NAME: "', trim(lux_param%source_element), '"'
  stop
elseif (n_loc > 1) then
  print *, 'MULTIPLE DETECTOR ELEMENTS FOUND MATCHING NAME: "', trim(lux_param%source_element), '"'
  stop
endif

lux_com%det_ele => eles(1)%ele
d_branch => lux_com%det_ele%branch

! Locate photon1 element

if (lux_param%photon1_element == '') then
  lux_com%photon1_ele => lux_com%det_ele
else
  call lat_ele_locator (lux_param%photon1_element, lat, eles, n_loc, err)
  if (n_loc == 0) then
    print *, 'NO PHOTON1 ELEMENT FOUND MATCHING NAME: "', trim(lux_param%source_element), '"'
    stop
  elseif (n_loc > 1) then
    print *, 'MULTIPLE PHOTON1 ELEMENTS FOUND MATCHING NAME: "', trim(lux_param%source_element), '"'
    stop
  endif
  lux_com%photon1_ele => eles(1)%ele 
endif

!

select case (lux_com%source_ele%key)
case (x_ray_source$)

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

case default
  print *, 'CANNOT SIMULATE PHOTONS GENERATED IN ELEMENT OF TYPE: ', trim(key_name(lux_com%source_ele%key))
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

if (.not. associated(detec_ele%photon)) then
  print *, 'DETECTOR ELEMENT: ', trim(detec_ele%name)
  print *, 'OF WRONG TYPE:    ', key_name(detec_ele%key)
  stop
endif
detector => detec_ele%photon%surface%grid

if (.not. allocated(detector%pt)) then
  print *, 'DETECTOR GRID NOT SET!'
  stop
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

if (lat%photon_type == coherent$) then
  if (photon_init%e_field_x == 0 .and. photon_init%e_field_y == 0) then
    print *, 'WARNING: INPUT E_FIELD IS ZERO SO RANDOM FILED WILL BE GENERATED WITH COHERENT PHOTONS!'
  endif
endif

n_live = 0
n_track_tot = 0

energy_loop: do ie = 1, 1   ! Energy loop not implemented yet.

  n_track_at_energy = 0
  intens_tot = 0
  intens_tot_x = 0
  intens_tot_y = 0

  do 
    if (lux_param%stop_num_photons > 0 .and. n_track_at_energy >= lux_param%stop_num_photons) exit
    n_track_at_energy = n_track_at_energy + 1
    n_track_tot = n_track_tot + 1
    photon%n_photon_generated = n_track_tot

    call lux_generate_photon (photon, lat, photon_init, lux_param)
    if (lux_param%debug) then
      call init_coord (photon%orb(0), lat%beam_start, d_branch%ele(0), downstream_end$, photon$, &
                                          1, d_branch%ele(0)%value(E_tot$) * (1 + lat%beam_start%vec(6)))
    endif

    call track_all (lat, photon%orb, d_branch%ix_branch, track_state)

    end_orb => photon%orb(nt)
    call offset_photon (detec_ele, end_orb, set$)

    intens_x = end_orb%field(1)**2
    intens_y = end_orb%field(2)**2
    intens = intens_x + intens_y

    ! Write results

    if (photon1_out_file /= '' .and. intens >= lux_param%intensity_min_photon1_cutoff) then
      this_orb => photon%orb(lux_com%photon1_ele%ix_ele)
      accept = .true.
      if (lux_param%reject_dead_at_det_photon1 .and. track_state /= moving_forward$) accept = .false.
      if (this_orb%state /= alive$) accept = .false.
      if (accept) then
        write (1, '(i6, 5f13.6, 3x, 5f13.6, 3x, f11.3, 2es13.4)') n_track_tot, 1d3*photon%orb(1)%vec(1:5), &
                        1d3*this_orb%vec(1:5), this_orb%p0c, end_orb%field(1)**2, end_orb%field(2)**2
        n_photon1_file = n_photon1_file + 1
      endif
    endif

    !

    if (track_state /= moving_forward$) then
      cycle
    endif

    ! Go to coordinates of the detector

    n_live = n_live + 1
    intens_tot_x = intens_tot_x + intens_x
    intens_tot_y = intens_tot_y + intens_y
    intens_tot   = intens_tot_x + intens_tot_y

    E_ref = detec_ele%value(e_tot$)
    x_sum  = x_sum  + intens * end_orb%vec(1)
    x2_sum = x2_sum + intens * end_orb%vec(1)**2
    y_sum  = y_sum  + intens * end_orb%vec(3)
    y2_sum = y2_sum + intens * end_orb%vec(3)**2
    e_ave  = e_ave  + intens * (end_orb%p0c - e_ref)
    e_rms = e_rms + intens * (end_orb%p0c - e_ref)**2

    nx = nint((end_orb%vec(1) - detector%r0(1)) / detector%dr(1))
    ny = nint((end_orb%vec(3) - detector%r0(2)) / detector%dr(2))

    if (nx_min <= nx .and. nx <= nx_max .and. ny_min <= ny .and. ny <= ny_max) then
      pix => detector%pt(nx,ny)
      pix%n_photon  = pix%n_photon + 1
      if (lat%photon_type == coherent$) then
        phase = end_orb%phase(1) 
        pix%E_x = pix%E_x + end_orb%field(1) * [cos(phase), sin(phase)]
        phase = end_orb%phase(2) 
        pix%E_y = pix%E_y + end_orb%field(2) * [cos(phase), sin(phase)]
      else
        pix%intensity_x = pix%intensity_x + intens_x
        pix%intensity_y = pix%intensity_y + intens_y
        pix%intensity   = pix%intensity   + intens
        pix%energy_ave  = pix%energy_ave  + intens * (end_orb%p0c - e_ref)
        pix%energy_rms  = pix%energy_rms  + intens * (end_orb%p0c - e_ref)**2
      endif

      if (photon_init%dE_relative_to_ref) then
        ix = nint((end_orb%p0c - lux_com%det_ele%value(p0c$) - lux_com%energy_bin(1)%energy_ave) / lux_com%dE_bin) + 1 
      else
        ix = nint((end_orb%p0c - lux_com%energy_bin(1)%energy_ave) / lux_com%dE_bin) + 1 
      endif
      if (ix < 1) ix = 1
      if (ix > ubound(lux_com%energy_bin, 1)) ix = ubound(lux_com%energy_bin, 1)
      lux_com%energy_bin(ix)%intensity   = intens
      lux_com%energy_bin(ix)%intensity_x = intens_x
      lux_com%energy_bin(ix)%intensity_y = intens_y
      lux_com%energy_bin(ix)%n_photon    = lux_com%energy_bin(ix)%n_photon + 1

    endif

  enddo

  ! Photon sets of different energies add incoherently

  if (lat%photon_type == coherent$) then
    do nx= nx_min, nx_max; do ny= ny_min, ny_max
      pix => detector%pt(nx,ny)
      intens_x = pix%E_x(1)**2 + pix%E_x(2)**2
      intens_y = pix%E_y(1)**2 + pix%E_y(2)**2
      intens = intens_x + intens_y
      pix%intensity_x = intens_x
      pix%intensity_y = intens_y
      pix%intensity   = pix%intensity + intens
      pix%energy_ave  = pix%energy_ave + intens * (end_orb%p0c - e_ref)
      pix%energy_rms  = pix%energy_rms + intens * (end_orb%p0c - e_ref)**2
    enddo; enddo
  endif

enddo energy_loop

area = fourpi
normalization = lux_param%intensity_normalization_coef * area / (n_track_tot * fourpi)

close(1)

!------------------------------------------
! det_pix_out_file

pix_in_file_intens = 0
intens_tot_x = 0
intens_tot_y = 0
intens_tot = 0

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

  write (3, '(3a)')        'master_input_file   = "', trim(param_file), '"'
  write (3, '(3a)')        'lattice_file        = "', trim(lattice_file), '"'
  write (3, '(a, es14.6)') 'normalization       =', normalization
  write (3, '(a, es16.5)') 'intensity_x_unnorm  =', sum(detector%pt(:,:)%intensity_x)
  write (3, '(a, es16.5)') 'intensity_x_norm    =', sum(detector%pt(:,:)%intensity_x) * normalization
  write (3, '(a, es16.5)') 'intensity_y_unnorm  =', sum(detector%pt(:,:)%intensity_y)
  write (3, '(a, es16.5)') 'intensity_y_norm    =', sum(detector%pt(:,:)%intensity_y) * normalization
  write (3, '(a, es16.5)') 'intensity_unnorm    =', sum(detector%pt(:,:)%intensity)
  write (3, '(a, es16.5)') 'intensity_norm      =', sum(detector%pt(:,:)%intensity) * normalization
  write (3, '(a, f10.6)')  'dx_pixel            =', detector%dr(1)
  write (3, '(a, f10.6)')  'dy_pixel            =', detector%dr(2)
  write (3, '(a, i8)')     'nx_active_min       =', nx_active_min
  write (3, '(a, i8)')     'nx_active_max       =', nx_active_max
  write (3, '(a, i8)')     'ny_active_min       =', ny_active_min
  write (3, '(a, i8)')     'ny_active_max       =', ny_active_max
  write (3, '(a)')         '#-----------------------------------------------------'
  write (3, '(a)')         '#     ix      iy        x_pix        y_pix   Intensity_x       Phase_x   Intensity_y       Phase_y     Intensity  N_photn     E_ave     E_rms'

  do i = nx_min, nx_max
  do j = ny_min, ny_max
    pix => detector%pt(i,j)
    intens_tot = intens_tot + pix%intensity 
    if (pix%intensity <= cut .or. pix%n_photon == 0) cycle
    pix%energy_ave = pix%energy_ave / pix%intensity
    pix%energy_rms = sqrt(max(0.0_rp, pix%energy_rms / pix%intensity - pix%energy_ave**2)) 
    phase_x = 0; phase_y = 0
    if (lat%photon_type == coherent$) then
      phase_x = atan2(pix%e_x(1), pix%e_x(2))
      phase_y = atan2(pix%e_y(1), pix%e_y(2))
    endif
    write (3, '(2i8, 2f13.8, 5es14.5, i8, 2f10.3)') i, j, [i,j]*detector%dr+detector%r0, &
           pix%intensity_x * normalization, phase_x, pix%intensity_y * normalization, phase_y, &
           pix%intensity * normalization, pix%n_photon, pix%energy_ave, pix%energy_rms
    pix_in_file_intens = pix_in_file_intens + pix%intensity 
  enddo
  enddo
  close(3)

  ! det_pix_out_file.x

  open (3, file = trim(det_pix_out_file) // '.x')
  write (3, '(a)')        '#-----------------------------------------------------'
  write (3, '(a)')        '#     ix        x_pix      ntensity_x     Intensity_y       Intensity  N_photn     E_ave     E_rms'
  do i = nx_min, nx_max
    pixel = surface_grid_pt_struct()
    pixel%intensity_x = sum(detector%pt(i,:)%intensity_x)
    pixel%intensity_y = sum(detector%pt(i,:)%intensity_y)
    pixel%intensity   = sum(detector%pt(i,:)%intensity)
    pixel%n_photon    = sum(detector%pt(i,:)%n_photon)
    if (pixel%intensity == 0) then
      pixel%energy_ave = 0
      pixel%energy_rms = 0
    else
      pixel%energy_ave = sum(detector%pt(i,:)%energy_ave) / pixel%intensity
      pixel%energy_rms = sqrt(max(0.0_rp, sum(detector%pt(i,:)%energy_rms) / pixel%intensity - pixel%energy_ave**2))
    endif
    write (3, '(i8, f13.8, 3es16.5, i9, 2f10.3)') i, i*detector%dr(1)+detector%r0(1), &
                       pixel%intensity_x * normalization, pixel%intensity_y * normalization, pixel%intensity * normalization, &
                       pixel%n_photon, pixel%energy_ave, pixel%energy_rms
  enddo
  close(3)

  ! det_pix_out_file.y

  open (3, file = trim(det_pix_out_file) // '.y')
  write (3, '(a)')        '#-----------------------------------------------------'
  write (3, '(a)')        '#     iy        y_pix      Intensity_x     Intensity_y       Intensity  N_photn     E_ave     E_rms'
  do j = ny_min, ny_max
    pixel = surface_grid_pt_struct()
    pixel%intensity_x = sum(detector%pt(:,j)%intensity_x)
    pixel%intensity_y = sum(detector%pt(:,j)%intensity_y)
    pixel%intensity   = sum(detector%pt(:,j)%intensity)
    pixel%n_photon    = sum(detector%pt(:,j)%n_photon)
    if (pixel%intensity == 0) then
      pixel%energy_ave = 0
      pixel%energy_rms = 0
    else
      pixel%energy_ave = sum(detector%pt(:,j)%energy_ave) / pixel%intensity
      pixel%energy_rms = sqrt(max(0.0_rp, sum(detector%pt(:,j)%energy_rms) / pixel%intensity - pixel%energy_ave**2)) 
    endif
    write (3, '(i8, f13.8, 3es16.5, i8, 2f10.3)') j, j*detector%dr(2)+detector%r0(2), &
                       pixel%intensity_x * normalization, pixel%intensity_y * normalization, pixel%intensity * normalization, &
                       pixel%n_photon, pixel%energy_ave, pixel%energy_rms
  enddo
  close(3)

  ! det_pix_out_file.energy

  open (3, file = trim(det_pix_out_file) // '.energy')
  write (3, '(a)')        '#-----------------------------------------------------'
  write (3, '(a)')        '#      Energy     Intensity_x     Intensity_y       Intensity  N_photn'
  do j = 1, ubound(lux_com%energy_bin, 1)
    pix => lux_com%energy_bin(j)
    write (3, '(f13.5, 3es16.5, i9)') pix%energy_ave, pix%intensity_x * normalization, &
                       pix%intensity_y * normalization, pix%intensity * normalization, pix%n_photon
  enddo
  close(3)

endif

!------------------------------------------

call run_timer ('READ', dtime)

if (lux_com%verbose) then
  print *, 'Photons Tracked:                       ', n_track_tot
  print *, 'Photons at detector:                   ', n_live
  print *, 'Normalization factor:                  ', normalization
  print *, 'Total intensity (unnormalized):        ', intens_tot
  print *, 'Total intensity (normalized):          ', intens_tot * normalization
  print *
  print *, 'Photons at detector making photon1 cut:', n_photon1_file
  if (intens_tot /= 0) then
    x_ave = x_sum / intens_tot; y_ave = y_sum / intens_tot
    e_ave = e_ave / intens_tot
    print '(a, f10.3)', &
            ' %Intensity at pixels not in det_pix file:', 100 * (intens_tot - pix_in_file_intens) / intens_tot
    print '(a, 2f12.5)', 'Average position at det (x, y) (mm):     ', 1000 * x_ave, 1000 * y_ave
    print '(a, 2f12.5)', 'RMS at det (x, y) (mm):                  ', 1000 * sqrt(x2_sum / intens_tot - x_ave**2), &
                                                                      1000 * sqrt(y2_sum / intens_tot - y_ave**2)
    print '(a, 2f12.5)', 'Average energy deviation at det (eV):    ', e_ave
    print '(a, 2f12.5)', 'RMS energy at det (eV):                  ', sqrt(max(0.0_rp, e_rms / intens_tot - e_ave**2))
  endif
  print '(a, f10.2)', &
          ' Simulation time (min):               ', dtime/60
  print *
  print *, 'Photon1 data file:        ', trim(photon1_out_file)
  print *, 'Detector pixel data file: ', trim(det_pix_out_file)
endif

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
