program spin_stroboscope

use bmad
use nr, only: indexx
use lmdif_mod
use mode3_mod

implicit none

! Structure for holding information for the N^th turn. 
! Beginning turn is turn 0.

type spin_track_point_struct
  real(rp) :: xfer_mat(3,3) = 0         ! 1-turn Spin transfer matrix starting at N^th turn 
  real(rp) :: invar_spin(3) = 0         ! Inavarient spin.
  real(rp) :: x_spin_axis(3) = 0        ! x-axis of coordinate system for calculating the spin
  real(rp) :: y_spin_axis(3) = 0        ! y-axis of coordinate system for calculating the spin
  real(rp) :: spin_track(3) = 0         ! Tracking for finding tune.
  real(rp) :: spin_track_angle = 0      ! phase angle in the transverse plane when tracking a spin
  real(rp) :: spin_track_dangle = 0     ! phase angle advance starting at the N^th turn.
  real(rp) :: orbit_vec(6) = 0          ! Phase space orbit vector
  real(rp) :: orbit_phase(3) = 0        ! Phase of orbital normal modes.
  real(rp) :: orbit_amp(3) = 0          ! Amplitude of orbital normal modes.
  integer :: ix_nearest(6) = -1         ! Nearest track point in phase space in directions (-x, +x, -y, +y, -z, +z).
  integer :: indexx(3)                  ! Sort order in phase.
  logical :: track_angle_corrected = .false.
end type

type spin_bin_struct
  real(rp) :: orbit_phase(3) = 0
  real(rp) :: invar_spin(3) = 0
  real(rp) :: x_spin_axis(3) = 0          ! x-axis of coordinate system for calculating the spin tune
  real(rp) :: y_spin_axis(3) = 0          ! y-axis of coordinate system for calculating the spin tune
end type

type z_global_struct 
  real(rp) :: z_axis(3) = 0     
  real(rp) :: merit = 1
end type

type fourier1_struct   ! A single Fourier term
  real(rp) :: cos = 0
  real(rp) :: sin = 0
  real(rp) :: amp = 0
  real(rp) :: phase = 0
end type

type fourier_component_struct
  type (fourier1_struct) :: component(3) = fourier1_struct()  ! Fourier transform for x, y, and z components of the spin axis
  real(rp) :: amp = 0
end type

type fourier_axis_struct
  type (fourier_component_struct) :: axis(2) = fourier_component_struct()  ! (x, y) spin axes
end type

integer, parameter :: n_term_max = 5
type fourier_term_struct
  type (fourier_axis_struct) :: term(0:n_term_max) = fourier_axis_struct()
end type

integer, parameter :: n_rot_max = 4
type fourier_rot_struct
  type (fourier_term_struct) :: rot(-n_rot_max:n_rot_max) = fourier_term_struct()
endtype

type fourier_mode_struct
  type (fourier_rot_struct) :: mode(3) = fourier_rot_struct() ! Fourier transform for the three orbital modes.
end type

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele0
type (coord_struct), allocatable, target :: closed_orb(:), orbit(:)
type (coord_struct) start_orb
type (spin_track_point_struct), allocatable :: s(:)
type (spin_bin_struct), allocatable :: bin(:,:,:)
type (z_global_struct) z_try(24)
type (fourier_mode_struct) fourier ! Fourier transform for the spin x-axis

real(rp) orbit_start(3), orbit_stop(3), delta(3), xfer_mat(3,3), closed_orb_invar_spin(3), dtune
real(rp) time, unit_mat(3,3), norm_max, norm, tune2, r, axis(3), ave_invar_spin(3), dphase, dphase2
real(rp) f(3), p_lim, angle, angle0, angle1, angle2, dangle, j_amp(3), old_relaxed_invar_spin(3), old_scatter_min_invar_spin(3)
real(rp) old_angle, spin_tune, phase(3), d, relaxation_tolerance, g, old_optimized_invar_spin(3), old_spin_tune
real(rp) x_global(3), y_global(3), z_global(3), rot_axis(3), merit_unnorm, merit_norm, merit_unnorm0, merit_norm0, d_invar_spin(3)
real(rp) invar_spin_tolerance, closed_orbit_spin_tune, n_in_bin, old_dangle, invar_mat(3,3), old_invar_spin(3)
real(rp) ct, st, x_spin_axis(3), y_spin_axis(3), z_global_axis(3), mat_1turn(6,6), N_mat_conj(6,6)
real(rp) orbit_tune(3), mode3_tune(3), hoff_invar_spin(3), old_hoff_invar_spin(3), a_spin(2), last_invar_spin(3)
real(rp) scatter_norm_cutoff, spin_tune_tolerance, methods_invar_spin0(3)
real(rp), allocatable :: merit_vec(:)

integer i, j, k, k2, ix, it, ir, j1, j2, jj, n, nn, ibx, iby, ibz, ix_x, ix_y, ix_z, ix_xyz(3), np
integer n_bin_max, i0, i1, i2, n_del, n_turns_min, n_turns_max, n0(3), n1(3), track_state, i_dependent, n_wind
integer i_turn, i_cycle, ix_lat_branch, status, ubound_bin(3), ia, iu_methods
integer calc_every, n_points(3)
integer, allocatable :: indx(:)

logical rf_on, first_time, linear_in_action, mode_is_oscillating(3), symmetric, err
logical scatter_minimization_calc, self_consistent_calc, is_indep(3), at_end, out_of_range, verbose, debug

character(200) bmad_lat, param_file, methods_data_file
character(6) num_str

namelist / strob_params / bmad_lat, rf_on, orbit_start, orbit_stop, n_points, n_turns_min, n_turns_max, calc_every, &
           n_in_bin, relaxation_tolerance, scatter_norm_cutoff, methods_data_file, methods_invar_spin0, &
           ix_lat_branch, linear_in_action, invar_spin_tolerance, spin_tune_tolerance, symmetric, &
           scatter_minimization_calc, self_consistent_calc, verbose, debug, z_global_axis

! Read parameter file

if (cesr_iargc() > 1) then
  print *, '??? MULTIPLE ARGUMENTS ON THE COMMAND LINE.'
  stop
endif

param_file = 'spin_stroboscope.init'
if (cesr_iargc() == 1) call cesr_getarg(1, param_file)

debug = .false.
ix_lat_branch = 0
linear_in_action = .true.
invar_spin_tolerance = 1d-3
calc_every = -1
spin_tune_tolerance = 1e-3
scatter_norm_cutoff = 1e-6
scatter_minimization_calc = .true.
self_consistent_calc = .true.
verbose = .true.
z_global_axis = 0
symmetric = .true.
methods_data_file = ''
methods_invar_spin0 = 0

open (1, file = param_file)
read (1, nml = strob_params)
close (1)

bmad_com%auto_bookkeeper = .false.
bmad_com%spin_tracking_on = .true.

allocate(s(0:n_turns_max), indx(0:n_turns_max))
allocate (merit_vec(10*n_turns_max+10))

call mat_make_unit(unit_mat)
first_time = .true.
old_dangle = 0

! Construct possible global coordinate z_global points which correspond to the vertices of two compond icosahedron.

g = (1 + sqrt(5.0_rp)) / 2
r = sqrt(1 + g*g)

z_try( 1)%z_axis = [0.0_rp, 1.0_rp, g] / r
z_try( 2)%z_axis = [g,      0.0_rp, 1.0_rp] / r
z_try( 3)%z_axis = [1.0_rp, g,      0.0_rp] / r

z_try( 4)%z_axis = [0.0_rp, 1.0_rp, -g] / r
z_try( 5)%z_axis = [g,      0.0_rp, -1.0_rp] / r
z_try( 6)%z_axis = [1.0_rp, g,      -0.0_rp] / r

z_try( 7)%z_axis = [0.0_rp, -1.0_rp, g] / r
z_try( 8)%z_axis = [g,      -0.0_rp, 1.0_rp] / r
z_try( 9)%z_axis = [1.0_rp, -g,      0.0_rp] / r

z_try(10)%z_axis = [-0.0_rp, 1.0_rp, g] / r
z_try(11)%z_axis = [-g,      0.0_rp, 1.0_rp] / r
z_try(12)%z_axis = [-1.0_rp, g,      0.0_rp] / r

do i = 1, 12
  z_try(i+12)%z_axis = -z_try(i)%z_axis
enddo

! Prepare the lattice and find the closed orbit

if (verbose) print *, 'Bmad lattice file: ', trim(bmad_lat)
call bmad_parser (bmad_lat, lat)
branch => lat%branch(ix_lat_branch)

if (.not. rf_on) then
  if (verbose) print *, 'Turning off RF...'
  call set_on_off (rfcavity$, lat, off$)
endif

call twiss_and_track (lat, closed_orb, status, ix_lat_branch)
call transfer_matrix_calc (lat, mat_1turn, ix_branch = ix_lat_branch)
call make_N (mat_1turn, N_mat_conj, err); if (err) stop
N_mat_conj = mat_symp_conj(N_mat_conj)

closed_orb_invar_spin = closed_orb(0)%spin
closed_orbit_spin_tune = branch%param%spin_tune/twopi

if (verbose) print '(a, 3f12.6)', 'Closed orbit invariant axis:', closed_orb_invar_spin
if (verbose) print '(a, 3f12.6)', 'Closed orbit spin tune:     ', closed_orbit_spin_tune
ele0 => branch%ele(0)

! Loop over all phase space amplitude pointspoints

call reallocate_coord (orbit, lat, ix_lat_branch)
call run_timer('START')

open (2, file = 'spin_stroboscope.dat', recl = 300)
write (2, '(a)') '# ix  iy  iy  n_turn         Initial position (x, y, pz)              Orbit Action (Jx, Jy, Jz)                   spin_tune      p_lim            Average Spin (x,y,z)             Invariant Spin at Start (x,y,z)'

if (methods_data_file /= '') then
  iu_methods = 3
  open (iu_methods, file = methods_data_file, recl = 250)
  write (iu_methods, '(a)') '# Turn                   Hoff_Invar_Spin                                    Self_Consistant_Invar_Spin                             Scatter_Min_Invar_Spin' 
endif

do ix_x = 1, n_points(1)
do ix_y = 1, n_points(2)
do ix_z = 1, n_points(3)
  ix_xyz = [ix_x, ix_y, ix_z]

  call init_coord (start_orb, closed_orb(0), branch%ele(0), upstream_end$)

  do i = 1, 3
    if (linear_in_action .and. i /= 3) then
      delta(i) = sqrt(orbit_start(i)**2 + (ix_xyz(i) - 1) * (orbit_stop(i)**2 - orbit_start(i)**2) / max(1, n_points(i) - 1))
    else
      delta(i) = orbit_start(i) + (ix_xyz(i) - 1) * (orbit_stop(i) - orbit_start(i)) / max(1, n_points(i) - 1)
    endif

    k = 2*i - 1
    if (i == 3) k = 6
    start_orb%vec(k) = start_orb%vec(k) + delta(i)
  enddo

  orbit(0) = start_orb

  mode_is_oscillating(1:3) = (delta /= 0)
  if (.not. rf_on) mode_is_oscillating(3) = .false.

  call run_timer('READ', time)
  if (verbose) then
    print *
    print '(a)', '!--------------------------------------------------------------'
    print '(a, 3i5, 6f12.8)', 'Orbit: ', ix_xyz, orbit(0)%vec - closed_orb(0)%vec
    print '(a, 3l4)',         'Oscillating Modes: ', mode_is_oscillating
    print '(a, f10.1)',       'Time From Start (min):', time/60
  endif

  ! Track turn-by-turn
  ! Note: When backtracking the spin, the spin is multiplied with the xfer matrix on the right.
  !       This works since the inverse of the xfer matrix is the transpose of the matrix.

  call calc_phase_space_phase_and_amp (s(0), orbit(0)%vec)
  old_spin_tune = 0

  old_hoff_invar_spin = 0
  old_relaxed_invar_spin = 0
  old_scatter_min_invar_spin = 0
  old_invar_spin = 0

  turns_loop: do i_turn = 0, n_turns_max
    do i = 1, 3
      orbit(0)%spin = 0
      orbit(0)%spin(i) = 1
      call track_all (lat, orbit, ix_lat_branch, track_state)
      xfer_mat(:,i) = orbit(branch%n_ele_track)%spin
      if (track_state /= moving_forward$) then
        print '(a, i0)', 'Particle lost in tracking at turn: ', i_turn
        print *, 'Stopping here.'
        stop
      endif
    enddo

    call calc_phase_space_phase_and_amp (s(i_turn), orbit(0)%vec, xfer_mat)
    orbit(0) = orbit(branch%n_ele_track)

    n = count(mode_is_oscillating)
    n_bin_max = nint((i_turn/n_in_bin)**(1.0/n))
    if (modulo(n_bin_max, 2) /= 0) n_bin_max = n_bin_max + 1 ! Make even

    !--------

    if (i_turn == 0) cycle
    if (calc_every < 0 .and. i_turn /= n_turns_max) cycle
    if (calc_every > 0 .and. i_turn /= n_turns_max .and. mod(i_turn, calc_every) /= 0) cycle
    if (i_turn < n_turns_min) cycle

    if (verbose) then
      print *
      print '(a)', '!----------------------------'
      print '(a, i8)', 'Turns Tracked: ', i_turn
    endif

    ! Calc orbit tune

    do j = 1, 3
      if (mode_is_oscillating(j)) then
        ubound_bin(j) = n_bin_max

        n_wind = 0
        do i = 1, i_turn
          if (s(i)%orbit_phase(j) < s(i-1)%orbit_phase(j)) n_wind = n_wind + 1
        enddo
        orbit_tune(j) = (n_wind + s(i_turn)%orbit_phase(j) - s(0)%orbit_phase(j)) / i_turn

      else
        ubound_bin(j) = 1
        orbit_tune(j) = 0
      endif
    enddo

    if (allocated(bin)) deallocate(bin)
    allocate(bin(0:ubound_bin(1), 0:ubound_bin(2), 0:ubound_bin(3)))

    !------------------------------------------------------
    ! For each track point: Find nearest neighbor track points in orbital phase space directions +/- x, +/- y, +/- z.
 
    do k = 1, 3
      s(:)%ix_nearest(2*k-1) = -1
      s(:)%ix_nearest(2*k) = -1
      if (.not. mode_is_oscillating(k)) cycle

      call indexx(s(0:i_turn)%orbit_phase(k), indx(0:i_turn))
      indx = indx - 1  ! Since s(:) is indexed from 0
      s(0:i_turn)%indexx(k) = indx

      ! Nearest in negative k^th direction

      do j = 0, i_turn
        i0 = indx(j)     ! Point to find nearest neighbors for 
        minus_dir_loop: do j2 = 1, i_turn   ! Loop over all possible nearest neighbors.
          jj = j - j2
          if (jj < 0) then
            jj = jj + i_turn + 1
            i2 = indx(jj)
            dphase = 1 + s(i0)%orbit_phase(k) - s(i2)%orbit_phase(k)
          else
            i2 = indx(jj)
            dphase = s(i0)%orbit_phase(k) - s(i2)%orbit_phase(k)
          endif
          ! Check if point is in the correct direction
          do k2 = 1, 3
            if (k2 == k) cycle
            if (.not. mode_is_oscillating(k2)) cycle
            dphase2 = modulo2 (s(i2)%orbit_phase(k2) - s(i0)%orbit_phase(k2), 0.5_rp)
            if (abs(dphase2) > dphase) cycle minus_dir_loop ! Reject since not in -k^th direction
          enddo
          s(i0)%ix_nearest(2*k-1) = i2  ! Found nearest neighbor in -k^th direction
          exit
        enddo minus_dir_loop
      enddo

      ! Nearest in positive k^th direction

      do j = 0, i_turn
        i0 = indx(j)     ! Point to find nearest neighbors for 
        plus_dir_loop: do j2 = 1, i_turn   ! Loop over all possible nearest neighbors.
          jj = j + j2
          if (jj > i_turn) then
            jj = jj - i_turn - 1
            i2 = indx(jj)
            dphase = 1 + s(i2)%orbit_phase(k) - s(i0)%orbit_phase(k)
          else
            i2 = indx(jj)
            dphase = s(i2)%orbit_phase(k) - s(i0)%orbit_phase(k)
          endif
          ! Check if point is in the correct direction
          do k2 = 1, 3
            if (k2 == k) cycle
            if (.not. mode_is_oscillating(k2)) cycle
            dphase2 = modulo2 (s(i2)%orbit_phase(k2) - s(i0)%orbit_phase(k2), 0.5_rp)
            if (abs(dphase2) > dphase) cycle plus_dir_loop ! Reject since not in k^th direction
          enddo
          s(i0)%ix_nearest(2*k) = i2  ! Found nearest neighbor in k^th direction
          exit
        enddo plus_dir_loop
      enddo

    enddo

    !------------------------------------------------------
    ! Hoffstaetter calc
    ! Make a first guess using the x, y, and z axis as a 0th guess.
    ! invar_mat(i,:) is the calculated invariant spin using the i^th axis as the guess.

    call mat_make_unit (invar_mat) 
    do nn = i_turn-1, 0, -1
      invar_mat = unit_mat + matmul(invar_mat, s(nn)%xfer_mat)
    enddo

    invar_mat = invar_mat / (i_turn + 1)
    norm_max = 0
    do i = 1, 3
      norm = norm2(invar_mat(i,:))
      if (norm > norm_max) then
        hoff_invar_spin = invar_mat(i,:) / norm
        norm_max = norm
      endif
    enddo

    if (verbose) then
      print '(a, 3f12.8)',        '  Hoff Calculation:'
      print '(a, 3f12.8)',        '    Invariant at Track pt 0: ', hoff_invar_spin
      print '(a, 4f12.8)',        '    Change in invariant:     ', hoff_invar_spin - old_hoff_invar_spin, norm2(hoff_invar_spin - old_hoff_invar_spin)
    endif

    old_hoff_invar_spin = hoff_invar_spin

    if (first_time) then
      first_time = .false.
      s(0)%invar_spin = hoff_invar_spin
    endif

    do nn = 1, i_turn
      s(nn)%invar_spin = matmul(s(nn-1)%xfer_mat, s(nn-1)%invar_spin)
    enddo

    call flip_spin_axis_if_needed_to_align()

    !------------------------------------------------------
    ! Relaxation to get the self-consistent solution.

    if (self_consistent_calc) then
      merit_norm0 = invariant_spin_merit_function(.true.)
      merit_unnorm0 = invariant_spin_merit_function(.false.)

      ! Relaxation calc...

      last_invar_spin = s(0)%invar_spin
      do i_cycle = 1, 100

        call average_invariant_spin_setup()

        ! Now compute the invarient spin turn-by-turn from the spin field calculated from the binned invariant.

        s(i_turn)%invar_spin = average_invariant_spin(s(i_turn), bin)  
        do nn = i_turn-1, 0, -1
          s(nn)%invar_spin = average_invariant_spin(s(nn), bin) + matmul(s(nn+1)%invar_spin, s(nn)%xfer_mat)
        enddo

        s(0)%invar_spin = s(0)%invar_spin / norm2(s(0)%invar_spin)
        do nn = 1, i_turn
          s(nn)%invar_spin = matmul(s(nn-1)%xfer_mat, s(nn-1)%invar_spin)
        enddo

        d_invar_spin = s(0)%invar_spin - last_invar_spin
        last_invar_spin = s(0)%invar_spin

        if (norm2(d_invar_spin) < relaxation_tolerance) exit
      enddo

      call flip_spin_axis_if_needed_to_align()

      if (verbose) then
        merit_norm = invariant_spin_merit_function(.true.)
        merit_unnorm = invariant_spin_merit_function(.false.)

        print '(a, i4)',              '  Self-consistant Calculation. Relaxation cycles:', i_cycle
        if (debug) then
          print '(a, es12.4, f12.8)', '    Initial figure of merit normalized/unnormalized:', merit_norm0, merit_unnorm0
          print '(a, es12.4, f12.8)', '    Final figure of merit normalized/unnormalized: m', merit_norm, merit_unnorm
        endif
        print '(a, 3f12.8)',          '    Invariant at Track pt 0: ', s(0)%invar_spin
        print '(a, 4f12.8)',          '    Change in invariant:     ', s(0)%invar_spin - old_relaxed_invar_spin, norm2(s(0)%invar_spin - old_relaxed_invar_spin)
      endif

      old_relaxed_invar_spin = s(0)%invar_spin
    endif

    !------------------------------------------------------
    ! Scatter minimization calc. 
    ! For maximum accuracy, the dependent variables are the two spin components with minimum amplitude.

    if (scatter_minimization_calc) then
      ! Minimize the merit function

      merit_norm0 = invariant_spin_merit_function(.true.)
      merit_unnorm0 = invariant_spin_merit_function(.false.)

      i_dependent = maxloc(abs(s(0)%invar_spin), 1)
      is_indep = .true.
      is_indep(i_dependent) = .false.

      a_spin = pack(s(0)%invar_spin, mask = is_indep)
      call initial_lmdif()
      at_end = .false.

      do i_cycle = 1, 100
        call scatter_min_merit_vec (a_spin, merit_vec, nn, out_of_range)
        if (out_of_range) then
          print '(i4, a)', i_cycle, ' Scatter minimization: Out of Range!'
        endif
        if (at_end) exit
        call suggest_lmdif(a_spin, merit_vec(1:nn), relaxation_tolerance, 100, at_end)
      enddo

      call flip_spin_axis_if_needed_to_align()

      if (verbose) then
        merit_norm = invariant_spin_merit_function(.true.)
        merit_unnorm = invariant_spin_merit_function(.false.)

        print '(a, i4, a)',           '  Scatter minimization. Optimization cycles:', i_cycle
        if (debug) then
          print '(a, es12.4, f12.8)', '    Initial figure of merit normalized/unnormalized:', merit_norm0, merit_unnorm0
          print '(a, es12.4, f12.8)', '    Final figure of merit normalized/unnormalized: m', merit_norm, merit_unnorm
        endif
        print '(a, 3f12.8)',          '    Invariant at Track pt 0: ', s(0)%invar_spin
        print '(a, 4f12.8)',          '    Change in invariant:     ', s(0)%invar_spin - old_scatter_min_invar_spin, norm2(s(0)%invar_spin - old_scatter_min_invar_spin)
      endif

      old_scatter_min_invar_spin = s(0)%invar_spin
    endif

    !------

    if (methods_data_file /= '') then
      write (iu_methods, '(i6, 3(4f13.8, 4x))') i_turn, old_hoff_invar_spin, norm2(old_hoff_invar_spin - methods_invar_spin0), &
                                                        old_relaxed_invar_spin, norm2(old_relaxed_invar_spin - methods_invar_spin0), &
                                                        old_scatter_min_invar_spin, norm2(old_scatter_min_invar_spin - methods_invar_spin0)
    endif

    !----------------------------
    ! Spin tune calculation...

    ! First find a direction in spin space that is not near being antiparallel to any invariant spin direction.

    z_try%merit = 1
    do i = 0, i_turn
      do j = 1, 24
        z_try(j)%merit = min(z_try(j)%merit, dot_product(z_try(j)%z_axis, s(i)%invar_spin))
      enddo
    enddo

    n = maxloc(z_try%merit, 1)
    z_global = z_try(n)%z_axis

    if (any(z_global_axis /= 0)) z_global = z_global_axis / norm2(z_global_axis)

    if (verbose) then
      print '(a, 3f12.6)', '  Z_global axis: ', z_global
    endif

    ! Now construct rather arbitrary transverse coords to this direction. 
    ! These are called the  global coordinants

    x_global = [0.0_rp, -z_global(3), z_global(2)]
    x_global = x_global / norm2(x_global)
    y_global = cross_product(z_global, x_global)

    ! At each track point construct transverse coordiantes that are close to the global transverse axes.

    do i = 0, i_turn
      rot_axis = cross_product(z_global, s(i)%invar_spin)
      norm = norm2(rot_axis)
      rot_axis = rot_axis / norm

      angle = atan2(norm, dot_product(z_global, s(i)%invar_spin))
      s(i)%x_spin_axis = rotate_vec_given_axis_angle(x_global, rot_axis, angle)
      s(i)%y_spin_axis = rotate_vec_given_axis_angle(y_global, rot_axis, angle)
    enddo

    ! Track a spin in the transverse plane.

    call transverse_plane_track(.true.)

    ! First pass at adjusting the axes to make the phase advance independent of angle.

    if (debug) call write_spin_track('before_track.dat')

    call rotate_axes_to_flatten_spin_phase_advance (10)
    call transverse_plane_track(.true.)

    if (debug) call write_spin_track('after_track.dat')

    !

    do i = 0, i_turn
      if (s(i)%track_angle_corrected) cycle
      do j = 0, i_turn
        print '(i6, 3f12.6, 3i6, l6, 6i6)', j, s(j)%orbit_phase, s(j)%ix_nearest, &
                                             s(j)%track_angle_corrected, s(j)%indexx
      enddo
      print *, 'ERROR IN ADJUSTING AXES. PLEASE GET HELP!', i
      stop
    enddo

    ! Fourier

    fourier = fourier_mode_struct()
    do j = 1, 3
      if (.not. mode_is_oscillating(j)) cycle

      do i = 0, i_turn
        do it = 0, n_term_max
          ct = cos(twopi * it * s(i)%orbit_phase(j))
          st = sin(twopi * it * s(i)%orbit_phase(j))

          angle = twopi * s(i)%orbit_phase(j)
          do ir = -n_rot_max, n_rot_max
            x_spin_axis = rotate_vec_given_axis_angle(s(i)%x_spin_axis, s(i)%invar_spin, ir*angle)
            y_spin_axis = rotate_vec_given_axis_angle(s(i)%y_spin_axis, s(i)%invar_spin, ir*angle)

            do ix = 1, 3
              fourier%mode(j)%rot(ir)%term(it)%axis(1)%component(ix)%cos = fourier%mode(j)%rot(ir)%term(it)%axis(1)%component(ix)%cos + ct * x_spin_axis(ix)
              fourier%mode(j)%rot(ir)%term(it)%axis(1)%component(ix)%sin = fourier%mode(j)%rot(ir)%term(it)%axis(1)%component(ix)%sin + st * x_spin_axis(ix)
              fourier%mode(j)%rot(ir)%term(it)%axis(2)%component(ix)%cos = fourier%mode(j)%rot(ir)%term(it)%axis(2)%component(ix)%cos + ct * y_spin_axis(ix)
              fourier%mode(j)%rot(ir)%term(it)%axis(2)%component(ix)%sin = fourier%mode(j)%rot(ir)%term(it)%axis(2)%component(ix)%sin + st * y_spin_axis(ix)
            enddo
          enddo
        enddo
      enddo

      do it = 0, n_term_max
        do ir = -n_rot_max, n_rot_max
          do ia = 1, 2
            do ix = 1, 3
              fourier%mode(j)%rot(ir)%term(it)%axis(ia)%component(ix)%cos = fourier%mode(j)%rot(ir)%term(it)%axis(ia)%component(ix)%cos / (i_turn + 1)
              fourier%mode(j)%rot(ir)%term(it)%axis(ia)%component(ix)%sin = fourier%mode(j)%rot(ir)%term(it)%axis(ia)%component(ix)%sin / (i_turn + 1)
              fourier%mode(j)%rot(ir)%term(it)%axis(ia)%component(ix)%amp = norm2([fourier%mode(j)%rot(ir)%term(it)%axis(ia)%component(ix)%cos, fourier%mode(j)%rot(ir)%term(it)%axis(ia)%component(ix)%sin])
              fourier%mode(j)%rot(ir)%term(it)%axis(ia)%component(ix)%phase = atan2(fourier%mode(j)%rot(ir)%term(it)%axis(ia)%component(ix)%sin, fourier%mode(j)%rot(ir)%term(it)%axis(ia)%component(ix)%cos) / twopi
            enddo
            fourier%mode(j)%rot(ir)%term(it)%axis(ia)%amp = norm2(fourier%mode(j)%rot(ir)%term(it)%axis(ia)%component(:)%amp)
          enddo
        enddo
      enddo
    enddo

    if (debug) then
      print '(a)', '# Mode  Axis_Rot |   Amp_x   Amp_y'
      do j = 1, 3
        if (.not. mode_is_oscillating(j)) cycle
        do ir = -n_rot_max, n_rot_max
          print '(i6, i10, 3x, 10(2f8.4, 4x))', j, ir, (fourier%mode(j)%rot(ir)%term(it)%axis(:)%amp, it = 0, n_term_max)
        enddo
      enddo
    endif

    ! Now calculate the tune

    do i = 1, i_turn
      s(i)%spin_track_angle = s(i-1)%spin_track_angle + s(i-1)%spin_track_dangle
    enddo

    spin_tune = (s(i_turn)%spin_track_angle - s(0)%spin_track_angle) / i_turn

    ! p_lim calc

    ave_invar_spin = [sum(s(0:i_turn)%invar_spin(1)), sum(s(0:i_turn)%invar_spin(2)), sum(s(0:i_turn)%invar_spin(3))] / (i_turn + 1)
    p_lim = norm2(ave_invar_spin)

    ! Print some results

    if (verbose) then
      print '(a, f10.5)',  '  Spin_Tune:  ', spin_tune
      print '(a, 3f10.5)', '  Average Spin:', ave_invar_spin
      print '(a, f10.5)',  '  P_lim:', p_lim
    endif

    d_invar_spin = s(0)%invar_spin - old_invar_spin
    dtune = abs(spin_tune - old_spin_tune)
    old_invar_spin = s(0)%invar_spin
    old_spin_tune = spin_tune

    if (norm2(d_invar_spin) < invar_spin_tolerance .and. dtune < spin_tune_tolerance) then
      if (verbose) then
        print '(a, f12.8, a)', '  Stopping tracking since:'
        print '(a, f12.8, a)', '    Chage in invariant spin, ', norm2(d_invar_spin), ', below invar_spin_tolerance. '
        print '(a, f12.8, a)', '    Chage in spin tune, ', dtune, ', below spin_tune_tolerance. '
      endif
      exit
    endif

    if (i_turn == n_turns_max) exit  ! Avoid i_turn = n_turns_max + 1
  enddo turns_loop

  !----------------------------
  ! Write some results

  write (num_str, '(3i2.2)') ix_xyz

  call write_spin_track('spin_track.dat' // trim(num_str))

  ! Sanity check

  if (debug) then
    open (1, file = 'sanity_check.dat')
    do i = 0, i_turn
      write (1, '(i6, 3f10.6, 3(5x, 3f10.6), 2x, f13.8)') i, s(i)%orbit_phase, norm2(s(i)%x_spin_axis), norm2(s(i)%y_spin_axis), norm2(s(i)%invar_spin), &
                   dot_product(s(i)%invar_spin, s(i)%y_spin_axis), dot_product(s(i)%invar_spin, s(i)%x_spin_axis), dot_product(s(i)%y_spin_axis, s(i)%x_spin_axis), dot_product(s(i)%spin_track, s(i)%invar_spin)
    enddo
    close (1)
  endif

  !

  open(1, file = 'orbit_track.dat' // trim(num_str))
  do i = 0, i_turn
    write (1, '(i6, 3f13.8, 5x, 6f13.8, 5x, 3f13.8)') i, s(i)%orbit_phase, s(i)%orbit_vec, s(i)%orbit_amp
  enddo
  close(1)

  !

  do i = 1, 3
    j_amp(i) = (ele0%value(e_tot$)/mass_of(branch%param%particle)) * delta(i)**2 / ele0%b%beta
  enddo

  write (2, '(3i4, i8, 2x, 3f13.8, 4x, 3es14.6, 5x f11.6, f11.6, 2(5x, 3f11.7))'), ix_xyz, i_turn, &
                            delta, j_amp, spin_tune,  p_lim, ave_invar_spin, s(0)%invar_spin

enddo
enddo
enddo

close (2)

! End stuff

call run_timer('READ', time)
if (verbose) print *
if (verbose) print '(a, f10.1)', 'Final Time (min):', time/60

!------------------------------------------------------------
contains

subroutine calc_phase_space_phase_and_amp (s_this, vec, xfer_mat)

type (spin_track_point_struct) s_this
real(rp) vec(6), aj(6)
real(rp), optional :: xfer_mat(3,3)

!

if (present(xfer_mat)) s_this%xfer_mat = xfer_mat
s_this%orbit_vec = vec

aj = matmul(N_mat_conj, vec)


s_this%orbit_phase(1) = -atan2(aj(2), aj(1)) / twopi
s_this%orbit_phase(2) = -atan2(aj(4), aj(3)) / twopi
s_this%orbit_phase(3) = -atan2(aj(6), aj(5)) / twopi

if (s_this%orbit_phase(1) < 0) s_this%orbit_phase(1) = 1 + s_this%orbit_phase(1)
if (s_this%orbit_phase(2) < 0) s_this%orbit_phase(2) = 1 + s_this%orbit_phase(2)
if (s_this%orbit_phase(3) < 0) s_this%orbit_phase(3) = 1 + s_this%orbit_phase(3)

s_this%orbit_amp(1) = sqrt(aj(1)**2 + aj(2)**2)
s_this%orbit_amp(2) = sqrt(aj(3)**2 + aj(4)**2)
s_this%orbit_amp(3) = sqrt(aj(5)**2 + aj(6)**2)


end subroutine calc_phase_space_phase_and_amp

!------------------------------------------------------------
! contains

function average_invariant_spin (s_this, bin) result (invar_spin)

type (spin_track_point_struct) s_this
type (spin_bin_struct) bin(0:,0:,0:)
real(rp) phase(3), invar_spin(3)
real(rp) f(3), r, dphase(3)
integer i, n2, kk, nb(3)

! Symmetric case

if (symmetric) then
  invar_spin = 0
  do kk = 1, 6
    n2 = s_this%ix_nearest(kk)
    if (n2 == -1) cycle
    dphase = phase_diff(s_this%orbit_phase - s(n2)%orbit_phase)
    r = 1 / max(scatter_norm_cutoff, norm2(dphase))
    invar_spin = invar_spin + r * s(n2)%invar_spin
  enddo

  invar_spin = invar_spin / norm2(invar_spin)

! Non-symmetric case

else

  do i = 1, 3
    if (mode_is_oscillating(i)) then
      nb(i) = int(s_this%orbit_phase(i) * n_bin_max)
      f(i) = s_this%orbit_phase(i) * n_bin_max - nb(i)
    else
      nb(i) = 0
      f(i) = 0
    endif
  enddo

  invar_spin = (1-f(1)) * (1-f(2)) * (1-f(3)) * bin(nb(1),   nb(2),   nb(3)  )%invar_spin + &
                  f(1)  * (1-f(2)) * (1-f(3)) * bin(nb(1)+1, nb(2),   nb(3)  )%invar_spin + &
               (1-f(1)) *    f(2)  * (1-f(3)) * bin(nb(1),   nb(2)+1, nb(3)  )%invar_spin + &
                  f(1)  *    f(2)  * (1-f(3)) * bin(nb(1)+1, nb(2)+1, nb(3)  )%invar_spin + &
               (1-f(1)) * (1-f(2)) *    f(3)  * bin(nb(1),   nb(2),   nb(3)+1)%invar_spin + &
                  f(1)  * (1-f(2)) *    f(3)  * bin(nb(1)+1, nb(2),   nb(3)+1)%invar_spin + &
               (1-f(1)) *    f(2)  *    f(3)  * bin(nb(1),   nb(2)+1, nb(3)+1)%invar_spin + &
                  f(1)  *    f(2)  *    f(3)  * bin(nb(1)+1, nb(2)+1, nb(3)+1)%invar_spin

  invar_spin = invar_spin / norm2(invar_spin)
endif

end function average_invariant_spin

!------------------------------------------------------------
! contains

function angle_of (vec1, vec2) result (angle)

real(rp) vec1(3), vec2(3), angle

angle = abs(atan2(norm2(cross_product(vec1, vec2)), dot_product(vec1, vec2)))

end function angle_of

!------------------------------------------------------------
! contains

subroutine average_invariant_spin_setup ()

real(rp) f(3)
integer i, j, k, n(3), np(3), jx

! symmetric calc does not use binning and there is no setup needed in this case.

if (symmetric) return

! Calculate the invariant spin at the bin points.

bin = spin_bin_struct()

do i = 0, i_turn
  do jx = 1, 3
    n(jx)  = int(s(i)%orbit_phase(jx) * ubound_bin(jx))
    f(jx)  = s(i)%orbit_phase(jx) * ubound_bin(jx) - n(jx)
    np(jx) = mod(n(jx)+1, ubound_bin(jx))
  enddo

  bin(n(1), n(2), n(3) )%invar_spin = bin(n(1), n(2), n(3) )%invar_spin + (1-f(1)) * (1-f(2)) * (1-f(3)) * s(i)%invar_spin
  bin(np(1),n(2), n(3) )%invar_spin = bin(np(1),n(2), n(3) )%invar_spin +    f(1)  * (1-f(2)) * (1-f(3)) * s(i)%invar_spin
  bin(n(1), np(2),n(3) )%invar_spin = bin(n(1), np(2),n(3) )%invar_spin + (1-f(1)) *    f(2)  * (1-f(3)) * s(i)%invar_spin
  bin(np(1),np(2),n(3) )%invar_spin = bin(np(1),np(2),n(3) )%invar_spin +    f(1)  *    f(2)  * (1-f(3)) * s(i)%invar_spin
  bin(n(1), n(2), np(3))%invar_spin = bin(n(1), n(2), np(3))%invar_spin + (1-f(1)) * (1-f(2)) *    f(3)  * s(i)%invar_spin
  bin(np(1),n(2), np(3))%invar_spin = bin(np(1),n(2), np(3))%invar_spin +    f(1)  * (1-f(2)) *    f(3)  * s(i)%invar_spin
  bin(n(1), np(2),np(3))%invar_spin = bin(n(1), np(2),np(3))%invar_spin + (1-f(1)) *    f(2)  *    f(3)  * s(i)%invar_spin
  bin(np(1),np(2),np(3))%invar_spin = bin(np(1),np(2),np(3))%invar_spin +    f(1)  *    f(2)  *    f(3)  * s(i)%invar_spin
enddo

do i = 1, 3
  bin(ubound_bin(1),:,:)%invar_spin(i) = bin(0,:,:)%invar_spin(i)
  bin(:,ubound_bin(2),:)%invar_spin(i) = bin(:,0,:)%invar_spin(i)
  bin(:,:,ubound_bin(3))%invar_spin(i) = bin(:,:,0)%invar_spin(i)
enddo

!

do i = 0, ubound_bin(1)
do j = 0, ubound_bin(2)
do k = 0, ubound_bin(3)
  bin(i,j,k)%orbit_phase(1) = real(i, rp) / ubound_bin(1)
  bin(i,j,k)%orbit_phase(2) = real(j, rp) / ubound_bin(2)
  bin(i,j,k)%orbit_phase(3) = real(k, rp) / ubound_bin(3)
  if (norm2(bin(i,j,k)%invar_spin) == 0) cycle
  bin(i,j,k)%invar_spin = bin(i,j,k)%invar_spin / norm2(bin(i,j,k)%invar_spin)
enddo
enddo
enddo

end subroutine average_invariant_spin_setup

!------------------------------------------------------------
! contains

function invariant_spin_merit_function (normalized, file_name) result (merit)

real(rp) merit, dspin, dphase(3), del, r(6), d_spin(6, 3)
integer j, kk, n, n2, num_merit_points
logical normalized
character(*), optional :: file_name

!

if (present(file_name)) then
  open(10, file = file_name)
endif

merit = 0
num_merit_points = 0

if (symmetric) then

   do n = 0, i_turn
    r = 0
    d_spin = 0
    do kk = 1, 6
      n2 = s(n)%ix_nearest(kk)
      if (n2 == -1) cycle
      dphase = phase_diff(s(n)%orbit_phase - s(n2)%orbit_phase)
      r(kk) = 1 / max(scatter_norm_cutoff, norm2(dphase))
      d_spin(kk,:) = s(n2)%invar_spin - s(n)%invar_spin
    enddo

    do j = 1, 3
      num_merit_points = num_merit_points + 1
      if (normalized) then
        merit = merit + dot_product(d_spin(:,j), r)**2
      else
        merit = merit + dot_product(d_spin(:,j), r)**2 / sum(r)**2
      endif
    enddo
  enddo

else
  do kk = 1, 6
    do n = 0, i_turn
      n2 = s(n)%ix_nearest(kk)
      if (n2 == -1) cycle
      num_merit_points = num_merit_points + 3
      dspin = norm2(s(n)%invar_spin - s(n2)%invar_spin)
      if (normalized) then
        dphase = phase_diff(s(n)%orbit_phase - s(n2)%orbit_phase)
        del = max(scatter_norm_cutoff, norm2(dphase))
        merit = merit + dspin**2 / del**2
        if (present(file_name)) then
          write (10, '(i6, 4es12.4)') num_merit_points, merit, dspin**2 / del**2, dspin, del
        endif
      else
        merit = merit + dspin
      endif
    enddo
  enddo
endif

merit = merit / num_merit_points

if (present(file_name)) then
  close(10)
endif


end function invariant_spin_merit_function

!------------------------------------------------------------
! contains

subroutine scatter_min_merit_vec (a_spin, merit_vec, n_pts, out_of_range)

real(rp) a_spin(2), del, dphase(3), amp2, r(6), dspin(6, 3)
real(rp), allocatable :: merit_vec(:)
integer n, kk, n2, j, n_pts
logical out_of_range

!

amp2 = dot_product(a_spin, a_spin)
out_of_range = (amp2 > 1)

if (.not. out_of_range) then
  s(0)%invar_spin = unpack (a_spin, is_indep, s(0)%invar_spin)
  s(0)%invar_spin(i_dependent) = sign_of(s(0)%invar_spin(i_dependent)) * sqrt(1 - amp2)

  do n = 1, i_turn
    s(n)%invar_spin = matmul(s(n-1)%xfer_mat, s(n-1)%invar_spin)
  enddo
endif


if (symmetric) then
  n_pts = 0
  do n = 0, i_turn
    r = 0
    dspin = 0
    do kk = 1, 6
      n2 = s(n)%ix_nearest(kk)
      if (n2 == -1) cycle
      dphase = phase_diff(s(n)%orbit_phase - s(n2)%orbit_phase)
      r(kk) = 1 / max(scatter_norm_cutoff, norm2(dphase))
      dspin(kk,:) = s(n2)%invar_spin - s(n)%invar_spin
    enddo

    do j = 1, 3
      n_pts = n_pts + 1
      if (out_of_range) then
        merit_vec(n_pts) = amp2 * sum(r)
      else
        merit_vec(n_pts) = dot_product (dspin(:,j), r)
      endif
    enddo
  enddo

else
  ! Only use +direction %ix_nearest(kk) which correspond to kk even.
  ! Do this since this halves the number of merit_vec points needed for a small penalty in precision. 

  n_pts = 0
  do kk = 2, 6, 2
    do n = 0, i_turn
      n2 = s(n)%ix_nearest(kk)
      if (n2 == -1) cycle
      dphase = phase_diff(s(n)%orbit_phase - s(n2)%orbit_phase)
      del = max(scatter_norm_cutoff, norm2(dphase))
      do j = 1, 3
        n_pts = n_pts + 1
        if (out_of_range) then
          merit_vec(n_pts) = amp2 / del
        else
          merit_vec(n_pts) = (s(n)%invar_spin(j) - s(n2)%invar_spin(j)) / del
        endif
      enddo
    enddo
  enddo
endif

end subroutine scatter_min_merit_vec

!------------------------------------------------------------
! contains

function phase_diff(dphase) result(diff)

real(rp) dphase(3), diff(3)
integer k

!

do k = 1, 3
  if (mode_is_oscillating(k)) then
    diff(k) = modulo2(dphase(k), 0.5_rp)
  else
    diff(k) = 0
  endif
enddo

end function phase_diff

!------------------------------------------------------------
! contains

subroutine transverse_plane_track(init)

integer i
logical init

! Track spin starting in transverse plane in order to calculate the tune.

if (init) s(0)%spin_track = s(0)%x_spin_axis
s(0)%spin_track_angle = atan2(dot_product(s(0)%spin_track, s(0)%y_spin_axis), dot_product(s(0)%spin_track, s(0)%x_spin_axis)) / twopi

do i = 1, i_turn
  s(i)%spin_track = matmul(s(i-1)%xfer_mat, s(i-1)%spin_track)
  s(i)%spin_track_angle = atan2(dot_product(s(i)%spin_track, s(i)%y_spin_axis), dot_product(s(i)%spin_track, s(i)%x_spin_axis)) / twopi
  s(i-1)%spin_track_dangle = s(i)%spin_track_angle - s(i-1)%spin_track_angle
enddo

s%track_angle_corrected = .false.
s(0)%track_angle_corrected = .true.
s(0)%spin_track_dangle = closed_orbit_spin_tune + modulo2(s(0)%spin_track_dangle - closed_orbit_spin_tune, 0.5_rp)

call adjust_track_dangle(s(0))

end subroutine transverse_plane_track

!------------------------------------------------------------
! contains

recursive subroutine adjust_track_dangle(s0)

type (spin_track_point_struct) s0
integer kk, i2

! Make sure that the spin phase advance is continuous as a function of orbit phase...

do kk = 1, 6
  i2 = s0%ix_nearest(kk)
  if (i2 == -1) cycle
  if (s(i2)%track_angle_corrected) cycle
  if (i2 == i_turn) then
    s(i2)%spin_track_dangle = s0%spin_track_dangle ! Dummy number to not throw off future adjustments.
  else
    s(i2)%spin_track_dangle = s(i2)%spin_track_dangle + nint(s0%spin_track_dangle - s(i2)%spin_track_dangle)
  endif
  s(i2)%track_angle_corrected = .true.
  call adjust_track_dangle(s(i2))
enddo

end subroutine adjust_track_dangle

!------------------------------------------------------------
! contains

subroutine write_spin_track(file_name)

character(*) file_name
integer i

!

open (1, file = file_name)

write (1, '(a)') '#   Ix     Orbit Phase (x, y, z)         spin_angle  spin_dangle               Spin_Track                        x_spin_axis                        y_spin_axis                      invariant_spin_axis'

do i = 0, i_turn
  write (1, '(i6, 3f10.6, 2x, 2f13.8, 4(5x, 3f10.6))') i, s(i)%orbit_phase, s(i)%spin_track_angle, s(i)%spin_track_dangle, &
                                                                           s(i)%spin_track, s(i)%x_spin_axis, s(i)%y_spin_axis, s(i)%invar_spin
enddo

close (1)

end subroutine write_spin_track

!------------------------------------------------------------
! contains

subroutine rotate_axes_to_flatten_spin_phase_advance (n_order)

real(rp) dangle_ave
real(rp) coef(440), phi_x, phi_y, phi_z, phi1, phi2
integer n_order, i, it, jx, jy, jz, nn, n_coef, n_dim
logical at_end

! Find best fit to given order of the phase advance

dangle_ave = sum(s(0:i_turn-1)%spin_track_dangle) / i_turn
n_dim = count(mode_is_oscillating)

n_coef = (1 + n_order)**n_dim - 1

call initial_lmdif()
coef = 0

do i = 1, 100

  do it = 0, i_turn - 1
    merit_vec(it+1) = s(it)%spin_track_dangle - dangle_ave

    n_coef = 0
    do jx = 0, n_order
      if (.not. mode_is_oscillating(1) .and. jx > 0) cycle
      phi_x = twopi * jx * s(it)%orbit_phase(1)

      do jy = 0, n_order
        if (.not. mode_is_oscillating(2) .and. jy > 0) cycle
        phi_y = twopi * jy * s(it)%orbit_phase(2)

        do jz = 0, n_order
          if (.not. mode_is_oscillating(3) .and. jz > 0) cycle
          if (jx + jy + jz > n_order) cycle
          if (jx + jy + jz == 0) cycle
          phi_z = twopi * jz * s(it)%orbit_phase(3)

          phi1 = phi_x + phi_y + phi_z
          phi2 = phi1 + twopi * (jx * orbit_tune(1) + jy * orbit_tune(2) + jz * orbit_tune(3))
          merit_vec(it+1) = merit_vec(it+1) - coef(n_coef+1) * (sin(phi2) - sin(phi1)) - coef(n_coef+2) * (cos(phi2) - cos(phi1))
          n_coef = n_coef + 2
        enddo
      enddo
    enddo
  enddo

  call suggest_lmdif(coef(1:n_coef), merit_vec(1:i_turn), 1e-8_rp, 100, at_end)
  if (at_end) exit
enddo

! Subtract off the fit

do it = 0, i_turn

  n_coef = 0
  do jx = 0, n_order
    if (.not. mode_is_oscillating(1) .and. jx > 0) cycle
    phi_x = twopi * jx * s(it)%orbit_phase(1)

    do jy = 0, n_order
      if (.not. mode_is_oscillating(2) .and. jy > 0) cycle
      phi_y = twopi * jy * s(it)%orbit_phase(2)

      do jz = 0, n_order
        if (.not. mode_is_oscillating(3) .and. jz > 0) cycle
        if (jx + jy + jz > n_order) cycle
        if (jx + jy + jz == 0) cycle
        phi_z = twopi * jz * s(it)%orbit_phase(3)

        phi1 = phi_x + phi_y + phi_z
        phi2 = twopi * (coef(n_coef+1) * sin(phi1) + coef(n_coef+2) * cos(phi1))
        s(it)%x_spin_axis = rotate_vec_given_axis_angle(s(it)%x_spin_axis, s(it)%invar_spin, phi2)
        s(it)%y_spin_axis = rotate_vec_given_axis_angle(s(it)%y_spin_axis, s(it)%invar_spin, phi2)
        n_coef = n_coef + 2
      enddo
    enddo
  enddo
enddo

end subroutine rotate_axes_to_flatten_spin_phase_advance

!------------------------------------------------------------
! contains

subroutine flip_spin_axis_if_needed_to_align()

! Flip axis if needed to align spins to be more nearly paralled to the closed orbit spin

ave_invar_spin = [sum(s(0:i_turn)%invar_spin(1)), sum(s(0:i_turn)%invar_spin(2)), sum(s(0:i_turn)%invar_spin(3))] / (i_turn + 1)
if (dot_product(ave_invar_spin, closed_orb_invar_spin)  < 0) then
  do i = 0, i_turn
    s(i)%invar_spin = -s(i)%invar_spin
  enddo
endif

end subroutine flip_spin_axis_if_needed_to_align

end program
