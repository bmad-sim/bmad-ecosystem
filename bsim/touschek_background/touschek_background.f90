!--------------------------------------------------------------------------
!+
! Program touschek_background
!
! This program produces a distribution of particles representing the rate
! of particles scattered into a particular momentum deviation range, and
! tracks those particles through the lattice to determine where they
! are lost.  This program takes the output file "aperture_by_slice.out" from 
! touschek_aperture as input.
!
! Modules needed:
!   use touschek_mod
!
!  Use: touschek_background <parameter file>
!
!-
PROGRAM touschek_background

USE bmad
use beam_mod
USE touschek_mod
USE ibs_mod
USE special_collimate_lattice_mod
USE count_lines_in_file_mod
USE slice_mod
USE spline_mod
USE touschek_background_mpi_mod
USE touschek_background_mod


IMPLICIT NONE

!-----------------------------------------------------------------------------------
! Define derived data types
!-----------------------------------------------------------------------------------

TYPE test_part_struct
  REAL(rp) delta_m      !Momentum deviation of test particle: deltap/p at creation
                        !Each test particle has initial phase space coords: (0,0,0,0,0,delta_m)
  REAL(rp) delta_eV     !Momentum deviation of test particle in eV
  INTEGER lost          !Was the test particle lost somewhere?
  INTEGER slix_lost     !If test particle was lost, this contains the slice where it was lost
  INTEGER plane_lost_at !If lost, this contains the dimension the test particle was lost in
                        ! x_plane$, y_plane$ (for apertures), or
                        ! z_plane$ (turned around in an lcavity)
END TYPE test_part_struct

TYPE aperture_by_s_struct
  REAL(rp) s                   !holds locations of slices
  REAL(rp) negative_aperture   !negative momentum aperture at this slice
  REAL(rp) positive_aperture   !positive momentum aperture at this slice
END TYPE aperture_by_s_struct

TYPE generation_tally_struct
  REAL(rp) pipe_this_slice     !number of particles per turn to strike beampipe at this slice
  REAL(rp) pipe_cumulative     !integrated from lattice start to this location number to strike pipe
  REAL(rp) stop_this_slice     !number of particles per turn to stop due to zero energy at this slice
  REAL(rp) stop_cumulative     !integrated from lattice start to this location number to stop
END TYPE generation_tally_struct

TYPE col_hash_struct
  INTEGER ele_ix      ! element number of collimator
  INTEGER slice_ix    ! last slice before end of collimator
END TYPE col_hash_struct

!-----------------------------------------------------------------------------------
! Declare variables
!-----------------------------------------------------------------------------------

!-----------------------------------------
!- Variables that manage program execution
!-----------------------------------------
INTEGER from_address
INTEGER templun
INTEGER part_offset
INTEGER n_sent
INTEGER n_received
INTEGER i, j, k, h, m
INTEGER outer_i
INTEGER slix
INTEGER eleix
INTEGER slix_a
INTEGER slix_b
INTEGER output
INTEGER n_loc
INTEGER funit
INTEGER last_good_traj
INTEGER max_job_size
INTEGER number_of_collimators
INTEGER n_snapshots
INTEGER nslave
INTEGER radcache
INTEGER, ALLOCATABLE :: col_funits(:)
CHARACTER(100) parameter_file
CHARACTER(5) col_num_str
CHARACTER(5) snap_file_str
CHARACTER(121) snap_file_line
REAL(rp) r_bucket
REAL(rp) term
LOGICAL err
LOGICAL stat
LOGICAL master
TYPE(JobDesc_struct),  ALLOCATABLE :: JobDescPackage(:)    !for sending jobs to multiple slaves
TYPE(JobDesc_struct)               :: JobDesc              !for sending jobs to one slave
TYPE(ele_struct) :: stubele
TYPE(spline_struct), ALLOCATABLE :: tou_spline(:)

!-----------------------------------------
!- Variables that manage simulation behavior
!-----------------------------------------
INTEGER tracking
INTEGER count_col_losses
INTEGER N_data_points
INTEGER N_test_part
INTEGER n_slices
integer n_turns
INTEGER number_of_test_particles
INTEGER start_stop_type
INTEGER num_cols
INTEGER slix_prod_start
INTEGER slix_prod_end
INTEGER slix_lost_start
INTEGER hist_bins
INTEGER slix_lost_end
INTEGER snapshot_start_slix
INTEGER snapshot_stop_slix
CHARACTER(100) lat_file
CHARACTER(200) aperture_file
CHARACTER(200) halo_aperture_file
CHARACTER(200) coll_file
CHARACTER(30) name_prod_start
CHARACTER(30) name_prod_end
CHARACTER(30) name_lost_start
CHARACTER(30) name_lost_end
REAL(rp) test_collimator
REAL(rp) bunch_charge
REAL(rp) pz_min
REAL(rp) pz_max
REAL(rp) multiplier
REAL(rp) distParam
REAL(rp) distXmax
REAL(rp) deltaX
REAL(rp) distX
REAL(rp) slice_len
REAL(rp) col_loc
REAL(rp) col_dia
REAL(rp) ignore_thresh
REAL(rp) s_prod_start
REAL(rp) s_prod_end
REAL(rp) s_lost_start
REAL(rp) s_lost_end
REAL(rp) hist_min
REAL(rp) hist_max
REAL(rp) hist_mid
REAL(rp) bin_size
REAL(rp) a_emittance
REAL(rp) b_emittance
REAL(rp) bunch_length
REAL(rp) energy_spread_eV
REAL(rp), PARAMETER :: c_q = 3.84e-13
REAL(rp), ALLOCATABLE :: slices(:)
LOGICAL halo
LOGICAL traj_snapshot
LOGICAL histogram_orbit
LOGICAL do_ibs
LOGICAL collimate
LOGICAL count_loss
LOGICAL coll_files

TYPE(ibs_sim_param_struct) :: ibs_sim_params
TYPE(aperture_by_s_struct), ALLOCATABLE :: aperture_by_slice(:)
TYPE(aperture_by_s_struct), ALLOCATABLE :: halo_aperture_by_slice(:)
TYPE(lat_struct) lat
TYPE(rad_int_all_ele_struct) :: rad_int_ele
TYPE(ele_pointer_struct), ALLOCATABLE :: eles(:)
TYPE(normal_modes_struct) mode
TYPE(normal_modes_struct) dummy_mode
TYPE(collimator_struct), ALLOCATABLE :: collimators(:)
TYPE(col_hash_struct), ALLOCATABLE :: col_hash(:)

!-----------------------------------------
!- Variables that contain data
!-----------------------------------------
INTEGER lost_to_col
INTEGER plane_lost_at
INTEGER cur_partnum
INTEGER cur_lost
INTEGER cur_slixlost
INTEGER cur_planelostat
INTEGER cur_lost_to_col
INTEGER beam_pipe_loss_count
INTEGER bmad_zero_energy_count
INTEGER slix_lost
INTEGER, ALLOCATABLE :: ntp(:)  ! total number of test particles passing through slice.  Used exclusively for snapshots
REAL(rp) piwirate     !Location on y-axis of piwinski touschek curve corresponding to delta_m
REAL(rp) gamma
REAL(rp) rate
REAL(rp) cur_deltam
REAL(rp) rad_delta_eV2
REAL(rp) ibs_delta_eV
REAL(rp) delta_r
REAL(rp) accumulator
REAL(rp) pipe_current
REAL(rp) noenergy_current
REAL(rp) total_current
REAL(rp) total_energy_deposited
REAL(rp) cur_upper_bnd
REAL(rp) cur_ne
REAL(rp) init_vec(6)
REAL(rp), ALLOCATABLE :: col_profile(:)
REAL(rp), ALLOCATABLE :: orbit_data_packed(:)
REAL(rp), ALLOCATABLE :: energy_spread(:)
REAL(rp), ALLOCATABLE :: data_points(:,:)
REAL(rp), ALLOCATABLE :: orbit_hist(:,:,:)
REAL(rp), ALLOCATABLE :: raw_scattering_rates(:,:)
REAL(rp), ALLOCATABLE :: Ndeposited_at_slice(:,:)
REAL(rp), ALLOCATABLE :: eV_at_slice(:)
REAL(rp), ALLOCATABLE :: ne(:)                ! total current of touschek particles passing through slice
REAL(rp), ALLOCATABLE :: sig_dbl_terms(:,:)   ! sig_dbl_terms(1:n_slices,{x_x,x_px,x_y,x_py,x_z,x_pz,px_px,px_y,px_py,px_z,px_pz,
                                              !                         y_y,y_py,y_z,y_pz,py_py,py_z,py_pz,z_z,z_pz,pz_pz})
REAL(rp), ALLOCATABLE :: sig_sgl_terms(:,:)   ! sig_sgl_terms(1:n_slices,{x,y,z,px,py,pz})
REAL(rp), ALLOCATABLE :: sig_avg_dbl_terms(:) ! sig_avg_dbl_terms({x_x,x_px,x_y,x_py,x_z,x_pz,px_px,px_y,px_py,px_z,px_pz,
                                              !                y_y,y_py,y_z,y_pz,py_py,py_z,py_pz,z_z,z_pz,pz_pz
REAL(rp), ALLOCATABLE :: sig_avg_sgl_terms(:) ! sig_avg_sgl_terms({x,px,y,py,z,pz})
REAL(rp), ALLOCATABLE :: energy_deposited_at_slice(:)  !contains energy deposited into each element by beam pipe collisions
TYPE(TrackResult_struct),   ALLOCATABLE :: TrackResults(:)
TYPE(JobResults_struct) :: JobResults
TYPE(orbit_data_struct),    ALLOCATABLE :: orbit_data(:)
TYPE(test_part_struct),     ALLOCATABLE :: test_particles(:)
TYPE(generation_tally_struct), ALLOCATABLE :: generation_at_slice(:)
TYPE(coord_struct), ALLOCATABLE :: slorbit(:)
TYPE(coord_struct), ALLOCATABLE :: orbit(:)
TYPE(coord_struct) initialized_coord_struct

type (beam_init_struct) :: beam_init
type (beam_struct), target :: beam
type (bunch_struct), pointer :: bunch
type (bunch_params_struct) :: bunch_params
logical :: use_beam, verbose

!-----------------------------------------------------------
! Namelise parameters contains variables set by the .in file
!-----------------------------------------------------------
NAMELIST /parameters/ lat_file, &              ! the lattice file
                      n_turns, &               ! for closed geometries, how many turns to track
                      aperture_file, &         ! Contains location of slices and pz momentum aperture. Format:
                                               !  s-position   negative_pz_aperture  positive_pz_aperture 
                                               ! This file can be generated by the aperture_by_tracking program.
                      halo_aperture_file, &    ! necessary if tracking for halo, rather than particle loss
                      halo, &                  ! if true, causes program to track particles with momentum deviations laying between
                                               ! the halo aperture and aperture for particle loss
                      a_emittance, &           ! a-mode emittance to use (normalized)
                      b_emittance, &           ! b-mode emittance to use (normalized)
                      bunch_charge, &          ! charge of each bunch
                      count_col_losses, &      ! 1=ignore losses to zero-length collimators, 2=never ignore, 3=ignore all losses to collimators
                      bunch_length, &          ! bunch length
                      energy_spread_eV, &      ! energy spread in eV
                      beam_init, &
                      use_beam, &
                      traj_snapshot, &         ! should simulation record the trajectory of particles destined to be lost?
                      snapshot_start_slix, &   !      -start recording trajectory here
                      snapshot_stop_slix, &    !      -stop recording trajectory here
                      test_collimator, &       ! radius of the test collimator.  Asks question at every slice:  If a collimator
                                               !      of radius test_colimator were here, how many particles would it catch that would
                                               !      otherwise be lost downstream?
                      collimate, &             ! should the collimators.in file be parsed?
                      coll_files, &            ! should coll_### files be generated?  contain coordinates of particles passing through, but
                                               ! not striking collimators.
                      start_stop_type, &       ! 1=use name_lost_start and name_lost_end to determine where to start and stop producing particles
                                               ! 2=use name_prod_start and name_prod_end to determine where to start and stop producing particles
                      s_prod_start, &          ! location in meters simulation is to start producing particles
                      s_prod_end, &            ! location in meters simulation is to stop producing particles
                      s_lost_start, &          ! only pay attention to particles lost between s_lost_start and s_lost_stop
                      s_lost_end, &            ! only pay attention to particles lost between s_lost_start and s_lost_stop
                      name_prod_start, &       ! first element matching this name is where to start producing particles
                      name_prod_end, &         ! first element matching this name is where to stop producing particles
                      name_lost_start, &       ! first element matching this name is where to start paying attention to losses
                      name_lost_end, &         ! first element matching this name is where to stop paying attention to losses
                      do_ibs, &                ! should IBS be calculated?  propagates sigma_p through linac
                      N_data_points, &         ! Number of delta_m at which to calculate Touschek rate when calculating distribution of
                                               !      scattered particles
                      N_test_part, &           ! Number of test particles Touschek curve is to be divided into
                      histogram_orbit, &       ! Should a histogram of the trajectories of the scattered particles be produced? 
                      ignore_thresh, &         ! Slices which produce scattered particles at a rate less than ignore_thresh are ignored
                      distParam, &             ! Determines how concentrated distribution of test particles is to the momentum aperture
                      verbose                  ! Verbose output to screen
!- Default namelist parameters
histogram_orbit = .false.

!- Default Parameters
distParam = 0.9999
ignore_thresh = 1.0E-5
coll_files = .false.
halo = .false.
histogram_orbit = .false.
do_ibs = .false.
traj_snapshot = .false.
use_beam = .false.
verbose = .false.

!- Read in namelist parameters
CALL GETARG(1,parameter_file)
templun = lunget()
OPEN(templun,FILE=parameter_file,STATUS='old')
READ(templun,NML=parameters)
CLOSE(templun)

!- Count slices in aperture file
CALL count_lines_in_file(aperture_file,n_slices)
ALLOCATE(slices(1:n_slices))
slices = 0.

!- mympi_initialize initializes mpi facilities
CALL mympi_initialize(master,nslave)  ! Calls wrapper to mpi routines which registers with MPI daemon and returns master=.true. if the process is
                                      ! the master.  The master process is designated by the MPI daemon.  All other process are slaves.
                                      ! nslave is the number of slaves in the cluster.

!- calcualte max_job_size
max_job_size = CEILING( (N_test_part*1.0) / (nslave*1.0) )   !the master sends packages of work units to the slaves.  max_job_size
                                                                   !is the largest number of work units to give a slave at one time.

!- mympi_bmad_parser is a simple wrapper for the bmad parser that ensures that the master 
!- process parses the lattice file before the slaves do.  This avoids a situation where all
!- the slaves and the master are trying to write the digested file simultaneously.
CALL mympi_bmad_parser(lat_file, lat)
if (lat%branch(0)%param%geometry .eq. closed$) then
  write(*,*) "WARNING: Support for non-open geometries (e.g. storage rings) is under development."
  write(*,*) "WARNING: Please consider results carefully."
endif
bmad_com%aperture_limit_on = .true.

!- Copy some parameters to there proper location
lat%param%particle = electron$
mode%a%emittance = a_emittance
mode%b%emittance = b_emittance
if (lat%branch(0)%param%geometry .eq. closed$) then
  mode%sigE_E = energy_spread_eV / lat%ele(0)%value(E_TOT$)
endif
lat%param%n_part = bunch_charge / e_charge

!----------------------------------------------------------------------
!- Allocate variables that will be needed by only the master
!----------------------------------------------------------------------
IF(master) THEN
  ALLOCATE(data_points(1:N_data_points,1:2))  ! delta_m and Touschek rate for values of delta_m spanning the momentum aperture to
                                              ! approximated infinity.
  ALLOCATE(test_particles(1:N_test_part))     ! distribution of test particles to track
  ALLOCATE(tou_spline(1:N_data_points))       ! used by spline_akima to get delta_m(rate) from rate(delta_m)
  ALLOCATE(TrackResults(1:N_test_part))       ! results of tracking each test particle
  ALLOCATE(orbit_data(1:N_test_part))         ! trajectory of each test particle
  DO i=1,N_test_part
    ALLOCATE(orbit_data(i)%orb(1:n_slices))   ! trajectory of each test particle
  ENDDO
  ALLOCATE(JobDescPackage(1:nslave))    ! Used for seeding.  Contains descriptions of which particles to track for all slaves
  DO i=1,nslave
    ALLOCATE(JobDescPackage(i)%job(1:max_job_size))  ! describes the test particles the slave is to track
  ENDDO
  ALLOCATE(aperture_by_slice(1:n_slices))
  ALLOCATE(halo_aperture_by_slice(1:n_slices))
  ALLOCATE(eV_at_slice(1:n_slices))
  ALLOCATE(raw_scattering_rates(1:n_slices,1:2))
  ALLOCATE(Ndeposited_at_slice(1:n_slices,0:2))
  ALLOCATE(energy_deposited_at_slice(1:n_slices))
  ALLOCATE(generation_at_slice(1:n_slices))
  ALLOCATE(sig_dbl_terms(1:n_slices,1:21))
  ALLOCATE(sig_sgl_terms(1:n_slices,1:6))
  ALLOCATE(sig_avg_dbl_terms(1:21))
  ALLOCATE(sig_avg_sgl_terms(1:6))
  ALLOCATE(ne(1:n_slices))
  IF(traj_snapshot) THEN
    ALLOCATE(ntp(1:n_slices))
    ntp = 0
  ENDIF

  !- Initialize the arrays where needed
  sig_dbl_terms = 0.
  sig_sgl_terms = 0.
  sig_avg_dbl_terms = 0.
  sig_avg_sgl_terms = 0.
  ne = 0.
  raw_scattering_rates = 0.0_rp
  Ndeposited_at_slice(1:n_slices,0:2) = 0.0_rp
  energy_deposited_at_slice(1:n_slices) = 0.0_rp
  generation_at_slice(1:n_slices)%pipe_this_slice = 0.0_rp
  generation_at_slice(1:n_slices)%pipe_cumulative = 0.0_rp
  generation_at_slice(1:n_slices)%stop_this_slice = 0.0_rp
  generation_at_slice(1:n_slices)%stop_cumulative = 0.0_rp

  !- Other, non-array variable initialization
  number_of_test_particles = 0
  beam_pipe_loss_count = 0
  bmad_zero_energy_count = 0
ENDIF

!- Allocate variables that will be needed by only the slaves
IF(.not. master) THEN
  ALLOCATE(slorbit(1:n_slices))
ENDIF

!- Allocate variables that will be needed by both the master and the slaves
ALLOCATE(JobDesc%job(1:max_job_size))    ! contains descriptions of several test particles for slave to track
ALLOCATE(JobResults%result(1:max_job_size))  ! contains results of slave tracking all the test particles that were sent to it
ALLOCATE(orbit_data_packed(1:(max_job_size*6*n_slices)))  ! contains the 6-dim trajectory data for all test particles a slave tracked.
                                                          ! Stored here as one long real array to facilitate passing by MPI.

!- impose a common tracking method for all non-accelerating cavity elements
!tracking = bmad_standard$ ! linear$, bmad_standard$, taylor$, symp_lie_ptc$
!DO i=1, lat%n_ele_track
!  IF(lat%ele(i)%key .ne. lcavity$) THEN
!    lat%ele(i)%tracking_method = tracking
!  ENDIF
!ENDDO

!- Calculate Twiss parameters for all elements in ring.  The formula for the Touschek scattering rate used by this program 
!- is based on Twiss parameters.
CALL twiss_and_track(lat,orbit)


! Beam tracking, for setting sigma_z only (TODO: Please generalize!!)
if (use_beam) then
  call init_beam_distribution (lat%ele(0), lat%param, beam_init, beam)
  bunch => beam%bunch(1)
  ! Set ele 0
  call calc_bunch_params (bunch, bunch_params, err, print_err = .true.)
  !lat%ele(0)%z%sigma = bunch_params%sigma(5,5)
  
  if (verbose .and. master) then
    write (*, '(a, f15.7, a)') ' charge            : ', 1e12_rp*bunch_params%charge_live, ' pC'
    write (*, '(a, f15.7, a)') ' norm_emit_a       : ', 1e6_rp*bunch_params%a%norm_emit, ' mm-mrad'
    write (*, '(a, f15.7, a)') ' norm_emit_b       : ', 1e6_rp*bunch_params%b%norm_emit, ' mm-mrad'
    write (*, '(a, f15.7, a)') ' sigma_z/c         : ', 1e12_rp*sqrt(bunch_params%sigma(5,5))/c_light, ' ps'
  endif
  
  do i = 1, lat%n_ele_track
    call track1_bunch (bunch, lat%ele(i), err)
    if (err) then
      print *, 'Bunch tracking error in ele: ', trim(lat%ele(i)%name)
      stop
    endif
    ! Don't worry about calc_bunc_params errors
    call calc_bunch_params (bunch, bunch_params, err, print_err = .false.)

    if (verbose .and. master) write(*, '(f12.3, a, f12.3, a)') lat%ele(i)%s, ' sigma_z/c: ', sqrt(bunch_params%sigma(5,5))*(1e15)/c_light, ' fs'

    ! Set element sigma z
    lat%ele(i)%z%sigma = sqrt(bunch_params%sigma(5,5))
  enddo
else
  !- Set bunch length same for every element, based on parameter obtained from .in file.
  DO i=0, lat%n_ele_track
    lat%ele(i)%z%sigma = bunch_length
  ENDDO
  if (lat%branch(0)%param%geometry .eq. closed$) then
    mode%sig_z = lat%ele(0)%z%sigma
  endif
endif

!-----------------------------------------------------------------
!- For high current, low emittance beams, 1-turn IBS growth can be significant (mostly only for energy spread).
!- This block of code calculates 1-turn IBS growth rates and propagates through linac.
!- Only the master process does the IBS calculation.  It then broadcasts the resulting ele-by-ele sigma_p to
!- the slaves.
!- Both IBS and SR cause energy spread to increase.  Here both are calculated and the squares are summed to obtain
!- the total energy spread increase.
!- A file 'sigma_p.out' is written containing the resulting energy spread.
!-----------------------------------------------------------------
IF(master) THEN
  dummy_mode=mode
  radcache = -1
  WRITE(*,*) "Starting radiation integrals..."
  CALL radiation_integrals(lat, orbit, dummy_mode, radcache, 0, rad_int_by_ele = rad_int_ele)
  WRITE(*,*) "Radiation integrals complete..."
  ALLOCATE(energy_spread(0:lat%n_ele_track))
  energy_spread = 0.
  energy_spread(0) = energy_spread_eV   ! parameter injected energy spread
  lat%ele(0)%z%sigma_p = energy_spread(0) / lat%ele(0)%value(E_TOT$)
  IF(do_ibs) WRITE(*,*) "Starting IBS energy spread calculations (this may take a minute)..."
  ibs_sim_params%formula = 'bjmt'
  ibs_sim_params%tau_a = 0
  ibs_sim_params%clog_to_use = 1
  
  DO i=1, lat%n_ele_track
    
    ! Take into account bunch compression by magnifying the previous energy spread by the compression ratio
    energy_spread_eV = energy_spread_eV * lat%ele(i-1)%z%sigma / lat%ele(i)%z%sigma
    
    !find change in energy spread due to radiation effects
    rad_delta_eV2 = (4./3.)*c_q*r_e*(m_electron**2)*rad_int_ele%branch(0)%ele(i)%lin_i3_E7
  
    !find change in energy spread due to ibs effects
    CALL convert_total_energy_to(lat%ele(i)%value(E_TOT$), lat%param%particle, gamma)
    lat%ele(i)%a%emit = mode%a%emittance / gamma
    lat%ele(i)%b%emit = mode%b%emittance / gamma
    lat%ele(i)%z%sigma_p = lat%ele(i-1)%z%sigma_p  !use sigma_p at start of ele to calc ibs blowup
    IF(do_ibs) THEN
      CALL ibs_delta_calc(lat, i, ibs_sim_params, delta_sigma_energy = ibs_delta_eV)
    ELSE
      ibs_delta_eV = 0.
    ENDIF
    energy_spread(i) = SQRT((energy_spread(i-1)+ibs_delta_eV)**2 + rad_delta_ev2)
    lat%ele(i)%z%sigma_p = energy_spread(i)/lat%ele(i)%value(E_TOT$)
  ENDDO
  DEALLOCATE(energy_spread)
  OPEN(82,FILE='sigma_p.out')
  WRITE(82,'(A41)') "#  location       sigma_p    Total energy"
  WRITE(82,'(A41)') "#       (s)       (dp/p0)            (eV)"
  DO i=0,lat%n_ele_track
    WRITE(82,'(F11.3,ES14.4, ES16.4)') lat%ele(i)%s, lat%ele(i)%z%sigma_p, lat%ele(i)%value(E_TOT$)
  ENDDO
  CLOSE(82)
  WRITE(*,*) "Energy spread calculations complete..."
ENDIF
! Broadcast energy spread calculation results from master to all slaves
CALL mympi_bcast_sigmap(lat)  !mpi wrapper that takes ele-by-ele energy spread from master process and
                              !broadcasts to all slaves.

!-----------------------------------------------------------------
!- The file aperture_by_s.out is generated by the companion program touschek_aperture.
!- It contains the momentum aperture for various locations s.  This program assumes
!- that aperture_by_s.out has already been generated and is present in the current
!- working directory.
!-----------------------------------------------------------------
IF(master) THEN
  !if master, then need to know slice locations and momentum aperture
  !read in slice locations and momentum aperture
  OPEN(79,FILE=TRIM(aperture_file),STATUS="OLD")
  DO i=1, n_slices
    READ(79,*) aperture_by_slice(i)%s, aperture_by_slice(i)%negative_aperture, aperture_by_slice(i)%positive_aperture
    slices(i)=aperture_by_slice(i)%s
    if (slices(i) > lat%param%total_length .or. slices(i) < 0) then
      print '(a, f13.3)', 'APERTURE S-POSITION FROM APERTURE_FILE OUT-OF-BOUNDS:', slices(i)
      print '(2a)',       'FROM FILE: ', trim(aperture_file)
      stop
    endif
  ENDDO
  slices(n_slices) = MIN(slices(n_slices),lat%param%total_length-0.000001)
  CLOSE(79)

  IF(halo) THEN
    !if running a halo simulation, also parse the halo aperture file
    !This program assumes that the halo aperture file and the aperture file for particle loss have the same slice indexing.
    !If they do not have the same indexing, things will break.
    OPEN(79,FILE=TRIM(halo_aperture_file),STATUS="OLD")
    DO i=1, n_slices
      READ(79, *) halo_aperture_by_slice(i)%s, halo_aperture_by_slice(i)%negative_aperture, &
                                                                 halo_aperture_by_slice(i)%positive_aperture
    ENDDO
    CLOSE(79)
  ENDIF

  !also calculate the beam energy at each slice location.
  eV_at_slice = 0.
  DO i=1,n_slices
    CALL twiss_and_track_at_s(lat,slices(i),stubele)
    eV_at_slice(i) = stubele%value(E_TOT$)
  ENDDO
ELSE
  !if slave, then all we need to know is the slice locations
  OPEN(79,FILE=TRIM(aperture_file),STATUS="OLD")
  DO i=1, n_slices
    READ(79, *) slices(i), r_bucket, r_bucket
  ENDDO
  CLOSE(79)
ENDIF
slices(n_slices) = slices(n_slices) - 0.000001_rp

!-----------------------------------------------------------------
!- Determine start and stop locations for paticle production.
!- Also determine start and stop locations for where to pay attention to losses.
!-----------------------------------------------------------------
IF(start_stop_type .eq. 1) THEN
  ! if locations are given in meters, this is easy and we are done.  Simply check
  ! for negative values, which indicate that the entire lattice is to be used.
  IF(s_prod_start .lt. 0.0) s_prod_start = 0.0
  IF(s_lost_start .lt. 0.0) s_lost_start = 0.0
  IF(s_prod_end .lt. 0.) s_prod_end = lat%param%total_length
  IF(s_lost_end .lt. 0.) s_lost_end = lat%param%total_length
ELSEIF(start_stop_type .eq. 2) THEN
  ! if locations given by element name, then use location of first element
  ! in lattice that matches that name.

  CALL lat_ele_locator(name_prod_start,lat,eles,n_loc,err) 
  IF(n_loc .gt. 0) THEN
    s_prod_start = lat%ele(eles(1)%ele%ix_ele)%s
  ELSE
    WRITE(*,*) "FATAL ERROR: name_prod_start not found"
    WRITE(*,*) "             name_prod_start = ", name_prod_start
    STOP
  ENDIF

  CALL lat_ele_locator(name_prod_end,lat,eles,n_loc,err) 
  IF(n_loc .gt. 0) THEN
    s_prod_end = lat%ele(eles(1)%ele%ix_ele)%s
  ELSE
    WRITE(*,*) "FATAL ERROR: name_prod_end not found"
    WRITE(*,*) "             name_prod_end = ", name_prod_end
    STOP
  ENDIF

  CALL lat_ele_locator(name_lost_start,lat,eles,n_loc,err) 
  IF(n_loc .gt. 0) THEN
    s_lost_start = lat%ele(eles(1)%ele%ix_ele)%s
  ELSE
    WRITE(*,*) "FATAL ERROR: name_lost_start not found"
    WRITE(*,*) "             name_lost_start = ", name_lost_start
    STOP
  ENDIF

  CALL lat_ele_locator(name_lost_end,lat,eles,n_loc,err) 
  IF(n_loc .gt. 0) THEN
    s_lost_end = lat%ele(eles(1)%ele%ix_ele)%s
  ELSE
    WRITE(*,*) "FATAL ERROR: name_lost_end not found"
    WRITE(*,*) "             name_lost_end = ", name_lost_end
    STOP
  ENDIF
ENDIF

!-----------------------------------------------------------------------
!- Find the slice indexes closest to the user-specified start and stop locations.
!- These indexes are the start and stop locations that will actually be used
!-----------------------------------------------------------------------
!- The following loop obtains the first slice before s_prod_start
!- slix_prod_start has a lower bound of 2 because we do not know the
!- length of slice 1. slice_len(i) = slices(i) - slices(i-1)
slix_prod_start = 2
DO i=3,n_slices
  IF(slices(i) .ge. s_prod_start) EXIT
  slix_prod_start = i
ENDDO
!- The following loop obtains the first slice after s_prod_end
DO i=slix_prod_start,n_slices
  slix_prod_end = i
  IF(slices(slix_prod_end) .ge. s_prod_end) EXIT
ENDDO
slix_prod_end = MIN(slix_prod_end,n_slices-1)
!- The following loop obtains the first slice before s_lost_start
slix_lost_start = 2
DO i=3,n_slices
  IF(slices(i) .ge. s_lost_start) EXIT
  slix_lost_start = i
ENDDO
!The following loop obtains the first slice after s_lost_end
DO i=slix_lost_start,n_slices
  slix_lost_end = i
  IF(slices(slix_lost_end) .ge. s_lost_end) EXIT
ENDDO
slix_lost_end = MIN(slix_lost_end,n_slices)
slix_prod_end = MIN(slix_prod_end,slix_lost_end-1)

!------------------------------------------------------------------
!- Allocate space for collimator profile file, which store information
!- about hypothetical collimator locations.
!------------------------------------------------------------------
IF(master) THEN
  ALLOCATE(col_profile(slix_prod_start:slix_lost_end))
  col_profile(slix_prod_start:slix_lost_end) = 0.0_rp
ENDIF

!-----------------------------------------------------------------
! The collimators.in file is a way to add collimators to the lattice
! via this simulation.  It does not change the lattice file.  It simply adds
! collimator elements to the in-memory lattice data structure.
! The format for the collimators.in file is: 
!         <collimator location in meters>   <collimator diameter>
!-----------------------------------------------------------------
IF(collimate) THEN
  coll_file = 'collimators.in'
  CALL count_lines_in_file(coll_file,num_cols)
  ALLOCATE(collimators(1:num_cols))
  OPEN(33,FILE=coll_file)
  DO i=1,num_cols
    READ(33,*) col_loc, col_dia
    collimators(i) = collimator_struct(col_loc,col_dia)
  ENDDO
  CLOSE(33)
  CALL special_collimate_lattice(lat,collimators)
  DEALLOCATE(collimators)
ENDIF

!-----------------------------------------------------------------
!- Note about file unit numbers.  Pretains to OPEN, CLOSE, WRITE, UNIT directives.
!- Collimator file units may range from 101 to 100+lat%n_ele_track.  
!- Snapshot file units may range from 101+lat%n_ele_track to 101+2*lat%n_ele_track.  
!- All other file units may range from 11 to 100.
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!- Open a log file for each collimator.  Filename is coll_<eleix>.dist
!- Also build hash table of collimator element indexes and slice locations.
!-----------------------------------------------------------------
IF(master) THEN
  IF(coll_files) THEN
    ! count number of collimators in lattice
    number_of_collimators = 0
    DO i=1, lat%n_ele_track
      IF( (lat%ele(i)%key .eq. rcollimator$) .or. (lat%ele(i)%key .eq. ecollimator$) ) THEN
        number_of_collimators = number_of_collimators + 1
      ENDIF
    ENDDO

    IF(number_of_collimators .gt. 0) THEN
      !-------------------------------------------------------------------------
      !- Build col_hash hash table.  This table has one entry for each collimator
      !- in the lattice.  Each entry contains the element index of the collimator
      !- along with the slice index of the first slice preceeding the collimator.
      !- Also open log file for each collimator.
      !-------------------------------------------------------------------------
      ALLOCATE(col_funits(1:number_of_collimators))
      ALLOCATE(col_hash(1:number_of_collimators))
      !- indexes in the loop that follows:
      !- k indexes collimators
      !- i indexes elements
      !- j indexes slices
      k=0
      DO i=1, lat%n_ele_track
        IF( (lat%ele(i)%key .eq. rcollimator$) .or. (lat%ele(i)%key .eq. ecollimator$) ) THEN
          k=k+1

          ! open log file for each collimator
          col_funits(k) = 100+k
          WRITE(col_num_str,'(I5.5)') i
          OPEN(col_funits(k), FILE='coll_'//col_num_str//'.dist')

          ! write header information
          WRITE(col_funits(k), '(6A14,A14,A15)') "x", "xp", "y", "yp", "z", "zp", &
                                                "charge in tp", "  Lost to Col?"
          !- populate col_hash with indexes of collimators
          col_hash(k)%ele_ix = i

          !- populate col_hash with index of slice following collimator entrance
          DO j=1, n_slices
            IF(slices(j) .gt. lat%ele(i)%s-lat%ele(i)%value(l$)) THEN
              col_hash(k)%slice_ix = j
              EXIT
            ENDIF
          ENDDO
        ENDIF
      ENDDO
    ENDIF
  ENDIF
ENDIF

IF(master) THEN
  WRITE(*,*)
  WRITE(*,*) "Touschek particle production set to start at: ", s_prod_start
  WRITE(*,*) "Touschek particle production set to stop at:  ", s_prod_end
  WRITE(*,*) "Touschek particle losses set to start at:     ", s_lost_start
  WRITE(*,*) "Touschek particle losses set to stop at:      ", s_lost_end
  WRITE(*,*)
ENDIF

!----------------------------------------------------------------
!- Open files for snapshots.  Snapshot files exist for every slice between
!- snapshot_start_slix and snapshot_stop_slix.  If a tracked particle passes
!- through a slice with a snapshot file, the particles coordinates are recorded.
!----------------------------------------------------------------
IF(master) THEN
  IF(traj_snapshot) THEN
    n_snapshots = snapshot_stop_slix - snapshot_start_slix + 1
    DO i=1,n_snapshots
      WRITE(snap_file_str,'(I5.5)') snapshot_start_slix+i-1
      funit = 101+lat%n_ele_track+snapshot_start_slix+i-1
      OPEN(funit,FILE='snapshot_'//snap_file_str//'.dist_temp')
      WRITE(funit,'(I8,A15)')  snapshot_start_slix+i-1,"  ! slice index"
      WRITE(funit,'(A12)')    "1  ! n_bunch"
      WRITE(funit,'(A16)')    "###  !n_particle"  !this line will be overwritten
      WRITE(funit,'(A13)')    "BEGIN_BUNCH"
      WRITE(funit,'(A34)')    "#.############E-###  ! bunch_charge"  !this line will be overwritten
      WRITE(funit,'(A40)')    "0.000000000E+000  ! z_center set to zero"
      WRITE(funit,'(A40)')    "0.000000000E+000  ! t_center set to zero"
    ENDDO
  ENDIF
ENDIF

!----------------------------------------------------------------------
!- Write out two log files.  One containing slice and name information. The
!- other containing physical apertures.
!----------------------------------------------------------------------
IF(master) THEN
  !- slice_index_from_track.out basically regurgitates slice information obtained
  !- from the aperture file, along with the index of the slice and name of ele it
  !- is inside.
  OPEN(59,FILE="slice_index_from_track.out")
  WRITE(59,'(A43)') "#     slice      location   element        "
  WRITE(59,'(A43)') "#     index           (s)   name           "
  DO i=1,n_slices
    eleix = element_at_s(lat,slices(i), .true.)
    WRITE(59,'(I11,F14.5,"   ",A)') i, slices(i), lat%ele(eleix)%name
  ENDDO
  CLOSE(59)

  !- physical_aperture_from_track.out contains the horizontal physical aperture 
  !- at each element.
  OPEN(59,FILE="physical_aperture_from_track.out")
  WRITE(59,'(A61)') "#       ele     location     aperture   element name         "
  WRITE(59,'(A61)') "#     index          (s)       radius                        "
  DO i=1, lat%n_ele_track
    IF( (lat%ele(i)%value(l$) .gt. 0.0) .or. \
        (lat%ele(i)%key .eq. ecollimator$) .or. \
        (lat%ele(i)%key .eq. rcollimator$) ) THEN
      WRITE(59,'(I11,F13.5,F13.5,"   ",A)') i, lat%ele(i)%s, lat%ele(i)%value(x1_limit$), lat%ele(i)%name
    ENDIF
  ENDDO
  CLOSE(59)
ENDIF

!--------------------------------------------------------------------------
!-  For slaves, override aperture_at property for every element, so that 
!-  particles lost in slices that are partway through an element are
!-  property flagged as lost.
!--------------------------------------------------------------------------
IF(.not. master) THEN
  DO i=1,lat%n_ele_track
    lat%ele(i)%aperture_at = continuous$
  ENDDO
ENDIF

!--------------------------------------------------------------------------
!- Hard-coded parameters for the orbit histogram.  The only reason for not
!- including these in the .in file is to reduce clutter.
!--------------------------------------------------------------------------
IF(master) THEN
  IF(histogram_orbit) THEN
    !orbit histogram parameters
    hist_min = -0.02_rp
    hist_max =  0.02_rp
    hist_bins = 41
    bin_size = (hist_max-hist_min)/hist_bins
    ALLOCATE(orbit_hist(slix_prod_start:slix_lost_end,1:hist_bins,1:2))
    orbit_hist = 0.
  ENDIF
ENDIF

!--------------------------------------------------------------------------
!- PREPARATION WORK COMPLETE, MAIN LOOP FOLLOWS
!-
!- What follows are the main loops of the program.  One loop is followed only
!- by the master process and loops over all slices.  The other is followed
!- by the slave processes and loops over work units received from the master.
!-
!- The master loop loops over all slices between slix_prod_start and slix_prod_end.
!- At each slice, Piwinski's formula for the Touschek scattering rate is evaluated
!- at a range of momentum deviations between the momentum aperture and approximate
!- infinity to produce a distribution of scattered particles generated at the slice.
!- These distributions are farmed out to the slaves, which track them through the
!- lattice to determine where they are lost.
!-
!- The slave loops block until a message arrives from the master.  The message
!- contains a slice location and a number of test particles to track.  These
!- test particles are tracked to determine where they are lost.  The results
!- of the tracking, including the slice-by-slice trajectory are returned to
!- the master
!--------------------------------------------------------------------------

IF(master) THEN
  DO slix=slix_prod_start, slix_prod_end
    slice_len = slices(slix) - slices(slix-1)
    CALL twiss_and_track_at_s(lat,slices(slix),stubele)   !get Twiss parameters at this slice
    ! k=1 for to evaluate positive momentum aperture.
    ! k=2 for to evaluate negative momentum aperture.
    DO k=1,2
      IF(halo) THEN
        !- pz_min is the halo aperture
        !- pz_max is the aperture for particle loss
        IF (k .eq. 1) THEN
          pz_min = halo_aperture_by_slice(slix)%positive_aperture
          pz_max = aperture_by_slice(slix)%positive_aperture
        ELSE
          pz_min = abs(halo_aperture_by_slice(slix)%negative_aperture)
          pz_max = abs(aperture_by_slice(slix)%negative_aperture)
        ENDIF
        IF(pz_min .eq. 0.0) CYCLE   !momentum aperture flagged by aperture program as to be ignored.
        IF(pz_min .gt. 0.09) CYCLE  !momentum aperture too large to have significant scattering rate
        IF(pz_max .eq. 0.0) THEN    !upper momentum aperture too large to have significant rate.  Search for 
                                    !new upper bound that does have a significant rate.
          pz_max=pz_min*8.0
        ENDIF

        mode%pz_aperture = pz_max
        CALL touschek_rate1(mode,rate,lat,s=slices(slix))
        DO WHILE(rate*slice_len/c_light .lt. 1.0E-30)
          !bring the max aperture 20% of the way towards the min aperture
          pz_max = pz_min + (pz_max-pz_min)*0.10
          mode%pz_aperture = pz_max
          CALL touschek_rate1(mode,rate,lat,s=slices(slix))
          rate = rate /2.0
          IF( (pz_max-pz_min)/pz_min .lt. 0.1 ) THEN
            ! upper bound has been brought too close to the lower bound.  results for this slice will be
            ! insignificant.
            EXIT
          ENDIF
        ENDDO

        ! this loop checks that the lower and upper momentum apertures are different enough to produce
        ! significant results.
        IF( (pz_max-pz_min)/pz_min .lt. 0.0100001 ) CYCLE
      ELSE
        !pz_min is the aperture for particle loss
        !pz_max is set as a multiple of pz_min
        IF (k .eq. 1) THEN
          pz_min = aperture_by_slice(slix)%positive_aperture
        ELSE
          pz_min = abs(aperture_by_slice(slix)%negative_aperture)
        ENDIF
        IF(pz_min .eq. 0.0) CYCLE   !momentum aperture flagged by aperture program as to be ignored.
        IF(pz_min .gt. 0.09) CYCLE  !momentum aperture too large to have significant scattering rate

        !find p_max as a multiple of p_min.
        multiplier = 12.0
        DO WHILE(multiplier .gt. 1.1)
          pz_max = pz_min * multiplier
          mode%pz_aperture = pz_max
          CALL touschek_rate1(mode,rate,lat,s=slices(slix))
          rate = rate /2.0
          IF( rate*slice_len/c_light .gt. 1.0E-30 ) THEN
            EXIT
          ELSE
            multiplier = multiplier - 1.0
          ENDIF
        ENDDO
        IF(multiplier .lt. 1.1) THEN
          CYCLE
        ENDIF
      ENDIF

      !- touschek_rate1 evaluated at the momentum aperture is the total rate of particles
      !- produced at this slice.
      !- Note: Touschek rate returned is in particles per second.
      mode%pz_aperture = pz_min
      CALL touschek_rate1(mode, rate, lat, s=slices(slix))
      rate = rate / 2.0 !touschek_rate1 assumes two particles are lost per scattering event.
      raw_scattering_rates(slix,k) = rate

      !Continue only of there is a non-negligible number of Touschek particles generated at slice slix
      IF (rate*slice_len/c_light .gt. ignore_thresh) THEN
        !-----------------------------------------------------------------------------------------
        !FIRST:  Generate Test Particle Distribution
        !-----------------------------------------------------------------------------------------

        !- This distribution is range of momentum deviations from pz_min to pz_max.
        !- In short, the distribition is the density of coils of a spring stretched
        !- between pz_min and pz_max, where the spring constant k varies linearly
        !- and is largest near pz_min.
        !- See Ehrlichman lab notebook 2 page 27 for the details.
        distXmax = 1.0_rp/distParam*LOG(1.0_rp-distParam)
        deltaX = distXmax/(N_data_points-1)
        DO i=1,N_data_points
          distX = deltaX*(i-1)
          data_points(N_data_points-i+1,1) = \
            pz_min+(pz_max-pz_min)*(1.0_rp+1.0_rp/distParam*(EXP(distParam*distX)-1.0_rp))
        ENDDO
        !- data_points(1,1) is pz_min  (closest to aperture)
        !- data_points(N_data_points,1) is pz_max  (furthest from aperture)

        !- For each point in the distribution, calculate the scattering rate
        DO i=1, N_data_points
          mode%pz_aperture = data_points(i,1)
          CALL touschek_rate1(mode, rate, lat, s=slices(slix))
          rate = rate / 2.0 !touschek_rate1 assumes two particles are lost per scattering event.
          data_points(i,2) = rate
        ENDDO

        !- Rearrange the data points so that we can invert them using an Akima spline.
        !- Spline is used to obtain delta_m(Rate) from Rate(delta_m)
        DO i=1,N_data_points
          tou_spline(N_data_points-i+1)%x0 = data_points(i,2)
          tou_spline(N_data_points-i+1)%y0 = data_points(i,1)
        ENDDO

        !- spline_akima populates tou_spline(:)%coef(0:3), which then
        !- allows spline_evaluate to be used to interpolate the data.
        CALL spline_akima(tou_spline,stat)
        IF(.not.stat) THEN
          WRITE(*,*) "WARNING: spline_akima error!"
        ENDIF

        !- calculate rate represented by each test particle
        !- Each test particle represents the same scattering rate, which is the total
        !- scattering rate at the slice divided by the number of test particles.
        delta_r = (data_points(1,2)-data_points(N_data_points,2)) / N_test_part   ! R_max / (number of test particles)
        cur_ne = delta_r*slice_len/c_light

        !- determine delta corresponding to each rate by interpolation of inverted spline
        DO i=1, N_test_part
          test_particles(i)%lost = 0
          piwirate = data_points(1,2) - (i-1+0.5)*delta_r
          CALL spline_evaluate(tou_spline, piwirate, stat, y=test_particles(i)%delta_m)
          test_particles(i)%delta_m = ((-1)**(k+1))*test_particles(i)%delta_m
          IF(.not.stat) THEN
            WRITE(*,*) "FATAL: spline_evaluate error!"
            STOP
          ENDIF
          test_particles(i)%delta_eV = test_particles(i)%delta_m * eV_at_slice(slix)
        ENDDO

        !-------------------------------------------------------------------------------------
        !SECOND:  Track test particles
        !-------------------------------------------------------------------------------------

        !Step 1 of 2: Seed slaves with test particles to track

        !Allocate jobs to slaves by round-robin
        n_sent=0
        JobDescPackage(1:nslave)%count = 0
        DO i=1,max_job_size
          DO j=1,nslave
            IF( n_sent .lt. N_test_part ) THEN
              n_sent=n_sent+1
              JobDescPackage(j)%count = JobDescPackage(j)%count + 1
              JobDescPackage(j)%job(JobDescPackage(j)%count)%partnum = n_sent
              JobDescPackage(j)%job(JobDescPackage(j)%count)%slix    = slix
              JobDescPackage(j)%job(JobDescPackage(j)%count)%deltam  = test_particles(n_sent)%delta_m
            ELSE
              EXIT
            ENDIF
          ENDDO
          IF( n_sent .eq. N_test_part ) THEN
            EXIT
          ENDIF
        ENDDO

        !- Send the job packages to the slaves for tracking
        DO i=1,nslave
          IF(JobDescPackage(i)%count .gt. 0) THEN
            CALL mympi_send_job_desc(JobDescPackage(i),i)
          ENDIF
        ENDDO

        !Step 2 of 2: Receive results from slaves and send remaining work units to slaves as work units are completed
 
        n_received = 0  ! holds number of results that have been received
        DO WHILE(.true.)
          CALL mympi_check_mailbox_with_timeout(from_address, JobResults%count) !Blocks until something arrives in mailbox.
                                                                                !from_address is id of slave the mail is from

          CALL mympi_receive_results(from_address, n_slices, JobResults, orbit_data_packed)

          !- Arguablly, the unpacking of data, which is done in the loop below, could be moved into the MPI module.  However, 
          !- it is left here because is useful to see that the slaves are returning.
          DO j=1,JobResults%count
            n_received = n_received + 1
            TrackResults(n_received)%partnum     = JobResults%result(j)%partnum
            TrackResults(n_received)%lost        = JobResults%result(j)%lost
            TrackResults(n_received)%slixlost    = JobResults%result(j)%slixlost
            TrackResults(n_received)%planelostat = JobResults%result(j)%planelostat
            TrackResults(n_received)%lost_to_col = JobResults%result(j)%lost_to_col

            part_offset = (j-1)*6*n_slices
            orbit_data(n_received)%orb(1:n_slices)%vec(1) = orbit_data_packed(part_offset+0*n_slices+1:part_offset+1*n_slices)
            orbit_data(n_received)%orb(1:n_slices)%vec(2) = orbit_data_packed(part_offset+1*n_slices+1:part_offset+2*n_slices)
            orbit_data(n_received)%orb(1:n_slices)%vec(3) = orbit_data_packed(part_offset+2*n_slices+1:part_offset+3*n_slices)
            orbit_data(n_received)%orb(1:n_slices)%vec(4) = orbit_data_packed(part_offset+3*n_slices+1:part_offset+4*n_slices)
            orbit_data(n_received)%orb(1:n_slices)%vec(5) = orbit_data_packed(part_offset+4*n_slices+1:part_offset+5*n_slices)
            orbit_data(n_received)%orb(1:n_slices)%vec(6) = orbit_data_packed(part_offset+5*n_slices+1:part_offset+6*n_slices)
          ENDDO

          IF(n_received .ge. N_test_part) THEN
            ! All test particles have been received
            EXIT
          ENDIF

          IF(n_sent .lt. N_test_part) THEN
            JobDesc%count = 0
            DO j=1, max_job_size
              n_sent=n_sent+1
              JobDesc%count = JobDesc%count + 1
              JobDesc%job(JobDesc%count)%partnum = n_sent
              JobDesc%job(JobDesc%count)%slix = slix
              JobDesc%job(JobDesc%count)%deltam = test_particles(n_sent)%delta_m

              IF(n_sent > N_test_part) EXIT
            ENDDO
            CALL mympi_send_job_desc(JobDesc,from_address)
          ENDIF
        ENDDO

        !------------------------------------------------------------------------------------
        !DO i=1, N_test_part
        !  WRITE(43,'(F14.4,ES14.4,I4)') slices(slix), orbit_data(i)%orb(slix)%vec(6), TrackResults(i)%lost
        !ENDDO
        !------------------------------------------------------------------------------------

        !------------------------------------
        CALL raw_touschek_analysis(stubele, cur_ne, TrackResults, orbit_data)
        !------------------------------------

        !-------------------------------------------------------------------------------------
        !THIRD:  Analyze trajectory of each test particle
        !-------------------------------------------------------------------------------------
        DO outer_i=1, N_test_part
          !- make short-hand variables for test particle being examined
          cur_partnum     = TrackResults(outer_i)%partnum
          cur_lost        = TrackResults(outer_i)%lost
          cur_slixlost    = TrackResults(outer_i)%slixlost
          cur_planelostat = TrackResults(outer_i)%planelostat
          cur_lost_to_col = TrackResults(outer_i)%lost_to_col

          !- If particle was lost, then the orbit data is valid from slix to cur_slixlost
          !- If the particle was not lost, then the orbit data is valid from slix to slix_lost_end
          IF(cur_lost .eq. 1) THEN
            last_good_traj = cur_slixlost
          ELSE
            last_good_traj = slix_lost_end
          ENDIF

          !- If we are taking snapshots of the slice-by-slice trajectory
          !- and if the particle trajectory intersects and
          !- the range over which we are taking snapshots, then
          !- write test particle data to the trajectory file for each slice.
          !- Data written consists 6D coords, current in amps represented by
          !- the test particle, and where it was lost.                       
          !- ntp is accumulated for use in snapshot file header written at end of simulation.
          IF(traj_snapshot) THEN                               
            IF(slix .lt. snapshot_stop_slix) THEN              
              IF(last_good_traj .ge. snapshot_start_slix) THEN 
                slix_a = MAX(snapshot_start_slix,slix)         
                slix_b = MIN(snapshot_stop_slix,last_good_traj)
                DO i=slix_a,slix_b
                  ntp(i-snapshot_start_slix+1) = ntp(i-snapshot_start_slix+1) + 1
                  WRITE(101+lat%n_ele_track+i,'(6ES14.5,ES14.4,I7,4F4.0)') \
                        orbit_data(outer_i)%orb(i)%vec(1:6), cur_ne*e_charge, cur_slixlost, 0.,0.,0.,0.
                ENDDO
              ENDIF
            ENDIF
          ENDIF

          !- Loop over all collimators.  If the trajectory of the test particle passes through
          !- the collimator, then record the 6D coords of the test particle, its current, and whether it was
          !- lost to a collimator.
          IF(coll_files) THEN
            DO i=1, number_of_collimators
              IF(lat%ele(col_hash(i)%ele_ix)%s .gt. slices(slix)) THEN
                IF( (lat%ele(col_hash(i)%ele_ix)%s-lat%ele(col_hash(i)%ele_ix)%value(l$)) .le. slices(last_good_traj)) THEN
                  WRITE(col_funits(i),'(6ES14.4,ES14.4,I10)') \
                          orbit_data(outer_i)%orb(col_hash(i)%slice_ix)%vec(1:6), cur_ne*e_charge, cur_lost_to_col
                ENDIF
              ENDIF
            ENDDO
          ENDIF

          IF (cur_lost .eq. 1) THEN  !- if test particle was lost

            !- First determine whether or not to count the particle loss.  Losses to non-collimator elements
            !- are always counted.  The .in parameter count_col_losses determines how losses to collimator
            !- elements are counted.
            IF(cur_lost_to_col .lt. 0) THEN  
              !- particle was not lost to a collimator
              count_loss = .true.
            ELSE 
              !- particle was lost to a collimator
              IF(count_col_losses .eq. 1) THEN 
                !- count as a loss only if the collimator is not zero length
                IF(lat%ele(cur_lost_to_col)%value(l$) .gt. 0.) THEN
                  count_loss = .true.
                ELSE
                  count_loss = .false.
                ENDIF
              ELSEIF(count_col_losses .eq. 2) THEN 
                !- always count as a loss
                count_loss = .true.
              ELSE !never count as a loss
                count_loss = .false.
              ENDIF
            ENDIF

            IF(count_loss) THEN
              !- The parameters lost_start and lost_end tell the program we are only interested in losses
              !- that occur between the elements numbered lost_start and lost_end, inclusive.
              IF((cur_slixlost .ge. slix_lost_start).AND.(cur_slixlost .le. slix_lost_end)) THEN
                test_particles(cur_partnum)%lost = 1
                test_particles(cur_partnum)%plane_lost_at = cur_planelostat
                test_particles(cur_partnum)%slix_lost = cur_slixlost

                !- Check is particle was lost due to collision with beampipe (as opposed to stopping
                !- with zero energy).
                IF((cur_planelostat .eq. x_plane$).OR.(cur_planelostat .eq. y_plane$)) THEN
                  DO i=slix, cur_slixlost-1
                    !- Loop over each slice along the test particles trajectory

                    !- sigma matrix calculations
                    m=1
                    DO j=1,6
                      DO h=j,6
                        sig_dbl_terms(i,m)=sig_dbl_terms(i,m)+orbit_data(outer_i)%orb(i)%vec(j)*orbit_data(outer_i)%orb(i)%vec(h)*cur_ne
                        m=m+1
                      ENDDO
                    ENDDO
                    DO j=1,6
                      sig_sgl_terms(i,j)=sig_sgl_terms(i,j)+orbit_data(outer_i)%orb(i)%vec(j)*cur_ne
                    ENDDO
                    ne(i) = ne(i) + cur_ne

                    !- Produce a profile for potential collimators
                    !- Asks question:  Would a collimator of radius test collimator places at this slice have stopped the particle?
                    IF(ABS(orbit_data(outer_i)%orb(i)%vec(1)) .ge. test_collimator) THEN
                      col_profile(i) = col_profile(i) + cur_ne
                    ENDIF

                    !- Add this test particle to the slice-by-slice histogram of the horizontal coordinate of test particles
                    !- passing through a slice.
                    IF(histogram_orbit) THEN
                      DO j=1, hist_bins
                        cur_upper_bnd = hist_min + bin_size*j
                        IF(orbit_data(outer_i)%orb(i)%vec(1) .lt. cur_upper_bnd) THEN
                          orbit_hist(i,j,1) = orbit_hist(i,j,1) + cur_ne
                          orbit_hist(i,j,2) = orbit_hist(i,j,2) + \
                                              (test_particles(cur_partnum)%delta_eV+eV_at_slice(i))*cur_ne
                          EXIT
                        ENDIF
                      ENDDO
                    ENDIF
                  ENDDO
                ENDIF

                !- Recall that slix is indexing the slice at which the simulation is currently generating test
                !- particles.  Add the current represented by this test particle to the total current of scattered
                !- particles generated at this slice.  Differentiate between particles scattered such that they collide
                !- with the beampipe and those that stop due to zero energy during deceleration.
                IF((cur_planelostat .eq. x_plane$).OR.(cur_planelostat .eq. y_plane$)) THEN
                  generation_at_slice(slix)%pipe_this_slice = generation_at_slice(slix)%pipe_this_slice + cur_ne
                ELSE
                  generation_at_slice(slix)%stop_this_slice = generation_at_slice(slix)%stop_this_slice + cur_ne
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDDO

        !--------------------------------------------------------------
        ! FOURTH: Produce loss profiles
        !--------------------------------------------------------------

        !- Calculate the rate of particles lost at each element, distinguishing between
        !- beampipe collisions and zero energy during deceleration
        !- these losses are outside the logic which veto's some collimator losses.  May want to
        !- move this into the preceeding loop.
        DO i=1, N_test_part
          number_of_test_particles = number_of_test_particles + 1

          IF(test_particles(i)%lost .eq. 1) THEN
            !- total losses
            Ndeposited_at_slice(test_particles(i)%slix_lost,0) = \
                 Ndeposited_at_slice(test_particles(i)%slix_lost,0) + \
                 delta_r*slice_len/c_light

            IF(test_particles(i)%plane_lost_at .eq. z_plane$) THEN
              bmad_zero_energy_count = bmad_zero_energy_count + 1

              !- particle stopped during deceleration
              Ndeposited_at_slice(test_particles(i)%slix_lost,2) = \
                   Ndeposited_at_slice(test_particles(i)%slix_lost,2) + \
                   delta_r*slice_len/c_light
            ELSE  !- must be x_plane$ or y_plane$
              beam_pipe_loss_count = beam_pipe_loss_count + 1

              !- beampipe collision
              Ndeposited_at_slice(test_particles(i)%slix_lost,1) = \
                   Ndeposited_at_slice(test_particles(i)%slix_lost,1) + \
                   delta_r*slice_len/c_light

              !- also record energy deposited into beampipe from energy of particle at time of collision
              energy_deposited_at_slice(test_particles(i)%slix_lost) = \
                  energy_deposited_at_slice(test_particles(i)%slix_lost) + \
                  (test_particles(i)%delta_eV + eV_at_slice(test_particles(i)%slix_lost))*delta_r * \
                  slice_len/c_light
            ENDIF
          ENDIF
        ENDDO
      ENDIF
    ENDDO
    IF(MOD(slix,10) .eq. 0) WRITE(*,*) "Slice ", slix, " of ", slix_prod_end, " complete."
  ENDDO

  !------------------------------------
  !- Send slaves terminate signal
  !------------------------------------
  JobDesc%count = 1
  JobDesc%job(1)%partnum = -1
  JobDesc%job(1)%slix = -1
  JobDesc%job(1)%deltam = -1.0
  DO i=1,nslave
    CALL mympi_send_job_desc(JobDesc,i)
  ENDDO

  !----------------------------------------------------------------------------------
  !- Compute the cumulative current generated, distinguished by pipe and stop losses
  !----------------------------------------------------------------------------------
  generation_at_slice(slix_prod_start)%pipe_cumulative = generation_at_slice(slix_prod_start)%pipe_this_slice
  generation_at_slice(slix_prod_start)%stop_cumulative = generation_at_slice(slix_prod_start)%stop_this_slice
  DO i=slix_prod_start+1, slix_prod_end
    generation_at_slice(i)%pipe_cumulative = generation_at_slice(i-1)%pipe_cumulative \
                                           + generation_at_slice(i)%pipe_this_slice
    generation_at_slice(i)%stop_cumulative = generation_at_slice(i-1)%stop_cumulative \
                                           + generation_at_slice(i)%stop_this_slice
  ENDDO

  !-----------------------------------------------------------------------------------
  !- Sum up total current of particles scattered, differentiating by how they are lost
  !-----------------------------------------------------------------------------------
  total_current = SUM(Ndeposited_at_slice(slix_lost_start:slix_lost_end,0))
  pipe_current = SUM(Ndeposited_at_slice(slix_lost_start:slix_lost_end,1))
  noenergy_current = SUM(Ndeposited_at_slice(slix_lost_start:slix_lost_end,2))
  total_energy_deposited = SUM(energy_deposited_at_slice(slix_lost_start:slix_lost_end))

  !-----------------------------------------------------------------------------------------------------------------------------
  ! Calculate the sigma matrix of the Touschek particles traveling through each element
  !-----------------------------------------------------------------------------------------------------------------------------
  OPEN(80, FILE="sigma_matrix.out")
  WRITE(80,'(A9,A11,21A16)') "# Slice  ","Location","Cov(x,x)","Cov(x,px)","Cov(x,y)","Cov(x,py)","Cov(x,z)","Cov(x,pz)", \
                                                   "Cov(px,px)","Cov(px,y)","Cov(px,py)","Cov(px,z)","Cov(px,pz)", \
                                                   "Cov(y,y)","Cov(y,py)","Cov(y,z)","Cov(y,pz)", \
                                                   "Cov(py,py)","Cov(py,z)","Cov(py,pz)", \
                                                   "Cov(z,z)","Cov(z,pz)", \
                                                   "Cov(pz,pz)"
  WRITE(80,'(A9,A11,21A16)') "# index  ","m","m^2","m","m^2","m","m^2","m","1","m","1","m","1","m^2","m","m^2","m","1","m","1","m^2","m","1"
  WRITE(80,'(A20,21A16)') "#    row number:    ","1","1","1","1","1","1","2","2","2","2","2","3","3","3","3","4","4","4","5","5","6"
  WRITE(80,'(A20,21A16)') "# column number:    ","1","2","3","4","5","6","2","3","4","5","6","3","4","5","6","4","5","6","5","6","6" 
  DO i=slix_prod_start,slix_lost_end
    IF(ne(i) .gt. 0.) THEN
      DO j=1,21
        sig_avg_dbl_terms(j) = sig_dbl_terms(i,j)/ne(i)
      ENDDO
      DO j=1,6
        sig_avg_sgl_terms(j) = sig_sgl_terms(i,j)/ne(i)
      ENDDO
      WRITE(80,'(I9,F11.3)',ADVANCE='no') i, slices(i)
      m=1
      DO j=1,6
        DO h=j,6
          term = sig_avg_dbl_terms(m)-sig_avg_sgl_terms(j)*sig_avg_sgl_terms(h)
          IF(ABS(term) .lt. 1.0Q-99) THEN
            WRITE(80,'(I16)',ADVANCE='no') 0
          ELSE
            WRITE(80,'(ES16.6)',ADVANCE='no') term
          ENDIF
          m=m+1
        ENDDO
      ENDDO
      WRITE(80,*) 
    ELSE
      WRITE(80,'(I9,F11.3,21I16)') i, slices(i), 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    ENDIF
  ENDDO
  CLOSE(80)

  !-------------------------------------------------------
  !- Write to stdout and touschek_track.details some information about how this simulation
  !- was run.
  !-------------------------------------------------------
  OPEN(80, FILE="touschek_track.details")
  DO output=6, 80, 74 !6 is standard output, 80 is the job_details.out file
    WRITE(output,*) "Lattice File: ", lat_file
!    WRITE(output,*) "Tracking: ", tracking
    WRITE(output,*) "Halo: ", halo
    IF(collimate) WRITE(output,*) "     ", collimators(1:num_cols)
    WRITE(output, '(A,I0)') "Number of data points: ", N_data_points
    WRITE(output, '(A,I0)') "Number of test particles per Touschek curve: ", N_test_part
    WRITE(output, '(A,F11.4)') "Total Integrated current (#e-): ", total_current
    WRITE(output, '(A,F11.4)') "  Beam Pipe Collisions: ", pipe_current
    WRITE(output, '(A,F11.4)') "  Zero Energy in Linac: ", noenergy_current
    WRITE(output, '(A,ES12.5)')  "Energy deposited into beam pipe (eV): ", total_energy_deposited

    WRITE(output,*) "Total number of test particles generated:      ", number_of_test_particles
    WRITE(output,*) "  Number lost to beam pipe:    ", beam_pipe_loss_count
    WRITE(output,*) "  Number lost to zero energy:  ", bmad_zero_energy_count
  ENDDO
  CLOSE(80)

  !---------------------------------------------------------------------------------------
  !- Write slice-by-slice data files showing where particles are lost and the amount
  !- of energy deposited into the beam pipe.
  !---------------------------------------------------------------------------------------
  OPEN(77, FILE="Npart_total_losses_by_s.out")
  accumulator = 0.
  WRITE(77,'(A61)') "#      slice   location     N deposited at       cumulative N"
  WRITE(77,'(A61)') "#      index        (s)    slice per bunch                   "
  DO i=slix_lost_start, slix_lost_end
    accumulator = accumulator + Ndeposited_at_slice(i,0)
    WRITE(77,'(I12,F11.3,2ES19.4)') i, slices(i), Ndeposited_at_slice(i,0), accumulator
  ENDDO
  CLOSE(77)
  OPEN(77, FILE="Npart_striking_pipe_by_s.out")
  accumulator = 0.
  WRITE(77,'(A61)') "#      slice   location     N deposited at       cumulative N"
  WRITE(77,'(A61)') "#      index        (s)    slice per bunch                   "
  DO i=slix_lost_start, slix_lost_end
    accumulator = accumulator + Ndeposited_at_slice(i,1)
    WRITE(77,'(I12,F11.3,2ES19.4)') i, slices(i), Ndeposited_at_slice(i,1), accumulator
  ENDDO
  CLOSE(77)
  OPEN(77, FILE="Npart_stopping_with_no_energy_by_s.out")
  accumulator = 0.
  WRITE(77,'(A61)') "#      slice   location     N deposited at       cumulative N"
  WRITE(77,'(A61)') "#      index        (s)    slice per bunch                   "
  DO i=slix_lost_start, slix_lost_end
    accumulator = accumulator + Ndeposited_at_slice(i,2)
    WRITE(77,'(I12,F11.3,2ES19.4)') i, slices(i), Ndeposited_at_slice(i,2), accumulator
  ENDDO
  CLOSE(77)
  OPEN(77, FILE="energy_deposited_at_s.out")
  accumulator = 0.
  WRITE(77,'(A61)') "#      slice   location   nrg deposited at         cumulative"
  WRITE(77,'(A61)') "#      index        (s)    slice per bunch                nrg"
  DO i=slix_lost_start, slix_lost_end
    accumulator = accumulator + energy_deposited_at_slice(i)
    WRITE(77,'(I12,F11.3,2ES19.4)') i, slices(i), energy_deposited_at_slice(i), accumulator
  ENDDO
  CLOSE(77)

  !----------------------------------------------------------------------------------------
  !- Write out the slice-by-slice histogram of the horizontal coordinate of scattered
  !- particles passing through each slice.
  !----------------------------------------------------------------------------------------
  IF(histogram_orbit) THEN
    OPEN(15,FILE="touschek_orbit_histogram_n.out")
      WRITE(15,'(A52)') "#     slice   location    midpoint of             N"
      WRITE(15,'(A52)') "#     index        (s)  histogram bin              "
      DO i=slix_prod_start, slix_lost_end
        DO j=1, hist_bins
          hist_mid = hist_min + bin_size*(j-0.5)
          WRITE(15,'(I10,F11.3,F15.5,ES14.4)') i, slices(i), hist_mid, orbit_hist(i,j,1)
        ENDDO
      ENDDO
    CLOSE(15)
    OPEN(16,FILE="touschek_orbit_histogram_energy.out")
      WRITE(16,'(A52)') "#     slice   location    midpoint of        energy"
      WRITE(16,'(A52)') "#     index        (s)  histogram bin              "
      DO i=slix_prod_start, slix_lost_end
        DO j=1, hist_bins
          hist_mid = hist_min + bin_size*(j-0.5)
          WRITE(16,'(I6,F11.3,F15.5,ES14.4)') i, slices(i), hist_mid, orbit_hist(i,j,2)
        ENDDO
      ENDDO
    CLOSE(16)
  ENDIF

  !------------------------------------------------------------------------
  !- Write results of potential collimator locations.
  !------------------------------------------------------------------------
  OPEN(17,FILE="collimator_profile.out")
    WRITE(17,'(A79)') "#     slice   location   potential N    element name                           "
    WRITE(17,'(A79)') "#     index        (s)        caught                                           "
    DO i=slix_prod_start, slix_lost_end
      eleix = element_at_s(lat,slices(i), .true.)
      WRITE(17,'(I11,F11.3,ES14.4,"    ",A40)') i, slices(i), col_profile(i), lat%ele(eleix)%name
    ENDDO
  CLOSE(17)

  !--------------------------------------------------------
  !- Close collimator log files
  !--------------------------------------------------------
  IF(coll_files) THEN
    DO i=1, number_of_collimators
      CLOSE(col_funits(i))
    ENDDO
  ENDIF

  !-----------------------------------------------------------------------
  !- Write generation profile files
  !-----------------------------------------------------------------------
  OPEN(18,FILE="generation_pipe.out")
    WRITE(18,'(A52)') "      slice   location    N generated   cumulative N"
    WRITE(18,'(A52)') "      index        (s)       at slice      generated"
    DO i=slix_prod_start, slix_prod_end
      WRITE(18,'(I11,F11.3,ES15.4,ES15.4)') i, slices(i), generation_at_slice(i)%pipe_this_slice, generation_at_slice(i)%pipe_cumulative
    ENDDO
  CLOSE(18)
  OPEN(19,FILE="generation_stop.out")
    WRITE(19,'(A52)') "      slice   location    N generated   cumulative N"
    WRITE(19,'(A52)') "      index        (s)       at slice      generated"
    DO i=slix_prod_start, slix_prod_end
      WRITE(19,'(I11,F11.3,ES15.4,ES15.4)') i, slices(i), generation_at_slice(i)%stop_this_slice, generation_at_slice(i)%stop_cumulative
    ENDDO
  CLOSE(19)
  OPEN(19,FILE="raw_rates.out")
    WRITE(19,'(A65)') "    slice    location     positive aperture     negative aperture"
    WRITE(19,'(A65)') "    index         (s)    rate (N/bunch/sec)    rate (N/bunch/sec)"
    DO i=slix_prod_start, slix_prod_end
      WRITE(19,'(I9,F12.3,ES22.4E3,ES22.4E3)') i, slices(i), raw_scattering_rates(i,1), raw_scattering_rates(i,2)
    ENDDO
  CLOSE(19)

  !-----------------------------------------------------------------
  !- Write header information for snapshot files, then append temporary snapshot files.
  !- Delete temporary files when done.
  !-----------------------------------------------------------------
  IF(traj_snapshot) THEN
    DO i=1,n_snapshots
      funit = 101+lat%n_ele_track+snapshot_start_slix+i-1
      WRITE(funit,'(A9)') "END_BUNCH"
      REWIND(funit)

      WRITE(snap_file_str,'(I5.5)') snapshot_start_slix+i-1
      OPEN(UNIT=79,FILE='snapshot_'//snap_file_str//'.dist')

      DO j=1,2
        READ(funit,'(A121)') snap_file_line
        WRITE(79,'(A121)') snap_file_line
      ENDDO
      READ(funit,'(A121)')
      WRITE(79,'(I8,A14)')    ntp(i), "  ! n_particle"
      DO j=1,1
        READ(funit,'(A121)') snap_file_line
        WRITE(79,'(A121)') snap_file_line
      ENDDO
      READ(funit,'(A121)')
      WRITE(79,'(ES14.4,A16)') ne(i)*e_charge,  "  ! bunch_charge"
      DO j=1,2+ntp(i)+1
        READ(funit,'(A121)') snap_file_line
        WRITE(79,'(A121)') snap_file_line
      ENDDO
      CLOSE(funit,STATUS='DELETE')
      CLOSE(79)
    ENDDO
  ENDIF
ELSE
  !slave
  DO WHILE(.true.)
    CALL mympi_check_mailbox_with_timeout(from_address, JobDesc%count) !Blocks until something arrives in mailbox.
                                                                       !from_address is id of slave the mail is from
    CALL mympi_recv_job_desc(JobDesc)

    IF( JobDesc%job(1)%slix .lt. 0 ) THEN
      !- negative slix is how master tells slave to terminate
      EXIT
    ELSE
      JobResults%count = JobDesc%count
      DO i=1, JobDesc%count
        cur_partnum = JobDesc%job(i)%partnum
        slix        = JobDesc%job(i)%slix
        cur_deltam  = JobDesc%job(i)%deltam

        !slorbit(slix) = initialized_coord_struct
        !slorbit(slix)%vec(:) = 0.0_rp
        !slorbit(slix)%vec(6) = cur_deltam

        init_vec(:) = 0.0_rp
        init_vec(6) = cur_deltam

        !------------------------------------------------------------------------
        !- Track from slice where test particle generated to slix_lost_end
        !------------------------------------------------------------------------
        CALL twiss_and_track_at_s(lat,slices(slix),stubele)   !get stubele
        CALL init_coord(slorbit(slix),init_vec,stubele, element_end = downstream_end$)
        if (lat%branch(0)%param%geometry .eq. open$) then
          CALL track_s_to_s(lat,slices,slix,slix_lost_end,slorbit,slix_lost,lost_to_col,plane_lost_at)
        else
          CALL track_s_to_s(lat,slices,slix,n_slices,slorbit,slix_lost,lost_to_col,plane_lost_at)
          if(slix_lost .ge. 0) then
            !particle was lost on first pass
            exit
          else
            !particle survived first past, track for n turns
            do j=1,n_turns
              slorbit(1) = slorbit(n_slices)
              call track_s_to_s(lat,slices,1,n_slices,slorbit,slix_lost,lost_to_col,plane_lost_at)
              if(slix_lost .ge. 0) then
                exit
              endif
            enddo
          endif
        endif

        JobResults%result(i)%partnum = cur_partnum
        IF(slix_lost .lt. 0) THEN
          JobResults%result(i)%lost = 0
        ELSE
          JobResults%result(i)%lost = 1
        ENDIF

        JobResults%result(i)%slixlost    = slix_lost
        JobResults%result(i)%planelostat = plane_lost_at
        JobResults%result(i)%lost_to_col = lost_to_col

        !----------------------------------------------------------------------------------------------
        !- Pack 6D trajectory for every test particle all into one long real array to facilitate
        !- passing by MPI.
        !----------------------------------------------------------------------------------------------
        part_offset = (i-1)*6*n_slices
        orbit_data_packed(part_offset+0*n_slices+1:part_offset+1*n_slices)=slorbit(1:n_slices)%vec(1)
        orbit_data_packed(part_offset+1*n_slices+1:part_offset+2*n_slices)=slorbit(1:n_slices)%vec(2)
        orbit_data_packed(part_offset+2*n_slices+1:part_offset+3*n_slices)=slorbit(1:n_slices)%vec(3)
        orbit_data_packed(part_offset+3*n_slices+1:part_offset+4*n_slices)=slorbit(1:n_slices)%vec(4)
        orbit_data_packed(part_offset+4*n_slices+1:part_offset+5*n_slices)=slorbit(1:n_slices)%vec(5)
        orbit_data_packed(part_offset+5*n_slices+1:part_offset+6*n_slices)=slorbit(1:n_slices)%vec(6)
      ENDDO

      !------------------------------------------------------------------------
      !- Now that all test particles received in job package have been tracked, send results back to master.
      !------------------------------------------------------------------------
      CALL mympi_send_results(JobResults,orbit_data_packed,n_slices)
    ENDIF
  ENDDO
ENDIF

!---------------------------------------------------------
!- Deallocate memory allocated by master process
!---------------------------------------------------------
IF(master) THEN
  DEALLOCATE(data_points)
  DEALLOCATE(test_particles)
  DEALLOCATE(tou_spline)
  DEALLOCATE(TrackResults)
  DO i=1,N_test_part
    DEALLOCATE(orbit_data(i)%orb)
  ENDDO
  DEALLOCATE(orbit_data)
  DO i=1,nslave
    DEALLOCATE(JobDescPackage(i)%job)
  ENDDO
  DEALLOCATE(JobDescPackage)
  DEALLOCATE(aperture_by_slice)
  DEALLOCATE(halo_aperture_by_slice)
  DEALLOCATE(eV_at_slice)
  DEALLOCATE(raw_scattering_rates)
  DEALLOCATE(Ndeposited_at_slice)
  DEALLOCATE(energy_deposited_at_slice)
  DEALLOCATE(generation_at_slice)
  DEALLOCATE(sig_dbl_terms)
  DEALLOCATE(sig_sgl_terms)
  DEALLOCATE(sig_avg_dbl_terms)
  DEALLOCATE(sig_avg_sgl_terms)
  DEALLOCATE(ne)
  IF(traj_snapshot) THEN
    DEALLOCATE(ntp)
  ENDIF
  DEALLOCATE(col_profile)
  IF(coll_files) THEN
    IF(number_of_collimators .gt. 0) THEN
      DEALLOCATE(col_funits)
      DEALLOCATE(col_hash)
    ENDIF
  ENDIF
  IF(histogram_orbit) THEN
    DEALLOCATE(orbit_hist)
  ENDIF
ENDIF


!------------------------------------------------------
!- Deallocate memory allocated by slaves
!------------------------------------------------------
IF(.not.master) THEN
  DEALLOCATE(slorbit)
ENDIF

!--------------------------------------------------------
!- Deallocate memory allocated by both master and slaves
!--------------------------------------------------------
DEALLOCATE(slices)
DEALLOCATE(JobDesc%job)        ! contains descriptions of several test particles for slave to track
DEALLOCATE(JobResults%result)  ! contains results of slave tracking all the test particles that were sent to it
DEALLOCATE(orbit_data_packed)  ! contains the 6-dim trajectory data for all test particles a slave tracked.

!----------------------------------------------------------------
!- Tells MPI daemon that this process will be terminating soon.
!----------------------------------------------------------------
CALL mympi_shutdown()

END PROGRAM touschek_background

