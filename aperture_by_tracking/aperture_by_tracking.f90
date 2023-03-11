!+
! Program aperture_by_tracking
!
! Use: aperture_by_tracking <lat_file>
!
! This program produces the location-by-location momentum aperture of a linear or circular
! accelerator.  The locations can be separated by a fixed distance, say, for example, 1 meter,
! or can be set to coincide with element ends.  A separate positive and negative
! momentum aperture is found at each location.  The aperture is found using a binary search
! for the largest in magnitude delta_p that can be applied to the closed orbit without
! the particle being lost.
!
! If the accelerator is linear, then the lattice is tracked through only once.  If the
! accelerator is circular, then the lattice is tracked through a number of times, set by
! the .in file.
!
! Input file format:
!     lat_file           --  Location of lattice file
!     slice_method       --  'byele' or 'byloc'.  If 'byloc', then slice_length must be set
!     slice length       --  length of element slices if slice_method is 'byloc'
!     target_accuracy    --  accuracy with which to determine the aperture 
!                            e.g.  0.001 would determine the aperture to 0.1% of its actual value
!
!     The following parameters are interesting only for circular accelerators
!     nturns             --  Number of turns to track particles

!     The following parametes interesting only for linear accelerators
!     start_stop_type    --  1 or 2.  1 to specify start and stop using start_s and end_s, and track_till
!                                     2 to specify element names using start_name and end_name, and track_end_name
!     start_s            --  first location (in meters) from start of lattice determine momentum aperture
!     end_s              --  last location (in meters) from start of lattice to determine momentum aperture
!     track_till         --  last location in meters in linear accelerator to track particles through
!     start_name         --  name of first element to determine momentum aperture.  first appearance of name is used.
!     end_name           --  name of last element to determine momentum aperture.  first appearance of name is used.
!     track_end_name     --  last element in linear accelerator to track particles through
!-
PROGRAM aperture_by_tracking

USE bmad
USE aperture_by_tracking_mod
USE aperture_by_tracking_mpi_mod
USE slice_mod

IMPLICIT NONE

!------------------------------Control variables
INTEGER i, k
INTEGER slix, eleix
INTEGER return_address
INTEGER templun
LOGICAL err
LOGICAL master
LOGICAL lost
INTEGER track_state
character(40), parameter :: r_name = 'aperture_by_tracking'

!------------------------------Variables Related to Binary Search
REAL(rp) delta_m
REAL(rp) accuracy
integer :: vec_ix   
real(rp) :: angle1, angle2
real(rp), allocatable :: angle(:)
integer :: n_apertures, n_angles
logical :: angle_scan

!------------------------------Hardcoded Parameters
REAL(rp), PARAMETER :: first_guess = 0.01

!------------------------------Parameters from .in
CHARACTER*100 parameter_file
CHARACTER lat_file*130
CHARACTER(20) start_name, end_name, track_end_name
CHARACTER slice_method*6
INTEGER n_loc, start_stop_type
INTEGER n_slices
INTEGER nturns
LOGICAL ring
LOGICAL halo
REAL(rp) halo_h_emittance
REAL(rp) halo_aperture
REAL(rp) target_accuracy
REAL(rp) track_till
REAL(rp) start_s, end_s, slice_length
REAL(rp) :: Qx = -1.0, Qy = -1.0, Qz = -999.0
character(5) :: aperture_type

!------------------------------Data
REAL(rp), ALLOCATABLE :: slices(:)
REAL(rp), ALLOCATABLE :: max_deltam(:,:)
REAL(rp), allocatable :: aperture(:)

!------------------------------Derived Data Types
TYPE(lat_struct) lat
TYPE(ApertureJob_struct) ApertureJob
TYPE(ApertureResult_struct) ApertureResult
TYPE(ele_pointer_struct), ALLOCATABLE :: eles(:)
TYPE(coord_struct), ALLOCATABLE :: co(:)
TYPE(coord_struct), ALLOCATABLE :: sco(:)
TYPE(coord_struct) vec_start
TYPE(coord_struct) vec_start0
TYPE(coord_struct) vec_end

!------------------------------Settings to be read in from .in
NAMELIST /parameters/ &
    lat_file, &         ! Latice file including possible path. 
    slice_method, &     ! 'bystep' or 'byelem': slice at element ends or take fixed length steps
    slice_length, &     ! If 'bystep', this is step size
    target_accuracy,&   ! Fractional accuracy.  eg. 0.001 finds aperture to within 0.1%
    start_stop_type, &  ! If 1, then use start_name, end_name, and track_end_name.
                        !     If 2, then use start_s, end_s, and track_till.
    start_name, &       ! Start generating momentum aperture here (first instance of start_name)
    end_name, &         ! Stop generating momentum aperture here (first instance of end_name)
    track_end_name, &   ! Stop tracking particles here (first instance of track_end_name)
    start_s, &          ! Start generating momentum aperture here (s location in meters)
    end_s, &            ! Stop generating momentum aperture here (s location in meters)
    track_till,&        ! Stop tracking particles here (s location in meters)
    nturns, &           ! For storage rings only: number of turns to track
    Qx, &               ! For storage rings only: qtune lattice
    Qy, &               ! For storage rings only: qtune lattice
    Qz, &               ! For storage rings only: qtune lattice
    halo, &             ! For linacs only: in addition to flagging particles as lost due to beam pipe collisions,
                        !     also count them as lost if they lay outside n-sigma of the beam phase-space
                        !     ellipse at the end of the linac, where n = halo_aperture.  The results from
                        !     halo mode can be extreme if the lattice ends with a dump.  In that case,
                        !     set track_end_name or track_till to some location before the dump.
    halo_aperture, &    ! For linacs only: see halo comment
    halo_h_emittance, & ! For linacs only: horizontal emittance to use when calculating beam phase space ellipse 
                        !     at end of lattice
    aperture_type, &    ! Momentum variable to scan, one of { 'px', 'py', 'pz' (default), 'angle'}
                        !     The output file with be prefixed with these characters
    angle1, angle2, &   ! Starting and ending angles, if aperture_type == 'angle'
    n_angles            ! Number of angles to scan

!------------------------------Read in parameter file
halo = .false. !default
aperture_type = 'pz'
angle1 = 0
angle2 = twopi
n_angles = 2

CALL GETARG(1,parameter_file)
templun = LUNGET()
OPEN(templun,FILE=parameter_file,STATUS='old')
READ(templun,NML=parameters)
CLOSE(templun)

!--- Set vec_ix and angle scan parameters
aperture_type = downcase(aperture_type)
angle_scan = .false.

select case (aperture_type)
  case('px')
  	vec_ix = 2
  	n_apertures = 2
  case('py')
  	vec_ix = 4  	
  	n_apertures = 2
  case('pz')
  	vec_ix = 6
  	n_apertures = 2
  case('angle')
    vec_ix = -1
    n_apertures = n_angles
    angle_scan = .true.
    if (n_apertures > aperture_by_tracking_mpi_mod_MAX_APERTURES) then
      print *, 'Number of angles > aperture_by_tracking_mpi_mod_MAX_APERTURES: ', aperture_by_tracking_mpi_mod_MAX_APERTURES
      stop
    endif
  case default
  	print *, 'Error: Bad aperture_type: '//aperture_type
  	stop
end select

allocate(aperture(n_apertures))
allocate(angle(n_apertures))

if (angle_scan) then
  ! Populate angle array
  if (n_apertures > 1) then
    do i=1, n_apertures
      angle(i) = angle1 + (i-1)*(angle2 - angle1)/(n_apertures -1)
    enddo
  else
    ! Just use first angle
     angle(1) = angle1  
  endif
endif


!------------------------------Initialize MPI and parse lattice
CALL mympi_initialize(master)
CALL mympi_bmad_parser(lat_file,lat)
bmad_com%aperture_limit_on = .true.
  
  
if (angle_scan .and. master) write (*, '(a, i0, a, f14.5, a, f14.5)') 'Scanning ', n_apertures, ' angles from ', angle1, ' to ', angle2

!------------------------------Detect lattice type
IF(lat%param%geometry == closed$) THEN
  ring = .true.
ELSE
  ring = .false.
ENDIF

!------------------------------Calculate twiss parameters and by-element closed orbit
IF(ring) THEN
  CALL prep_lat_ring(lat,co,Qx,Qy,Qz,master)
ELSE
  CALL twiss_and_track(lat,co)
ENDIF

!------------------------------Override certain parameters if lattice is a ring
IF(ring) THEN
  !If the accelerator is a ring, then override the start and stop parameters to ensure
  !the entire ring is tracked through.
  start_stop_type = 1
  start_s = 0.0_rp
  end_s = -1.0
  track_till = -1.0
ENDIF

!--------------------------------Determine start and stop locations from .in parameters
IF(start_stop_type .eq. 1) THEN
  !Start and stop locations in parameter file are given by location in meters
  IF(end_s .lt. 0.) end_s = lat%param%total_length
  IF(track_till .lt. 0) THEN
    track_till = lat%param%total_length
  ENDIF
ELSEIF(start_stop_type .eq. 2) THEN
  !Start and stop locations in parameter file are given by element name
  CALL lat_ele_locator(start_name,lat,eles,n_loc,err)
  start_s = lat%ele(eles(1)%ele%ix_ele-1)%s
   
  if (n_loc == 0) call out_io (s_fatal$, r_name, 'Start element not found: '//trim(start_name) )

  CALL lat_ele_locator(end_name,lat,eles,n_loc,err)
  if (n_loc == 0) call out_io (s_fatal$, r_name, 'End element not found: '//trim(end_name) )
  end_s = lat%ele(eles(1)%ele%ix_ele)%s

  CALL lat_ele_locator(track_end_name,lat,eles,n_loc,err)
  if (n_loc == 0) call out_io (s_fatal$, r_name, 'Track end element not found: '//trim(track_end_name) )
  track_till = eles(1)%ele%s
  
ENDIF

!--------------------------------Set up slices
IF(master) THEN
  CALL make_slices(lat,slice_method,end_s,start_s,slice_length,slices,n_slices)
  ALLOCATE(max_deltam(1:n_slices,1:n_apertures))
ENDIF

!-----------------------------Obtain by-slice closed orbit
IF(MASTER) THEN
  ALLOCATE(sco(1:n_slices))
  sco(1)%vec = co(0)%vec
  CALL track_s_to_s(lat,slices,1,n_slices,sco)
ENDIF

!-----------------------------Log some data to files
IF(master) THEN
  templun=lunget()
  OPEN(templun,FILE="physical_aperture.out")
  WRITE(templun,'(A71)') "#       ele     location     aperture   element   element name         "
  WRITE(templun,'(A71)') "#     index          (s)       radius       key                        "
  DO i=1, lat%n_ele_track
    IF( (lat%ele(i)%value(l$) .gt. 0.0) .or. \
        (lat%ele(i)%key .eq. ecollimator$) .or. \
        (lat%ele(i)%key .eq. rcollimator$) ) THEN
      WRITE(templun,'(I11,F13.5,F13.5,"  ",A,"  ",A)') i, lat%ele(i)%s, lat%ele(i)%value(x1_limit$), \
                                                     key_name(lat%ele(i)%key), lat%ele(i)%name
    ENDIF
  ENDDO
  CLOSE(templun)

  templun=lunget()
  OPEN(templun,FILE="slice_index.out")
  WRITE(templun,'(A43)') "#     slice      location   element        "
  WRITE(templun,'(A43)') "#     index           (s)   name           "
  DO i=1,n_slices
    eleix = element_at_s(lat,slices(i),.true.)
    WRITE(templun,'(I11,F14.5,"   ",A)') i, slices(i), lat%ele(eleix)%name
  ENDDO
  CLOSE(templun)
ENDIF

!-----------------------------Preparation is complete.  Begin computing momentum apertures.
IF(master) THEN
  WRITE(*,'(A)') "Preparation complete.  Beginning processing..."

  ! This first loop seeds each slave with a tracking job.
  slix=1
  DO i=1, min(nslave,n_slices)
    ApertureJob%slix = i
    ApertureJob%s = slices(i)
    ApertureJob%sco(1:6) = sco(i)%vec(1:6)
    CALL mympi_send_job(i, ApertureJob)
  ENDDO
  slix = i 

  ! This second loop receives each completed job as it comes in, and reseeds the 
  ! slave if there are more tracking jobs left to do.
  DO i=1, n_slices
    CALL mympi_check_mailbox_with_timeout(return_address)        !Blocks until something arrives in mailbox.
                                                                 !return_address is id of slave the mail is from
    CALL mympi_receive_result(return_address, ApertureResult)    !Gets the mail, stores it in ApertureResult

    !If there are work units left to complete, the following IF block sends them to the slave at return_address
    !If there are no work units left, send terminate signal to the slave.
    IF(slix .le. n_slices) THEN
      !send new job
      ApertureJob%slix = slix
      ApertureJob%s = slices(slix)
      ApertureJob%sco(1:6) = sco(slix)%vec(1:6)
      CALL mympi_send_job(return_address, ApertureJob)
      slix = slix + 1
    ELSE
      !send terminate signal
      ApertureJob%slix = -1 
      CALL mympi_send_job(return_address, ApertureJob)
    ENDIF

    max_deltam(ApertureResult%slix,1:n_apertures) = ApertureResult%aperture(1:n_apertures)

    CALL progress_indicator(i,n_slices) !Prints percent done
  ENDDO
ELSE !slave
  vec_start0 = vec_start
  DO WHILE(.true.)
    CALL mympi_check_mailbox_with_timeout(return_address)  !Blocks until something arrives in mailbox.
    CALL mympi_receive_job(ApertureJob)

    IF(ApertureJob%slix .lt. 0) THEN
      ! slix < 0 is how the master signals to the slaves that the program is to terminate.
      EXIT
    ELSE
      !k=1 determines the positive aperture, k=2 determines the negative aperture
      DO k= 1, n_apertures
        !The following loop conducts a binary search for the largest energy kick that can be given to a particle
        !such that it is lost.
        !This loop first brackets the largest energy kick, then uses a binary search to obtain the
        !desired precision.

        !Initialization
        CALL binary_search(.false.,0._rp,accuracy,.true.) !calling binary_search with last arg .true. resets its state variables
        delta_m = first_guess  

        DO WHILE (accuracy .gt. target_accuracy)
          vec_start = vec_start0 
          vec_start%vec = ApertureJob%sco                 ! start from the closed orbit
          
          if (angle_scan) then
            vec_start%vec(2) = vec_start%vec(2) + cos(angle(k))*delta_m  ! Kick at an angle
            vec_start%vec(4) = vec_start%vec(4) + sin(angle(k))*delta_m
          else
            ! Simple positive or negative kick
            ! if k=1, then pos aperture.  if k=2, then neg aperture.
            vec_start%vec(vec_ix) = vec_start%vec(vec_ix) + ((-1)**(k+1))*delta_m   ! add the momentum kick
          endif
          
          IF(ring) THEN
            CALL check_if_lost_ring(lat,ApertureJob%s,vec_start,vec_end,nturns, track_state)
          ELSE
            CALL check_if_lost_linac(lat,ApertureJob%s,vec_start,track_state,track_till,halo,halo_aperture,halo_h_emittance)
          ENDIF

          lost = (track_state .ne. moving_forward$)

          CALL binary_search(lost, delta_m, accuracy, .false.) !updates delta_m according to one iteration of binary search
        ENDDO

        
          aperture(k) = delta_m
          ! Include sign when not doing angle scan
          if (.not. angle_scan) aperture(k) = ((-1)**(k+1))*aperture(k) 
        
      ENDDO

      !Package the results and send them back to the master.
      ApertureResult%slix = ApertureJob%slix
      ApertureResult%aperture(1:n_apertures) = aperture(1:n_apertures)
      CALL mympi_send_result(ApertureResult)
    ENDIF
  ENDDO
ENDIF

!--------------------------------------Write results to file
IF(master) THEN
  templun = lunget()
  IF(halo) THEN
    OPEN(templun, FILE=trim(aperture_type)//"_halo_aperture_by_s.out")
  ELSE
    OPEN(templun, FILE=trim(aperture_type)//"_aperture_by_s.out")
  ENDIF
  !This file contains the element-by-element momentum aperture
  
  if (angle_scan) then
    write(templun,'(a15)', advance = 'NO') 's(m),angle:'
    do k = 1, n_apertures -1
      write(templun,'(es15.7)', advance = 'NO')  angle(k)
    enddo
    write(templun,'(es15.7)', advance = 'YES')  angle(n_apertures) 
  endif
  
  DO i=1,n_slices
    if (angle_scan) then
      ! Write a long line
      write(templun,'(es15.7)', advance = 'NO') slices(i)
      do k = 1, n_apertures -1
        write(templun,'(es15.7)', advance = 'NO')  max_deltam(i,k)
      enddo
      write(templun,'(es15.7)', advance = 'YES')  max_deltam(i,n_apertures)
    else
      WRITE(templun,'(3es15.7)') slices(i), max_deltam(i,2), max_deltam(i,1)
    endif
  ENDDO
  CLOSE(templun)
  WRITE(*,*) "Done!"
ENDIF

!--------------------------------------Shutdown program
CALL mympi_shutdown()  !Disengage from MPI daemon
IF(master) THEN
  DEALLOCATE(slices)
  DEALLOCATE(max_deltam)
  DEALLOCATE(sco)
ENDIF

END PROGRAM aperture_by_tracking




