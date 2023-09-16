!+
! Module touschek_aperture_mpi_mod
!
! This module consolidates all of the native MPI calls and code in the touschek_aperture
! program.
!-
MODULE touschek_background_mpi_mod

USE precision_def
USE bmad !needed for err_exit
USE sim_utils !needed for milli_sleep
use mpi

IMPLICIT NONE

INTEGER status(MPI_STATUS_SIZE)           ! stores return codes from MPI daemon
INTEGER mpistatus(MPI_STATUS_SIZE)        ! also stores MPI daemon codes, used when overwriting status is not practical
INTEGER myrank    ! holds rank of process.  rank 0 is master, all others are slaves
INTEGER mpierr

!--------------------------------
!- TrackJob_struct is used to send job information to the slaves.  It contains:
!- partnum: ID number of the test particle
!- slix: slice index where tracking starts
!- deltam: momentum offset of particle
!--------------------------------
TYPE TrackJob_struct
  SEQUENCE
  INTEGER partnum
  INTEGER slix
  REAL(rp) deltam
END TYPE TrackJob_struct
INTEGER TrackJob_type

!---------------------------------
!- JobDesc_struct is a container for multiple TrackJob_structs.  It contains:
!- count: number of jobs
!- job(:): instances of TrackJob_struct
!---------------------------------
TYPE JobDesc_struct
  SEQUENCE
  INTEGER count
  INTEGER padding
  TYPE(TrackJob_struct), ALLOCATABLE :: job(:)
END TYPE JobDesc_struct

!---------------------------------
!- Used by slave to send results of tracking jobs back to the master.  Contains:
!- partnum: ID of test particle
!- lost: 0 if particle not lost. 1 if lost
!- slixlost: slice index where particle was lost
!- planelostat: plane where particle was lost.  x_plane$, y_plant$, or z_plane$
!- lost_to_col: set to true if particle was lost to a collimator.
!---------------------------------
TYPE TrackResult_struct
  SEQUENCE
  INTEGER partnum
  INTEGER lost
  INTEGER slixlost
  INTEGER planelostat
  INTEGER lost_to_col
END TYPE TrackResult_struct
INTEGER TrackResult_type

!---------------------------------
!- Container for multiple TrackResult_structs.  Contains:
!- count: number of TrackResults
!- result(:): instances of TrackResult_struct.  One for each tracking job.
!---------------------------------
TYPE JobResults_struct
  SEQUENCE
  INTEGER count
  INTEGER padding
  TYPE(TrackResult_struct), ALLOCATABLE :: result(:)
END TYPE JobResults_struct

!Data
PRIVATE status
PRIVATE mpistatus
PRIVATE myrank
PRIVATE mpierr

!Data structures
PUBLIC TrackJob_struct
PUBLIC JobDesc_struct
PUBLIC TrackResult_struct
PUBLIC JobResults_struct

!Functions
PUBLIC mympi_initialize
PRIVATE mympi_register_derived_types

CONTAINS

!+
! Subroutine mympi_initialize(master)
!
! First mpi call made by a program.  This subroutine calls mpi_init, which announces
! the executible to the MPI daemon.  mpi_comm_rank asks the MPI daemon to tell the executible
! what its rank number is.  One node is given rank 0 and becomes the master.  All other nodes
! receive numbers counting up from 1 and become workers.
! mympi_register_derived_types tells the MPI daemon about the structures we will be sending.
! nslave is returned by this program because the main program needs to know the number of slaves
! in order to decide how to divvy out jobs.
!
! Input:
!   None
! Output:
!   master:  LOGICAL, INTENT(OUT): Set to true if master, false otherwise.
!   nslave:  INTEGER, INTENT(OUT): number of slaves in the cluster.
!-
SUBROUTINE mympi_initialize(master,nslave)
  LOGICAL, INTENT(OUT) :: master
  INTEGER, INTENT(OUT) :: nslave

  INTEGER i
  INTEGER cluster_size

  CALL mpi_init(mpierr)                             ! Introduce yourself to the MPI daemon
  CALL mpi_comm_rank(MPI_COMM_WORLD,myrank,mpierr)  ! Get your rank number, store in myrank.  Master is rank 0.
  CALL mympi_register_derived_types()               ! Tell the daemon about derived types that will be
                                                    ! shared among the cluster
  DO i=1,MPI_STATUS_SIZE
    status(i) = 1
    mpistatus(i) = 1
  ENDDO

  !Check that cluster has at least two nodes
  CALL mpi_comm_size(MPI_COMM_WORLD,cluster_size,mpierr)
  nslave=cluster_size-1
  IF(myrank .eq. 0) THEN
    IF(nslave .eq. 0) THEN
      WRITE(*,*) "ERROR: no slaves found in cluster.  At least two nodes"
      WRITE(*,*) "must be available to run this program."
      STOP
    ENDIF
  ENDIF

  !Return master=TRUE if I'm master, else false
  IF(myrank .eq. 0) THEN
    master=.true.
  ELSE
    master=.false.
  ENDIF

END SUBROUTINE mympi_initialize

!+
! Subroutine mympi_register_derived_types()
!
! Tells the MPI daemon about the derived data types that will be passed
! between master and slaves.
! Sets module data TrackJob_type and TrackResult_type which are IDs used when calling
! mpi_send and mpi_recv.
!
! Input:
!   None
! Output:
!   None
!-
SUBROUTINE mympi_register_derived_types()
  INTEGER mpierr
  INTEGER oldtypes(0:1),blockcounts(0:1)
  integer(mpi_address_kind) lb, extent, offsets(0:1)

  !Setup mpi type for TrackJob_struct
  offsets(0)=0
  oldtypes(0)=MPI_INTEGER
  blockcounts(0)=2
  CALL mpi_type_get_extent(MPI_INTEGER,lb,extent,mpierr)
  offsets(1)=2*extent
  oldtypes(1)=MPI_DOUBLE_PRECISION
  blockcounts(1)=1
  CALL mpi_type_create_struct(2,blockcounts,offsets,oldtypes,TrackJob_type,mpierr)
  CALL mpi_type_commit(TrackJob_type,mpierr)

  !Setup mpi type for TrackResult_struct
  offsets(0)=0
  oldtypes(0)=MPI_INTEGER
  blockcounts(0)=5
  CALL mpi_type_create_struct(1,blockcounts,offsets,oldtypes,TrackResult_type,mpierr)
  CALL mpi_type_commit(TrackResult_type,mpierr)
END SUBROUTINE mympi_register_derived_types

!+
! Subroutine mympi_bmad_parser(lat_file,lat)
!
! Front-end to the bmad_parser subroutine.  If master, then bmad_parser is called immediately.
! If slave, then wait until master call to bmad_parser is complete.
!
! Having the slaves wait ensures that the digested file will be available when they call bmad_parser.
! If the digested file were not available, then the slaves would simultaneously begin recreating
! the same digested file, which would result in file corruption.
!
! Input:
!   lat_file    -- CHARACTER(130), INTENT(IN): filename of lattice
! Output:
!   lat         -- TYPE(lat_struct): parsed lattice
!-
SUBROUTINE mympi_bmad_parser(lat_file, lat)
  USE bmad
  CHARACTER(*), INTENT(IN) :: lat_file
  TYPE(lat_struct) lat
  INTEGER mpierr

  ! Master calls bmad parser first, which will write the digested file if it does not exist.
  ! After completing bmad parser, it reaches the barrier, signaling the slaves to call bmad parser.
  IF(myrank .eq. 0) THEN
    CALL fullfilename(lat_file, lat_file)
    CALL bmad_parser(lat_file, lat)
    CALL mpi_barrier(MPI_COMM_WORLD,mpierr)  ! Each process waits at barrier till all processes have reached the barrier
  ELSE
    CALL mpi_barrier(MPI_COMM_WORLD,mpierr)  ! Each process waits at barrier till all processes have reached the barrier
    CALL fullfilename(lat_file, lat_file)
    CALL bmad_parser(lat_file, lat)
  ENDIF
END SUBROUTINE mympi_bmad_parser

!+
! Subroutine mympi_bcast_sigmap(lat)
!
! Broadcasts lat%param%sigmap from the master to the slaves.
!
! Input:
!     lat%ele(:)%z%sigma_p: If called by master, this is input
! Output:
!     lat%ele(:)%z%sigma_p: If called by slave, this is output
!-
SUBROUTINE mympi_bcast_sigmap(lat)
  USE bmad

  TYPE(lat_struct) lat
  REAL(rp), ALLOCATABLE :: temp_array(:)
  INTEGER mpierr

  ALLOCATE(temp_array(1:lat%n_ele_track+1))

  IF(myrank .eq. 0) THEN
    temp_array(1:lat%n_ele_track+1) = lat%ele(0:lat%n_ele_track)%z%sigma_p
  ENDIF
  CALL mpi_bcast(temp_array, lat%n_ele_track+1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
  IF(myrank .ne. 0) THEN
    lat%ele(0:lat%n_ele_track)%z%sigma_p = temp_array(1:lat%n_ele_track+1)
  ENDIF
  DEALLOCATE(temp_array)
END SUBROUTINE mympi_bcast_sigmap

!+
! Subroutine mympi_check_mailbox_with_timeout(return_address)
!
! This subroutine blocks until a message arrives in the mailbox.  When something arrives,
! it returns with return_address set to the ID of the node the message is from.  This subroutine
! does not open the message, it simply says that something arrived and who it is from.
!
! wait_ms = 100 means check the mailbox every 0.1 seconds.
! max_wait_ix = 6000 means check the mailbox 6000 times before timing out.
!               Timing out terminates the program.
!-
SUBROUTINE mympi_check_mailbox_with_timeout(return_address, package_size)
  INTEGER, INTENT(OUT) :: return_address
  INTEGER, INTENT(OUT) :: package_size
  INTEGER wait_ms, wait_ix, max_wait_ix
  INTEGER mpierr
  LOGICAL flag !set to true by mpi_iprobe when there is something in the mailbox

  wait_ms=10
  wait_ix=0
  max_wait_ix=3000
  DO WHILE (.true.)
    CALL mpi_iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,flag,status,mpierr)  !check the mailbox
    IF(flag) THEN   !Did somethine arrive?
      package_size = status(MPI_TAG)
      EXIT ! Good! Stop checking mailbox.
    ELSE
      CALL milli_sleep(wait_ms)  !wait for wait_ms milliseconds
      wait_ix = wait_ix + 1
      IF(wait_ix .gt. max_wait_ix) THEN
        WRITE(*,*) "Notice: Communication Taking Longer Than Expected..."
        wait_ix = 0
!        CALL mpi_finalize(mpierr)  !terminate the program
!        CALL err_exit
      ENDIF
    ENDIF
  ENDDO

  return_address = status(MPI_SOURCE)  !return address is which node the message is from.
END SUBROUTINE mympi_check_mailbox_with_timeout

!+
! Subroutine mympi_send_job_desc(JobDesc,slave_id)
!
! Called by master to send JobDesc, which contains multiple tracking jobs, to slave slave_id.
! Counterpart to mympi_recv_job_desc.
!
! Input:
!   JobDesc%count   : INTEGER: number of jobs
!   JobDesc%job(:)  : TYPE(TrackJob_struct): job descriptions
!   slave_id        : ID of slave to send JobDesc to
!-
SUBROUTINE mympi_send_job_desc(JobDesc,slave_id)
  TYPE(JobDesc_struct), INTENT(IN) :: JobDesc
  INTEGER, INTENT(IN) :: slave_id

  CALL mpi_send(JobDesc%job(1:JobDesc%count), JobDesc%count, TrackJob_type, &
                  slave_id, JobDesc%count, MPI_COMM_WORLD, mpierr)
END SUBROUTINE mympi_send_job_desc

!+
! Subroutine mympi_recv_job_desc(JobDesc)
!
! Called by slave to receive JobDesc from master.
! Counterpart to mympi_send_job_desc.
!
! Input:
!   None
! Output:
!   JobDesc:   Type(JobDesc_struct), INTENT(INOUT)
!-
SUBROUTINE mympi_recv_job_desc(JobDesc)
  TYPE(JobDesc_struct), INTENT(INOUT) :: JobDesc
  INTEGER package_size

  package_size = JobDesc%count

  CALL mpi_recv(JobDesc%job(1:package_size), package_size, TrackJob_type, &
                  0, MPI_ANY_TAG, MPI_COMM_WORLD, mpistatus, mpierr)
END SUBROUTINE mympi_recv_job_desc

!+
! Subroutine mympi_receive_results(from_address, n_slices, JobResults, orbit_data_package)
!
! Called by master to receive results of tracking jobs from slaves.
! Counterpart to mympi_send_results
! 
! Input:
!   from_address  :  INTEGER, INTENT(IN): ID of slave from which to receive the results
!   n_slices      :  INTEGER, INTENT(IN): number of slices in lattice.  needed for packing orbit data
!   JobResults%count : INTEGER, INTENT(IN): number of results to receive.
! Output:
!   JobResults%result(:)    :  TYPE(TrackResult_struct), INTENT(OUT): results of tracking job
!   orbit_data_packed(:)    :  REAL(rp), INTENT(OUT): 6D trajectory data for all particles
!-
SUBROUTINE mympi_receive_results(from_address, n_slices, JobResults, orbit_data_packed)
  INTEGER, INTENT(IN) :: from_address
  INTEGER, INTENT(IN) :: n_slices
  TYPE(JobResults_struct), INTENT(INOUT) :: JobResults
  REAL(rp), INTENT(INOUT) :: orbit_data_packed(:)
  INTEGER package_size

  package_size = JobResults%count

  CALL mpi_recv(JobResults%result(1:package_size), package_size, TrackResult_type, from_address, MPI_ANY_TAG, MPI_COMM_WORLD, mpistatus, mpierr)

  CALL mpi_recv(orbit_data_packed(1:(package_size*6*n_slices)), package_size*6*n_slices, \
                MPI_DOUBLE_PRECISION, from_address, MPI_ANY_TAG, MPI_COMM_WORLD, mpistatus, mpierr)

END SUBROUTINE mympi_receive_results

!+
! Subroutine mympi_send_results(from_address, n_slices, JobResults, orbit_data_package)
!
!  Called by slave to send results of tracking jobs to master.
!  Counterpart to mympi_receive_results
!
!  Input:
!    JobResults%count   : INTEGER, INTENT(IN): number of results
!    orbit_data_packed(:)  : REAL(rp), INTENT(IN): 6D trajectory data for all jobs packed into one array
!    n_slices : INTEGER, INTENT(IN): number of slices in lattice.  needed for packing orbit data.
!
!-
SUBROUTINE mympi_send_results(JobResults,orbit_data_packed,n_slices)
  TYPE(JobResults_struct), INTENT(IN) :: JobResults
  REAL(rp), INTENT(IN) :: orbit_data_packed(:)
  INTEGER, INTENT(IN) :: n_slices

  INTEGER package_size

  package_size = JobResults%count

  CALL mpi_send(JobResults%result(1:package_size),package_size,TrackResult_type,0,package_size,MPI_COMM_WORLD,mpierr)
  CALL mpi_send(orbit_data_packed(1:package_size*6*n_slices),package_size*6*n_slices,MPI_DOUBLE_PRECISION,0,package_size,MPI_COMM_WORLD,mpierr)
END SUBROUTINE mympi_send_results

!+
! Subroutine mympi_shutdown()
!
! This routine prepares for program termination by calling mpi_finalize, which disconnects
! the process from the MPI daemon.
!
! Input:
!   None
! Output:
!   None
!-
SUBROUTINE mympi_shutdown()
  INTEGER mpierr
  CALL mpi_finalize(mpierr)
END SUBROUTINE mympi_shutdown

END MODULE touschek_background_mpi_mod






