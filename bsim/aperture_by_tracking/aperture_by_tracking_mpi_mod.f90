!+
! Module touschek_aperture_mpi_mod
!
! This module consolidates all of the native MPI calls and code in the aperture_by_tracking
! program.  
!-
MODULE aperture_by_tracking_mpi_mod

USE precision_def
USE bmad !needed for err_exit
USE sim_utils !needed for milli_sleep
use mpi

IMPLICIT NONE

integer, parameter :: aperture_by_tracking_mpi_mod_MAX_APERTURES = 360

!+
! TYPE ApertureJob_struct
!
! This structure is used for sending job instructions from the master to the workers.
!  slix - tells the worker the index number of the slice it is working.  This field
!         effectively is a label for the job.
!  s    - tells the worker the s location (in meters) at which to determine the
!         the momentum aperture.  i.e.  %s = 10.123 tells the worker to determine
!         the momentum aperture 10.123 meters from the start of the accelerator.
!  sco(1:6) - The closed orbit at the location s.
!- 
TYPE ApertureJob_struct
  SEQUENCE !Tells FORTRAN to store this data in memory in order.  Needed so that the MPI daemon
           !can ship the structre around.
  INTEGER ipad !for alignment
  INTEGER slix
  REAL(rp) s
  REAL(rp) sco(1:6)
END TYPE ApertureJob_struct
INTEGER ApertureJob_mpi_id

!+
! TYPE ApertureResult_struct
!
! This structure is used for sending job results from a worker back to the master.
!  slix - tells the master the index number of the slice it was working on.  This field
!         is effectively a label fo the job.
!  p_aperture - The positive momentum aperture at the s location associated with index slix.
!  n_aperture - The negative momentum aperture at the s location associated with index slix.
!-
TYPE ApertureResult_struct
  SEQUENCE !Tells FORTRAN to store this data in memory in order.  Needed so that the MPI daemon
           !can ship the structre around.
  INTEGER ipad !for alignment
  INTEGER slix
  REAL(rp) aperture(aperture_by_tracking_mpi_mod_MAX_APERTURES)   ! Hard coded maximum. 
END TYPE ApertureResult_struct
INTEGER ApertureResult_mpi_id

!status and mpistatus are used to send and receive codes from the MPI daemon.
INTEGER status(MPI_STATUS_SIZE),mpistatus(MPI_STATUS_SIZE)

!myrank is populated by the mympi_initialize and contains the rank number.
!nslave contains the number of workers in the cluster.  Typically #nodes-1
INTEGER myrank, nslave

!Subroutines
PUBLIC mympi_initialize         !how the executible announces itself to the MPI daemon
PUBLIC mympi_bmad_parser        !front end to bmad parser
PUBLIC mympi_check_mailbox_with_timeout ! block with timeout waiting for something to arrive in mailbox
PUBLIC mympi_receive_result             ! called by master to receive job result from a worker
PUBLIC mympi_send_result
PUBLIC mympi_shutdown
PRIVATE mympi_register_derived_types

!Data
PUBLIC ApertureJob_struct
PUBLIC ApertureResult_struct
PRIVATE myrank
PRIVATE status, mpistatus

CONTAINS

!+
! Subroutine mympi_initialize(master)
!
! First mpi call made by a program.  This subroutine calls mpi_init, which announces
! the executible to the MPI daemon.  mpi_comm_rank asks the MPI daemon to tell the executible
! what its rank number is.  One node is given rank 0 and becomes the master.  All other nodes
! receive numbers counting up from 1 and become workers.
! mympi_register_derived_types tells the MPI daemon about the structures we will be sending.
!
! Input:
!   None
! Output:
!   master:  LOGICAL, INTENT(OUT): Set to true if master, false otherwise.
!-
SUBROUTINE mympi_initialize(master)
  LOGICAL, INTENT(OUT) :: master
  INTEGER i
  INTEGER mpierr
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
  IF(myrank .eq. 0) THEN
    CALL mpi_comm_size(MPI_COMM_WORLD,cluster_size,mpierr)
    nslave=cluster_size-1
    IF(nslave .eq. 0) THEN
      WRITE(*,*) "ERROR: this code must be run with mpi with at least 2 nodes."
      WRITE(*,*) "try the following: mpirun -n 2 aperture_by_tracking aperture_by_tracking.in"
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
  CHARACTER(130), INTENT(IN) :: lat_file
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
! Subroutine mympi_register_derived_types()
!
! Because this program uses structures, rather than native data types, for communication,
! the MPI daemon needs to be told about form of the structure.
! Later when mpi_send and mpi_recv are called, the integer parameters ApertureJob_mpi_id
! and ApertureResult_mpi_id are used to tell MPI which data type is being sent.
!
! Input:
!   None
! Output:
!   None
!-
SUBROUTINE mympi_register_derived_types()
  INTEGER mpierr
  INTEGER oldtypes(0:1),blockcounts(0:1)
  integer(mpi_address_kind) lb, extent,offsets(0:1)

  !Setup mpi type for ApertureJob_struct
  offsets(0)=0
  oldtypes(0)=MPI_INTEGER
  blockcounts(0)=2
  CALL mpi_type_get_extent(MPI_INTEGER,lb,extent,mpierr)
  offsets(1)=2*extent
  oldtypes(1)=MPI_DOUBLE_PRECISION
  blockcounts(1)=7
  CALL mpi_type_create_struct(2,blockcounts,offsets,oldtypes,ApertureJob_mpi_id,mpierr)
  CALL mpi_type_commit(ApertureJob_mpi_id,mpierr)

  !Setup mpi type for ApertureResult_struct
  offsets(0)=0
  oldtypes(0)=MPI_INTEGER
  blockcounts(0)=2
  CALL mpi_type_get_extent(MPI_INTEGER,lb,extent,mpierr)
  offsets(1)=2*extent
  oldtypes(1)=MPI_DOUBLE_PRECISION
  blockcounts(1)=360
  CALL mpi_type_create_struct(2,blockcounts,offsets,oldtypes,ApertureResult_mpi_id,mpierr)
  CALL mpi_type_commit(ApertureResult_mpi_id,mpierr)
END SUBROUTINE mympi_register_derived_types

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
SUBROUTINE mympi_check_mailbox_with_timeout(return_address)
  INTEGER, INTENT(OUT) :: return_address
  INTEGER wait_ms, wait_ix, max_wait_ix
  INTEGER mpierr
  LOGICAL flag !set to true by mpi_iprobe when there is something in the mailbox

  wait_ms=10
  wait_ix=0
  max_wait_ix=60000
  DO WHILE (.true.)
    CALL mpi_iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,flag,status,mpierr)  !check the mailbox
    IF(flag) THEN   !Did somethine arrive?
      EXIT ! Good! Stop checking mailbox.
    ELSE
      CALL milli_sleep(wait_ms)  !wait for wait_ms milliseconds
      wait_ix = wait_ix + 1
      IF(wait_ix .gt. max_wait_ix) THEN
        wait_ix = 0
        WRITE(*,*) "Notice: Communication Taking Longer Than Expected..."
!        CALL mpi_finalize(mpierr)  !terminate the program
!        CALL err_exit
      ENDIF
    ENDIF
  ENDDO

  return_address = status(MPI_SOURCE)  !return address is which node the message is from.
END SUBROUTINE mympi_check_mailbox_with_timeout

!+
! Subroutine mympi_receive_result(address, ApertureResult)
!
! Receive momentum aperture calculation results from a worker.
! This routine is complementary to mympi_send_result.  The worker calls mympi_send_result,
! while the master calls mympi_receive_result.
!
! This simple subroutine is a front-end for the native mpi call mpi_recv.
! Its purpose is to remove from the main program arguments which are irrelevant or obscure.
!
! Input:
!   address        -- INTEGER, INTENT(IN): address of slave from whom to receive the result.
! Output:
!   ApertureResult -- TYPE(ApertureResult_struct), INTENT(OUT): job results received from slave.
!-
SUBROUTINE mympi_receive_result(address, ApertureResult)
  INTEGER, INTENT(IN) :: address
  TYPE(ApertureResult_struct), INTENT(OUT) :: ApertureResult
  INTEGER mpierr

  CALL mpi_recv(ApertureResult,1,ApertureResult_mpi_id,address,1,MPI_COMM_WORLD,mpistatus,mpierr)
END SUBROUTINE mympi_receive_result

!+
! Subroutine mympi_send_result(ApertureResult)
!
! Send job result to a master.
! This routine is complementary to mympi_receive_result.  The worker calls mympi_send_result,
! while the master calls mympi_receive_result.
!
! This simple subroutine is a front-end for the native mpi call mpi_send.
! Its purpose is to remove from the main program arguments which are irrelevant or obscure.
!
! Note that unlike mympi_receive_result, this subroutine takes no address.  This is because
! the only node the workers ever talk to is the master, which always has address zero.
!
! Input:
!   ApertureResult -- TYPE(ApertureResult_struct), INTENT(OUT): job results received from slave.
! Output:
!   None
!-
SUBROUTINE mympi_send_result(ApertureResult)
  TYPE(ApertureResult_struct), INTENT(IN) :: ApertureResult
  INTEGER mpierr

  CALL mpi_send(ApertureResult,1,ApertureResult_mpi_id,0,1,MPI_COMM_WORLD,mpierr)
END SUBROUTINE mympi_send_result

!+
! Subroutine mympi_receive_job(ApertureJob)
!
! Receive job specification from the master.
! This routine is complementary to mympi_send_job.  The worker calls mympi_receive_job,
! while the master calls mympi_send_job.
!
! This simple subroutine is a front-end for the native mpi call mpi_recv.
! Its purpose is to remove from the main program arguments which are irrelevant or obscure.

! Note that unlike mympi_send_job, this subroutine takes no address.  This is because
! the only node the workers ever talk to is the master, which always has address zero.
!
! Input:
!   None
! Output:
!   ApertureJob -- TYPE(ApertureJob_struct), INTENT(OUT): job specification from master
!-
SUBROUTINE mympi_receive_job(ApertureJob)
  TYPE(ApertureJob_struct), INTENT(OUT) :: ApertureJob
  INTEGER mpierr

  CALL mpi_recv(ApertureJob,1,ApertureJob_mpi_id,0,1,MPI_COMM_WORLD,mpistatus,mpierr)
END SUBROUTINE mympi_receive_job

!+
! Subroutine mympi_send_job(address,ApertureJob)
!
! Send job specification to a slave.
! This routine is complementary to mympi_receive_job.  The worker calls mympi_receive_job,
! while the master calls mympi_send_job.
!
! This simple subroutine is a front-end for the native mpi call mpi_send.
! Its purpose is to remove from the main program arguments which are irrelevant or obscure.

! Input:
!   address     -- INTEGER, INTENT(IN) :: address of worker to send job to
!   ApertureJob -- TYPE(ApertureJob_struct), INTENT(IN): job specification to send to worker
! Output:
!   None
!-
SUBROUTINE mympi_send_job(address, ApertureJob)
  INTEGER, INTENT(IN) :: address
  TYPE(ApertureJob_struct) ApertureJob
  INTEGER mpierr

  CALL mpi_send(ApertureJob,1,ApertureJob_mpi_id,address,1,MPI_COMM_WORLD,mpierr)
END SUBROUTINE mympi_send_job

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

END MODULE aperture_by_tracking_mpi_mod






