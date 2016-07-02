!+
! Module tune_tracker_mod
!
! This module simulates a digital tune tracker.  In short, it consists of a PLL
! that locks onto bpm data.  The PLL output is used to modulate a kicker that fires
! on every turn to excite the beam at resonance.  The result is a kicker modulator
! that excites the beam at its natural fractional tune.
!
! This module is coded as an object and supports multiple instances.  Each instance
! is initiated by calling a constructor which returns an ID that identifies the instance.
! An outline of how to call the tune tracker is:
!   id = init_dTT(tt_param_struct)
!   DO
!     <track a turn, store in orb>
!     bpm_meas = orb(somewhere)%vec(1)
!     sinphi = TT_update(bpm_meas, id)
!     kick = sinphi * kickamplitude
!     <set element kick attribute>
!   ENDDO
!   CALL dest_dTT(id)
!
! init_dTT is the constructor and is passed the tune tracker parameters in tt_param_struct
! The contents of tt_param_struct are below
! TT_update is passed one BPM measurement per call.  
! It returns the kick amplitude (normalized to one) to be used on the next turn.
! dest_dTT is the destructor.  It is passed the id of the tune tracker to destroy.
!-

module tune_tracker_mod

use bmad
use windowls_mod

implicit none

integer, parameter :: max_tt = 3 !number of tune trackers to allocate memory for,
                                 !effectively maximum number of tune trackers

! Structure to hold parameters of tune tracker.
TYPE tt_param_struct
  LOGICAL :: useSaveState = .false. ! If true then save TT state vars to file on exit,
                                    ! and load TT state vars on initialization.
  REAL(rp) phi_to_kicker     ! phase advance from x phase at bpm to xdot phase at kicker
  REAL(rp) Dt                ! time between data points (usually set to ring period)
  REAL(rp) LPinertia         ! inertia of lowpass filter that follows mixer. Try 2^15
  REAL(rp) Ki                ! Integrator gain
  REAL(rp) Kp                ! Proportional gain
  REAL(rp) Kvco              ! VCO gain
  REAL(rp) modw0             ! VCO base frequency (initial guess of beam frac tune)
  REAL(rp) offset            ! Closed orbit at bpm.
  CHARACTER(3) mixmode       ! To mix BPM data with square wave or sine wave
  INTEGER cyc_per_turn       ! Number of times tune tracker cycles per turn.
                             ! Usually an integer fraction of the harmonic number.
                             ! Try 183
  INTEGER :: log_period = 1  ! How often to write to tt_log files.  Zero means do not log, one
                             ! means log every call to TT_update, three means every third call, etc.

  !Parameters specific to D channel
  !Note:  The D channel is not compatible with save states.
  LOGICAL :: use_D_chan = .FALSE.      ! Whether to use the D channel.
                                       ! Disable D channel if not in use
  REAL(rp) Kd             ! Differential gain
  INTEGER wls_N           ! Number of data poins for window LS
  INTEGER wls_order       ! Order of fit polynomial

  ! The following parameters are overwritten by the init_dTT
  REAL(rp) LPalpha        ! low-pass filter constant. Calculated from LPinertia
  INTEGER wls_id          ! Instance ID for window LS module for D channel
  LOGICAL :: allocated = .FALSE.   ! Keeps track of whic tt_ids have been instantiated by user.

  ! The following are parameters that are likely to be needed by any program
  ! that makes use of the tune tracker module.
  ! This module does not actually use these parameters.
  ! For examples, see the tune tracker driver program.
  INTEGER :: bpm_loc = -1        ! Location of BPM
  INTEGER :: kck_loc = -1        ! Location of kicker
  CHARACTER(40) :: bpm_name = ''    ! Name of BPM
  CHARACTER(40) :: kck_name = ''    ! Name of Kicker
  REAL(rp) modTfrac0      ! initial fractional tune of kicker modulator
  REAL(rp) kickAmplitude  ! amplitude of kicks
  CHARACTER orientation   ! 'h', 'v', or 'l'
  INTEGER Onum            ! Element of coord_struct used by bpm.
                          ! 1 for horiz, 3 for vert, 1 for longitudinal
END TYPE tt_param_struct

! Structures to hold the state variables of the tune tracker.
TYPE tt_state_struct
  REAL(rp) deltaw               ! w_vco = w0 + deltaw
  REAL(rp) intDphi              ! integrated Dphi
  REAL(rp) bpm_msmt_last        ! stashes previous bpm measurement
  REAL(rp) psi                  ! keeps track of modulator angle
  REAL(rp) Dphi                 ! Delta Phi - angle that comes out of filters mixer signal
  REAL(rp) gain                 ! bpm gain
  INTEGER  counter              ! counts number of times TT_update called.  Needed for log file.
END TYPE tt_state_struct

! Variables related to multiple TT instances
INTEGER, PRIVATE, SAVE :: tts_instantiated = 0  !number of tune trackers instantiated
TYPE(tt_param_struct), PRIVATE, SAVE :: tt_param(max_tt)
TYPE(tt_state_struct), PRIVATE, SAVE :: tt_state(max_tt)
INTEGER, PRIVATE, SAVE :: log_luns(max_tt)

PUBLIC init_dTT
PUBLIC TT_update
PUBLIC dest_dTT
PUBLIC get_dTT
PRIVATE modulator
PRIVATE check_id

CONTAINS

!+
! Function id = init_dTT(incoming_tt_param)
!
! Constructor for TT_update module.  This function creates a new instance of the tune tracker by assigning
! an id and copying the tune tracker parameters to the module saved data.  It also initializes the tune tracker's
! state variables.
!
! This function opens one log file "tt_log_<id>.out", which needs to be closed by calling the destructor.
!
! Modules needed:
!   use tune_tracker_mod
!
! Input:
!   incoming_tt_param    -- TYPE(tt_param_struct): Tune tracker parameters.  See structure definition for details.
!
! Output:
!   <return value>       -- INTEGER: id of the tune tracker instance created.
!-
FUNCTION init_dTT(incoming_tt_param,saved_coords) RESULT(id)
  TYPE(tt_param_struct) :: incoming_tt_param
  TYPE(coord_struct), OPTIONAL :: saved_coords
  INTEGER i
  INTEGER id
  INTEGER ios
  CHARACTER id_str
  CHARACTER(12) tt_log_name
  CHARACTER(10) state_name
  LOGICAL state_read

  IF(tts_instantiated + 1 .gt. max_tt) THEN
    WRITE(*,*) "Maximum number of tune trackers hard-coded to max_tt = ", max_tt
    STOP
  ENDIF

  tts_instantiated = tts_instantiated + 1

  DO i=1,max_tt
    IF( .not. tt_param(i)%allocated ) THEN
      id = i 
      tt_param(id) = incoming_tt_param
      tt_param(i)%allocated = .true.
      EXIT
    ENDIF
  ENDDO

  WRITE(id_str,'(I1)') id

  ! Initialize state variables
  state_read = .false.
  IF(tt_param(id)%useSaveState) THEN
    IF(PRESENT(saved_coords)) THEN
      state_name = "tt_state."//id_str
      OPEN(UNIT=9999,FILE=state_name,STATUS='OLD',ACTION='READ',IOSTAT=ios)
      IF(ios == 0) THEN
        READ(9999,*) saved_coords, tt_state(id)
        state_read = .true.
        CLOSE(9999)
        WRITE(*,'(A,I2,A)') "NOTICE: Save state found for TT #", id, ":"
        WRITE(*,'(A,F10.7)') "        VCO frequency set to ", (tt_state(id)%deltaw+tt_param(id)%modw0)* &
                                                    tt_param(id)%Dt / 2.0_rp / pi
      ELSE
        WRITE(*,*) "NOTICE: Tune tracker warning for TT #", id, ":"
        WRITE(*,*) "        Unable to open file: ", state_name
        WRITE(*,*) "        Setting initial TT state variables to defaults."
        WRITE(*,*) "        Leaving initial coordinates unchanged."
      ENDIF
    ENDIF
  ENDIF
  IF(.not. state_read) THEN
    tt_state(id)%deltaw = 0.0_rp
    tt_state(id)%intDphi = 0.0_rp
    tt_state(id)%bpm_msmt_last = 0.0_rp
    tt_state(id)%psi = 0.0_rp
    tt_state(id)%Dphi = 0.0_rp
    tt_state(id)%gain = 0.001_rp
    tt_state(id)%counter = 0
  ENDIF

  !calculate low-pass filter parameter
  tt_param(id)%LPalpha = 1.0_rp / ( 1.0_rp + (tt_param(id)%LPinertia/2.0_rp/pi) )

  tt_log_name = "tt_log_"//id_str//".out"
  log_luns(id) = lunget()
  OPEN(UNIT=log_luns(id),FILE=tt_log_name,STATUS='REPLACE')

  IF(tt_param(id)%use_D_chan) THEN
    IF( .not. tt_param(id)%useSaveState ) THEN
      tt_param(id)%wls_id = initFixedWindowLS(tt_param(id)%wls_N,tt_param(id)%Dt,tt_param(id)%wls_order,1)
    ELSE
      WRITE(*,*) "WARNING: D channel calculations are not compatible with save states."
      WRITE(*,*) "         D channel will be disabled."
      tt_param(id)%use_D_chan = .FALSE.
    ENDIF
  ENDIF
END FUNCTION init_dTT

!+
! Subroutine reset_dTT(id,tt_param)
!
! Reset a tune tracker and optionally change its parameters.
!
! Modules Needed:
!   use tune_tracker_mod
!
! Input:
!   id        -- INTEGER, INTENT(IN): Tune tracker instance ID obtained from constructor.
!   tt_param  -- 
! Output:
!   none
!-
SUBROUTINE reset_dTT(id,tt_param_update)
  INTEGER, INTENT(IN) :: id
  TYPE(tt_param_struct), OPTIONAL :: tt_param_update

  CALL check_id(id)

  ! Reset state variables
  tt_state(id)%deltaw = 0.0_rp
  tt_state(id)%intDphi = 0.0_rp
  tt_state(id)%bpm_msmt_last = 0.0_rp
  tt_state(id)%psi = 0.0_rp
  tt_state(id)%Dphi = 0.0_rp
  tt_state(id)%gain = 0.001_rp
  tt_state(id)%counter = 0

  ! Update tt_param if tt_param_update is present.  Only certain fields can be updated.
  IF(PRESENT(tt_param_update)) THEN
    tt_param(id)%phi_to_kicker = tt_param_update%phi_to_kicker
    tt_param(id)%Dt = tt_param_update%Dt
    tt_param(id)%LPinertia = tt_param_update%LPinertia
    tt_param(id)%Ki = tt_param_update%Ki
    tt_param(id)%Kp = tt_param_update%Kp
    tt_param(id)%Kvco = tt_param_update%Kvco
    tt_param(id)%modw0 = tt_param_update%modw0
    tt_param(id)%offset = tt_param_update%offset
    tt_param(id)%mixmode = tt_param_update%mixmode
    tt_param(id)%cyc_per_turn = tt_param_update%cyc_per_turn

    tt_param(id)%LPalpha = 1.0_rp / ( 1.0_rp + (tt_param(id)%LPinertia/2.0_rp/pi) )
  ENDIF

  ! Write two blank lines to the log file to indicate that a reset has occurred
  WRITE(log_luns(id),*)
  WRITE(log_luns(id),*)

END SUBROUTINE reset_dTT

!+
! Subroutine dest_dTT(id)
!
! Destructor for tune tracker.  Should be called end of program to ensure all files are properly closed.
!
! Modules Needed:
!   use tune_tracker_mod
!
! Input:
!   id     -- INTEGER, INTENT(IN): Tune tracker instance ID obtained from constructor.
! Output:
!   none
!-
SUBROUTINE dest_dTT(id,coords)
  INTEGER, INTENT(IN) :: id
  CHARACTER id_str
  CHARACTER(10) state_name
  TYPE(coord_struct), OPTIONAL :: coords

  CALL check_id(id)

  tts_instantiated = tts_instantiated - 1

  tt_param(id)%allocated = .false.

  IF(tt_param(id)%useSaveState) THEN
    IF(PRESENT(coords)) THEN
      WRITE(id_str,'(I1)') id
      state_name = "tt_state."//id_str
      OPEN(UNIT=9999,FILE=state_name,STATUS='REPLACE',ACTION='WRITE')
      WRITE(9999,*) coords, tt_state(id)
      CLOSE(9999)
    ENDIF
  ENDIF

  CLOSE(log_luns(id))
END SUBROUTINE dest_dTT

!+
! Function z = TT_update(bpm_msmt,id)
!
! Main function of the tune tracker module.  This funcion is given one new data point each turn,
! and it returns the phase of the kicker to be used on the next turn.  The incoming data is mixed with
! the modulator signal and passed to a low-pass filter which returns the phase difference between the
! bpm data and the modulator.  The phase difference is passed to a PID controller which adjusts the VCO frequency.
!
! This function returns the phi of the VCO plus an offset set during initialization.  The offset is the phase
! advance from the bpm to the kicker
!
! Modules needed:
!   use tune_tracker_mod
!
! Input:
!   bpm_msmt     -- REAL(rp): new bpm measurement. passed by value.
!   id           -- INTEGER, INTENT(IN): instance id of tune tracker
! Output:
!   z            -- Double precision: sin(modulator_angle + phase_advance_to_kicker)
!-
FUNCTION TT_update(bpm_msmt,id) RESULT(z)
  IMPLICIT NONE

  REAL(rp), PARAMETER :: ga = 0.05_rp     !bpm gain time constant

  REAL(rp) :: bpm_msmt
  REAL(rp) bpm_msmt_nrml
  INTEGER, INTENT(IN) :: id
  REAL(rp) :: z

  REAL(rp) AB  ! mixed signal
  REAL(rp) sinout, sqrout  ! modulator output
  REAL(rp) proDphi  !proportional to Dphi
  REAL(rp) dirDphi  !first derivative of Dphi
  REAL(rp) PIDout ! P + I + D
  REAL(rp) bpm
  REAL(rp) fastPeriod

  INTEGER i

  CALL check_id(id)

  tt_state(id)%counter = tt_state(id)%counter + 1    !needed for log file

  ! Adjust bpm data for closed orbit offset and apply gain
  bpm_msmt_nrml = bpm_msmt - tt_param(id)%offset
  ! Update gain.  Gain adjusted to keep bpm average at unity.
  tt_state(id)%gain = ga*ABS(bpm_msmt_nrml) + (1.0_rp-ga)*tt_state(id)%gain
  ! Apply gain.
  bpm_msmt_nrml = bpm_msmt_nrml/tt_state(id)%gain

  fastPeriod = tt_param(id)%Dt / tt_param(id)%cyc_per_turn

  !This loop represents the tune tracker fast elecronics (typically 14ns)
  DO i=1, tt_param(id)%cyc_per_turn
    IF( i <= FLOOR(tt_param(id)%cyc_per_turn / 2.0) ) THEN
      bpm = tt_state(id)%bpm_msmt_last
    ELSE
      bpm = bpm_msmt_nrml
    ENDIF
    CALL modulator(tt_state(id)%psi,sinout,sqrout)

    ! Phase Detector: mixer followed by low-pass filter
    IF( tt_param(id)%mixmode == 'sqr' ) THEN
      AB = bpm*sqrout
    ELSEIF( tt_param(id)%mixmode == 'sin' ) THEN
      AB = bpm*sinout
    ELSE
      WRITE(*,*) "Unknown mixmode specified.  Must be 'sin' or 'sqr'.  Terminating."
      STOP
    ENDIF
    tt_state(id)%Dphi = AB*tt_param(id)%LPalpha + (1.0_rp-tt_param(id)%LPalpha)*tt_state(id)%Dphi

    ! Calculate Integral Channel
    tt_state(id)%intDphi = tt_state(id)%intDphi + &
                           (tt_param(id)%Ki/tt_param(id)%cyc_per_turn)*tt_state(id)%Dphi

    ! Update modulator angle
    tt_state(id)%psi = tt_state(id)%psi + (tt_param(id)%modw0 + tt_state(id)%deltaw)*fastPeriod
    tt_state(id)%psi = MOD(tt_state(id)%psi,2.0_rp*pi)

  ENDDO
  tt_state(id)%bpm_msmt_last = bpm_msmt_nrml

  !The following calculates the proportional and derivative channels
  proDphi = tt_param(id)%Kp * tt_state(id)%Dphi
  IF( tt_param(id)%use_D_chan ) THEN
    dirDphi = fixedWindowLS(proDphi,tt_param(id)%wls_id)
    PIDout = tt_state(id)%intDphi + proDphi - dirDphi*tt_param(id)%Kd
  ELSE
    PIDout = tt_state(id)%intDphi + proDphi
  ENDIF

  !Update modulator frequency
  tt_state(id)%deltaw = tt_param(id)%Kvco * PIDout

  !Write to log file
  IF( tt_param(id)%log_period .gt. 0 ) THEN
    IF( MOD( (tt_state(id)%counter+tt_param(id)%log_period-1), tt_param(id)%log_period) == 0 ) THEN
      IF( tt_param(id)%use_D_chan ) THEN
        !This statement writes the state of each PID channel to the tt_log_n.out file.
        WRITE(log_luns(id),'(I8,3ES14.4)') tt_state(id)%counter, tt_param(id)%Ki*tt_state(id)%intDphi, &
                                                                 tt_param(id)%Kp*proDphi, &
                                                                 -tt_param(id)%Kd*dirDphi
      ELSE
        WRITE(log_luns(id),'(I8,2ES14.4)') tt_state(id)%counter, tt_param(id)%Ki*tt_state(id)%intDphi, &
                                                                 tt_param(id)%Kp*proDphi
      ENDIF
    ENDIF
  ENDIF

  CALL modulator(tt_state(id)%psi + tt_param(id)%phi_to_kicker, sinout, sqrout)
  z = sinout ! z is result, the value returned by this function

END FUNCTION TT_update

!+
! Function z = get_dTT(name,id)
!
! A get function to obtain the tune tracker VCO frequency.
!
! Modules needed:
!   use tune_tracker_mod
!
! Input:
!   name     -- CHARACTER(2), INTENT(IN):  'dw' to get VCO trim, 'wf' to get VCO frequency
!   id       -- INTEGER, INTENT(IN): Instance id
! Output:
!   <return value>  -- REAL(rp): value of requested data
!-
FUNCTION get_dTT(name,id) RESULT(z)
  CHARACTER(2), INTENT(IN) :: name
  INTEGER, INTENT(IN) :: id
  REAL(rp) z

  CALL check_id(id)

  IF(name == 'dw') THEN
    ! return VCO trim 
    z = tt_state(id)%deltaw
  ELSEIF(name == 'wf') THEN
    ! return VCO frequency
    z = tt_state(id)%deltaw + tt_param(id)%modw0
  ELSE
    z = -999.
  ENDIF
END FUNCTION get_dTT

!+
! Subroutine modulator(psi,sinout,sqrout)
!
! Tune tracker modulator.  Returns sin(psi) and corresponding square wave.
!
! Modules needed:
!   use tune_tracker_mod
!
! Input:
!   psi     -- REAL(rp), INTENT(IN):  modulator angle
! Output:
!   sinout  -- REAL(rp), INTENT(OUT):  sin(psi)
!   sqrout  -- REAL(rp), INTENT(OUT):  +1 if sin(psi)>=0, -1 else
!-
SUBROUTINE modulator(psi,sinout,sqrout)
  REAL(rp), INTENT(IN) :: psi
  REAL(rp), INTENT(OUT) :: sinout
  REAL(rp), INTENT(OUT) :: sqrout

  sinout = sin(psi)
  
  IF(sinout >= 0._rp) THEN
    sqrout = 1.0_rp
  ELSE
    sqrout = -1.0_rp
  ENDIF

END SUBROUTINE modulator

SUBROUTINE check_id(id)
  INTEGER, INTENT(IN) :: id
  
  IF( (id .LT. 1) .OR. id .GT. max_tt) THEN
    !Absurd TT id passed to TT module
    WRITE(*,'(A,I3)') "FATAL IN CALL TO TUNE TRACKER: ID received is ", id
    WRITE(*,'(A,I3)') "                  ID should be between 1 and ", max_tt
    WRITE(*,'(A)')    "                  Check that init_dTT has been called."
    WRITE(*,'(A)')    "                  Check that ID is an integer."
    CALL err_exit
  ENDIF

  IF( .not. tt_param(id)%allocated ) THEN
    !Tune tracker not initialized
    WRITE(*,'(A,I3)') "FATAL IN CALL TO TUNE TRACKER: ID received is ", id
    WRITE(*,'(A)')    "                A tune tracker with this ID has not been instantiated."
    WRITE(*,'(A)')    "                Check that init_dTT has been called and dest_dTT has not been called."
    WRITE(*,'(A)')    "                Check that ID is an integer."
    CALL err_exit
  ENDIF
END SUBROUTINE check_id

END MODULE tune_tracker_mod
