!+
! Module tune_tracker_mod
!
! This module simulates the CESR Digital tune tracker.  In short, it consists of a PLL
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
!     sinphi = cesr_dTT(bpm_meas, id)
!     kick = sinphi * kickamplitude
!     <set element kick attribute>
!   ENDDO
!   CALL dest_dTT(id)
!
! init_dTT is the constructor and is passed the tune tracker parameters in tt_param_struct
! The contents of tt_param_struct are below
! cesr_dTT is passed one BPM measurement per call.  It returns the kick amplitude (normalized to one) to be used on the next turn.
! dest_dTT is the destructor.  It is passed the id of the tune tracker to destroy.
!-
MODULE tune_tracker_mod

USE bmad
USE windowLS_mod

IMPLICIT NONE

INTEGER, PARAMETER :: max_tt = 3 !number of tune trackers to allocate memory for,
                                 !effectively maximum number of tune trackers

! Structure to hold parameters of tune tracker.
TYPE tt_param_struct
  REAL(rp) phi_to_kicker  ! phase advance from x phase at bpm to xdot phase at kicker
  REAL(rp) Dt             ! time between data points (usually set to ring period)
  REAL(rp) fastPeriod     ! period of PLL electronics.  e.g. 14 ns for CESR 
  REAL(rp) LPalpha        ! low-pass filter constant
  REAL(rp) alpha          ! Twiss alpha.  Needed by elliptical_oscillator
  REAL(rp) beta           ! Twiss alpha.  Needed by elliptical_oscillator
  REAL(rp) Ki             ! Integrator gain
  REAL(rp) Kp             ! Proportional gain
  REAL(rp) Kvco           ! Redundant, because total gain looks like Kvco*(Ki*x + Kp*y)
  REAL(rp) modw0          ! Modulator base frequency (initial guess of beam frac tune)
  REAL(rp) offset         ! Closed orbit at bpm.
  CHARACTER(3) mixmode    ! To mix BPM data with square wave or sine wave

  !Parameters specific to D channel
  LOGICAL :: use_D_chan = .FALSE.      ! Whether to use the D channel.  Disable D channel if not in use
  REAL(rp) Kd             ! Differential gain
  INTEGER wls_id          ! Instance ID for window LS module for differential channel
  INTEGER wls_N           ! Number of data poins for window LS
  INTEGER wls_order       ! Order of fit polynomial

  ! The following entries are to facilitate namelist input in whatever program is used to drive this module.
  ! This module does not use these parameters.
  INTEGER bpm_loc         ! Location of BPM
  INTEGER kck_loc         ! Location of kicker
  REAL(rp) modTfrac0      ! initial fractional tune of kicker modulator
  REAL(rp) LPinertia      ! inertial of lowpass filter that follows mixer
  REAL(rp) kickAmplitude  ! amplitude of kicks
  CHARACTER orientation   ! 'h', 'v', or 'l'
  INTEGER Onum            ! 1 for horiz, 3 for vert, 5 for longitudinal
END TYPE tt_param_struct

! Structures to hold the state variables of the tune tracker.
TYPE tt_state_struct
  REAL(rp) deltaw               ! w_vco = w0 + deltaw
  REAL(rp) intDphi              ! integrated Dphi
  REAL(rp) bpm_msmt_last        ! stashes previous bpm measurement
  REAL(rp) t                    ! time
  REAL(rp) psi                  ! modulator angle
  REAL(rp) Dphi                 ! Delta Phi - angle that comes out of filters mixer signal
  REAL(rp) gain                 ! bpm gain
END TYPE tt_state_struct

! Variables related to multiple TT instances
INTEGER, PRIVATE, SAVE :: tt_ids = 0  !number of tune trackers instantiated
TYPE(tt_param_struct), PRIVATE, SAVE :: tt_param(max_tt)
TYPE(tt_state_struct), PRIVATE, SAVE :: tt_state(max_tt)

PUBLIC init_dTT
PUBLIC cesr_dTT
PUBLIC dest_dTT
PUBLIC get_dTT
PRIVATE modulator
PRIVATE elliptical_modulator

CONTAINS

!+
! Function id = init_dTT(incoming_tt_param)
! 
! Constructor for cesr_dTT module.  This function creates a new instance of the tune tracker by assigning
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
FUNCTION init_dTT(incoming_tt_param) RESULT(id)
  TYPE(tt_param_struct) :: incoming_tt_param
  INTEGER id
  CHARACTER id_str
  CHARACTER(12) tt_log_name

  tt_ids = tt_ids + 1
  IF(tt_ids .gt. max_tt) THEN
    WRITE(*,*) "Maximum number of tune trackers hard-coded to max_tt = ", max_tt
    STOP
  ENDIF
  id = tt_ids

  tt_param(id) = incoming_tt_param

  ! Initialize state variables
  tt_state(id)%deltaw = 0.0
  tt_state(id)%intDphi = 0.0
  tt_state(id)%bpm_msmt_last = 0.0
  tt_state(id)%gain = 0.001
  tt_state(id)%psi = 0.0

  !calculate lp filter parameter
  tt_param(id)%LPalpha = 1.0_rp / ( 1.0_rp + (tt_param(id)%LPinertia/2.0/pi) )

  WRITE(id_str,'(I1)') id
  tt_log_name = "tt_log_"//id_str//".out"
  OPEN(UNIT=500+id,NAME=tt_log_name,STATUS='REPLACE')

  IF(tt_param(id)%use_D_chan) THEN
    tt_param(id)%wls_id = initFixedWindowLS(tt_param(id)%wls_N,tt_param(id)%Dt,tt_param(id)%wls_order,2)
  ENDIF
END FUNCTION init_dTT

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
SUBROUTINE dest_dTT(id)
  INTEGER, INTENT(IN) :: id

  IF(tt_param(id)%use_D_chan) THEN
    CLOSE(500+id)
  ENDIF

END SUBROUTINE dest_dTT

!+
! Function z = cesr_dTT(bpm_msmt,id)
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
!   bpm_msmt     -- REAL(rp) [VALUE]: new bpm measurement. passed by value.
!   id           -- INTEGER, INTENT(IN): instance id of tune tracker
! Output:
!   z            -- Double precision: sin(modulator_angle + phase_advance_to_kicker)
!-
FUNCTION cesr_dTT(bpm_msmt,id) RESULT(z)
  IMPLICIT NONE

  REAL(rp), PARAMETER :: ga = 0.05     !bpm gain time constant

  REAL(rp) :: bpm_msmt [VALUE]
  INTEGER, INTENT(IN) :: id
  REAL(rp) :: z

  REAL(rp) t0  ! set to time at beginning of current cycle
  REAL(rp) AB  ! mixed signal
  REAL(rp) sinout, sqrout  ! modulator output
  REAL(rp) proDphi  !proportional to Dphi
  REAL(rp) dirDphi  !first derivative of Dphi
  REAL(rp) PIout ! P + I
  REAL(rp) bpm

  ! Adjust bpm data for closed orbit offset and apply gain
  bpm_msmt = bpm_msmt - tt_param(id)%offset
  ! Gain adjusted to keep bpm average at unity.
  tt_state(id)%gain = ga*ABS(bpm_msmt) + (1.0-ga)*tt_state(id)%gain
  bpm_msmt = bpm_msmt/tt_state(id)%gain

  t0=tt_state(id)%t
  !This loop represents the tune tracker fast elecronics (typically 183 per cesr period).
  DO WHILE(tt_state(id)%t < t0+tt_param(id)%Dt)
    IF( tt_state(id)%t < t0+tt_param(id)%Dt/2.0 ) THEN
      bpm = tt_state(id)%bpm_msmt_last
    ELSE
      bpm = bpm_msmt
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
                           tt_param(id)%Ki*tt_state(id)%Dphi*tt_param(id)%fastPeriod

    ! Update modulator angle
    tt_state(id)%psi = tt_state(id)%psi + &
                       (tt_param(id)%modw0 + tt_state(id)%deltaw)*tt_param(id)%fastPeriod
    tt_state(id)%psi = MOD(tt_state(id)%psi,2.0_rp*pi)

    tt_state(id)%t = tt_state(id)%t + tt_param(id)%fastPeriod
  ENDDO
  tt_state(id)%bpm_msmt_last = bpm_msmt

  !The following calculates the proportional channel and derivative channel obtained from an LS fit.
  proDphi = tt_param(id)%Kp * tt_state(id)%Dphi
  IF( tt_param(id)%use_D_chan ) THEN
    dirDphi = fixedWindowLS(tt_state(id)%intDphi,tt_param(id)%wls_id)
    PIout = tt_state(id)%intDphi + proDphi + dirDphi*tt_param(id)%Kd
    !This statement writes the state of each PID channel to the tt_log_n.out file.
    WRITE(500+id,'(4ES14.4)') tt_state(id)%t, tt_state(id)%intDphi, proDphi, dirDphi
  ELSE
    PIout = tt_state(id)%intDphi + proDphi
    WRITE(500+id,'(3ES14.4)') tt_state(id)%t, tt_state(id)%intDphi, proDphi
  ENDIF
  tt_state(id)%deltaw = tt_param(id)%Kvco * PIout

  CALL modulator(tt_state(id)%psi+tt_param(id)%phi_to_kicker,sinout,sqrout)
  !CALL elliptical_modulator(id,sinout)
  z = sinout

END FUNCTION cesr_dTT

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

!+
! Subroutine elliptical_modulator(id,waveout)
!
! Experimental modulator.  Not working and not in use.
!-
SUBROUTINE elliptical_modulator(id,waveout)
  INTEGER, INTENT(IN) :: id
  REAL(rp), INTENT(OUT) :: waveout
  REAL(rp) alpha, beta, x, psi

  psi = tt_state(id)%psi + tt_param(id)%phi_to_kicker
  alpha = tt_param(id)%alpha
  beta = tt_param(id)%beta

  x = (-alpha*SIN(psi)+COS(psi)) / -ABS(COS(psi)) / beta
  waveout = x/SQRT(1+x*x)
END SUBROUTINE elliptical_modulator

END MODULE tune_tracker_mod










