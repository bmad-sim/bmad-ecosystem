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
!   do
!     <track a turn, store in orb>
!     bpm_meas = orb(somewhere)%vec(1)
!     sinphi = TT_update(bpm_meas, id)
!     kick = sinphi * kickamplitude
!     <set element kick attribute>
!   enddo
!   call dest_dTT(id)
!
! init_dTT is the constructor and is passed the tune tracker parameters in tt_param_struct
! The contents of tt_param_struct are below
! TT_update is passed one BPM measurement per call.  
! It returns the kick amplitude (normalized to one) to be used on the next turn.
! dest_dTT is the destructor.  It is passed the id of the tune tracker to destroy.
!-

module tune_tracker_mod

use bmad_routine_interface
use windowls_mod

implicit none

integer, parameter :: max_tt = 3 !number of tune trackers to allocate memory for,
                                 !effectively maximum number of tune trackers

! Structure to hold parameters of tune tracker.
type tt_param_struct
  logical :: useSaveState = .false. ! If true then save TT state vars to file on exit,
                                    ! and load TT state vars on initialization.
  real(rp) phi_to_kicker     ! phase advance from x phase at bpm to xdot phase at kicker
  real(rp) Dt                ! time between data points (usually set to ring period)
  real(rp) LPinertia         ! inertia of lowpass filter that follows mixer. Try 2^15
  real(rp) Ki                ! Integrator gain
  real(rp) Kp                ! Proportional gain
  real(rp) Kvco              ! VCO gain
  real(rp) modw0             ! VCO base frequency (initial guess of beam frac tune)
  real(rp) offset            ! Closed orbit at bpm.
  character(3) mixmode       ! To mix BPM data with square wave or sine wave
  integer cyc_per_turn       ! Number of times tune tracker cycles per turn.
                             ! Usually an integer fraction of the harmonic number.
                             ! Try 183
  integer :: log_period = 1  ! How often to write to tt_log files.  Zero means do not log, one
                             ! means log every call to TT_update, three means every third call, etc.

  !Parameters specific to D channel
  !Note:  The D channel is not compatible with save states.
  logical :: use_D_chan = .FALSE.      ! Whether to use the D channel.
                                       ! Disable D channel if not in use
  real(rp) Kd             ! Differential gain
  integer wls_N           ! Number of data poins for window LS
  integer wls_order       ! Order of fit polynomial

  ! The following parameters are overwritten by the init_dTT
  real(rp) LPalpha        ! low-pass filter constant. Calculated from LPinertia
  integer wls_id          ! Instance ID for window LS module for D channel
  logical :: allocated = .FALSE.   ! Keeps track of whic tt_ids have been instantiated by user.

  ! The following are parameters that are likely to be needed by any program
  ! that makes use of the tune tracker module.
  ! This module does not actually use these parameters.
  ! For examples, see the tune tracker driver program.
  integer :: bpm_loc = -1        ! Location of BPM
  integer :: kck_loc = -1        ! Location of kicker
  character(40) :: bpm_name = ''    ! Name of BPM
  character(40) :: kck_name = ''    ! Name of Kicker
  real(rp) modTfrac0      ! initial fractional tune of kicker modulator
  real(rp) kickAmplitude  ! amplitude of kicks
  character orientation   ! 'h', 'v', or 'l'
  integer onum            ! Element of coord_struct used by bpm.
                          ! 1 for horiz, 3 for vert, 1 for longitudinal
end type tt_param_struct

! Structures to hold the state variables of the tune tracker.
type tt_state_struct
  real(rp) deltaw               ! w_vco = w0 + deltaw
  real(rp) intDphi              ! integrated Dphi
  real(rp) bpm_msmt_last        ! stashes previous bpm measurement
  real(rp) psi                  ! keeps track of modulator angle
  real(rp) Dphi                 ! Delta Phi - angle that comes out of filters mixer signal
  real(rp) gain                 ! bpm gain
  integer  counter              ! counts number of times TT_update called.  Needed for log file.
end type tt_state_struct

! Variables related to multiple TT instances
integer, private, save :: tts_instantiated = 0  !number of tune trackers instantiated
type(tt_param_struct), private, save :: tt_param(max_tt)
type(tt_state_struct), private, save :: tt_state(max_tt)
integer, private, save :: log_luns(max_tt)

public init_dTT
public TT_update
public dest_dTT
public get_dTT
private modulator
private check_id

contains

!+
! Function id = init_dTT(incoming_tt_param)
!
! Constructor for TT_update module.  This function creates a new instance of the tune tracker by assigning
! an id and copying the tune tracker parameters to the module saved data.  It also initializes the tune tracker's
! state variables.
!
! This function opens one log file "tt_log_<id>.out", which needs to be closed by calling the destructor.
!
! Input:
!   incoming_tt_param    -- type(tt_param_struct): Tune tracker parameters.  See structure definition for details.
!
! Output:
!   <return value>       -- integer: id of the tune tracker instance created.
!-
function init_dTT(incoming_tt_param,saved_coords) result(id)
  type(tt_param_struct) :: incoming_tt_param
  type(coord_struct), optional :: saved_coords
  integer i
  integer id
  integer ios
  character id_str
  character(12) tt_log_name
  character(10) state_name
  logical state_read

  if(tts_instantiated + 1 .gt. max_tt) then
    write(*,*) "Maximum number of tune trackers hard-coded to max_tt = ", max_tt
    stop
  endif

  tts_instantiated = tts_instantiated + 1

  do i=1,max_tt
    if( .not. tt_param(i)%allocated ) then
      id = i 
      tt_param(id) = incoming_tt_param
      tt_param(i)%allocated = .true.
      exit
    endif
  enddo

  write(id_str,'(I1)') id

  ! Initialize state variables
  state_read = .false.
  if(tt_param(id)%useSaveState) then
    if(present(saved_coords)) then
      state_name = "tt_state."//id_str
      open(unit=9999,FILE=state_name,status='OLD',ACTION='READ',IOSTAT=ios)
      if(ios == 0) then
        READ(9999,*) saved_coords, tt_state(id)
        state_read = .true.
        close(9999)
        write(*,'(A,I2,A)') "NOTICE: Save state found for TT #", id, ":"
        write(*,'(A,F10.7)') "        VCO frequency set to ", (tt_state(id)%deltaw+tt_param(id)%modw0)* &
                                                    tt_param(id)%Dt / 2.0_rp / pi
      else
        write(*,*) "NOTICE: Tune tracker warning for TT #", id, ":"
        write(*,*) "        Unable to open file: ", state_name
        write(*,*) "        Setting initial TT state variables to defaults."
        write(*,*) "        Leaving initial coordinates unchanged."
      endif
    endif
  endif
  if(.not. state_read) then
    tt_state(id)%deltaw = 0.0_rp
    tt_state(id)%intDphi = 0.0_rp
    tt_state(id)%bpm_msmt_last = 0.0_rp
    tt_state(id)%psi = 0.0_rp
    tt_state(id)%Dphi = 0.0_rp
    tt_state(id)%gain = 0.001_rp
    tt_state(id)%counter = 0
  endif

  !calculate low-pass filter parameter
  tt_param(id)%LPalpha = 1.0_rp / ( 1.0_rp + (tt_param(id)%LPinertia/2.0_rp/pi) )

  tt_log_name = "tt_log_"//id_str//".out"
  log_luns(id) = lunget()
  open(unit=log_luns(id),FILE=tt_log_name,status='REPLACE')

  write(log_luns(id),'(a)') "# Tune tracker channels"
  write(log_luns(id),'(a1,a7,3a14)') "#", "turn", "I", "P", "D"

  if(tt_param(id)%use_D_chan) then
    if( .not. tt_param(id)%useSaveState ) then
      tt_param(id)%wls_id = initFixedWindowLS(tt_param(id)%wls_N,tt_param(id)%Dt,tt_param(id)%wls_order,1)
    else
      write(*,*) "WARNING: D channel calculations are not compatible with save states."
      write(*,*) "         D channel will be disabled."
      tt_param(id)%use_D_chan = .FALSE.
    endif
  endif
end function init_dTT

!+
! Subroutine reset_dTT(id,tt_param)
!
! Reset a tune tracker and optionally change its parameters.
!
! Input:
!   id        -- integer, intent(in): Tune tracker instance ID obtained from constructor.
!   tt_param  -- 
! Output:
!   none
!-
subroutine reset_dTT(id,tt_param_update)
  integer, intent(in) :: id
  type(tt_param_struct), optional :: tt_param_update

  call check_id(id)

  ! Reset state variables
  tt_state(id)%deltaw = 0.0_rp
  tt_state(id)%intDphi = 0.0_rp
  tt_state(id)%bpm_msmt_last = 0.0_rp
  tt_state(id)%psi = 0.0_rp
  tt_state(id)%Dphi = 0.0_rp
  tt_state(id)%gain = 0.001_rp
  tt_state(id)%counter = 0

  ! Update tt_param if tt_param_update is present.  Only certain fields can be updated.
  if(present(tt_param_update)) then
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
  endif

  ! Write two blank lines to the log file to indicate that a reset has occurred
  write(log_luns(id),*)
  write(log_luns(id),*)

end subroutine reset_dTT

!+
! Subroutine dest_dTT(id)
!
! Destructor for tune tracker.  Should be called end of program to ensure all files are properly closed.
!
! Input:
!   id     -- integer, intent(in): Tune tracker instance ID obtained from constructor.
! Output:
!   none
!-
subroutine dest_dTT(id,coords)
  integer, intent(in) :: id
  character id_str
  character(10) state_name
  type(coord_struct), optional :: coords

  call check_id(id)

  tts_instantiated = tts_instantiated - 1

  tt_param(id)%allocated = .false.

  if(tt_param(id)%useSaveState) then
    if(present(coords)) then
      write(id_str,'(I1)') id
      state_name = "tt_state."//id_str
      open(unit=9999,FILE=state_name,status='REPLACE',ACTION='write')
      write(9999,*) coords, tt_state(id)
      close(9999)
    endif
  endif

  close(log_luns(id))
end subroutine dest_dTT

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
! Input:
!   bpm_msmt     -- real(rp): new bpm measurement. passed by value.
!   id           -- integer, intent(in): instance id of tune tracker
! Output:
!   z            -- Double precision: sin(modulator_angle + phase_advance_to_kicker)
!-
function TT_update(bpm_msmt,id) result(z)
  implicit none

  real(rp), parameter :: ga = 0.05_rp     !bpm gain time constant

  real(rp) :: bpm_msmt
  real(rp) bpm_msmt_nrml
  integer, intent(in) :: id
  real(rp) :: z

  real(rp) AB  ! mixed signal
  real(rp) sinout, sqrout  ! modulator output
  real(rp) proDphi  !proportional to Dphi
  real(rp) dirDphi  !first derivative of Dphi
  real(rp) PIDout ! P + I + D
  real(rp) bpm
  real(rp) fastPeriod

  integer i

  call check_id(id)

  tt_state(id)%counter = tt_state(id)%counter + 1    !needed for log file

  ! Adjust bpm data for closed orbit offset and apply gain
  bpm_msmt_nrml = bpm_msmt - tt_param(id)%offset
  ! Update gain.  Gain adjusted to keep bpm average at unity.
  tt_state(id)%gain = ga*ABS(bpm_msmt_nrml) + (1.0_rp-ga)*tt_state(id)%gain
  ! Apply gain.
  bpm_msmt_nrml = bpm_msmt_nrml/tt_state(id)%gain

  fastPeriod = tt_param(id)%Dt / tt_param(id)%cyc_per_turn

  !This loop represents the tune tracker fast elecronics (typically 14ns)
  do i=1, tt_param(id)%cyc_per_turn
    if( i <= FLOOR(tt_param(id)%cyc_per_turn / 2.0) ) then
      bpm = tt_state(id)%bpm_msmt_last
    else
      bpm = bpm_msmt_nrml
    endif
    call modulator(tt_state(id)%psi,sinout,sqrout)

    ! Phase Detector: mixer followed by low-pass filter
    if( tt_param(id)%mixmode == 'sqr' ) then
      AB = bpm*sqrout
    elseif( tt_param(id)%mixmode == 'sin' ) then
      AB = bpm*sinout
    else
      write(*,*) "Unknown mixmode specified.  Must be 'sin' or 'sqr'.  Terminating."
      stop
    endif
    tt_state(id)%Dphi = AB*tt_param(id)%LPalpha + (1.0_rp-tt_param(id)%LPalpha)*tt_state(id)%Dphi

    ! Calculate Integral Channel
    tt_state(id)%intDphi = tt_state(id)%intDphi + &
                           (tt_param(id)%Ki/tt_param(id)%cyc_per_turn)*tt_state(id)%Dphi

    ! Update modulator angle
    tt_state(id)%psi = tt_state(id)%psi + (tt_param(id)%modw0 + tt_state(id)%deltaw)*fastPeriod
    tt_state(id)%psi = mod(tt_state(id)%psi,2.0_rp*pi)

  enddo
  tt_state(id)%bpm_msmt_last = bpm_msmt_nrml

  !The following calculates the proportional and derivative channels
  proDphi = tt_param(id)%Kp * tt_state(id)%Dphi
  if( tt_param(id)%use_D_chan ) then
    dirDphi = fixedWindowLS(proDphi,tt_param(id)%wls_id)
    PIDout = tt_state(id)%intDphi + proDphi - dirDphi*tt_param(id)%Kd
  else
    PIDout = tt_state(id)%intDphi + proDphi
  endif

  !Update modulator frequency
  tt_state(id)%deltaw = tt_param(id)%Kvco * PIDout

  !Write to log file
  if( tt_param(id)%log_period .gt. 0 ) then
    if( mod( (tt_state(id)%counter+tt_param(id)%log_period-1), tt_param(id)%log_period) == 0 ) then
      if( tt_param(id)%use_D_chan ) then
        !This statement writes the state of each PID channel to the tt_log_n.out file.
        write(log_luns(id),'(I8,3ES14.4)') tt_state(id)%counter, tt_param(id)%Ki*tt_state(id)%intDphi, &
                                                                 tt_param(id)%Kp*proDphi, &
                                                                 -tt_param(id)%Kd*dirDphi
      else
        write(log_luns(id),'(I8,2ES14.4)') tt_state(id)%counter, tt_param(id)%Ki*tt_state(id)%intDphi, &
                                                                 tt_param(id)%Kp*proDphi
      endif
    endif
  endif

  call modulator(tt_state(id)%psi + tt_param(id)%phi_to_kicker, sinout, sqrout)
  z = sinout ! z is result, the value returned by this function

end function TT_update

!+
! Function z = get_dTT(name,id)
!
! A get function to obtain the tune tracker VCO frequency.
!
! Input:
!   name     -- character(2), intent(in):  'dw' to get VCO trim, 'wf' to get VCO frequency
!   id       -- integer, intent(in): Instance id
! Output:
!   <return value>  -- real(rp): value of requested data
!-
function get_dTT(name,id) result(z)
  character(2), intent(in) :: name
  integer, intent(in) :: id
  real(rp) z

  call check_id(id)

  if(name == 'dw') then
    ! return VCO trim 
    z = tt_state(id)%deltaw
  elseif(name == 'wf') then
    ! return VCO frequency
    z = tt_state(id)%deltaw + tt_param(id)%modw0
  elseif(name == 'ps') then
    ! return VCO frequency
    z = tt_state(id)%psi
  else
    z = -999.
  endif
end function get_dTT

!+
! Subroutine modulator(psi,sinout,sqrout)
!
! Tune tracker modulator.  Returns sin(psi) and corresponding square wave.
!
! Input:
!   psi     -- real(rp), intent(in):  modulator angle
! Output:
!   sinout  -- real(rp), intent(out):  sin(psi)
!   sqrout  -- real(rp), intent(out):  +1 if sin(psi)>=0, -1 else
!-
subroutine modulator(psi,sinout,sqrout)
  real(rp), intent(in) :: psi
  real(rp), intent(out) :: sinout
  real(rp), intent(out) :: sqrout

  sinout = sin(psi)
  
  if(sinout >= 0._rp) then
    sqrout = 1.0_rp
  else
    sqrout = -1.0_rp
  endif

end subroutine modulator

subroutine check_id(id)
  integer, intent(in) :: id
  
  if( (id .lt. 1) .or. id .gt. max_tt) then
    !Absurd TT id passed to TT module
    write(*,'(A,I3)') "FATAL in call TO TUNE TRACKER: ID received is ", id
    write(*,'(A,I3)') "                  ID should be between 1 and ", max_tt
    write(*,'(A)')    "                  Check that init_dTT has been called."
    write(*,'(A)')    "                  Check that ID is an integer."
    call err_exit
  endif

  if( .not. tt_param(id)%allocated ) then
    !Tune tracker not initialized
    write(*,'(A,I3)') "FATAL in call TO TUNE TRACKER: ID received is ", id
    write(*,'(A)')    "                A tune tracker with this ID has not been instantiated."
    write(*,'(A)')    "                Check that init_dTT has been called and dest_dTT has not been called."
    write(*,'(A)')    "                Check that ID is an integer."
    call err_exit
  endif
end subroutine check_id

end module tune_tracker_mod
