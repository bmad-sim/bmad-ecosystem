!+
! Program tune_tracker
!
! Driver program for tune_tracker_mod.  This module creates tune trackers based on the settings in tune_tracker.in
! Tracking is done with damping on and fluctuations off.
! BPM data at element 1 is stored and FFT analysis is done to determine horizontal and vertical tune.
!
! The following log files are opened for each tune tracker:
!   bpm_meas_<id>.out  -- turn number, x at bpm, xdot at bpm
!   kck_meas_<id>.out  -- turn number, x at kicker, xdot at kicker
!   vco_stat_<id>.out  -- turn number, vco frequency in fractional tune units
!   mod_stat_<id>.out  -- tuen number, modulator amplitude
!
! The FFT is stored in the following two files.
! The units of the x-axis of the FFT output files is fractional tune.
!   xfft.out    -- FFT of horizontal data at element 1 over last half of turns
!   yfft.out    -- FFT of vertical data at element 1 over last half of turns
!-
program tune_tracker_driver

  use bmad
  use mode3_mod
  use tune_tracker_mod
  implicit none

  type(lat_struct) ring
  type(coord_struct), allocatable :: orb(:)
  type(coord_struct), allocatable :: co(:)
  type(tt_param_struct) :: tt_params(max_tt) !max_tt set on tune_tracker module

  character(20) inputfile
  character(200) lat_file,bpm_out_file,vco_out_file,kck_out_file,mod_out_file
  character(200) xfft_out_file, yfft_out_file

  integer nTTs      ! number of tune trackers
  integer i,j       ! loop counter
  integer nturns    ! number of turns
  integer nfft
  integer n_ele_track
  integer id(max_tt)

  real(rp) Tring  ! period of ring
  real(rp) bpmdata
  real(rp) sinphi
  real(rp) aFracTune, bFracTune, zFracTune, kick
  real(rp) Deltaphi
  real(rp) TTw
  real(rp), allocatable :: xdata(:), ydata(:)
  complex(rp), allocatable :: xfft(:), yfft(:)
  character(1) i_str

  logical error

  namelist /TT_in/ lat_file, nturns, nTTs, tt_params

  inputfile = "tune_tracker.in"
  open(unit=2,file=inputfile)
  read(2,nml=TT_in)
  close(2)

  ! Parse the ring and set it up for tracking
  call bmad_parser(lat_file,ring)
  n_ele_track = ring%n_ele_track
  call set_on_off(rfcavity$, ring, on$)
  bmad_com%radiation_damping_on = .true.
  bmad_com%radiation_fluctuations_on = .false.
  call twiss3_at_start(ring,error)
  call twiss3_propagate_all(ring)
  call closed_orbit_calc(ring,orb,6)

  ! Save Closed Orbit
  allocate(co(0:n_ele_track))
  co = orb

  ! Check if the element for each tune tracker's kicker has kick attributes.
  do i=1,nTTs
    if( ring%ele(tt_params(i)%kck_loc)%key == hkicker$ ) then
      if( tt_params(i)%orientation /= 'h' ) then
        write(*,'(A,I2,A)') "ERROR: Orientation of tune tracker ", i," does not match kicker orientation"
        write(*,'(A)') "TERMINATING"
        stop
      endif
      if( .not. attribute_free(tt_params(i)%kck_loc,'KICK',ring,.false.) ) then
        write(*,'(A,I2,A)') "ERROR: kick attribute for TT ", i, " kicker is not free."
        write(*,'(A)') "TERMINATING"
        stop
      endif
    elseif( ring%ele(tt_params(i)%kck_loc)%key == vkicker$ ) then
      if( tt_params(i)%orientation /= 'v' ) then
        write(*,'(A,I2,A)') "ERROR: Orientation of tune tracker ", i," does not match kicker orientation"
        write(*,'(A)') "TERMINATING"
        stop
      endif
      if( .not. attribute_free(tt_params(i)%kck_loc,'KICK',ring,.false.) ) then
        write(*,'(A,I2,A)') "ERROR: kick attribute for TT ", i, " kicker is not free."
        write(*,'(A)') "TERMINATING"
        stop
      endif
    else
      if( .not. has_hkick_attributes(ring%ele(tt_params(i)%kck_loc)%key) ) then
        write(*,'(A,I2,A)') "ERROR: Element for TT ", i, " kicker does not have kick attributes."
        write(*,'(A)') "TERMINATING"
        stop
      endif
      if( tt_params(i)%orientation == 'h' ) then
        if( .not. attribute_free(tt_params(i)%kck_loc,'HKICK',ring,.false.) ) then
          write(*,'(A,I2,A)') "ERROR: hkick attribute for TT ", i, " kicker is not free."
          write(*,'(A)') "TERMINATING"
          stop
        endif
      elseif( tt_params(i)%orientation == 'v' ) then
        if( .not. attribute_free(tt_params(i)%kck_loc,'VKICK',ring,.false.) ) then
          write(*,'(A,I2,A)') "ERROR: vkick attribute for TT ", i, " kicker is not free."
          write(*,'(A)') "TERMINATING"
          stop
        endif
      elseif( tt_params(i)%orientation == 'z' ) then
        if( .not. attribute_free(tt_params(i)%kck_loc,'PHI0',ring,.false.) ) then
          write(*,'(A,I2,A)') "ERROR: phi0 attribute for TT ", i, " kicker is not free."
          write(*,'(A)') "TERMINATING"
          stop
        endif
      endif
    endif
  enddo

  !calculate ring period
  Tring = ring%ele(n_ele_track)%s / c_light

  aFracTune = mod(ring%ele(n_ele_track)%mode3%a%phi,2.0_rp*pi)/2.0_rp / pi
  bFracTune = mod(ring%ele(n_ele_track)%mode3%b%phi,2.0_rp*pi)/2.0_rp / pi
  zFracTune = mod(ring%ele(n_ele_track)%mode3%c%phi,2.0_rp*pi)/2.0_rp / pi
  write(*,'(A50,F10.7)') "Calculated (twiss_and_track) Fractional Tune (a): ", aFracTune
  write(*,'(A50,F10.7)') "Calculated (twiss_and_track) Fractional Tune (b): ", bFracTune
  write(*,'(A50,F10.7)') "Calculated (twiss_and_track) Fractional Tune (z): ", zFracTune

  !initialize each tune tracker
  do i=1,nTTs
    !calculate initial frequency of modulator
    write(*,'(A14,I3,A21,F8.5,A22)') &
         "Tune tracker #",i, " VCO base frequency: ", tt_params(i)%modTfrac0, " oscillations per turn"
    tt_params(i)%modw0 = 2.0_rp*pi/(Tring/tt_params(i)%modTfrac0)   ! 2pi/(initial period of modulator)

    ! Determine phase between BPM and kicker
    if(tt_params(i)%orientation == 'h') then
      !The following calculation works because the PLL leads the BPM data by pi/2
      Deltaphi = ring%ele(n_ele_track)%mode3%a%phi - &
                 ring%ele(tt_params(i)%bpm_loc)%mode3%a%phi + ring%ele(tt_params(i)%kck_loc)%mode3%a%phi
      tt_params(i)%Onum = 1  !element of coord_struct for BPM to observe
      !Normalize kick amplitude by beta at kicker
      tt_params(i)%kickAmplitude = tt_params(i)%kickAmplitude / SQRT(ring%ele(tt_params(i)%kck_loc)%mode3%a%beta)
    elseif(tt_params(i)%orientation == 'v') then
      !The following calculation works because the PLL leads the BPM data by pi/2
      Deltaphi = ring%ele(n_ele_track)%mode3%b%phi - &
                 ring%ele(tt_params(i)%bpm_loc)%mode3%b%phi + ring%ele(tt_params(i)%kck_loc)%mode3%b%phi
      tt_params(i)%Onum = 3  !element of coord_struct for BPM to observe
      !Normalize kick amplitude by beta at kicker
      tt_params(i)%kickAmplitude = tt_params(i)%kickAmplitude / SQRT(ring%ele(tt_params(i)%kck_loc)%mode3%b%beta)
    elseif(tt_params(i)%orientation == 'z') then
      Deltaphi = -pi/2.0_rp + ring%ele(n_ele_track)%mode3%c%phi - &
                 ring%ele(tt_params(i)%bpm_loc)%mode3%c%phi + ring%ele(tt_params(i)%kck_loc)%mode3%c%phi
      tt_params(i)%Onum = 1  !element of coord_struct for BPM to observe
    endif
    tt_params(i)%phi_to_kicker = Deltaphi
    tt_params(i)%phi_to_kicker = mod(tt_params(i)%phi_to_kicker,2.0*pi) 

    tt_params(i)%Dt = Tring
    tt_params(i)%offset = orb(tt_params(i)%bpm_loc)%vec(tt_params(i)%Onum)

    id(i) = init_dTT(tt_params(i), orb(0))
  enddo

  ! Allocate data arrays for FFT
  allocate(xdata(nturns))
  allocate(ydata(nturns))
  allocate(xfft(nturns/2))
  allocate(yfft(nturns/2))

  ! Open output files for each kicker
  do i=1,nTTs
    write(i_str,'(I1)') i
    bpm_out_file = "bpm_meas_"//i_str//".out"
    kck_out_file = "kck_meas_"//i_str//".out"
    vco_out_file = "vco_stat_"//i_str//".out"
    mod_out_file = "mod_stat_"//i_str//".out"
    open(unit=(i*100+20),file=bpm_out_file)
    open(unit=(i*100+21),file=kck_out_file)
    open(unit=(i*100+22),file=vco_out_file)
    open(unit=(i*100+23),file=mod_out_file)

    write(i*100+20,'(a1,a,i6,"   ",a)') "#", "Measurements at bpm ", tt_params(i)%bpm_loc, ring%ele(tt_params(i)%bpm_loc)%name
    write(i*100+20,'(a1,a9,2a14)') "#", "turn", tt_params(i)%orientation, tt_params(i)%orientation//"'"

    write(i*100+21,'(a1,a,i6,"   ",a)') "#", "Measurements at kicker ", tt_params(i)%kck_loc, ring%ele(tt_params(i)%kck_loc)%name
    write(i*100+21,'(a1,a9,2a14)') "#", "turn", tt_params(i)%orientation, tt_params(i)%orientation//"'"

    write(i*100+22,'(a1,a9,a)') "#", "turn", "   VCO Frequency (fractional, unit)"

    write(i*100+23,'(a1,a9,a)') "#", "turn", "   Kicker modulator sinphi"
  enddo

  ! Print information about the location of each TT's kicker and BPM
  do i=1,nTTs
    write(*,'(A,I2,A,A,A,A)') "Kicker #", i, " at element: ", ring%ele(tt_params(i)%kck_loc)%name, &
                                  " Type: ", key_name(ring%ele(tt_params(i)%kck_loc)%key)
    write(*,'(A,I2,A,A,A,A)') "   BPM #", i, " at element: ", ring%ele(tt_params(i)%bpm_loc)%name, &
                                  " Type: ", key_name(ring%ele(tt_params(i)%bpm_loc)%key)
  enddo


  ! Main Loop
  sinphi = 0.0_rp
  do i=1, nturns
    call track_all(ring,orb)
    orb(0) = orb(n_ele_track)

    do j=1,nTTs
      bpmdata = orb(tt_params(j)%bpm_loc)%vec(tt_params(j)%Onum)
      write((j*100+23),'(I10,3ES14.6)') i, sinphi
      sinphi = TT_update(bpmdata,id(j))
      kick = sinphi*tt_params(j)%kickAmplitude
      if( ring%ele(tt_params(j)%kck_loc)%key == hkicker$ ) then
        ring%ele(tt_params(j)%kck_loc)%value(kick$) = kick
      elseif( ring%ele(tt_params(j)%kck_loc)%key == vkicker$ ) then
        ring%ele(tt_params(j)%kck_loc)%value(kick$) = kick
      else
        if(tt_params(j)%orientation == 'h') then
          ring%ele(tt_params(j)%kck_loc)%value(hkick$) = kick
        elseif(tt_params(j)%orientation == 'v') then
          ring%ele(tt_params(j)%kck_loc)%value(vkick$) = kick
        elseif(tt_params(j)%orientation == 'z') then
          ring%ele(tt_params(j)%kck_loc)%value(phi0$) = kick
        endif
      endif

      write((j*100+20),'(I10,2ES14.6)') i, orb(tt_params(j)%bpm_loc)%vec(tt_params(j)%Onum), &
                                           orb(tt_params(j)%bpm_loc)%vec(tt_params(j)%Onum+1)
      write((j*100+21),'(I10,2ES14.6)') i, orb(tt_params(j)%kck_loc)%vec(tt_params(j)%Onum), &
                 orb(tt_params(j)%kck_loc)%vec(tt_params(j)%Onum+1)-co(tt_params(j)%kck_loc)%vec(tt_params(j)%Onum+1)
      TTw = get_dTT('wf',id(j))
      write((j*100+22),'(I10,ES14.6)') i, TTw*Tring/2.0_rp/pi
    enddo

    ! FFT data gathered at element 1
    xdata(i) = orb(1)%vec(1) - co(1)%vec(1)
    ydata(i) = orb(1)%vec(3) - co(1)%vec(3)

    if( mod(i,1000) == 0 ) write(*,*) "Progress: ", i, "/", nturns
  enddo

  do i=1,nTTs
    close(i*100+20)
    close(i*100+21)
    close(i*100+22)
    close(i*100+23)
    write(*,'(A,I2,A,F10.7,A)') "Tune tracker ",i," VCO final frequency: ", &
                            get_dTT('wf',id(i))*Tring/2.0_rp/pi, " oscillations per turn"
    call dest_dTT(id(i),orb(n_ele_track))
  enddo

  ! FFT Calculation
  nfft = nturns / 2
  do i=1,nfft
    xfft(i) = cmplx(xdata(nturns-nfft+i),0.)
    yfft(i) = cmplx(ydata(nturns-nfft+i),0.)
  enddo
  call fft_1d(xfft,1)
  call fft_1d(yfft,1)
  xfft_out_file = "xfft.out"
  yfft_out_file = "yfft.out"
  open(unit=24,file=xfft_out_file)
  open(unit=25,file=yfft_out_file)
  do i=1,nfft
    write(24,'(F11.7,ES14.5)') (i-1)/real(nfft),abs(xfft(i))
    write(25,'(F11.7,ES14.5)') (i-1)/real(nfft),abs(yfft(i))
  enddo
  close(24)
  close(25)
  deallocate(xdata)
  deallocate(ydata)
  deallocate(xfft)
  deallocate(yfft)

end program tune_tracker_driver





