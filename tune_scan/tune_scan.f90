!+
! Program       : TUNE_SCAN
!
! Description : Program to simulate a tunescan.
!
! Arguments  :
!
! Mod/Commons:
!
! Calls      :
!
! Author     : 
!
! Modified   :
!-
!........................................................................
!
! $Id$
!
! $Log$
!
!Revision 2007/06/08 Jim Shanks
!Added support for space charge module
!
!
! Revision 1.2  2007/01/30 16:14:32  dcs
! merged with branch_bmad_1.
!
! Revision 1.1.1.1.2.1  2006/12/22 20:30:43  dcs
! conversion compiles.
!
! Revision 1.1.1.1  2005/06/14 14:59:02  cesrulib
! Beam Simulation Code
!
!
!........................................................................
!
#include "CESR_platform.inc"

program tune_scan

  use bmad_struct
  use bmad_interface
  use bmadz_interface
  use bsim_interface
  use scan_parameters
  use trans_space_charge_mod
 
  implicit none

  type (lat_struct) ring, ring_ok, ring_in
  type (coord_struct), allocatable :: orb(:)
  type (coord_struct) init(100)
  type (scan_params_struct) scan_params
!mode declared for radiation integrals
  type (normal_modes_struct) mode


!Local variables for scan_params removed
  integer i, j, ie
  integer n_x, n_y, i_x, i_y, i_ele, n_z
  integer i_z
  integer int_Q_x, int_Q_y, ix, ixx,ixxx
  integer kx, n_arg
  integer ios
  integer version
  integer ix_cache

  real(rp), allocatable :: dk1(:) 
  real(rp) phi_x, phi_y
  real(rp) a(6,6), wr(6), wi(6), tune(6)
  real(rp) test
  real(rp)  Q_x0, Q_x1, dQ_x, Q_y0, Q_y1, dQ_y
  real(rp) Q_z
  real(rp) Q_x, Q_y, Q_z0, Q_z1, dQ_z
  real(rp) Q_x_init, Q_y_init
!New real(rp) allowing the user to define the a- and b-mode (horizontal and vertical) emittance
  real(rp) a_emittance, b_emittance

  character(200) :: file_name = '', out_file, new_lat_file
  character(10) date, time
  character(2) wordx, wordy

  logical ok, found
!New logical for user to switch space charge module on
  logical trans_space_charge_on
  logical fft/.false./


!Now setting scan_params here, instead of local variables to be input into scan_params
  namelist / parameters / scan_params, Q_x0, Q_x1, dQ_x, n_x, Q_y0, Q_y1, dQ_y, n_y, &
                Q_z0, Q_z1, dQ_z, Q_x_init, Q_y_init, trans_space_charge_on, &
                a_emittance, b_emittance, fft
                

! read in the parameters

  n_arg = cesr_iargc()
  if (n_arg > 1) then
    print *, 'Usage: tune_scan <input_file>'
    print *, 'Default: <input_file> = tune_scan.in'
    stop
  endif

  if (n_arg == 1) call cesr_getarg(1, file_name)

  if (n_arg == 0) then
    print '(a, $)', ' Input command file <CR=tune_scan.in>: '
    read (*, '(a)') file_name
  endif

  if (file_name == '')   file_name = 'tune_scan.in'

  print *, 'Opening: ', trim(file_name)
  open (unit= 1, file = file_name, status = 'old')

  scan_params%particle = positron$
  scan_params%beambeam_ip = .false.
  scan_params%close_pretz = .false.
  scan_params%close_vert = .false.
  scan_params%rec_taylor = .false.
  scan_params%slices = 1
  scan_params%min_sig = 0.
  scan_params%coupling_sb = 0.01
  n_x =1
  n_y=1
  n_z = 1
  scan_params%radiation = .false.
!Initialize the use of space charge module to 'false'
  trans_space_charge_on = .false.
  a_emittance = 0
  b_emittance = 0
  ix_cache = 0
  

  read(1, nml = parameters)



  if (n_x == 1 .and. dQ_x /= 0.)n_x = nint(abs((Q_x1 - Q_x0) / dQ_x)) + 1
  if (n_y == 1 .and. dQ_y /= 0.)n_y = nint(abs((Q_y1 - Q_y0) / dQ_y)) + 1
  if (n_z == 1 .and. dQ_z /= 0.)n_z = nint(abs((Q_z1 - Q_z0) / dQ_z)) + 1 

  if(n_x > 1)dQ_x = (Q_x1 - Q_x0) / (n_x - 1)
  if(n_y > 1)dQ_y = (Q_y1 - Q_y0) / (n_y - 1)


  if(n_x /= 1)then
    dQ_x = (Q_x1 - Q_x0) / (n_x - 1)
  else
    dQ_x = 0.
  endif
  if(n_y /= 1)then
    dQ_y = (Q_y1 - Q_y0) / (n_y - 1)
  else
   dQ_y = 0.
  endif
  if(n_z /= 1)then
   dQ_z = (Q_z1 - Q_z0) / (n_z - 1)
  else
   dQ_z = 0.
  endif
  if(Q_z0 == 0.)Q_z0 = Q_z
! init lattice

  if(scan_params%lat_file(1:8) == 'digested')then
    call read_digested_bmad_file (scan_params%lat_file, ring, version)
   else
    call bmad_parser (scan_params%lat_file, ring)
  endif

!Set the bunch size from the user-defined ring current

  ring%param%n_part   = scan_params%current*0.001*(ring%param%total_length/c_light)/e_charge

! find how many starting points there are

  print *, 'n_x, n_y:', n_x, n_y
  print *, 'dQ_x, dQ_y:', dQ_x, dQ_y


  ix = index(file_name,'.')
  out_file = file_name(1:ix)//'hed'

  ixx = max(index(scan_params%lat_file,'/'),1)
  ixxx = len_trim(scan_params%lat_file)
  kx = ixxx-ixx+1
  new_lat_file = scan_params%lat_file(ixx:ixxx)
  new_lat_file(kx+1:kx+1) = "'"

  call date_and_time(date, time)


  open(unit=21, file = out_file)
  write (21, '(a13,a200)') 'Lattice   = `', new_lat_file
  write (21, *) 'File_name = `', trim(file_name)// "'"
  write (21, '(a13,1x,a4,a1,a2,a1,a2,1x,a2,a1,a2,a1,a2,a1)') 'data_date = `', &
                   date(1:4),'-',date(5:6),'-',date(7:8),time(1:2),':',time(3:4),':',time(5:6), "'"
  write (21,'(a11,i5)' )'  n_part = ',scan_params%n_part
  write (21,'(a11, i7)')'  n_turn = ',scan_params%n_turn
  write (21, *) 'Return'
  close( unit=21)

  out_file = file_name(1:ix)//'out'
  open(unit=21, file = out_file, status='replace')
  write (21, '(a13,a200)') 'Lattice   = `', new_lat_file
  write (21, *) 'File_name = `', trim(file_name)// "'"
  write (21, '(a13,1x,a4,a1,a2,a1,a2,1x,a2,a1,a2,a1,a2,a1)') 'data_date = `', &
                   date(1:4),'-',date(5:6),'-',date(7:8),time(1:2),':',time(3:4),':',time(5:6), "'"
  write (21, *) 'Return'
  write (21, *)

  write(21, '(5(a, f5.3))') '!  Q_x0 = ', Q_x0, '  Q_y0 = ', Q_y0, &
                            '  Q_x1 = ', Q_x1, '  Q_y1 = ', Q_y1, '  Q_z = ', ring%z%tune/twopi

  write(21,*)
  if(.not. fft)write(21,'(a78,a10)')'!   Q_x       Q_y       Q_z      x_end       y_end      amp_x_max   amp_y_max ',&
                   '  Turns   '
  if(fft)write(21,'(a78,11x,a5,14x,a8,14x,8x,a20)')'!   Q_x       Q_y       Q_z      x_end       y_end      amp_x_max   amp_y_max ',&
                   'Turns', 'Q(1:512)', 'Q(n_turns-512+1:512)'

  scan_params%file_name = out_file

  close(unit=21)

  do i = 1, ring%n_ele_max
    if(ring%ele(i)%value(x1_limit$) == 0) ring%ele(i)%value(x1_limit$) = 0.05
    if(ring%ele(i)%value(x2_limit$) == 0) ring%ele(i)%value(x2_limit$) = 0.05
    if(ring%ele(i)%value(y1_limit$) == 0) ring%ele(i)%value(y1_limit$) = 0.05
    if(ring%ele(i)%value(y2_limit$) == 0) ring%ele(i)%value(y2_limit$) = 0.05
  enddo

  ring.param.aperture_limit_on = .true.

  allocate(dk1(ring%n_ele_max))
  call reallocate_coord (orb, ring%n_ele_max)

  if(Q_x_init == 0)Q_x_init = Q_x0
  if(Q_y_init == 0)Q_y_init = Q_y0

! find which quads to change

  call twiss_and_track(ring, orb)

  call choose_quads(ring, dk1)
  int_Q_x = int(ring%ele(ring%n_ele_track)%a%phi / twopi)
  int_Q_y = int(ring%ele(ring%n_ele_track)%b%phi / twopi)
  phi_x = (int_Q_x + Q_x_init) * twopi
  phi_y = (int_Q_y + Q_y_init) * twopi


! Call radiation integrals -- ix_cache initialized to zero
! to account for wiggler-dominated emittance
    call radiation_integrals(ring, orb, mode,ix_cache)

 !If a-mode emittance is defined by the user,
 !otherwise use user-defined a-mode emittance value

  if(a_emittance .ne. 0)then
     mode%a%emittance = a_emittance
  endif

! Set b-mode emittance from user-defined value
! (will initialize to 1% coupling to horizontal if no user input)

  if(b_emittance == 0)then
     mode%b%emittance = 0.01 * mode%a%emittance
    else
     mode%b%emittance = b_emittance
  endif

write(*,*) "A-MODE EMITTANCE: ", mode%a%emittance
write(*,*) "B-MODE EMITTANCE: ", mode%b%emittance


! Call setup_trans_space_charge_calc

  call setup_trans_space_charge_calc(trans_space_charge_on, ring, ring%param%n_part, mode)

  call custom_set_tune(phi_x, phi_y, dk1, ring, orb, ok)    

! scan the tune

  ring_ok = ring

  print *, 'Number of starting positions:', scan_params%n_part
  print *, 'n_x, n_y, n_z:', n_x, n_y, n_z
  print *, 'dQ_x, dQ_y, dQ_z:', dQ_x, dQ_y, dQ_z

  if (scan_params%n_part == 0) then
    print *, 'ERROR: NO INITIAL POSITIONS SPECIFIED'
    stop
  endif

  ring%param%ixx = 0
  if(fft)ring%param%ixx = 1

  do i_z = 1, n_z

    Q_z = Q_z0 + dQ_z * i_z
    ring%z%tune = Q_z * twopi
    call set_z_tune(ring)

    do i_x = 0, n_x - 1
      Q_x = Q_x0 + dQ_x * i_x
      phi_x = (int_Q_x + Q_x) * twopi

      do i_y = 0, n_y - 1
        Q_y = Q_y0 + dQ_y * i_y
        phi_y = (int_Q_y + Q_y) * twopi

       ring_in = ring
       call beambeam_initialize(ring_in, scan_params, phi_x, phi_y)

      enddo

    enddo

  enddo

  close(unit = 21)
  deallocate(orb)

end program












