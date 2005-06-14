!+
! Program       : SCAN_DRIVER
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
! Revision 1.1  2005/06/14 14:59:02  cesrulib
! Initial revision
!
!
!........................................................................
!
#include "CESR_platform.h"

program scan_driver

  use bmad_struct
  use bmad_interface
  use bmadz_interface
  use bsim_interface
  use scan_parameters
 
  implicit none

  type (ring_struct) ring, ring_ok, ring_in
  type (coord_struct), allocatable :: orb_(:)
  type (coord_struct) init(100)
  type (scan_params_struct) scan_params
  

  integer i, j, ie
  integer n_x, n_y, i_x, i_y, i_ele, n_z
  integer i_z
  integer int_Q_x, int_Q_y, ix, ixx,ixxx
  integer kx
  integer ios
  integer version
  integer n_part, n_turn, particle, i_train, j_car, n_trains_tot, n_cars, slices 

  real(rdef), allocatable :: dk1(:) 
  real(rdef) phi_x, phi_y
  real(rdef) a(6,6), wr(6), wi(6), tune(6)
  real(rdef) test
  real(rdef)  Q_x0, Q_x1, dQ_x, Q_y0, Q_y1, dQ_y
  real(rdef) Q_z, current
  real(rdef) Q_x, Q_y, Q_z0, Q_z1, dQ_z
  real(rdef) min_sig
  real(rdef) Q_x_init, Q_y_init
  real(rdef) coupling_sb, coupling_wb

  character lat_file*80, file_name*80, out_file*80, new_lat_file*80
  character date*10, time*10
  character*2 wordx, wordy

  logical ok, found
  logical lrbbi, beambeam_ip, close_pretz, close_vert, rec_taylor, radiation

  namelist / parameters / Q_x0, Q_x1, dQ_x, n_x, Q_y0, Q_y1, dQ_y, n_y, &
                Q_z, Q_z0, Q_z1, dQ_z, n_part, lat_file, n_turn, particle, &
                i_train, j_car, n_trains_tot, n_cars, current, &
                lrbbi, beambeam_ip, close_pretz,  &
                close_vert, slices, rec_taylor, n_part, min_sig, Q_x_init, Q_y_init, &
                coupling_sb, coupling_wb, radiation, init


! read in the parameters


  do 

    print '(a, $)', ' Input command file <CR=tunescan.in>: '
    read(*, '(a)') file_name
    call string_trim (file_name, file_name, ix)
    if (ix .eq. 0) file_name = 'tunescan.in'

    open (unit= 1, file = file_name, status = 'old', iostat = ios)

    if (ios == 0) then
      exit
    else
      type *
      type *, 'ERROR: CANNOT OPEN FILE: ', trim(file_name)
    endif

  enddo

  particle = positron$
  beambeam_ip = .false.
  close_pretz = .false.
  close_vert = .false.
  rec_taylor = .false.
  slices = 1
  min_sig = 0.
  coupling_sb = 0.01
  n_x =1
  n_y=1
  n_z = 1
  radiation = .false.

  read(1, nml = parameters)

  scan_params%Q_z    =        Q_z 
  scan_params%n_part =        n_part 
  scan_params%lat_file =      lat_file 
  scan_params%n_turn =        n_turn 
  scan_params%particle =      particle
  scan_params%i_train =       i_train 
  scan_params%j_car =         j_car 
  scan_params%n_trains_tot =  n_trains_tot 
  scan_params%n_cars =        n_cars 
  scan_params%current =       current
  scan_params%lrbbi =         lrbbi
  scan_params%beambeam_ip =   beambeam_ip 
  scan_params%close_pretz =   close_pretz
  scan_params%close_vert =    close_vert 
  scan_params%slices =        slices 
  scan_params%rec_taylor =    rec_taylor
  scan_params%min_sig    =    min_sig
  scan_params%coupling_sb =   coupling_sb
  scan_params%coupling_wb =   coupling_wb
  scan_params%radiation  =    radiation
  scan_params%init(1:n_part)  =    init(1:n_part)
  scan_params%particle    = particle

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

  if(lat_file(1:8) == 'digested')then
    call read_digested_bmad_file (lat_file, ring, version)
   else
    call bmad_parser (lat_file, ring)
  endif



! find how many starting points there are

  type *, 'n_x, n_y:', n_x, n_y
  type *, 'dQ_x, dQ_y:', dQ_x, dQ_y


  ix = index(file_name,'.')
  out_file = file_name(1:ix)//'hed'

  ixx = max(index(lat_file,'/'),1)
  ixxx = len_trim(lat_file)
  kx = ixxx-ixx+1
  new_lat_file = lat_file(ixx:ixxx)
  new_lat_file(kx+1:kx+1) = "'"

  call date_and_time(date, time)


  open(unit=21, file = out_file)
  write (21, '(a13,a80)') 'Lattice   = `', new_lat_file
  write (21, *) 'File_name = `', trim(file_name)// "'"
  write (21, '(a13,1x,a4,a1,a2,a1,a2,1x,a2,a1,a2,a1,a2,a1)') 'data_date = `', &
                   date(1:4),'-',date(5:6),'-',date(7:8),time(1:2),':',time(3:4),':',time(5:6), "'"
  write (21,'(a11,i5)' )'  n_part = ',n_part
  write (21,'(a11, i7)')'  n_turn = ',n_turn
  write (21, *) 'Return'
  close( unit=21)

  out_file = file_name(1:ix)//'out'
  open(unit=21, file = out_file, status='replace')
  write (21, '(a13,a80)') 'Lattice   = `', new_lat_file
  write (21, *) 'File_name = `', trim(file_name)// "'"
  write (21, '(a13,1x,a4,a1,a2,a1,a2,1x,a2,a1,a2,a1,a2,a1)') 'data_date = `', &
                   date(1:4),'-',date(5:6),'-',date(7:8),time(1:2),':',time(3:4),':',time(5:6), "'"
  write (21, *) 'Return'
  write (21, *)
  write(21,1) Q_x0, Q_y0, Q_x1, Q_y1, ring%z%tune/twopi
1    format(1x,'!','  Q_x0 = ',f5.3, '  Q_y0 = ',f5.3, '  Q_x1 = ',f5.3, &
                           '  Q_y1 = ', f5.3,' Q_z = ',f5.3)
2  format(1x,'! Traj',i3,' start ',2x,6e10.2)

  write(21,*)
  write(21,'(a78,a10)')'!   Q_x       Q_y       Q_z      x_end       y_end      amp_x_max   amp_y_max ',&
                   '  Turns   '

  scan_params%file_name = out_file

  close(unit=21)

  do i = 1, ring.n_ele_max
    if(ring.ele_(i).value(x_limit$) == 0.)ring.ele_(i).value(x_limit$) = 0.05
    if(ring.ele_(i).value(y_limit$) == 0.)ring.ele_(i).value(y_limit$) = 0.05
  enddo

  ring.param.aperture_limit_on = .true.

  allocate(dk1(ring%n_ele_max))
  allocate(orb_(0:ring%n_ele_max))

  if(Q_x_init == 0)Q_x_init = Q_x0
  if(Q_y_init == 0)Q_y_init = Q_y0

! find which quads to change

  call twiss_at_start(ring)
  orb_(0)%vec=0
  call closed_orbit_at_start(ring, orb_(0), 4, .true.)
  call track_all(ring, orb_)
  call ring_make_mat6(ring, -1, orb_)
  call twiss_at_start(ring)
  call twiss_propagate_all(ring)
  call choose_quads(ring, dk1)
  int_Q_x = int(ring%ele_(ring%n_ele_ring)%x%phi / twopi)
  int_Q_y = int(ring%ele_(ring%n_ele_ring)%y%phi / twopi)
  phi_x = (int_Q_x + Q_x_init) * twopi
  phi_y = (int_Q_y + Q_y_init) * twopi

  call custom_set_tune(phi_x, phi_y, dk1, ring, orb_, ok)    

! scan the tune
  ring_ok = ring

! find how many starting points there are

!  do n_part = 0, size(init)
!    if (all(init(size(init))%vec == 0)) exit
!  enddo

  type *, 'Number of starting positions:', n_part
  type *, 'n_x, n_y, n_z:', n_x, n_y, n_z
  type *, 'dQ_x, dQ_y, dQ_z:', dQ_x, dQ_y, dQ_z

  if (n_part == 0) then
    type *, 'ERROR: NO INITIAL POSITIONS SPECIFIED'
    stop
  endif

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


!       write(wordx,'(i2.2)')i_x
!       write(wordy,'(i2.2)')i_y
!       scan_params%file_name = file_name(1:ix-1)//'_'//wordx//'_'//wordy
!       type *,' scan_params%file_name ',scan_params%file_name

       ring_in = ring
       call beambeam_initialize(ring_in, scan_params, phi_x, phi_y)

!
!       open(unit=21, file = out_file,status='old', access='append' )
!       write(21,'(2f10.5,2x,2f10.5,3e12.4,2x,3e12.4,2x,i5,2x,e12.4)')Q_x,Q_y, &
!            ring_in%x%tune/twopi, ring_in%y%tune/twopi, &
!                     scan_params%sig_in(1:3), &
!                                    scan_params%sig_out(1:3), scan_params%n_part_out, scan_params%lum
!       close(unit=21)
!       print *,Q_x, Q_y, scan_params%lum
      enddo

    enddo

  end do

  close(unit = 21)
  deallocate(orb_)

end program












