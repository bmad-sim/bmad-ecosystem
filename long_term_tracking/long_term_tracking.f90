program long_term_tracking

use beam_mod
use twiss_and_track_mod

implicit none

type (lat_struct) lat
type (beam_init_struct) beam_init
type (bunch_struct) bunch, bunch_init
type (coord_struct), allocatable :: closed_orb(:)

real(rp) time0, time

integer ix_ele_start, n_turns
integer i, n, ie
logical err_flag, rfcavity_on

character(200) lat_file, init_file, end_dat_file, turn_by_turn_dat_file

namelist / long_term_tracking_params / lat_file, n_turns, ix_ele_start, bmad_com, end_dat_file, &
                    turn_by_turn_dat_file, beam_init, rfcavity_on 

! Parse command line

init_file = 'long_term_tracking.init'

if (cesr_iargc() == 1) then
  call cesr_getarg(1, init_file)
elseif (cesr_iargc() > 1) then
  print *, 'Extra stuff on the command line!? Stopping here.'
  stop
endif

! Read parameters

ix_ele_start = 0
rfcavity_on = .true.
end_dat_file = ''
turn_by_turn_dat_file = ''

print '(2a)', 'Initialization file: ', trim(init_file)

open (1, file = init_file, status = 'old')
read (1, nml = long_term_tracking_params)
close (1)

! Lattice init

call bmad_parser (lat_file, lat)
if (.not. rfcavity_on) call set_on_off (rfcavity$, lat, off$)
call twiss_and_track (lat, closed_orb)

! Bunch init

call init_bunch_distribution (lat%ele(ix_ele_start), lat%param, beam_init, 0, bunch, err_flag)
if (err_flag) stop

do n = 1, size(bunch%particle)
  bunch%particle(n)%vec = bunch%particle(n)%vec + closed_orb(ix_ele_start)%vec
enddo

allocate(bunch_init%particle(size(bunch%particle)))
bunch_init%particle = bunch%particle

! Track

call run_timer('START')
call run_timer('READ', time0)

if (turn_by_turn_dat_file /= '') then
  open (1, file = turn_by_turn_dat_file, recl = 200)
  write (1, '(a)') '#  Turn           x           px            y           py            z           pz             spin_x       spin_y       spin_z'
  write (1, '(i9, 6f13.8, 4x, 3f13.8)'), 0, bunch%particle(1)%vec, bunch%particle(1)%spin
endif

do n = 1, n_turns
  do ie = ix_ele_start+1, lat%n_ele_track
    call track1_bunch(bunch, lat, lat%ele(ie), bunch, err_flag)
  enddo
  do ie = 1, ix_ele_start
    call track1_bunch(bunch, lat, lat%ele(ie), bunch, err_flag)
  enddo
  call run_timer('READ', time)

  if (time-time0 > 60) then
    print '(a, f10.2)', 'Ellapsed time (min): ', time/60
    time0 = time
  endif

  if (turn_by_turn_dat_file /= '') then
    write (1, '(i9, 6f13.8, 4x, 3f13.8)'), n, bunch%particle(1)%vec, bunch%particle(1)%spin
  endif
enddo

if (turn_by_turn_dat_file /= '') close (1)

! Output results

if (end_dat_file /= '') then
  open (1, file = end_dat_file, recl = 300)
  write (1, '(2a)') '#                                                    Start                                                                        |', &
                        '                                                 End'
  write (1, '(2a)') '#  Ix          x           px            y           py            z           pz             spin_x       spin_y       spin_z    |', &
                        '           x           px            y           py            z           pz             spin_x       spin_y       spin_z'

  do i = 1, size(bunch%particle)
    write (1, '(i6, 6f13.8, 4x, 3f13.8, 6x, 6f13.8, 4x, 3f13.8)'), i, &
                      bunch_init%particle(i)%vec, bunch_init%particle(i)%spin, &
                      bunch%particle(i)%vec, bunch%particle(i)%spin
  enddo

  close (1)
endif

end program
