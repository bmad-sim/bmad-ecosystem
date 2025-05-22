program backwards_time_track_test

use bmad
use tpsa

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele, ele2
type (ele_struct) ele0
type (coord_struct) start_orb, end_orb, start2_orb, d

character(100) :: lat_file  = 'backwards_time_track_test.bmad', slat_file = 's_to_s.bmad'
character(60) str

real(rp) mat6(6,6), vec0(6), m_unit(6,6), beta, merit, global_merit, ele_merit, s1, s2
integer j, ib, ie, nargs
integer, parameter :: n_methods = ubound(tracking_method_name, 1)

logical debug_mode, loc_equal, global_loc_equal, ele_loc_equal
 
!

global_com%exit_on_error = .false.
call mat_make_unit(m_unit)

debug_mode = .false.
nargs = command_argument_count()

if (nargs > 0) then
  call get_command_argument(1, lat_file)
  print *, 'Using ', trim(lat_file)
  slat_file = lat_file
  debug_mode = .true.
endif

open (1, file = 'output.now')

! To Do:
!   * edge fields
!   * ele%orientation = -1
!   * absolute and rel time tracking
!   * radiation
!   * track_from_s_to_s

if (debug_mode) then
  print '(a, t36, 7es18.10)', 'Start:', lat%particle_start%vec
  print '(a, t36, 7es18.10)', 'Spin:', lat%particle_start%spin
  print *
endif

!

call bmad_parser (slat_file, lat, .false.)
do ie = 1, lat%n_ele_track-1, 2
  ele => lat%ele(ie); ele2 => lat%ele(ie+1)
  s1 = 0.5_rp * (ele%s_start + ele%s) + 0.1_rp
  s2 = 0.5_rp * (ele2%s_start + ele2%s) + 0.1_rp
  call init_coord (start_orb, lat%particle_start, ele, inside$, s_pos = s1)
  call track_from_s_to_s (lat, s1, s2, start_orb, end_orb)
  end_orb%time_dir = -1
  call track_from_s_to_s (lat, s2, s1, end_orb, start2_orb)

  beta = start_orb%beta
  d%vec  =  start2_orb%vec  - start_orb%vec
  d%spin =  start2_orb%spin - start_orb%spin
  d%t    = (start2_orb%t    - start_orb%t) * c_light * beta
  d%s    =  start2_orb%s    - start_orb%s
  d%p0c  =  start2_orb%p0c  - start_orb%p0c
  d%beta =  start2_orb%beta - start_orb%beta
  loc_equal = (start2_orb%location == start_orb%location)
  merit = maxval([real(rp):: abs(d%vec), abs(d%t), abs(d%s), abs(d%spin), abs(d%beta), abs(d%p0c)])

  str = trim(key_name(ele%key)) // ':' // trim(key_name(ele2%key))

  write (1, '(2a, 7es18.10)')    quote(trim(str) // '-end'), '                ABS 1e-9', end_orb%vec, 1e-3*c_light*beta*end_orb%t
  write (1, '(2a, 7es18.10)')    quote(trim(str) // '-dendSpin'), '           ABS 1e-9', end_orb%spin - start_orb%spin
  write (1, '(2a, 6es18.10)')    quote(trim(str) // '-dOrb'), '               ABS 1e-9', d%vec
  write (1, '(2a, 6es18.10)')    quote(trim(str) // '-dSpin'), '              ABS 1e-9', d%spin
  write (1, '(2a, 4es18.10)')    quote(trim(str) // '-c*dt,dp0c,ds,dbeta'), ' ABS 1e-9', d%t, d%p0c, d%s, d%beta
  write (1, '(2a, es18.10, l4)') quote(trim(str) // '-Merit'),  '             ABS 1e-9', merit

  if (debug_mode) then
    print *
    print '(2a, 7es18.10)',    quote(trim(str) // '-end'), '                ABS 1e-9', end_orb%vec, 1e-3*c_light*beta*end_orb%t
    print '(2a, 7es18.10)',    quote(trim(str) // '-dendSpin'), '           ABS 1e-9', end_orb%spin - start_orb%spin
    print '(a)', '------------------------------------------------------------------------------------'
    print '(2a, 6es18.10)',    quote(trim(str) // '-dOrb'), '               ABS 1e-9', d%vec
    print '(2a, 6es18.10)',    quote(trim(str) // '-dSpin'), '              ABS 1e-9', d%spin
    print '(2a, 4es18.10)',    quote(trim(str) // '-c*dt,dp0c,ds,dbeta'), ' ABS 1e-9', d%t, d%p0c, d%s, d%beta
    print '(2a, es18.10, l4)', quote(trim(str) // '-Merit'), '              ABS 1e-9', merit, loc_equal
  endif
enddo

!

call bmad_parser (lat_file, lat, .false.)

if (any(lat%particle_start%spin /= 0)) then
  bmad_com%spin_tracking_on = .true.
endif

global_loc_equal = .true.
global_merit = 0

do ib = 0, ubound(lat%branch,1)
  branch => lat%branch(ib)

  do ie = 1, branch%n_ele_max - 1
    ele => branch%ele(ie)
    ele_loc_equal = .true.
    ele_merit = 0


    do j = 1, n_methods
      if (.not. valid_tracking_method(ele, branch%param%particle, j)) cycle
      select case (j)
      case (bmad_standard$, runge_kutta$, time_runge_kutta$, linear$, taylor$, symp_lie_bmad$)
      case default;   cycle
      end select
      ele%tracking_method = j

      str = trim(ele%name) // ': ' // tracking_method_name(j)
      call init_coord (start_orb, lat%particle_start, ele, upstream_end$)

      call track1 (start_orb, ele, branch%param, end_orb)
      if (end_orb%state /= alive$) then
        print *, '!!!! ', trim(str), '  Forward tracking particle lost'
        cycle
      endif

      end_orb%time_dir = -1
      call track1 (end_orb, ele, lat%param, start2_orb)
      if (start2_orb%state /= alive$) then
        print *, '!!!! ', trim(str), '  Backwards tracking particle lost'
        cycle
      endif

      beta = start_orb%beta
      d%vec  =  start2_orb%vec  - start_orb%vec
      d%spin =  start2_orb%spin - start_orb%spin
      d%t    = (start2_orb%t    - start_orb%t) * c_light * beta
      d%s    =  start2_orb%s    - start_orb%s
      d%p0c  =  start2_orb%p0c  - start_orb%p0c
      d%beta =  start2_orb%beta - start_orb%beta
      loc_equal = (start2_orb%location == start_orb%location)
      merit = maxval([real(rp):: abs(d%vec), abs(d%t), abs(d%s), abs(d%spin), abs(d%beta), abs(d%p0c)])

      write (1, '(2a, 7es18.10)')    quote(trim(str) // '-end'), '                ABS 1e-9', end_orb%vec, 1d-3*c_light*beta*end_orb%t
      write (1, '(2a, 7es18.10)')    quote(trim(str) // '-dendSpin'), '           ABS 1e-9', end_orb%spin - start_orb%spin
      write (1, '(2a, 6es18.10)')    quote(trim(str) // '-dOrb'), '               ABS 1e-9', d%vec
      write (1, '(2a, 6es18.10)')    quote(trim(str) // '-dSpin'), '              ABS 1e-9', d%spin
      write (1, '(2a, 4es18.10)')    quote(trim(str) // '-c*dt,dp0c,ds,dbeta'), ' ABS 1e-9', d%t, d%p0c, d%s, d%beta
      write (1, '(2a, es18.10, l4)') quote(trim(str) // '-Merit'),  '             ABS 1e-9', merit

      if (debug_mode) then
        print *
        print '(2a, 7es18.10)',    quote(trim(str) // '-end'), '                ABS 1e-9', end_orb%vec, 1d-3*c_light*beta*end_orb%t
        print '(2a, 7es18.10)',    quote(trim(str) // '-dendSpin'), '           ABS 1e-9', end_orb%spin - start_orb%spin
        print '(a)', '------------------------------------------------------------------------------------'
        print '(2a, 6es18.10)',    quote(trim(str) // '-dOrb'), '               ABS 1e-9', d%vec
        print '(2a, 6es18.10)',    quote(trim(str) // '-dSpin'), '              ABS 1e-9', d%spin
        print '(2a, 4es18.10)',    quote(trim(str) // '-c*dt,dp0c,ds,dbeta'), ' ABS 1e-9', d%t, d%p0c, d%s, d%beta
        print '(2a, es18.10, 2l4)', quote(trim(str) // '-Merit'), '             ABS 1e-9 ', merit, loc_equal
        ele_merit = max(ele_merit, merit)
        ele_loc_equal = (ele_loc_equal .and. loc_equal)
      endif
    enddo

    if (debug_mode .and. branch%n_ele_max > 2) then
      print '(a, es18.10, l4)', '*************** "Ele-Merit"  ', ele_merit, ele_loc_equal
      print *
    endif

    global_merit = max(ele_merit, global_merit)
    global_loc_equal = (ele_loc_equal .and. global_loc_equal)

  end do
  if (debug_mode .and. ie > 1) then
    print '(a, es18.10, l4)', '**************** "Global-Merit"  ', global_merit, global_loc_equal
  endif
enddo

end program
