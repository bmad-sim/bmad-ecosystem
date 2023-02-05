program backwards_time_track_test

use bmad
use tpsa

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele
type (ele_struct) ele0
type (coord_struct) start_orb, end_orb, start2_orb, d

character(100) :: lat_file  = 'backwards_time_track_test.bmad'
character(60) str

real(rp) mat6(6,6), vec0(6), m_unit(6,6), beta
integer j, ib, ie, nargs

logical debug_mode, loc_equal
 
!

global_com%exit_on_error = .false.
call mat_make_unit(m_unit)

debug_mode = .false.
nargs = command_argument_count()

if (nargs > 0) then
  call get_command_argument(1, lat_file)
  print *, 'Using ', trim(lat_file)
  debug_mode = .true.
endif

call bmad_parser (lat_file, lat, .false.)

if (any(lat%particle_start%spin /= 0)) then
  bmad_com%spin_tracking_on = .true.
endif

open (1, file = 'output.now')

! To Do:
!   * edge fields
!   * ele%orientation = -1
!   * absolute and rel time tracking
!   * radiation

if (debug_mode) then
  print '(a, t36, 7es18.10)', 'Start:', lat%particle_start%vec
  print '(a, t36, 7es18.10)', 'Spin:', lat%particle_start%spin
  print *
endif

!

do ib = 0, 0
  branch => lat%branch(ib)

  do ie = 1, branch%n_ele_max - 1
    ele => branch%ele(ie)

    do j = 1, n_methods$
      if (.not. valid_tracking_method(ele, branch%param%particle, j)) cycle
      select case (j)
      case (bmad_standard$, runge_kutta$, time_runge_kutta$, linear$, taylor$)
      case default;   cycle
      end select
      ele%tracking_method = j

      str = trim(ele%name) // ': ' // tracking_method_name(j)
      call init_coord (start_orb, lat%particle_start, ele, upstream_end$)

      call track1 (start_orb, ele, branch%param, end_orb)
      end_orb%time_dir = -1
      call track1 (end_orb, ele, lat%param, start2_orb)

      beta = start_orb%beta
      d%vec  =  start2_orb%vec  - start_orb%vec
      d%spin =  start2_orb%spin - start_orb%spin
      d%t    = (start2_orb%t    - start_orb%t) * c_light * beta
      d%s    =  start2_orb%s    - start_orb%s
      d%p0c  =  start2_orb%p0c  - start_orb%p0c
      d%beta =  start2_orb%beta - start_orb%beta
      loc_equal = (start2_orb%location == start_orb%location)

      write (1, '(2a, 7es18.10)')    quote(trim(str) // '-end'), '              ABS 1E-13', end_orb%vec, c_light*beta*end_orb%t
      write (1, '(2a, 7es18.10)')    quote(trim(str) // '-dendSpin'), '         ABS 1E-13', end_orb%spin - start_orb%spin
      write (1, '(2a, 6es18.10)')    quote(trim(str) // '-dOrb'), '             ABS 1E-13', d%vec
      write (1, '(2a, 6es18.10)')    quote(trim(str) // '-dSpin'), '            ABS 1E-13', d%spin
      write (1, '(2a, 4es18.10)')    quote(trim(str) // '-dt,dp0c,ds,dbeta'), ' ABS 1E-13', d%t, d%p0c, d%s, d%beta
      write (1, '(2a, es18.10, l4)') quote(trim(str) // '-Merit'),  '           ABS 1E-13', &
                              maxval([abs(d%vec), abs(d%t), abs(d%s), abs(d%spin), abs(d%beta), abs(d%p0c)]), loc_equal

      if (debug_mode) then
        print '(2a, 7es18.10)',    quote(trim(str) // '-end'), '              ABS 1E-13', end_orb%vec, c_light*beta*end_orb%t
        print '(2a, 7es18.10)',    quote(trim(str) // '-dendSpin'), '         ABS 1E-13', end_orb%spin - start_orb%spin
        print '(a)', '------------------------------------------------------------------------------------'
        print '(2a, 6es18.10)',    quote(trim(str) // '-dOrb'), '             ABS 1E-13', d%vec
        print '(2a, 6es18.10)',    quote(trim(str) // '-dSpin'), '            ABS 1E-13', d%spin
        print '(2a, 4es18.10)',    quote(trim(str) // '-dt,dp0c,ds,dbeta'), ' ABS 1E-13', d%t, d%p0c, d%s, d%beta
        print '(2a, es18.10, l4)', quote(trim(str) // '-Merit'), '            ABS 1E-13', &
                            maxval([abs(d%vec), abs(d%t), abs(d%s), abs(d%spin), abs(d%beta), abs(d%p0c)]), loc_equal
        print *
      endif
    enddo
  end do
enddo

end program
