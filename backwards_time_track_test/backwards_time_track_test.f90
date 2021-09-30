program backwards_time_track_test

use bmad
use tpsa

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele
type (ele_struct) ele0
type (coord_struct) start_orb, end_orb, start2_orb

character(100) :: lat_file  = 'backwards_time_test.bmad'

real(rp) mat6(6,6), vec0(6), m_unit(6,6)
integer j, ib, ie, nargs

logical debug_mode
 
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


!

do ib = 0, 0
  branch => lat%branch(ib)

  do ie = 1, branch%n_ele_max - 1
    ele => branch%ele(ie)
    ele%spin_tracking_method = tracking$

!    do j = 1, n_methods$

    call init_coord (start_orb, lat%particle_start, ele, downstream_end$)

    bmad_com%backwards_time_tracking_on = .false.
    call make_mat6 (ele, branch%param, start_orb, end_orb)
    ele0 = ele
    bmad_com%backwards_time_tracking_on = .true.
    ele%orientation = -1
    call make_mat6 (ele, branch%param, end_orb, start2_orb)
    mat6 = matmul(ele0%mat6, ele%mat6) - m_unit
    vec0 = matmul(ele%mat6, ele0%vec0) + ele%vec0
    print '(3a, 6f12.8)', '"', trim(ele%name), ':dOrb"', start2_orb%vec - start_orb%vec
    print '(3a, 2f12.8)', '"', trim(ele%name), ':mat6"', maxval(abs(mat6)), maxval(abs(vec0))
    print '(3a, 2f12.8)', '"', trim(ele%name), ':dt"  ', c_light*(start2_orb%t - start_orb%t), (start2_orb%p0c - start_orb%p0c)/1d9


!    enddo
  end do
enddo

end program
