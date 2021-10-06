program spin_method_test

use bmad
use mad_mod
use s_def_kind

implicit none

type (lat_struct), target :: lat
type (coord_struct) start_orb, end_orb
type (ele_struct), pointer :: ele
type (branch_struct), pointer :: branch

character(40) :: lat_file  = 'spin_method_test.bmad'
character(100) line

integer :: i, j, k, ib, nargs, ns
logical custom_test, err

!

!old_thick_bend = .true.

global_com%exit_on_error = .false.

custom_test = .false.
nargs = command_argument_count()
if (nargs == 1) then
  call get_command_argument(1, lat_file)
  print *, 'Using ', trim(lat_file)
  custom_test = .true.
elseif (nargs > 1) then
  print *, 'Only one command line arg permitted.'
  call err_exit
endif

!

bmad_com%spin_tracking_on = .true.
call bmad_parser (lat_file, lat, make_mats6 = .false.)

call lattice_bookkeeper (lat)

open (1, file = 'output.now', recl = 200)

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  do i = 1, branch%n_ele_max
    ele => branch%ele(i)
    if (i == branch%n_ele_track .and. ele%name == 'END') cycle

    do j = 1, n_methods$
      if (.not. valid_mat6_calc_method(ele, branch%param%particle, j) .or. j == static$ .or. j == custom$ .or. j == mad$) cycle
      if (ele%key /= taylor$) call kill_taylor(ele%taylor)
      ele%mat6_calc_method = j
      call init_coord (start_orb, lat%particle_start, ele, upstream_end$, branch%param%particle)
      call make_mat6 (ele, branch%param, start_orb, end_orb, err_flag = err)
      if (custom_test .and. ele%mat6_calc_method == bmad_standard$) then
        write (1, '(a, 6es16.8)') 'Start track:', start_orb%vec
        write (1, '(a, 6es16.8)') 'End track:  ', end_orb%vec 
        write (1, *)
      endif
    enddo


  end do  ! ele
enddo   ! branch

close(1)

end program

