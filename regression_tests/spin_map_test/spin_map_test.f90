program spin_map_test

use bmad

implicit none

type (lat_struct), target :: lat
type (coord_struct) start_orb, end_orb
type (ele_struct), pointer :: ele
type (branch_struct), pointer :: branch

character(40) :: lat_file  = 'spin_map_test.bmad'
character(100) line

integer, parameter :: n_methods = ubound(tracking_method_name, 1)
integer :: i, j, k, ib, nargs, ns, iq
logical custom_test, err

type spin_q_struct
  character(16) :: method = '?'
  real(rp) q(0:3,0:6)
end type

type (spin_q_struct) spin_q(10)

!

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
    spin_q%method = '?'
    iq = 0

    if (i == branch%n_ele_track .and. ele%name == 'END') cycle

    do j = 1, n_methods
      ele%spin_q(0,0) = real_garbage$
      if (.not. valid_mat6_calc_method(ele, branch%param%particle, j) .or. j == custom$) cycle
      if (ele%key /= taylor$) then
        call kill_taylor(ele%taylor)
        call kill_taylor(ele%spin_taylor)
      endif

      ele%mat6_calc_method = j
      ele%spin_tracking_method = tracking$
      if (j == bmad_standard$) then
        ele%spin_tracking_method = sprint$
        spin_q(iq+1)%method = 'Sprint'
      else
        spin_q(iq+1)%method = 'PTC'
      endif

      call init_coord (start_orb, lat%particle_start, ele, upstream_end$, branch%param%particle)

      call make_mat6 (ele, branch%param, start_orb, end_orb, err)
      if (.not. associated(ele%spin_taylor(1)%term)) cycle

      ele%spin_q = spin_taylor_to_linear(ele%spin_taylor, .true., start_orb%vec-ele%spin_taylor_ref_orb_in, ele%is_on)
      iq = iq + 1
      spin_q(iq)%q = ele%spin_q

    enddo

    if (custom_test) then
      print '(a, 6f13.8)', 'Start orbit:', start_orb%vec
      print '(a, 6f13.8)', 'End orbit:  ', end_orb%vec
    endif

    do k = 0, 3
      do j = 1, iq
        write (1, '(5a, i0, a, t35, f15.10, 4x, 6f15.10)') '"', trim(ele%name), '-', trim(spin_q(j)%method), '-', k, '" ABS 1E-10', spin_q(j)%q(k,:)
      enddo
      write (1, *)
    enddo

  end do  ! ele
enddo   ! branch

close(1)

end program

