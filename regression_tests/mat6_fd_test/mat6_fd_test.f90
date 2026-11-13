!+
! Program mat6_fd_test
!
! Cross-check the transfer matrix (mat6) produced by the analytic/map calculation of each
! element against a matrix constructed by finite-differencing the *same* tracking routine
! used to propagate the orbit.
!
! For a given element:
!   * mat6_analytic -- Computed by make_mat6 using the element's mat6_calc_method. For a
!                      bmad_standard element this is the differentiated-formula matrix that
!                      lives inside the low-level track_a_* routines.
!   * mat6_fd       -- Computed by make_mat6_tracking which perturbs the starting orbit
!                      (step = bmad_com%d_orb) and tracks with the element's tracking_method.
!                      The tracking_method is forced to match the mat6_calc_method so that
!                      the finite-difference matrix exercises the exact same physics routine.
!
! The per-row residual (mat6_analytic - mat6_fd) is written to output.now. A correct matrix
! calculation gives residuals at the level of the finite-difference truncation/roundoff error.
! A mistake in a differentiated formula shows up as a large residual in the corresponding row.
!
! For custom testing:
!   mat6_fd_test <lattice_file>
!-

program mat6_fd_test

use bmad

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele
type (coord_struct) start_orb, end_orb, end_orb_fd

real(rp) mat6_analytic(6,6), mat6_fd(6,6), dmat(6,6)

integer ib, i, k, method, nargs
logical err, err_fd

character(200) :: lat_file = 'mat6_fd_test.bmad'
character(60) :: final_str
character(20) :: method_name

!

global_com%exit_on_error = .false.
bmad_com%auto_bookkeeper = .false.

nargs = command_argument_count()
if (nargs > 0) then
  call get_command_argument(1, lat_file)
  print *, 'Using ', trim(lat_file)
endif

call bmad_parser (lat_file, lat, make_mats6 = .false.)
call lattice_bookkeeper (lat)

open (1, file = 'output.now', recl = 300)

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)

  do i = 1, branch%n_ele_max
    ele => branch%ele(i)
    if (i == branch%n_ele_track .and. ele%name == 'END') cycle

    method = ele%mat6_calc_method

    ! The mat6_calc_method defaults to "auto", which make_mat6 resolves to a concrete method
    ! based on the tracking_method. Resolve it here the same way.

    if (method == auto$) then
      select case (ele%tracking_method)
      case (bmad_standard$, linear$);  method = bmad_standard$
      case (symp_lie_bmad$);           method = symp_lie_bmad$
      case (symp_lie_ptc$);            method = symp_lie_ptc$
      case (taylor$);                  method = taylor$
      case default;                    cycle   ! runge_kutta, custom, mad: nothing to cross-check.
      end select
    endif

    ! Only elements whose matrix is built analytically (or from a map) give a meaningful
    ! cross-check. When the mat6_calc_method is itself finite-difference tracking there is
    ! nothing to compare against.

    select case (method)
    case (bmad_standard$, symp_lie_bmad$, symp_lie_ptc$, taylor$)
    case default
      cycle
    end select

    if (.not. valid_mat6_calc_method (ele, branch%param%particle, method)) cycle
    if (.not. valid_tracking_method (ele, branch%param%particle, method)) cycle

    method_name = mat6_calc_method_name(method)

    ! Force the tracking_method to match the mat6_calc_method so the finite-difference
    ! matrix is built from the same physics routine as the analytic matrix.

    ele%mat6_calc_method = method
    ele%tracking_method  = method

    call init_coord (start_orb, lat%particle_start, ele, upstream_end$, branch%param%particle)

    ! Analytic / map matrix.

    call make_mat6 (ele, branch%param, start_orb, end_orb, err_flag = err)
    if (err .or. end_orb%state /= alive$) then
      write (1, '(a)') '"' // trim(ele%name) // ':' // trim(method_name) // ':Status"  STR  Analytic_Track_Lost'
      cycle
    endif
    mat6_analytic = ele%mat6

    ! Finite-difference matrix through the same tracking routine.

    call make_mat6_tracking (ele, branch%param, start_orb, end_orb_fd, err_fd)
    if (err_fd .or. end_orb_fd%state /= alive$) then
      write (1, '(a)') '"' // trim(ele%name) // ':' // trim(method_name) // ':Status"  STR  FD_Track_Lost'
      cycle
    endif
    mat6_fd = ele%mat6

    dmat = mat6_analytic - mat6_fd

    do k = 1, 6
      final_str = '"' // trim(ele%name) // ':' // trim(method_name) // ':dMat6_Row' // int_str(k) // '"'
      write (1, '(a, a, 6es16.7)') trim(final_str), '  ABS  1E-8', dmat(k,:)
    enddo

    write (1, *)
  enddo
enddo

close (1)

end program
