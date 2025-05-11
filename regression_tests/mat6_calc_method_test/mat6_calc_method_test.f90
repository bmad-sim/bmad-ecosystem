program mat6_calc_method_test

use bmad
use mad_mod
use s_def_kind

implicit none

type (lat_struct), target :: lat
type (ele_struct), target, allocatable :: eles(:)
type (coord_struct) start_orb, end_orb
type (ele_struct), pointer :: ele, ele2
type (lat_param_struct) param
type (branch_struct), pointer :: branch
type (taylor_struct) t_map(6)

real(rp) err_mat(6,6)

character(40) :: lat_file  = 'mat6_calc_method_test.bmad'
character(46) :: final_str
character(20)  :: fmt1 = '(a,a,6es22.13)'
character(20)  :: fmt2 = '(a,a,es22.13)'
character(100) line

integer, parameter :: n_methods = ubound(tracking_method_name, 1)
integer :: i, j, k, ib, nargs, ns, tracking_method
logical custom_test, err, abs_time

!

!old_thick_bend = .true.

global_com%exit_on_error = .false.

custom_test = .false.
nargs = command_argument_count()
if (nargs == 1) then
  call get_command_argument(1, lat_file)
  print *, 'Using ', trim(lat_file)
  custom_test = .true.
  fmt1 = '(a, 6es15.7)' ! Don't need as much precison for test purposes
  fmt2 = '(a, a, es18.9)'
elseif (nargs > 1) then
  print *, 'Only one command line arg permitted.'
  call err_exit
endif

!

bmad_com%auto_bookkeeper = .false.
bmad_com%spin_tracking_on = .true.
call bmad_parser (lat_file, lat, make_mats6 = .false.)
abs_time = bmad_com%absolute_time_tracking

if (custom_test) then
  print '(a, 6f12.6)', 'Init orb: ', lat%particle_start%vec
endif

if (custom_test) then
  if (lat%param%geometry == open$) then
    bmad_com%convert_to_kinetic_momentum = .false.
    print *, '*** Note: wiggler end kicks not cancelled (so like PTC tracking).'
  else
    bmad_com%convert_to_kinetic_momentum = .true.
    print *, '*** Note: wiggler end kicks cancelled (so like RUNGE_KUTTA tracking).'
  endif
endif

call lattice_bookkeeper (lat)

open (1, file = 'output.now', recl = 200)

allocate (eles(n_methods)) 

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  do i = 1, branch%n_ele_max
    ele => branch%ele(i)
    tracking_method = ele%tracking_method

    if (i == branch%n_ele_track .and. ele%name == 'END') cycle
    ns = len_trim(ele%name) + 28

    if (index(ele%name, 'ABS_TIME') /= 0) bmad_com%absolute_time_tracking = .true.

    do j = 1, n_methods
      if (.not. valid_mat6_calc_method(ele, branch%param%particle, j) .or. j == custom$ .or. j == mad$) cycle
      if (j == auto$) cycle
      if (ele%key /= taylor$) call kill_taylor(ele%taylor)

      ele%mat6_calc_method = j
      if (j == bmad_standard$) then
        ele%tracking_method = bmad_standard$
      else
        ele%tracking_method = tracking_method
      endif

      call init_coord (start_orb, lat%particle_start, ele, upstream_end$, branch%param%particle)
      call make_mat6 (ele, branch%param, start_orb, end_orb, err_flag = err)
      call transfer_ele(ele, eles(j), .true.)
      if (custom_test .and. ele%mat6_calc_method == bmad_standard$) then
        print '(a, 6es16.8)', 'Start track:', start_orb%vec
        print '(a, 6es16.8)', 'End track:  ', end_orb%vec 
        print *
      endif
    enddo

    if (index(ele%name, 'ABS_TIME') /= 0) bmad_com%absolute_time_tracking = abs_time

    do k = 1, 8  ! Output line index
      do j = 1, n_methods
        ! if (j == mad$ .and. custom_test) cycle
        if (j == fixed_step_runge_kutta$ .or. j == fixed_step_time_runge_kutta$) cycle
        if (j == auto$) cycle
        if(.not. valid_mat6_calc_method(ele, branch%param%particle, j) .or. j == custom$ .or. j == mad$) cycle
        ele2 => eles(j)
        if (k < 7) then
          if (custom_test) then
            final_str = '"' // trim(ele2%name) // ':' // trim(mat6_calc_method_name(j)) // ':MatrixRow' // int_str(k) // '"' 
            if (err) then
              print '(2a)', final_str(1:ns), '  -------------------------------------- LOST -------------------------------'
            else
              print fmt1, final_str(1:ns), ele2%mat6(k,:)
            endif
          else
            final_str = '"' // trim(ele2%name) // ':' // trim(mat6_calc_method_name(j)) // ':MatrixRow' // int_str(k) // '"' 
            write (1, fmt1) final_str, tolerance(final_str), ele2%mat6(k,:)
          endif
        else if (k == 7) then
          if (custom_test) then
            final_str = '"' // trim(ele2%name) // ':' // trim(mat6_calc_method_name(j)) // ':Vector"' 
            print fmt1, final_str(1:ns), ele2%vec0
          else
            final_str = '"' // trim(ele2%name) // ':' // trim(mat6_calc_method_name(j)) // ':Vector"' 
            write (1, fmt1) final_str, tolerance(final_str), ele2%vec0
          endif
        else if (k == 8) then
          final_str = '"' // trim(ele2%name) // ':' // trim(mat6_calc_method_name(j)) // ':Symp_Err"' 
          if (custom_test) then
            print fmt2, final_str, tolerance(final_str), mat_symp_error(ele2%mat6, ele2%value(p0c$)/ele2%value(p0c_start$), err_mat)
          else
            write (1, fmt2) final_str, tolerance(final_str), mat_symp_error(ele2%mat6, ele2%value(p0c$)/ele2%value(p0c_start$), err_mat)
          endif
        end if
      end do

      if (custom_test .and. k == 8) then
        if (valid_mat6_calc_method(ele, branch%param%particle, bmad_standard$) .and. &
            valid_mat6_calc_method(ele, branch%param%particle, tracking$)) then
          err_mat = abs(eles(bmad_standard$)%mat6 - eles(tracking$)%mat6)
          print '(a, 2i4, es12.2)',   'Max diff |BS - track|:   ', maxloc(err_mat), maxval(err_mat)
        endif

        if (valid_mat6_calc_method(ele, branch%param%particle, bmad_standard$) .and. &
            valid_mat6_calc_method(ele, branch%param%particle, symp_lie_ptc$)) then
          err_mat = abs(eles(bmad_standard$)%mat6 - eles(symp_lie_ptc$)%mat6)
          print '(a, 2i4, es12.2)',   'Max diff |BS - PTC|:     ', maxloc(err_mat), maxval(err_mat)
        endif

        if (valid_mat6_calc_method(ele, branch%param%particle, symp_lie_bmad$) .and. &
            valid_mat6_calc_method(ele, branch%param%particle, symp_lie_ptc$)) then
          err_mat = abs(eles(symp_lie_bmad$)%mat6 - eles(symp_lie_ptc$)%mat6)
          print '(a, 2i4, es12.2)',   'Max diff |SLBmad - PTC|:     ', maxloc(err_mat), maxval(err_mat)
        endif

      endif

      if (custom_test) then
        print *
      else
        write (1,*)
      endif
    end do  ! k

  end do  ! ele
enddo   ! branch

close(1)

!----------------------------------------------------------------------
contains

character(10) function tolerance(instr)
character(44) :: instr

select case (instr)

case ('"E_GUN1:Tracking:MatrixRow1"')              ; tolerance = 'ABS 5e-10'
case ('"E_GUN1:Tracking:MatrixRow3"')              ; tolerance = 'ABS 7e-10'
case ('"E_GUN1:Tracking:MatrixRow5"')              ; tolerance = 'ABS 4e-07'
case ('"E_GUN1:Tracking:MatrixRow6"')              ; tolerance = 'ABS 6e-08'
case ('"E_GUN1:Tracking:Symp_Err"')                ; tolerance = 'ABS 4e-07'
case ('"E_GUN1:Tracking:Vector"')                  ; tolerance = 'ABS 6e-10'

case ('"ELSEPARATOR1:Tracking:MatrixRow5"')        ; tolerance = 'ABS 2e-10'

case ('"LCAVITY1:Tracking:MatrixRow5"')            ; tolerance = 'ABS 2e-10'
case ('"LCAVITY1:Tracking:MatrixRow6"')            ; tolerance = 'ABS 6e-10'
case ('"LCAVITY1:Tracking:Symp_Err"')              ; tolerance = 'ABS 3e-9'
case ('"LCAVITY2:Tracking:MatrixRow5"')            ; tolerance = 'ABS 2e-10'
case ('"LCAVITY2:Tracking:Symp_Err"')              ; tolerance = 'ABS 1e-7'
case ('"LCAVITY3:Tracking:MatrixRow5"')            ; tolerance = 'ABS 2e-10'
case ('"LCAVITY3:Tracking:MatrixRow6"')            ; tolerance = 'ABS 4e-10'
case ('"LCAVITY3:Tracking:Symp_Err"')              ; tolerance = 'ABS 1e-9'


case ('"LCAVITY1_ABS_TIME:Tracking:MatrixRow5"')   ; tolerance = 'ABS 2e-10'
case ('"LCAVITY1_ABS_TIME:Tracking:MatrixRow6"')   ; tolerance = 'ABS 6e-10'
case ('"LCAVITY1_ABS_TIME:Tracking:Symp_Err"')     ; tolerance = 'ABS 3e-9'
case ('"LCAVITY2_ABS_TIME:Tracking:MatrixRow5"')   ; tolerance = 'ABS 2e-10'
case ('"LCAVITY3_ABS_TIME:Tracking:MatrixRow5"')   ; tolerance = 'ABS 2e-10'
case ('"LCAVITY3_ABS_TIME:Tracking:MatrixRow6"')   ; tolerance = 'ABS 5e-10'
case ('"LCAVITY3_ABS_TIME:Tracking:Symp_Err"')     ; tolerance = 'ABS 3e-9'

case ('"RFCAVITY2:Tracking:MatrixRow5"')           ; tolerance = 'ABS 2e-10'
case ('"RFCAVITY2:Tracking:MatrixRow6"')           ; tolerance = 'ABS 1e-09'

case ('"RBEND4:Symp_Lie_PTC:MatrixRow1"')          ; tolerance = 'ABS 4e-09'
case ('"RBEND4:Taylor:MatrixRow1"')                ; tolerance = 'ABS 4e-09'
case ('"RBEND4:Symp_Lie_PTC:MatrixRow3"')          ; tolerance = 'ABS 4e-08'
case ('"RBEND4:Taylor:MatrixRow3"')                ; tolerance = 'ABS 4e-08'
case ('"RBEND4:Symp_Lie_PTC:MatrixRow5"')          ; tolerance = 'ABS 3e-07'
case ('"RBEND4:Taylor:MatrixRow5"')                ; tolerance = 'ABS 3e-07'
case ('"RBEND4:Symp_Lie_PTC:Vector"')              ; tolerance = 'ABS 4e-08'
case ('"RBEND4:Taylor:Vector"')                    ; tolerance = 'ABS 4e-08'
case ('"RBEND4:Symp_Lie_PTC:Symp_Err"')            ; tolerance = 'ABS 4e-07'
case ('"RBEND4:Taylor:Symp_Err"')                  ; tolerance = 'ABS 4e-07'
case ('"RBEND4:Bmad_Standard:MatrixRow3"')         ; tolerance = 'ABS 2e-10'
case ('"RBEND4:Bmad_Standard:MatrixRow5"')         ; tolerance = 'ABS 1e-07'
case ('"RBEND4:Bmad_Standard:Vector"')             ; tolerance = 'ABS 4e-10'
case ('"RBEND4:Tracking:MatrixRow1"')              ; tolerance = 'ABS 2e-08'
case ('"RBEND4:Tracking:MatrixRow2"')              ; tolerance = 'ABS 2e-10'
case ('"RBEND4:Tracking:MatrixRow3"')              ; tolerance = 'ABS 1e-07'
case ('"RBEND4:Tracking:MatrixRow5"')              ; tolerance = 'ABS 1e-05'
case ('"RBEND4:Tracking:Vector"')                  ; tolerance = 'ABS 1e-07'
case ('"RBEND4:Tracking:Symp_Err"')                ; tolerance = 'ABS 2e-06'
case ('"SBEND5:Tracking:MatrixRow1"')              ; tolerance = 'ABS 8e-09'
case ('"SBEND5:Tracking:MatrixRow3"')              ; tolerance = 'ABS 8e-10'
case ('"SBEND5:Tracking:MatrixRow4"')              ; tolerance = 'ABS 5e-10'
case ('"SBEND5:Tracking:MatrixRow5"')              ; tolerance = 'ABS 1e-09'
case ('"SBEND5:Tracking:Symp_Err"')                ; tolerance = 'ABS 4e-09'
case ('"RBEND6:Symp_Lie_PTC:MatrixRow1"')          ; tolerance = 'ABS 3e-09'
case ('"RBEND6:Taylor:MatrixRow1"')                ; tolerance = 'ABS 3e-09'
case ('"RBEND6:Symp_Lie_PTC:MatrixRow3"')          ; tolerance = 'ABS 3e-08'
case ('"RBEND6:Taylor:MatrixRow3"')                ; tolerance = 'ABS 3e-08'
case ('"RBEND6:Symp_Lie_PTC:MatrixRow5"')          ; tolerance = 'ABS 2e-07'
case ('"RBEND6:Taylor:MatrixRow5"')                ; tolerance = 'ABS 2e-07'
case ('"RBEND6:Symp_Lie_PTC:Vector"')              ; tolerance = 'ABS 2e-07'
case ('"RBEND6:Taylor:Vector"')                    ; tolerance = 'ABS 4e-08'
case ('"RBEND6:Symp_Lie_PTC:Symp_Err"')            ; tolerance = 'ABS 2e-07'
case ('"RBEND6:Taylor:Symp_Err"')                  ; tolerance = 'ABS 2e-07'
case ('"SBEND7:Tracking:MatrixRow1"')              ; tolerance = 'ABS 1e-08'
case ('"SBEND7:Tracking:MatrixRow3"')              ; tolerance = 'ABS 1e-09'
case ('"SBEND7:Tracking:MatrixRow4"')              ; tolerance = 'ABS 4e-10'
case ('"SBEND7:Tracking:Symp_Err"')                ; tolerance = 'ABS 5e-09'

case ('"SOL_QUAD2:Tracking:MatrixRow1"')           ; tolerance = 'ABS 8e-10'
case ('"SOL_QUAD2:Tracking:MatrixRow4"')           ; tolerance = 'ABS 2e-10'

case default                                       ; tolerance = 'ABS 1E-10'
end select

end function tolerance

end program
