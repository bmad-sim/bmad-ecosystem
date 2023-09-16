program test

use bmad

implicit none

type (lat_struct), target :: lat
type (ele_pointer_struct), allocatable :: eles(:)
type (ele_struct), pointer :: ele

character(40) mad_name, bmad_name
character(200) lat_file, mad_err_file, bmad_err_file
character(5000) str

real(rp) :: dx = 0, dy = 0, ds = 0, dphi = 0, dtheta = 0, dpsi = 0
real(rp) :: mrex = 0, mrey = 0, mscalx = 0, mscaly = 0, arex = 0, arey = 0

integer i, ie, ix, n_col, col_max, n_loc
integer :: ix_dx = 0, ix_dy = 0, ix_ds = 0, ix_dphi = 0, ix_dtheta = 0, ix_dpsi = 0, n_inst_err = 0
integer :: ix_mrex = 0, ix_mrey = 0, ix_mscalx = 0, ix_mscaly = 0, ix_arex = 0, ix_arey = 0

logical err

! Get file names.

call get_command_argument(1, lat_file)
call get_command_argument(2, mad_err_file)
bmad_err_file = trim(mad_err_file) // '.bmad'

if (lat_file == '' .or. mad_err_file == '') then
  print *, 'Syntax: mad_errors_to_bmad <bmad-lattice-file-name> <error-file-name>'
  stop
endif

print *, 'Lattice file:    ', trim(lat_file)
print *, 'MAD Error file:  ', trim(mad_err_file)
print *, 'Bmad Error file: ', trim(bmad_err_file)

call bmad_parser(lat_file, lat)

! Parse error header.

open (1, file = mad_err_file, status = 'old', action = 'read')
open (2, file = bmad_err_file)

do
  read(1, '(a)') str
  if (str(1:1) == '$') exit

  if (str(1:1) == '*') then
    call string_trim(str(2:), str, ix)

    do n_col = 1, 1000
      select case (str(:ix))
      case('DX');      ix_dx     = n_col
      case('DY');      ix_dy     = n_col
      case('DS');      ix_ds     = n_col
      case('DPHI');    ix_dphi   = n_col
      case('DTHETA');  ix_dtheta = n_col
      case('DPSI');    ix_dpsi   = n_col
      case('MREX');    ix_mrex   = n_col
      case('MREY');    ix_mrey   = n_col
      case('MSCALX');  ix_mscalx = n_col
      case('MSCALY');  ix_mscaly = n_col
      case('AREX');    ix_arex   = n_col
      case('AREY');    ix_arey   = n_col
      end select

      call string_trim(str(ix+1:), str, ix)
      if (ix == 0) exit
    enddo
  endif
enddo

! Parse error table.

col_max = max(ix_dx, ix_dy, ix_ds, ix_dphi, ix_dtheta, ix_dpsi, &
              ix_mrex, ix_mrey, ix_mscalx, ix_mscaly, ix_arex, ix_arey)


do ie = 1, 1000000
  dx = 0; dy = 0; ds = 0; dphi = 0; dtheta = 0; dpsi = 0
  mrex = 0; mrey = 0; mscalx = 0; mscaly = 0; arex = 0; arey = 0

  read(1, '(a)', end = 9000) str
  if (str == '') exit
  call string_trim(str, str, ix)
  mad_name = unquote(str(:ix))

  do n_col = 2, col_max
    call string_trim(str(ix+1:), str, ix)
    if (n_col == ix_dx)      dx     = real_val(str(:ix))
    if (n_col == ix_dy)      dy     = real_val(str(:ix))
    if (n_col == ix_ds)      ds     = real_val(str(:ix))
    if (n_col == ix_dphi)    dphi   = real_val(str(:ix))
    if (n_col == ix_dtheta)  dtheta = real_val(str(:ix))
    if (n_col == ix_dpsi)    dpsi   = real_val(str(:ix))
    if (n_col == ix_mrex)    mrex   = real_val(str(:ix))
    if (n_col == ix_mrey)    mrey   = real_val(str(:ix))
    if (n_col == ix_mscalx)  mscalx = real_val(str(:ix))
    if (n_col == ix_mscaly)  mscaly = real_val(str(:ix))
    if (n_col == ix_arex)    arex   = real_val(str(:ix))
    if (n_col == ix_arey)    arey   = real_val(str(:ix))
  enddo

  if (index(mad_name, '$') /= 0) cycle

  call lat_ele_locator(mad_name, lat, eles, n_loc, err)

  if (err) then
    print *, 'ERROR IN ELEMENT LOCATION CODE FOR: ' // trim(mad_name)
    cycle
  endif

  if (n_loc == 0) then
    print *, 'ERROR LOCATING ELEMENT: ' // trim(mad_name)
    cycle
  endif

  ele => eles(1)%ele
  ele%ix_pointer = ele%ix_pointer+1

  if (n_loc == 1) then
    bmad_name = mad_name
  else
    bmad_name = trim(mad_name) // '##' // int_str(ele%ix_pointer)
  endif

  if (mrex /= 0 .or. mrey /= 0 .or. mscalx /= 0 .or. mscaly /= 0) then
    if (.not. has_attribute(ele, 'X_GAIN_CALIB')) then
      n_inst_err = n_inst_err + 1
      if (n_inst_err < 6) then
        print *, 'Note: Nonzero measurement error (MREX, MREY, MSCALX, OR MSALY) but Bmad element does not have instrumental measurement attributes: ' // mad_name
      elseif (n_inst_err == 6) then
        print *, 'Will not print instrumentation messages for other elements.'
      endif

    else
      call value_out (bmad_name, '[x_offset]', -mrex)
      call value_out (bmad_name, '[y_offset]', -mrey)
      call value_out (bmad_name, '[x_gain_err]', mscalx)
      call value_out (bmad_name, '[y_gain_err]', mscaly)
    endif
  endif

  if (arex /= 0) then
    if (ele%value(x1_limit$) == 0 .or. ele%value(x2_limit$) == 0) then
      print *, 'ERROR: "AREX" SET BUT APERTURE NOT SET FOR ELEMENT: ' // trim(mad_name)
    else
      call value_out (bmad_name, '[x1_limit]', ele%value(x1_limit$) - arex)
      call value_out (bmad_name, '[x2_limit]', ele%value(x2_limit$) + arex)
    endif
  endif

  if (arey /= 0) then
    if (ele%value(y1_limit$) == 0 .or. ele%value(y2_limit$) == 0) then
      print *, 'ERROR: "AREY" SET BUT APERTURE NOT SET FOR ELEMENT: ' // trim(mad_name)
    else
      call value_out (bmad_name, '[y1_limit]', ele%value(y1_limit$) - arey)
      call value_out (bmad_name, '[y2_limit]', ele%value(y2_limit$) + arey)
    endif
  endif

  dx = dx - 0.5_rp * dtheta * ele%value(l$)
  dy = dy - 0.5_rp * dphi   * ele%value(l$)

  call value_out (bmad_name, 'x_offset', dx)
  call value_out (bmad_name, 'y_offset', dy)
  call value_out (bmad_name, 'z_offset', ds)
  call value_out (bmad_name, 'x_pitch', dtheta)
  call value_out (bmad_name, 'y_pitch', dphi)
  if (ele%key == sbend$) then
    call value_out (bmad_name, 'roll', dpsi)
  else
    call value_out (bmad_name, 'tilt', dpsi)
  endif
enddo

9000 continue

!---------------------------------------------------------
contains

subroutine value_out(name, who, value)

character(*) name, who
real(rp) value

!

if (value == 0) return
write (2, '(5a)') trim(name), '[', trim(who), '] = ', real_str(value, 14)

end subroutine value_out

!---------------------------------------------------------
! contains

function real_val(str) result (re)

character(*) str
real(rp) re

read (str, *) re

end function real_val

end program 
