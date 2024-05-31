program twiss_test

use bmad

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct) ele0
type (ele_struct), pointer :: elen
type (coord_struct), allocatable :: orbit(:)
type (twiss_struct) dt, tt
type (xy_disp_struct) dxy, xy

real(rp) dc(2,2), max_diff
integer n, ib, status
logical err

character(200) lat_file

! Idea is to check set_twiss

open (1, file = 'output.now')

lat_file = 'twiss_test.bmad'
call bmad_parser(lat_file, lat)

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  call reallocate_coord(orbit, branch%n_ele_max)
  call twiss_and_track(lat, orbit, status, ib, .true.)
  n = branch%n_ele_track
  elen => branch%ele(n-1)
  call transfer_twiss(branch%ele(n), ele0)
  call set_twiss(branch, ele0, n-1, err)

  max_diff = 0
  if (ib /= 0) write (1, '(a)')

  tt = ele0%a
  dt = twiss_diff(ele0%a, elen%a)
  max_diff = max(max_diff, dt%emit)
  write (1, '(2a, 3es14.6)') quote('br' // int_str(ib) // '-a-Twiss'), ' REL 1e-7', tt%beta, tt%alpha, tt%gamma
  write (1, '(2a, 3es14.6)') quote('br' // int_str(ib) // '-a-dTwiss'), ' REL 1e-7', dt%beta, dt%alpha, dt%gamma
  write (1, '(2a, 3es14.6)') quote('br' // int_str(ib) // '-a-Disp'),  ' REL 1e-7', tt%eta, tt%etap, tt%deta_ds
  write (1, '(2a, 3es14.6)') quote('br' // int_str(ib) // '-a-dDisp'),  ' REL 1e-7', dt%eta, dt%etap, dt%deta_ds

  tt = ele0%b
  dt = twiss_diff(ele0%b, elen%b)
  max_diff = max(max_diff, dt%emit)
  write (1, '(2a, 3es14.6)') quote('br' // int_str(ib) // '-b-Twiss'), ' REL 1e-7', tt%beta, tt%alpha, tt%gamma
  write (1, '(2a, 3es14.6)') quote('br' // int_str(ib) // '-b-dTwiss'), ' REL 1e-7', dt%beta, dt%alpha, dt%gamma
  write (1, '(2a, 3es14.6)') quote('br' // int_str(ib) // '-b-Disp'),  ' REL 1e-7', tt%eta, tt%etap, tt%deta_ds
  write (1, '(2a, 3es14.6)') quote('br' // int_str(ib) // '-b-dDisp'),  ' REL 1e-7', dt%eta, dt%etap, dt%deta_ds

  tt = ele0%z
  dt = twiss_diff(ele0%z, elen%z)
  max_diff = max(max_diff, dt%emit)
  write (1, '(2a, 3es14.6)') quote('br' // int_str(ib) // '-z-Twiss'), ' REL 1e-7', tt%beta, tt%alpha, tt%gamma
  write (1, '(2a, 3es14.6)') quote('br' // int_str(ib) // '-z-dTwiss'), ' REL 1e-7', dt%beta, dt%alpha, dt%gamma
  write (1, '(2a, 3es14.6)') quote('br' // int_str(ib) // '-z-Disp'),  ' REL 1e-7', tt%eta, tt%etap, tt%deta_ds
  write (1, '(2a, 3es14.6)') quote('br' // int_str(ib) // '-z-dDisp'),  ' REL 1e-7', dt%eta, dt%etap, dt%deta_ds

  dc = ele0%c_mat - elen%c_mat
  max_diff = max(max_diff, maxval(abs(dc)))
  write (1, '(2a, 4es14.6)') quote('br' // int_str(ib) // '-C_mat'), ' ABS 1e-7', ele0%c_mat
  write (1, '(2a, 4es14.6)') quote('br' // int_str(ib) // '-dC_mat'), ' ABS 1e-7', dc

  xy = ele0%x
  dxy = xy_disp_diff(ele0%x, elen%x)
  write (1, '(2a, 3es14.6)') quote('br' // int_str(ib) // '-x-Disp'), ' REL 1e-7', xy%eta, xy%etap, xy%deta_ds
  write (1, '(2a, 3es14.6)') quote('br' // int_str(ib) // '-x-dDisp'), ' REL 1e-7', dxy%eta, dxy%etap, dxy%deta_ds

  xy = ele0%y
  dxy = xy_disp_diff(ele0%y, elen%y)
  write (1, '(2a, 3es14.6)') quote('br' // int_str(ib) // '-y-Disp'), ' REL 1e-7', xy%eta, xy%etap, xy%deta_ds
  write (1, '(2a, 3es14.6)') quote('br' // int_str(ib) // '-y-dDisp'), ' REL 1e-7', dxy%eta, dxy%etap, dxy%deta_ds

  write (1, '(2a, 2l3, a)') quote('br' // int_str(ib) // '-Flib'), ' STR "', ele0%mode_flip, ele0%mode_flip .eqv. elen%mode_flip, '"'
  write (1, '(2a, es14.6)') quote('br' // int_str(ib) // '-dMax'), ' ABS 1e-7', max_diff

enddo

close (1)


!--------------------------------------------------------------------------------------------
contains

function twiss_diff(twiss1, twiss2) result (twiss_d)

type (twiss_struct) twiss1, twiss2, twiss_d

!

twiss_d%beta      = twiss1%beta - twiss2%beta
twiss_d%alpha     = twiss1%alpha - twiss2%alpha
twiss_d%gamma     = twiss1%gamma - twiss2%gamma
twiss_d%phi       = twiss1%phi - twiss2%phi
twiss_d%eta       = twiss1%eta - twiss2%eta
twiss_d%etap      = twiss1%etap - twiss2%etap
twiss_d%deta_ds   = twiss1%deta_ds - twiss2%deta_ds
twiss_d%sigma     = twiss1%sigma - twiss2%sigma
twiss_d%sigma_p   = twiss1%sigma_p - twiss2%sigma_p
twiss_d%emit      = twiss1%emit - twiss2%emit
twiss_d%norm_emit = twiss1%norm_emit - twiss2%norm_emit

twiss_d%emit = max(abs(twiss_d%beta), abs(twiss_d%alpha), abs(twiss_d%gamma), &
               abs(twiss_d%eta), abs(twiss_d%etap), abs(twiss_d%deta_ds))

end function twiss_diff

!--------------------------------------------------------------------------------------------
! contains

function xy_disp_diff(xy_disp1, xy_disp2) result (xy_disp_d)

type (xy_disp_struct) xy_disp1, xy_disp2, xy_disp_d

xy_disp_d%eta       = xy_disp1%eta - xy_disp2%eta
xy_disp_d%etap      = xy_disp1%etap - xy_disp2%etap
xy_disp_d%deta_ds   = xy_disp1%deta_ds - xy_disp2%deta_ds
xy_disp_d%sigma     = xy_disp1%sigma - xy_disp2%sigma

end function xy_disp_diff

end program
