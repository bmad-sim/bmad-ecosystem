program twiss_test

use bmad

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct) t_ele, ds_ele
type (ele_struct), pointer :: ele0, elen, ele
type (coord_struct), allocatable :: orb0(:), orb1(:)
type (twiss_struct) dt, tt
type (xy_disp_struct) dxy, xy

real(rp) dc(2,2), max_diff, dpz, eta_vec(6), dorb(6)
integer n, ib, status
logical err

character(200) lat_file
character(3) b_str

! Check twiss_propagate vs true dispersion calculated via tracking.

open (1, file = 'output.now')

dpz = 1d-8
lat_file = 'twiss_test.bmad'
call bmad_parser(lat_file, lat)

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  b_str = 'br' // int_str(ib)
  ele0 => branch%ele(0)
  call reallocate_coord(orb0, branch%n_ele_max)
  call twiss_and_track(lat, orb0, status, ib, .true.)
  n = branch%n_ele_track
  eta_vec = [ele0%x%eta, ele0%x%etap, ele0%y%eta, ele0%y%etap, ele0%z%eta, 1.0_rp] 

  call reallocate_coord(orb1, branch%n_ele_max)
  call init_coord(orb1(0), orb0(0)%vec + dpz * eta_vec, ele0, downstream_end$)
  call track_all(lat, orb1, ib)

  dorb = (orb1(n)%vec - orb0(n)%vec) / (orb1(n)%vec(6) - orb0(n)%vec(6))
  elen => branch%ele(n)
  eta_vec = [elen%x%eta, elen%x%etap, elen%y%eta, elen%y%etap, elen%z%eta, 1.0_rp] 

  if (ib /= 0) write (1, '(a)')
  write (1, '(2a, 6es14.6)') quote(b_str // '-Fwd-Orb'), ' ABS 1e-7', dorb
  write (1, '(2a, 6es14.6)') quote(b_str // '-Fwd-dOrb'), ' ABS 1e-14', dorb - eta_vec
  max_diff = maxval(abs(dorb - eta_vec))
  write (1, '(2a, es14.6)') quote(b_str // '-dMax-Fwd'), ' ABS 1e-7', max_diff

  ! Check of set_twiss by setting Twiss at ele n-2 with the Twiss at ele n-1.
  ! Note: (a,b) %deta_ds are not matched to.

  call transfer_twiss(branch%ele(n-1), t_ele)
  call set_twiss(branch, t_ele, n-2, .true., err)
  elen => branch%ele(n-2)
  call transfer_twiss(elen, ds_ele)
  call set_twiss(branch, t_ele, n-2, .false., err)

  max_diff = 0
  tt = t_ele%a
  dt = twiss_diff(t_ele%a, ds_ele%a, elen%a, max_diff)
  write (1, '(a)')
  write (1, '(2a, 3es14.6)') quote(b_str // '-a-Twiss'), ' REL 1e-7', tt%beta, tt%alpha, tt%gamma
  write (1, '(2a, 3es14.6)') quote(b_str // '-a-dTwiss'), ' ABS 1E-14', dt%beta, dt%alpha, dt%gamma
  write (1, '(2a, 3es14.6)') quote(b_str // '-a-Disp'),  ' REL 1e-7', tt%eta, tt%etap, tt%deta_ds
  write (1, '(2a, 3es14.6)') quote(b_str // '-a-dDisp'),  ' ABS 1E-14', dt%eta, dt%etap ! , dt%deta_ds

  tt = t_ele%b
  dt = twiss_diff(t_ele%b, ds_ele%b, elen%b, max_diff)
  write (1, '(a)')
  write (1, '(2a, 3es14.6)') quote(b_str // '-b-Twiss'), ' REL 1e-7', tt%beta, tt%alpha, tt%gamma
  write (1, '(2a, 3es14.6)') quote(b_str // '-b-dTwiss'), ' ABS 1E-14', dt%beta, dt%alpha, dt%gamma
  write (1, '(2a, 3es14.6)') quote(b_str // '-b-Disp'),  ' REL 1e-7', tt%eta, tt%etap, tt%deta_ds
  write (1, '(2a, 3es14.6)') quote(b_str // '-b-dDisp'),  ' ABS 1E-14', dt%eta, dt%etap ! , dt%deta_ds

  tt = t_ele%z
  dt = twiss_diff(t_ele%z, ds_ele%z, elen%z, max_diff)
  write (1, '(a)')
  write (1, '(2a, 3es14.6)') quote(b_str // '-z-Twiss'), ' REL 1e-7', tt%beta, tt%alpha, tt%gamma
  write (1, '(2a, 3es14.6)') quote(b_str // '-z-dTwiss'), ' ABS 1E-14', dt%beta, dt%alpha, dt%gamma
  write (1, '(2a, 3es14.6)') quote(b_str // '-z-Disp'),  ' REL 1e-7', tt%eta, tt%etap, tt%deta_ds
  write (1, '(2a, 3es14.6)') quote(b_str // '-z-dDisp'),  ' ABS 1E-14', dt%eta, dt%etap, dt%deta_ds

  xy = t_ele%x
  dxy = xy_disp_diff(t_ele%x, ds_ele%x, elen%x, max_diff)
  write (1, '(a)')
  write (1, '(2a, 3es14.6)') quote(b_str // '-x-Disp'), ' REL 1e-7', xy%eta, xy%etap, xy%deta_ds
  write (1, '(2a, 3es14.6)') quote(b_str // '-x-dDisp'), ' ABS 1E-14', dxy%eta, dxy%etap, dxy%deta_ds

  xy = t_ele%y
  dxy = xy_disp_diff(t_ele%y, ds_ele%y, elen%y, max_diff)
  write (1, '(2a, 3es14.6)') quote(b_str // '-y-Disp'), ' REL 1e-7', xy%eta, xy%etap, xy%deta_ds
  write (1, '(2a, 3es14.6)') quote(b_str // '-y-dDisp'), ' ABS 1E-14', dxy%eta, dxy%etap, dxy%deta_ds

  dc = t_ele%c_mat - elen%c_mat
  max_diff = max(max_diff, maxval(abs(dc)))
  write (1, '(a)')
  write (1, '(2a, 4es14.6)') quote(b_str // '-C_mat'), ' ABS 1e-7', t_ele%c_mat
  write (1, '(2a, 4es14.6)') quote(b_str // '-dC_mat'), ' ABS 1e-14', dc

  write (1, '(2a, 2l1, a)') quote(b_str // '-flip'), ' STR "', t_ele%mode_flip, t_ele%mode_flip .eqv. elen%mode_flip, '"'
  write (1, '(2a, es14.6)') quote(b_str // '-dMax-Rev'), ' ABS 1e-7', max_diff
enddo

!----------------------------------------------------------

branch => lat%branch(0)
ele => branch%ele(0)

write (1, '(a, 3es16.8)') quote('dispersion0-x') // ' ABS 1E-10', ele%x%eta, ele%x%etap, ele%x%deta_ds
write (1, '(a, 3es16.8)') quote('dispersion0-y') // ' ABS 1E-10', ele%y%eta, ele%y%etap, ele%y%deta_ds
write (1, '(a, 3es16.8)') quote('dispersion0-a') // ' ABS 1E-10', ele%a%eta, ele%a%etap, ele%a%deta_ds
write (1, '(a, 3es16.8)') quote('dispersion0-b') // ' ABS 1E-10', ele%b%eta, ele%b%etap, ele%b%deta_ds



!

close (1)

!--------------------------------------------------------------------------------------------
contains

function twiss_diff(twiss1, ds_twiss, twiss2, max_diff) result (twiss_d)

type (twiss_struct) twiss1, ds_twiss, twiss2, twiss_d
real(rp) max_diff

!

twiss_d%beta      = twiss1%beta - twiss2%beta
twiss_d%alpha     = twiss1%alpha - twiss2%alpha
twiss_d%gamma     = twiss1%gamma - twiss2%gamma
twiss_d%phi       = twiss1%phi - twiss2%phi
twiss_d%eta       = twiss1%eta - twiss2%eta
twiss_d%etap      = twiss1%etap - twiss2%etap
twiss_d%deta_ds   = twiss1%deta_ds - ds_twiss%deta_ds
twiss_d%sigma     = twiss1%sigma - twiss2%sigma
twiss_d%sigma_p   = twiss1%sigma_p - twiss2%sigma_p
twiss_d%emit      = twiss1%emit - twiss2%emit
twiss_d%norm_emit = twiss1%norm_emit - twiss2%norm_emit

max_diff = max(max_diff, abs(twiss_d%beta), abs(twiss_d%alpha), abs(twiss_d%gamma), &
                    abs(twiss_d%eta), abs(twiss_d%etap))  ! , abs(twiss_d%deta_ds))

end function twiss_diff

!--------------------------------------------------------------------------------------------
! contains

function xy_disp_diff(xy_disp1, ds_xy, xy_disp2, max_diff) result (xy_disp_d)

type (xy_disp_struct) xy_disp1, ds_xy, xy_disp2, xy_disp_d
real(rp) max_diff

xy_disp_d%eta       = xy_disp1%eta - xy_disp2%eta
xy_disp_d%etap      = xy_disp1%etap - xy_disp2%etap
xy_disp_d%deta_ds   = xy_disp1%deta_ds - ds_xy%deta_ds
xy_disp_d%sigma     = xy_disp1%sigma - xy_disp2%sigma

max_diff = max(max_diff, abs(xy_disp_d%eta), abs(xy_disp_d%etap), abs(xy_disp_d%deta_ds))

end function xy_disp_diff

end program
