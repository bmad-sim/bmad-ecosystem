program fixer_test

use bmad
use fixer_mod

implicit none

type (lat_struct), target :: lat, lat2
type (ele_pointer_struct), allocatable :: eles(:)
type (coord_struct), allocatable :: orbit(:), orbit2(:)

real(rp) :: chrom_a, chrom_b, delta_e = 0
integer n_loc
logical ok
character(2) :: run_str(3) = ['R1', 'R2', 'R3']

!

open (1, file = 'output.now')
bmad_com%auto_bookkeeper = .false.

!

call bmad_parser ('fixer_test.bmad', lat2)

call twiss_and_track(lat2, orbit2, orb_start = lat2%particle_start)
call chrom_calc(lat2, delta_e, chrom_a, chrom_b)
lat = lat2

call lat_ele_locator('fixer::*', lat, eles, n_loc)

call this_fixer_test(1, eles(1)%ele)
call this_fixer_test(2, eles(2)%ele)
call this_fixer_test(3, lat%ele(0))


!------------------------------------------------------------------------------------------
contains

subroutine this_fixer_test(indx, fixer)

type (ele_struct) fixer, ele0, eleN
type (coord_struct), allocatable :: orbit(:)

integer indx, nn

!

nn = lat%n_ele_track
fixer = lat2%ele(fixer%ix_ele)
ok = transfer_fixer_params(fixer, .true., orbit2(fixer%ix_ele))
call set_active_fixer(fixer)

if (fixer%ix_ele /= 0) call init_this_ele(lat%ele(0))
call init_this_ele(lat%ele(nn))

call twiss_and_track(lat, orbit, orb_start = orbit2(fixer%ix_ele))
call chrom_calc(lat, delta_e, chrom_a, chrom_b, orb0 = orbit2(fixer%ix_ele))

call orbit_delta_write(indx, 0, orbit(0), orbit2(0))
call twiss_delta_write(indx, 0, 'A', lat%ele(0)%a, lat2%ele(0)%a)
call twiss_delta_write(indx, 0, 'B', lat%ele(0)%b, lat2%ele(0)%b)
call twiss_delta_write(indx, 0, 'Z', lat%ele(0)%z, lat2%ele(0)%z)
call disp_delta_write(indx, 0, 'X', lat%ele(0)%x, lat2%ele(0)%x)
call disp_delta_write(indx, 0, 'Y', lat%ele(0)%y, lat2%ele(0)%y)

call orbit_delta_write(indx, nn, orbit(nn), orbit2(nn))
call twiss_delta_write(indx, nn, 'A', lat%ele(nn)%a, lat2%ele(nn)%a)
call twiss_delta_write(indx, nn, 'B', lat%ele(nn)%b, lat2%ele(nn)%b)
call twiss_delta_write(indx, nn, 'Z', lat%ele(nn)%z, lat2%ele(nn)%z)
call disp_delta_write(indx, nn, 'X', lat%ele(nn)%x, lat2%ele(nn)%x)
call disp_delta_write(indx, nn, 'Y', lat%ele(nn)%y, lat2%ele(nn)%y)

deallocate(orbit)

end subroutine this_fixer_test

!------------------------------------------------------------------------------------------
! contains

subroutine init_this_ele(ele)

type (ele_struct) ele
type (coord_struct) orbit

ele%a = twiss_struct()
ele%b = twiss_struct()
ele%z = twiss_struct()
ele%x = xy_disp_struct()
ele%y = xy_disp_struct()

end subroutine init_this_ele

!------------------------------------------------------------------------------------------
! contains

subroutine twiss_delta_write(indx, ix_ele, mode, t1, t2)

type (twiss_struct) t1, t2
integer indx, ix_ele
character(*) mode
character(20) str

!

str = run_str(indx) // '-E' // int_str(ix_ele) // '-' // mode // '-Twiss-'

write (1, '(2a, es16.8)') quote(trim(str) // 'beta'), ' ABS 1e-5', t1%beta - t2%beta
write (1, '(2a, es16.8)') quote(trim(str) // 'alpha'), ' ABS 1e-5', t1%alpha - t2%alpha
write (1, '(2a, es16.8)') quote(trim(str) // 'gamma'), ' ABS 1e-5', t1%gamma - t2%gamma
write (1, '(2a, es16.8)') quote(trim(str) // 'phi'), ' ABS 1e-5', t1%phi - t2%phi
write (1, '(2a, es16.8)') quote(trim(str) // 'eta'), ' ABS 1e-5', t1%eta - t2%eta
write (1, '(2a, es16.8)') quote(trim(str) // 'etap'), ' ABS 1e-5', t1%etap - t2%etap
write (1, '(2a, es16.8)') quote(trim(str) // 'deta_ds'), ' ABS 1e-5', t1%deta_ds - t2%deta_ds
write (1, '(2a, es16.8)') quote(trim(str) // 'dbeta_dpz'), ' ABS 1e-5', t1%dbeta_dpz - t2%dbeta_dpz
write (1, '(2a, es16.8)') quote(trim(str) // 'dalpha_dpz'), ' ABS 1e-5', t1%dalpha_dpz - t2%dalpha_dpz
write (1, '(2a, es16.8)') quote(trim(str) // 'deta_dpz'), ' ABS 1e-5', t1%deta_dpz - t2%deta_dpz
write (1, '(2a, es16.8)') quote(trim(str) // 'detap_dpz'), ' ABS 1e-5', t1%detap_dpz - t2%detap_dpz

end subroutine twiss_delta_write

!------------------------------------------------------------------------------------------
! contains

subroutine disp_delta_write(indx, ix_ele, mode, t1, t2)

type (xy_disp_struct) t1, t2
integer indx, ix_ele
character(*) mode
character(20) str

!

str = run_str(indx) // '-E' // int_str(ix_ele) // '-' // mode // '-Disp-'

write (1, '(2a, es16.8)') quote(trim(str) // 'eta'), ' ABS 1e-5', t1%eta - t2%eta
write (1, '(2a, es16.8)') quote(trim(str) // 'etap'), ' ABS 1e-5', t1%etap - t2%etap
write (1, '(2a, es16.8)') quote(trim(str) // 'deta_ds'), ' ABS 1e-5', t1%deta_ds - t2%deta_ds
write (1, '(2a, es16.8)') quote(trim(str) // 'deta_dpz'), ' ABS 1e-5', t1%deta_dpz - t2%deta_dpz
write (1, '(2a, es16.8)') quote(trim(str) // 'detap_dpz'), ' ABS 1e-5', t1%detap_dpz - t2%detap_dpz

end subroutine disp_delta_write

!------------------------------------------------------------------------------------------
! contains

subroutine orbit_delta_write(indx, ix_ele, orb1, orb2)

type (coord_struct) orb1, orb2
integer indx, ix_ele
character(20) str

!

str = run_str(indx) // '-E' // int_str(ix_ele) // '-Orbit-'

write (1, '(2a, es16.8)') quote(trim(str) // 'X'),  ' ABS 1E-10', orb1%vec(1) - orb2%vec(1)
write (1, '(2a, es16.8)') quote(trim(str) // 'PX'), ' ABS 1E-10', orb1%vec(2) - orb2%vec(2)
write (1, '(2a, es16.8)') quote(trim(str) // 'Y'), ' ABS 1E-10', orb1%vec(3) - orb2%vec(3)
write (1, '(2a, es16.8)') quote(trim(str) // 'PY'), ' ABS 1E-10', orb1%vec(4) - orb2%vec(4)
write (1, '(2a, es16.8)') quote(trim(str) // 'Z'), ' ABS 1E-10', orb1%vec(5) - orb2%vec(5)
write (1, '(2a, es16.8)') quote(trim(str) // 'PZ'), ' ABS 1E-10', orb1%vec(6) - orb2%vec(6)
write (1, '(2a, es16.8)') quote(trim(str) // 'T'), ' ABS 1E-10', orb1%t - orb2%t
write (1, '(2a, es16.8)') quote(trim(str) // 'P0C'), ' ABS 1E-10', orb1%p0c - orb2%p0c

end subroutine orbit_delta_write

end program
