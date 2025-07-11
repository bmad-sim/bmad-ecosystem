program taylor_test

use bmad
use transfer_map_mod

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele1, ele2
type (taylor_struct) t_map(6), t2_map(6), spin_map(0:3), t, t2

real(rp) max_diff, diff
integer i, j, nargs, tdiff
logical err, debug_mode

character(100) lat_file

!

debug_mode = .false.
nargs = command_argument_count()
lat_file = 'taylor_test.bmad'

if (nargs > 0) then
  debug_mode = .true.
  call get_command_argument(1, lat_file)
endif

!

call bmad_parser (lat_file, lat)
open (1, file = 'output.now', recl = 200)

!call concat_taylor (lat%ele(3)%taylor, lat%ele(4)%taylor, t_map)
!call type_taylors(t_map)
!
!
!ele1 => lat%ele(1)
!ele1%tracking_method = taylor$
!call attribute_bookkeeper (ele1, lat%param, .true.)  ! So taylor map will not be killed in tracking
!call ele_to_taylor (ele1, lat%particle_start)
!
!call concat_ele_taylor(ele1%taylor, ele2, ele2%taylor)

call transfer_map_calc (lat, t_map, err, 0, lat%n_ele_track, lat%particle_start, 0, .true.)
!call transfer_map_calc (lat, t2_map, err, 0, lat%n_ele_track, lat%particle_start, 0, .true., spin_map = spin_map)
max_diff = 0

do i = 1, 6
  call sort_taylor_terms (t_map(i), t, 1e-20_rp)
!  call sort_taylor_terms (t2_map(i), t2, 1e-20_rp)
  do j = 1, size(t%term)
!    tdiff = maxval(abs(t%term(j)%expn-t2%term(j)%expn))
!    diff = abs(t%term(j)%coef - t2%term(j)%coef)
    write (1, '(a, i1, a, 6i1, a, es18.10)')     '"TM', i, '-', t%term(j)%expn, '"  REL 1E-7', t%term(j)%coef
!    write (1, '(a, i1, a, 6i1, a, i6, es12.4)') '"dTM', i, '-', t%term(j)%expn, '"  ABS 1E-10', tdiff, diff
                        
!    max_diff = max(max_diff, diff, real(tdiff))
  enddo
  write (1, *)
enddo

if (debug_mode) then
!  print '(a, es12.4)', 'max_diff: ', max_diff
  stop
endif

!

call transfer_map_from_s_to_s (lat, t_map, 0.1_rp, lat%param%total_length-2.0_rp, lat%particle_start)

do i = 1, 6
  call sort_taylor_terms (t_map(i), t, 1e-20_rp)
  do j = 1, size(t%term)
    write (1, '(a, i1, a, 6i1, a, es18.10)') '"TMSS', i, '-', t%term(j)%expn, '"  REL 4E-7', t%term(j)%coef
  enddo
  write (1, *)
enddo


end program
