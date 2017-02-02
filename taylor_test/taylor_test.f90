program taylor_test

use bmad
use transfer_map_mod

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele1, ele2
type (taylor_struct) t_map(6), t

integer i, j
logical err

!

call bmad_parser ('taylor_test.bmad', lat)
open (1, file = 'output.now', recl = 200)

!

!call concat_taylor (lat%ele(3)%taylor, lat%ele(4)%taylor, t_map)
!call type_taylors(t_map)
!
!
!ele1 => lat%ele(1)
!ele1%tracking_method = taylor$
!call attribute_bookkeeper (ele1, lat%param, .true.)  ! So taylor map will not be killed in tracking
!call ele_to_taylor (ele1, lat%param, ele1%taylor, lat%beam_start)
!
!call concat_ele_taylor(ele1%taylor, ele2, ele2%taylor)

!

call transfer_map_calc (lat, t_map, err, 0, lat%n_ele_track, lat%beam_start, 0, .true.)

do i = 1, 6
  call sort_taylor_terms (t_map(i), t, 1e-20_rp)
  do j = 1, size(t%term)
    write (1, '(a, i1, a, 6i1, a, es18.10)') '"TM', i, '-', t%term(j)%expn, '"  REL 1E-7', t%term(j)%coef
  enddo
  write (1, *)
enddo

call transfer_map_from_s_to_s (lat, t_map, 0.1_rp, lat%param%total_length-2.0_rp, lat%beam_start)

do i = 1, 6
  call sort_taylor_terms (t_map(i), t, 1e-20_rp)
  do j = 1, size(t%term)
    write (1, '(a, i1, a, 6i1, a, es18.10)') '"TMSS', i, '-', t%term(j)%expn, '"  REL 1E-7', t%term(j)%coef
  enddo
  write (1, *)
enddo


end program
