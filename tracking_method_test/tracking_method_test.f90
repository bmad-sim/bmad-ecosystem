program tracking_method_test

use bmad

implicit none

type (lat_struct), target :: lat

character(40) :: input_file  = 'tracking_method_test.bmad'
character(32) :: final_str
integer :: i, j
 
type (coord_struct) end_orb

call bmad_parser ('tracking_method_test.bmad', lat)

open (1, file = 'output.now')

DO i = 1, lat%n_ele_max
   DO j = 1, n_methods$
      if(.not. valid_tracking_method(lat%ele(i),j) .or. j == symp_map$ .or. j == custom$) cycle
      lat%ele(i)%tracking_method = j
      if (j == linear$) then
        lat%ele(i)%tracking_method = bmad_standard$
        call make_mat6 (lat%ele(i), lat%param, lat%beam_start)
        lat%ele(i)%tracking_method = j
      endif
      call init_coord (lat%beam_start, lat%beam_start, ele = lat%ele(i), at_exit_end = .false.)
      call track1 (lat%beam_start, lat%ele(i), lat%param, end_orb)
      final_str = '"' // trim(key_name(lat%ele(i)%key)) // ':' // trim(calc_method_name(j)) // '"'  
      write (1,'(a,a,es24.15,es24.15,es24.15,es24.15,es24.15,es24.15,es24.15)',advance='no') final_str, 'REL  1E-10', end_orb%vec(1), end_orb%vec(2), end_orb%vec(3), end_orb%vec(4), end_orb%vec(5), end_orb%vec(6)
      write (1,*)
   END DO
   write (1,*)
END DO

close(1)

end program
