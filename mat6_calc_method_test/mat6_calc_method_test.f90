program mat6_calc_method_test

use bmad

implicit none

type (lat_struct), target :: lat

character(40) :: input_file  = 'mat6_calc_method_test.bmad'
character(44) :: final_str
character(*), PARAMETER  :: fmt = '(a,a,es24.15,es24.15,es24.15,es24.15,es24.15,es24.15)'
integer :: i, j, k
 
type (coord_struct) end_orb

call bmad_parser ('mat6_calc_method_test.bmad', lat)

open (1, file = 'output.now')

DO i = 1, lat%n_ele_max - 1
   DO j = 1, n_methods$
      if(.not. valid_mat6_calc_method(lat%ele(i),j) .or. j == static$ .or. j == tracking$ .or. j == custom$) cycle     
      lat%ele(i)%mat6_calc_method = j  
      call make_mat6 (lat%ele(i), lat%param, lat%beam_start)
      DO k = 1, 6
         final_str = '"' // trim(lat%ele(i)%name) // ':' // trim(calc_method_name(j)) // ':MatrixRow' // trim(convert_to_string(k)) // '"' 
         write (1,fmt,advance='no') final_str, 'REL  1E-10', lat%ele(i)%mat6(k,1), lat%ele(i)%mat6(k,2), lat%ele(i)%mat6(k,3), lat%ele(i)%mat6(k,4), lat%ele(i)%mat6(k,5), lat%ele(i)%mat6(k,6)
         write (1,*)
      END DO
      final_str = '"' // trim(lat%ele(i)%name) // ':' // trim(calc_method_name(j)) // ':Vector"' 
      write (1,fmt,advance='no') final_str, 'REL  1E-10', lat%ele(i)%vec0(1), lat%ele(i)%vec0(2), lat%ele(i)%vec0(3), lat%ele(i)%vec0(4), lat%ele(i)%vec0(5), lat%ele(i)%vec0(6)
      write (1,*)
   END DO
   write (1,*)
END DO

close(1)

contains
character(8) function convert_to_string(a)
  integer :: a
  write(convert_to_string, '(I1.1)') a
end function convert_to_string

end program
