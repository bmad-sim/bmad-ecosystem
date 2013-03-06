program mat6_calc_method_test

use bmad

use mad_mod

implicit none

type (lat_struct), target :: lat
type (ele_struct), dimension(:,:), allocatable :: temp_ele

character(40) :: lat_file  = 'mat6_calc_method_test.bmad'
character(44) :: final_str
character(*), PARAMETER  :: fmt1 = '(a,a,es24.15,es24.15,es24.15,es24.15,es24.15,es24.15)'
character(*), PARAMETER  :: fmt2 = '(a,a,es24.15)'
integer :: i, j, k, nargs
 
type (coord_struct) end_orb

nargs = cesr_iargc()
if (nargs == 1)then
   call cesr_getarg(1, lat_file)
   print *, 'Using ', trim(lat_file)
elseif (nargs > 1) then
  print *, 'Only one command line arg permitted.'
  call err_exit
endif

call bmad_parser (lat_file, lat)

open (1, file = 'output.now')

mad_print_if_misaligned = .false.

allocate (temp_ele(lat%n_ele_max - 1, n_methods$)) 

DO i = 1, lat%n_ele_max - 1
   DO j = 1, n_methods$
      if(.not. valid_mat6_calc_method(lat%ele(i),j) .or. j == static$ .or. j == custom$) cycle   
      call kill_taylor(lat%ele(i)%taylor)
      lat%ele(i)%mat6_calc_method = j  
      call make_mat6 (lat%ele(i), lat%param, lat%beam_start)
      call transfer_ele(lat%ele(i), temp_ele(i,j), .true.)
   END DO
END DO

DO k = 1, 8
   DO i = 1, lat%n_ele_max - 1
      DO j = 1, n_methods$
         if(.not. valid_mat6_calc_method(lat%ele(i),j) .or. j == static$ .or. j == custom$) cycle
         if (k < 7) then
            final_str = '"' // trim(temp_ele(i,j)%name) // ':' // trim(mat6_calc_method_name(j)) // ':MatrixRow' // trim(convert_to_string(k)) // '"' 
            write (1,fmt1,advance='no') final_str, 'REL  1E-10', temp_ele(i,j)%mat6(k,1), temp_ele(i,j)%mat6(k,2), temp_ele(i,j)%mat6(k,3), temp_ele(i,j)%mat6(k,4), temp_ele(i,j)%mat6(k,5), temp_ele(i,j)%mat6(k,6)
         else if (k == 7) then
            final_str = '"' // trim(temp_ele(i,j)%name) // ':' // trim(mat6_calc_method_name(j)) // ':Vector"' 
            write (1,fmt1,advance='no') final_str, 'REL  1E-10', temp_ele(i,j)%vec0(1), temp_ele(i,j)%vec0(2), temp_ele(i,j)%vec0(3), temp_ele(i,j)%vec0(4), temp_ele(i,j)%vec0(5), temp_ele(i,j)%vec0(6)
         else if (k == 8) then
            final_str = '"' // trim(temp_ele(i,j)%name) // ':' // trim(mat6_calc_method_name(j)) // ':Symp_Err"' 
            write (1,fmt2,advance='no') final_str, 'REL  1E-10', mat_symp_error(temp_ele(i,j)%mat6)
         end if
         write (1,*)
      END DO
      write (1,*)
   END DO
END DO

deallocate (temp_ele)

close(1)

contains
character(8) function convert_to_string(a)
  integer :: a
  write(convert_to_string, '(I1.1)') a
end function convert_to_string

end program
