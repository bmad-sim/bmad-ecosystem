program mat6_calc_method_test

use bmad

use mad_mod

implicit none

type (lat_struct), target :: lat
type (ele_struct), target, allocatable :: temp_ele(:,:)
type (coord_struct) end_orb
type (ele_struct), pointer :: ele

character(40) :: lat_file  = 'mat6_calc_method_test.bmad'
character(44) :: final_str
character(*), PARAMETER  :: fmt1 = '(a,a,6es22.13)'
character(*), PARAMETER  :: fmt2 = '(a,a,es22.13)'

integer :: i, j, k, nargs
logical print_extra

!

print_extra = .false.
nargs = cesr_iargc()
if (nargs == 1)then
   call cesr_getarg(1, lat_file)
   print *, 'Using ', trim(lat_file)
   print_extra = .true.
elseif (nargs > 1) then
  print *, 'Only one command line arg permitted.'
  call err_exit
endif

call bmad_parser (lat_file, lat)

open (1, file = 'output.now', recl = 200)

allocate (temp_ele(lat%n_ele_max - 1, n_methods$)) 

DO i = 1, lat%n_ele_max - 1
   DO j = 1, n_methods$
      if(.not. valid_mat6_calc_method(lat%ele(i),j) .or. j == static$ .or. j == custom$) cycle   
      call kill_taylor(lat%ele(i)%taylor)
      lat%ele(i)%mat6_calc_method = j  
      call make_mat6 (lat%ele(i), lat%param, lat%beam_start, end_orb)
      if (print_extra .and. lat%ele(i)%mat6_calc_method == bmad_standard$) then
        write (1, '(a, 6es22.13)'), 'Start track:', lat%beam_start%vec
        write (1, '(a, 6es22.13)'), 'End track:  ', end_orb%vec 
        write (1, *)
      endif
      call transfer_ele(lat%ele(i), temp_ele(i,j), .true.)
   END DO
END DO

DO i = 1, lat%n_ele_max - 1
  DO k = 1, 8
    DO j = 1, n_methods$
      if(.not. valid_mat6_calc_method(lat%ele(i),j) .or. j == static$ .or. j == custom$) cycle
        ele => temp_ele(i,j)
      if (k < 7) then
        final_str = '"' // trim(ele%name) // ':' // trim(mat6_calc_method_name(j)) // ':MatrixRow' // trim(convert_to_string(k)) // '"' 
        write (1, fmt1) final_str, tolerance(final_str), ele%mat6(k,:)
      else if (k == 7) then
        final_str = '"' // trim(ele%name) // ':' // trim(mat6_calc_method_name(j)) // ':Vector"' 
        write (1, fmt1) final_str, tolerance(final_str), ele%vec0
      else if (k == 8) then
        final_str = '"' // trim(ele%name) // ':' // trim(mat6_calc_method_name(j)) // ':Symp_Err"' 
        write (1, fmt2) final_str, tolerance(final_str), mat_symp_error(ele%mat6, ele%value(p0c$)/ele%value(p0c_start$))
      end if
    END DO
    write (1,*)
  END DO
END DO

deallocate (temp_ele)

close(1)

!----------------------------------------------------------------------
contains

character(8) function convert_to_string(a)
integer :: a
write(convert_to_string, '(I1.1)') a
end function convert_to_string

!----------------------------------------------------------------------
! contains

character(10) function tolerance(instr)
character(44) :: instr

select case (instr)
case ('"E_GUN1:Tracking:MatrixRow1"')              ; tolerance = 'ABS 5e-10'
case ('"E_GUN1:Tracking:MatrixRow2"')              ; tolerance = 'ABS 5e-11'
case ('"E_GUN1:Tracking:MatrixRow3"')              ; tolerance = 'ABS 7e-10'
case ('"E_GUN1:Tracking:MatrixRow4"')              ; tolerance = 'ABS 1e-10'
case ('"E_GUN1:Tracking:MatrixRow5"')              ; tolerance = 'ABS 4e-07'
case ('"E_GUN1:Tracking:MatrixRow6"')              ; tolerance = 'ABS 6e-08'
case ('"E_GUN1:Tracking:Symp_Err"')                ; tolerance = 'ABS 4e-07'
case ('"E_GUN1:Tracking:Vector"')                  ; tolerance = 'ABS 6e-10'
case ('"RFCAVITY1:Tracking:MatrixRow5"')           ; tolerance = 'ABS 2e-11'
case ('"RFCAVITY1:Tracking:MatrixRow6"')           ; tolerance = 'ABS 1e-10'
case ('"RFCAVITY1:Tracking:Symp_Err"')             ; tolerance = 'ABS 9e-11'
case ('"RFCAVITY2:Tracking:MatrixRow1"')           ; tolerance = 'ABS 4e-12'
case ('"RFCAVITY2:Tracking:MatrixRow2"')           ; tolerance = 'ABS 6e-12'
case ('"RFCAVITY2:Tracking:MatrixRow3"')           ; tolerance = 'ABS 2e-12'
case ('"RFCAVITY2:Tracking:MatrixRow5"')           ; tolerance = 'ABS 8e-11'
case ('"RFCAVITY2:Tracking:MatrixRow6"')           ; tolerance = 'ABS 7e-10'
case ('"RFCAVITY2:Tracking:Symp_Err"')             ; tolerance = 'ABS 6e-11'
case ('"RFCAVITY2:Tracking:Vector"')               ; tolerance = 'ABS 4e-12'
case ('"SBEND2:Tracking:MatrixRow5"')              ; tolerance = 'ABS 2e-12'
case ('"SBEND4:Bmad_Standard:MatrixRow3"')         ; tolerance = 'ABS 9e-12'
case ('"SBEND4:Bmad_Standard:MatrixRow5"')         ; tolerance = 'ABS 4e-09'
case ('"SBEND4:Bmad_Standard:Symp_Err"')           ; tolerance = 'ABS 1e-11'
case ('"SBEND4:Bmad_Standard:Vector"')             ; tolerance = 'ABS 2e-11'
case ('"SBEND4:Tracking:MatrixRow2"')              ; tolerance = 'ABS 2e-12'
case ('"SBEND5:Tracking:MatrixRow2"')              ; tolerance = 'ABS 2e-12'
case ('"SBEND5:Tracking:MatrixRow5"')              ; tolerance = 'ABS 4e-10'
case ('"SBEND5:Tracking:Symp_Err"')                ; tolerance = 'ABS 4e-12'
case ('"SBEND5:Tracking:Vector"')                  ; tolerance = 'ABS 2e-12'
case ('"SBEND7:Tracking:MatrixRow1"')              ; tolerance = 'ABS 3e-09'
case ('"SBEND7:Tracking:MatrixRow2"')              ; tolerance = 'ABS 4e-12'
case ('"SBEND7:Tracking:MatrixRow3"')              ; tolerance = 'ABS 2e-10'
case ('"SBEND7:Tracking:Symp_Err"')                ; tolerance = 'ABS 2e-10'
case ('"SBEND7:Tracking:Vector"')                  ; tolerance = 'ABS 4e-12'
case ('"SOL_QUAD2:Bmad_Standard:MatrixRow1"')      ; tolerance = 'ABS 3e-12'
case ('"SOL_QUAD2:Bmad_Standard:MatrixRow4"')      ; tolerance = 'ABS 2e-12'
case ('"SOL_QUAD2:Tracking:MatrixRow1"')           ; tolerance = 'ABS 8e-10'
case ('"SOL_QUAD2:Tracking:MatrixRow2"')           ; tolerance = 'ABS 2e-11'
case ('"SOL_QUAD2:Tracking:MatrixRow3"')           ; tolerance = 'ABS 8e-11'
case ('"SOL_QUAD2:Tracking:MatrixRow4"')           ; tolerance = 'ABS 2e-10'
case ('"SOL_QUAD2:Tracking:MatrixRow5"')           ; tolerance = 'ABS 5e-12'
case ('"SOL_QUAD2:Tracking:Symp_Err"')             ; tolerance = 'ABS 2e-10'
case ('"SOL_QUAD2:Tracking:Vector"')               ; tolerance = 'ABS 6e-12'
case ('"WIGGLER_PERIODIC1:Tracking:MatrixRow4"')   ; tolerance = 'ABS 6e-12'
case ('"WIGGLER_PERIODIC1:Tracking:Symp_Err"')     ; tolerance = 'ABS 7e-12'
case default                                       ; tolerance = 'ABS 1E-12'
end select

end function tolerance

end program
