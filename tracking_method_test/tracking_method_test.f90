program tracking_method_test

use bmad

use mad_mod

implicit none

type (lat_struct), target :: lat

character(40) :: lat_file  = 'tracking_method_test.bmad'
character(38) :: final_str
integer :: i, j, nargs
 
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

DO i = 1, lat%n_ele_max - 1
   DO j = 1, n_methods$
      if(.not. valid_tracking_method(lat%ele(i),j) .or. j == symp_map$ .or. j == custom$) cycle
      if(lat%ele(i)%key == elseparator$ .and. (j == runge_kutta$ .or. j == boris$ .or. j == time_runge_kutta$)) cycle
      call kill_taylor(lat%ele(i)%taylor)
      lat%ele(i)%tracking_method = j
      if (j == linear$) then
        lat%ele(i)%tracking_method = symp_lie_ptc$
        if(lat%ele(i)%key == beambeam$) lat%ele(i)%tracking_method = bmad_standard$
        call make_mat6 (lat%ele(i), lat%param, lat%beam_start)
        lat%ele(i)%tracking_method = j
      endif
      call init_coord (lat%beam_start, lat%beam_start, ele = lat%ele(i), at_downstream_end = .false.)
      call track1 (lat%beam_start, lat%ele(i), lat%param, end_orb)
      final_str = '"' // trim(lat%ele(i)%name) // ':' // trim(tracking_method_name(j)) // '"' 
      write (1,'(2a,7es22.13)') final_str, tolerance(final_str), end_orb%vec, c_light * (end_orb%t - lat%beam_start%t)
   END DO
   write (1,*)
END DO

close(1)

contains
  
character(10) function tolerance(instr)
  character(38) :: instr

  select case (instr)
     case('"RFCAVITY1:Time_Runge_Kutta"          ') ; tolerance = 'REL  1E-08'
     case('"RFCAVITY2:Linear"                    ') ; tolerance = 'REL  1E-07'
     case('"RFCAVITY2:Time_Runge_Kutta"          ') ; tolerance = 'REL  1E-08'
     case('"SBEND4:Symp_Lie_PTC"                 ') ; tolerance = 'REL  1E-03'
     case('"SBEND4:Linear"                       ') ; tolerance = 'REL  1E-03'
     case('"SBEND6:Symp_Lie_PTC"                 ') ; tolerance = 'REL  1E-03'
     case('"SBEND6:Linear"                       ') ; tolerance = 'REL  1E-03'
     case('"LCAVITY1:Time_Runge_Kutta"           ') ; tolerance = 'REL  1E-08'
     case('"LCAVITY3:Time_Runge_Kutta"           ') ; tolerance = 'REL  1E-08'
     case default ; tolerance = 'REL  1E-10'
  end select

end function tolerance

end program
