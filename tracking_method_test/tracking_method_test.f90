program tracking_method_test

use bmad

use mad_mod

implicit none

type (lat_struct), target :: lat
type (coord_struct) start_orb, end_orb
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele

character(40) :: lat_file  = 'tracking_method_test.bmad'
character(38) :: final_str
integer :: i, j, ib, nargs

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

open (1, file = 'output.now')

if (print_extra) then
  print '(a, 42x, 7es18.10)', 'Start:', lat%beam_start%vec
  print *
endif

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  DO i = 1, branch%n_ele_max - 1
     ele => branch%ele(i)
     DO j = 1, n_methods$
        if(.not. valid_tracking_method(ele, branch%param%particle, j) .or. j == symp_map$ .or. j == custom$) cycle
        if(ele%key == elseparator$ .and. (j == runge_kutta$ .or. j == boris$ .or. j == time_runge_kutta$)) cycle
        call kill_taylor(ele%taylor)
        ele%tracking_method = j
        if (j == linear$) then
          ele%tracking_method = symp_lie_ptc$
          if(ele%key == beambeam$) ele%tracking_method = bmad_standard$
          call make_mat6 (ele, branch%param, lat%beam_start)
          ele%tracking_method = j
        endif
        start_orb = lat%beam_start
        call init_coord (start_orb, start_orb, ele, .false., branch%param%particle, E_photon = ele%value(p0c$) * 1.006)
        start_orb%field = [1, 2]
        call track1 (start_orb, ele, branch%param, end_orb)
        final_str = '"' // trim(ele%name) // ':' // trim(tracking_method_name(j)) // '"' 
        write (1,'(2a,7es18.10)') final_str, tolerance(final_str), end_orb%vec, c_light * (end_orb%t - start_orb%t)
        if (branch%param%particle == photon$) then
          write (1, '(4a, 2es18.10)') '"', trim(ele%name), ':E_Field', '"              REL 5E-08', end_orb%field
        endif
     END DO
     write (1,*)
  END DO
enddo

close(1)

!--------------------------------------------------------------------------------------
contains
  
character(10) function tolerance(instr)
character(38) :: instr

  select case (instr)
    case('"MIRROR1:Bmad_Standard"')      ; tolerance = 'ABS  1E-09'
    case('"RFCAVITY1:Time_Runge_Kutta"') ; tolerance = 'REL  1E-08'
    case('"RFCAVITY2:Linear"')           ; tolerance = 'REL  1E-07'
    case('"RFCAVITY2:Time_Runge_Kutta"') ; tolerance = 'REL  1E-08'
    case('"SBEND4:Symp_Lie_PTC"')        ; tolerance = 'REL  2E-05'
    case('"SBEND4:RUNGE_Kutta"')         ; tolerance = 'REL  3E-07'
    case('"SBEND4:Linear"')              ; tolerance = 'REL  2E-05'
    case('"SBEND4:Taylor"')              ; tolerance = 'REL  2E-07'
    case('"SBEND6:Symp_Lie_PTC"')        ; tolerance = 'REL  2E-05'
    case('"SBEND6:Linear"')              ; tolerance = 'REL  2E-05'
    case('"SBEND6:Taylor"')              ; tolerance = 'REL  2E-07'
    case('"LCAVITY1:Time_Runge_Kutta"')  ; tolerance = 'REL  1E-08'
    case('"LCAVITY3:Time_Runge_Kutta"')  ; tolerance = 'REL  1E-08'
    case default                         ; tolerance = 'REL  1E-10'
  end select

end function tolerance
end program
