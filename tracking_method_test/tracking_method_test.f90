program tracking_method_test

use bmad
use mad_mod
use spin_mod

implicit none

type (lat_struct), target :: lat
type (coord_struct) start_orb, end_orb, end_bs, end_ptc
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele
real(rp) start_p(3), end_p(3)

character(200) :: line(4)
character(40) :: lat_file  = 'tracking_method_test.bmad'
character(38) :: final_str, fmt
integer :: i, j, ib, nargs, isn

logical print_extra
 
!

global_com%exit_on_error = .false.

fmt = '(a,t42, a, 7es18.10)'

print_extra = .false.
nargs = cesr_iargc()
if (nargs > 1) then
  print *, 'Only one command line arg permitted.'
  call err_exit

elseif (nargs > 0)then
  call cesr_getarg(1, lat_file)
  print *, 'Using ', trim(lat_file)
  print_extra = .true.
  fmt = '(a, t42, a, 7es14.6)'
endif

call bmad_parser (lat_file, lat)

if (print_extra) then
  if (lat%param%geometry == open$) then
    bmad_com%convert_to_kinetic_momentum = .false.
    print *, '*** Note: wiggler end kicks not cancelled (so like PTC tracking).'
  else
    bmad_com%convert_to_kinetic_momentum = .true.
    print *, '*** Note: wiggler end kicks cancelled (so like RUNGE_KUTTA tracking).'
  endif
endif

if (lat%beam_start%spin(1) /= 0 .or. lat%beam_start%spin(2) /= 0) then
  bmad_com%spin_tracking_on = .true.
endif

open (1, file = 'output.now')

if (print_extra) then
  print '(a, t36, 7es18.10)', 'Start:', lat%beam_start%vec
  print *
endif

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  do i = 1, branch%n_ele_max - 1
    ele => branch%ele(i)
    ele%spin_tracking_method = tracking$
    isn = 0
    do j = 1, n_methods$
      if(.not. valid_tracking_method(ele, branch%param%particle, j) .or. j == symp_map$ .or. j == custom$) cycle
      if (ele%key /= taylor$) call kill_taylor(ele%taylor)
!      if ((ele%key == rfcavity$ .or. ele%key == lcavity$) .and. (j == runge_kutta$ .or. j == time_runge_kutta$ .or. j == boris$)) cycle
      ele%tracking_method = j
      if (ele%tracking_method == symp_lie_ptc$) then
        ele%spin_tracking_method = symp_lie_ptc$
      else
        ele%spin_tracking_method = tracking$
      endif
      if (j == linear$) then
        ele%tracking_method = symp_lie_ptc$
        if(ele%key == beambeam$) ele%tracking_method = bmad_standard$
        call make_mat6 (ele, branch%param, lat%beam_start)
        ele%tracking_method = j
      endif
      start_orb = lat%beam_start
      call init_coord (start_orb, start_orb, ele, upstream_end$, branch%param%particle, E_photon = ele%value(p0c$) * 1.006)
      start_orb%field = [1, 2]
      call track1 (start_orb, ele, branch%param, end_orb)
      final_str = trim(ele%name) // ':' // trim(tracking_method_name(j))
      write (1,fmt) '"' // trim(final_str) // '"' , tolerance(final_str), end_orb%vec, c_light * (end_orb%t - start_orb%t)

      if (ele%key == wiggler$) then
        if (j == symp_lie_bmad$) end_bs = end_orb
      else
        if (j == bmad_standard$) end_bs = end_orb
      endif
      if (j == symp_lie_ptc$)  end_ptc = end_orb

      if (j == bmad_standard$ .or. j == runge_kutta$ .or. j == symp_lie_ptc$ .or. j == time_runge_kutta$) then
        isn = isn + 1
        start_p = spinor_to_vec(start_orb%spin)
        end_p = spinor_to_vec(end_orb%spin)
        final_str = trim(final_str) // ' dSpin'
        write (line(isn), '(a, t42, a,  4f14.9, 4x, f14.9)') '"' // trim(final_str) // '"', tolerance_spin(final_str), &
              end_p-start_p, norm2(end_p) - norm2(start_p)
      endif

      if (branch%param%particle == photon$) then
        write (1, '(3a, t42, a, 2es18.10)') '"', trim(ele%name), ':E_Field"', 'REL 5E-08', end_orb%field
      endif
    end do

    do j = 1, isn
      write (1, '(a)') trim(line(j))
    enddo

    if (print_extra) then
      print '(a, t36, 7es18.10)', 'Diff PTC - BS:', end_ptc%vec - end_bs%vec
      print *
    endif

    write (1,*)
  end do
enddo

close(1)

!--------------------------------------------------------------------------------------
contains
  
character(10) function tolerance(instr)
character(38) :: instr

  select case (instr)
    case('RFCAVITY1:Time_Runge_Kutta')           ; tolerance = 'REL 1E-09'
    case('RFCAVITY2:Time_Runge_Kutta')           ; tolerance = 'REL 1E-09'
    case('SBEND4:Symp_Lie_PTC')                  ; tolerance = 'REL 1E-09'
    case('SBEND4:Bmad_Standard')                 ; tolerance = 'REL 1E-09'
    case('SBEND4:Linear')                        ; tolerance = 'REL 1E-09'
    case('SBEND6:Symp_Lie_PTC')                  ; tolerance = 'REL 4E-05'
    case('SBEND6:Linear')                        ; tolerance = 'REL 4E-05'
    case('SBEND6:Taylor')                        ; tolerance = 'REL 2E-08'
    case('LCAVITY1:Time_Runge_Kutta')            ; tolerance = 'REL 1E-09'
    case('LCAVITY3:Time_Runge_Kutta')            ; tolerance = 'REL 2E-09'
    case default                                 ; tolerance = 'REL 1E-10'
  end select

end function tolerance

!--------------------------------------------------------------------------------------
! contains
  
character(10) function tolerance_spin(instr)
character(38) :: instr

  select case (instr)
    case('WIGGLER_PERIODIC1:Runge_Kutta dSpin')  ; tolerance_spin = 'ABS 2E-7'
    case default                                 ; tolerance_spin = 'ABS 1E-8'
  end select

end function tolerance_spin
end program
