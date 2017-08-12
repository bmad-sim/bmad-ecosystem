program tracking_method_test

use bmad
use mad_mod
use spin_mod

implicit none

type (lat_struct), target :: lat
type (coord_struct) start_orb, end_orb, end_bs, end_ptc
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele
type (track_struct) track

character(200) :: line(10)
character(40) :: lat_file  = 'tracking_method_test.bmad'
character(38) :: out_str, fmt, tol_str, track_method
integer :: i, j, ib, nargs, isn

logical print_extra
 
!

global_com%exit_on_error = .false.

fmt = '(a,t47, a, 7es18.10)'
track_method = ''

print_extra = .false.
nargs = cesr_iargc()

if (nargs > 0) then
  call cesr_getarg(1, lat_file)
  call cesr_getarg(2, track_method)
  print *, 'Using ', trim(lat_file)
  print_extra = .true.
  fmt = '(a, t47, a, 7es14.6)'
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

if (any(lat%beam_start%spin /= 0)) then
  bmad_com%spin_tracking_on = .true.
endif

open (1, file = 'output.now')

if (print_extra) then
  print '(a, t36, 7es18.10)', 'Start:', lat%beam_start%vec
  print *
endif

call track_it (1, 1)
if (.not. print_extra) call track_it (1, -1)

close(1)

!------------------------------------------------
contains

subroutine track_it(d_sign, p_sign)

integer d_sign, p_sign

!

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  if (branch%param%particle == photon$ .and. p_sign /= 1) cycle

  do i = 1, branch%n_ele_max - 1
    ele => branch%ele(i)
    if (p_sign /= 1 .and. (ele%key == sbend$ .or. ele%key == e_gun$)) cycle
    ele%spin_tracking_method = tracking$

    isn = 0
    do j = 1, n_methods$
      if (track_method /= '' .and. upcase(tracking_method_name(j)) /= upcase(track_method)) cycle
      if (.not. valid_tracking_method(ele, branch%param%particle, j)) cycle
      if (j == symp_map$ .or. j == custom$) cycle
      if (j == mad$) cycle   ! Ignore MAD
      if (j == taylor$ .and. lat%beam_start%direction == -1) cycle
      if (p_sign /= 1 .and. (j == taylor$ .or. j == linear$)) cycle
      ele%tracking_method = j

      if (ele%key /= taylor$) call kill_taylor(ele%taylor)

      if (ele%tracking_method == symp_lie_ptc$) then
        ele%spin_tracking_method = symp_lie_ptc$
      else
        ele%spin_tracking_method = tracking$
      endif

      if (j == linear$) then
        ele%tracking_method = symp_lie_ptc$
        if(ele%key == beambeam$ .or. ele%key == ac_kicker$) ele%tracking_method = bmad_standard$
        call make_mat6 (ele, branch%param, lat%beam_start)
        ele%tracking_method = j
      endif

      start_orb = lat%beam_start
      call init_coord (start_orb, start_orb, ele, upstream_end$, &
                                    default_tracking_species(branch%param), E_photon = ele%value(p0c$) * 1.006)
      if (p_sign == -1) start_orb%species = antiparticle(start_orb%species)
      start_orb%field = [1, 2]

      if (print_extra) then
        track%n_pt = -1  ! Reset
        call track1 (start_orb, ele, branch%param, end_orb, track = track)
      else
        call track1 (start_orb, ele, branch%param, end_orb)
      endif

      tol_str = trim(ele%name) // ':' // trim(tracking_method_name(j))
      if (p_sign == 1) then
        out_str = trim(ele%name) // ':' // trim(tracking_method_name(j))
      else
        out_str = trim(ele%name) // '-Anti:' // trim(tracking_method_name(j))
      endif

      if (ele%key == e_gun$) then
        write (1,fmt) '"' // trim(out_str) // '"' , tolerance(tol_str), end_orb%vec, c_light * (end_orb%t - start_orb%t)
      else
        write (1,fmt) '"' // trim(out_str) // '"' , tolerance(tol_str), end_orb%vec, (end_orb%vec(5) - start_orb%vec(5)) - &
                c_light * (end_orb%beta * (ele%ref_time - end_orb%t) - start_orb%beta * (ele%ref_time - ele%value(delta_ref_time$) - start_orb%t))
      endif

      if (ele%key == wiggler$) then
        if (j == symp_lie_bmad$) end_bs = end_orb
      else
        if (j == bmad_standard$) end_bs = end_orb
      endif
      if (j == symp_lie_ptc$)  end_ptc = end_orb

      if (j == bmad_standard$ .or. j == runge_kutta$ .or. j == symp_lie_ptc$ .or. j == time_runge_kutta$ .or. j == taylor$) then
        isn = isn + 1
        out_str = trim(out_str) // ' dSpin'
        write (line(isn), '(a, t47, a,  4f14.9, 4x, f14.9)') '"' // trim(out_str) // '"', tolerance_spin(tol_str), &
              end_orb%spin-start_orb%spin, norm2(end_orb%spin) - norm2(start_orb%spin)
      endif

      if (branch%param%particle == photon$) then
        write (1, '(3a, t47, a, 2es18.10)') '"', trim(ele%name), ':E_Field"', 'REL 5E-08', end_orb%field
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

end subroutine track_it

!--------------------------------------------------------------------------------------
! contains
  
character(10) function tolerance(instr)
character(38) :: instr

  select case (instr)
    case('AC_KICKER2:Time_Runge_Kutta')          ; tolerance = 'ABS 1e-11'
    case('QUADRUPOLE1:Time_Runge_Kutta')         ; tolerance = 'ABS 5e-12'
    case('QUADRUPOLE2:Time_Runge_Kutta')         ; tolerance = 'ABS 2e-11'
    case('QUADRUPOLE4:Time_Runge_Kutta')         ; tolerance = 'ABS 2e-11'
    case('QUADRUPOLE5:Time_Runge_Kutta')         ; tolerance = 'ABS 1e-11'
    case('RFCAVITY1:Runge_Kutta')                ; tolerance = 'ABS 1e-13'
    case('RFCAVITY1:Time_Runge_Kutta')           ; tolerance = 'ABS 2e-11'
    case('RFCAVITY2:Runge_Kutta')                ; tolerance = 'ABS 1e-13'
    case('RFCAVITY2:Time_Runge_Kutta')           ; tolerance = 'ABS 2e-11'
    case('SBEND2:Time_Runge_Kutta')              ; tolerance = 'ABS 1e-11'
    case('SBEND4:Bmad_Standard')                 ; tolerance = 'ABS 1e-11'
    case('SBEND4:Linear')                        ; tolerance = 'ABS 1e-11'
    case('SBEND4:Time_Runge_Kutta')              ; tolerance = 'ABS 4e-13'
    case('SBEND5:Bmad_Standard')                 ; tolerance = 'ABS 5e-13'
    case('SBEND5:Linear')                        ; tolerance = 'ABS 5e-13'
    case('SBEND6:Taylor')                        ; tolerance = 'ABS 2e-11'
    case('SBEND7:Bmad_Standard')                 ; tolerance = 'ABS 2e-13'
    case('SBEND7:Linear')                        ; tolerance = 'ABS 2e-13'
    case('SOLENOID1:Time_Runge_Kutta')           ; tolerance = 'ABS 1e-11'
    case('SOLENOID2:Symp_Lie_Bmad')              ; tolerance = 'ABS 2e-14'
    case('SOL_QUAD1:Time_Runge_Kutta')           ; tolerance = 'ABS 1e-11'
    case('SOL_QUAD2:Time_Runge_Kutta')           ; tolerance = 'ABS 2e-10'
    case('LCAVITY1:Bmad_Standard')               ; tolerance = 'ABS 2e-12'
    case('LCAVITY1:Time_Runge_Kutta')            ; tolerance = 'ABS 3e-11'
    case('LCAVITY1:Runge_Kutta')                 ; tolerance = 'ABS 1e-13'
    case('LCAVITY2:Time_Runge_Kutta')            ; tolerance = 'ABS 2e-13'
    case('LCAVITY2:Runge_Kutta')                 ; tolerance = 'ABS 1e-13'
    case('LCAVITY3:Runge_Kutta')                 ; tolerance = 'ABS 1e-13'
    case('LCAVITY3:Time_Runge_Kutta')            ; tolerance = 'ABS 2e-11'
    case('WIGGLER_MAP1:Time_Runge_Kutta')        ; tolerance = 'ABS 2e-13'
    case('WIGGLER_MAP1:Runge_Kutta')             ; tolerance = 'ABS 1e-13'
    case('WIGGLER_PERIODIC1:Runge_Kutta')        ; tolerance = 'ABS 5e-13'
    case('WIGGLER_PERIODIC1:Time_Runge_Kutta')   ; tolerance = 'ABS 2e-13'
    case default                                 ; tolerance = 'ABS 1e-14'
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
