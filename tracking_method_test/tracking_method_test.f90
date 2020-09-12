program tracking_method_test

use bmad
use tpsa

implicit none

type (lat_struct), target :: lat

character(200) :: line(10), line_debug(10)
character(100) :: lat_file  = 'tracking_method_test.bmad'
character(46) :: out_str, fmt, track_method
integer :: i, j, ib, nargs, isn

logical debug_mode
 
!
!switch_bessel = .false.
global_com%exit_on_error = .false.

fmt = '(a, t49, a, 7es18.10)'
track_method = ''

debug_mode = .false.
nargs = cesr_iargc()

if (nargs > 0) then
  call cesr_getarg(1, lat_file)
  call cesr_getarg(2, track_method)
  print *, 'Using ', trim(lat_file)
  debug_mode = .true.
  fmt = '(a, t49, a, 7es14.6)'
endif

call bmad_parser (lat_file, lat, .false.)

if (debug_mode) then
  if (lat%param%geometry == open$) then
    bmad_com%convert_to_kinetic_momentum = .false.
    print *, '*** Note: wiggler end kicks not cancelled (so like PTC tracking).'
  else
    bmad_com%convert_to_kinetic_momentum = .true.
    print *, '*** Note: wiggler end kicks cancelled (so like RUNGE_KUTTA tracking).'
  endif
endif

if (any(lat%particle_start%spin /= 0)) then
  bmad_com%spin_tracking_on = .true.
endif

open (1, file = 'output.now')

if (debug_mode) then
  print '(a, t36, 7es18.10)', 'Start:', lat%particle_start%vec
  print *
  print '(a, t46, a, t64, a, t82, a, t100, a, t118, a, t136, a, t143, a)', &
                            'Name: Tracking_Method', 'x', 'px', 'y', 'py', 'z', 'pz', 'dz-d(v*(t_ref-t))'
endif

call track_it (lat, 1, 1)
if (.not. debug_mode) call track_it (lat, 1, -1)

close(1)

!------------------------------------------------
contains

subroutine track_it(lat, d_sign, p_sign)

type (lat_struct), target :: lat
type (coord_struct) start_orb, end_orb, end_bs, end_ptc
type (ele_struct), pointer :: ele
type (branch_struct), pointer :: branch
type (track_struct) track
integer d_sign, p_sign

!

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  if (branch%param%particle == photon$ .and. p_sign == -1) cycle

  do i = 1, branch%n_ele_max - 1
    ele => branch%ele(i)
    if (p_sign == -1 .and. ele%key == e_gun$) cycle
    ele%spin_tracking_method = tracking$

    isn = 0
    do j = 1, n_methods$
      if ((j == fixed_step_runge_kutta$ .or. j == fixed_step_time_runge_kutta$)) cycle
      if (track_method /= '' .and. upcase(tracking_method_name(j)) /= upcase(track_method)) cycle
      if (.not. valid_tracking_method(ele, branch%param%particle, j)) cycle
      if (j == symp_map$ .or. j == custom$) cycle
      if (j == mad$) cycle   ! Ignore MAD
      if (j == taylor$ .and. lat%particle_start%direction == -1) cycle
      if (p_sign == -1 .and. (j == taylor$ .or. j == linear$)) cycle
      ele%tracking_method = j

      if (ele%key == e_gun$ .and. (j == runge_kutta$ .or. j == fixed_step_runge_kutta$)) cycle

      if (ele%key /= taylor$) call kill_taylor(ele%taylor)

      if (ele%tracking_method == symp_lie_ptc$) then
        ele%spin_tracking_method = symp_lie_ptc$
      else
        ele%spin_tracking_method = tracking$
      endif

      if (j == linear$) then
        ele%tracking_method = symp_lie_ptc$
        if (ele%key == ac_kicker$) ele%tracking_method = bmad_standard$
        if (lat%particle_start%direction == 1) then
          call make_mat6 (ele, branch%param, lat%particle_start)
        else  ! Can happen with a test lattice file
          call make_mat6 (ele, branch%param)
        endif
        ele%tracking_method = j
      endif

      start_orb = lat%particle_start

      if (p_sign == -1) then
        start_orb%direction = -1
        start_orb%species = antiparticle(start_orb%species)
        lat%absolute_time_tracking = .true.
      else
        start_orb%species = default_tracking_species(branch%param)
      endif

      call init_coord (start_orb, start_orb, ele, start_end$, start_orb%species, start_orb%direction, E_photon = ele%value(p0c$) * 1.006)

      start_orb%field = [1, 2]

      if (debug_mode) then
        track%n_pt = -1  ! Reset
        call track1 (start_orb, ele, branch%param, end_orb, track = track)
      else
        call track1 (start_orb, ele, branch%param, end_orb)
      endif

      if (p_sign == 1) then
        out_str = trim(ele%name) // ':' // trim(tracking_method_name(j))
        if (debug_mode) out_str = trim(ele%name) // ': ' // trim(tracking_method_name(j))

      else
        out_str = trim(ele%name) // '-Anti:' // trim(tracking_method_name(j))
      endif

      if (ele%key == e_gun$) then
        write (1,fmt) '"' // trim(out_str) // '"' , tolerance(out_str), end_orb%vec, c_light * (end_orb%t - start_orb%t)
        if (debug_mode) print '(a30, 3x, 7es18.10)', out_str,  end_orb%vec, c_light * (end_orb%t - start_orb%t)
      else
        write (1,fmt) '"' // trim(out_str) // '"' , tolerance(out_str), end_orb%vec, (end_orb%vec(5) - start_orb%vec(5)) - &
                c_light * (end_orb%beta * (ele%ref_time - end_orb%t) - start_orb%beta * (ele%ref_time - ele%value(delta_ref_time$) - start_orb%t))
        if (debug_mode) print '(a30, 3x, 7es18.10)', out_str,  end_orb%vec, (end_orb%vec(5) - start_orb%vec(5)) - &
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
        write (line(isn), '(a, t49, a,  3f14.9, 4x, f14.9)') '"' // trim(out_str) // '"', tolerance_spin(out_str), &
              end_orb%spin-start_orb%spin, norm2(end_orb%spin) - norm2(start_orb%spin)
        if (debug_mode) write(line_debug(isn), '(a40, 3f14.9, 4x, f14.9)') out_str, end_orb%spin-start_orb%spin, norm2(end_orb%spin) - norm2(start_orb%spin)
      endif

      if (branch%param%particle == photon$) then
        write (1, '(3a, t49, a, 2es18.10)') '"', trim(ele%name), ':E_Field"', 'REL 1E-07', end_orb%field
      endif
    end do

    if (debug_mode) print '(t46, a, t60, a, t74, a, t91, a)', 'dSpin_x', 'dSpin_y', 'dSpin_z', 'dSpin_amp'
    do j = 1, isn
      write (1, '(a)') trim(line(j))
      if (debug_mode) print '(a)', trim(line_debug(j))
    enddo

    if (debug_mode) then
      print *
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
character(*) :: instr

  select case (instr)
    case('CRYSTAL2:Bmad_Standard')                     ; tolerance = 'ABS 1e-13'
    case('AC_KICKER2:Time_Runge_Kutta')                ; tolerance = 'ABS 1e-11'
    case('QUADRUPOLE1:Time_Runge_Kutta')               ; tolerance = 'ABS 5e-12'
    case('QUADRUPOLE2:Time_Runge_Kutta')               ; tolerance = 'ABS 2e-11'
    case('QUADRUPOLE4:Time_Runge_Kutta')               ; tolerance = 'ABS 2e-11'
    case('QUADRUPOLE5:Time_Runge_Kutta')               ; tolerance = 'ABS 1e-11'
    case('RFCAVITY1:Runge_Kutta')                      ; tolerance = 'ABS 1e-13'
    case('RFCAVITY1:Time_Runge_Kutta')                 ; tolerance = 'ABS 2e-11'
    case('RFCAVITY2:Runge_Kutta')                      ; tolerance = 'ABS 1e-13'
    case('RFCAVITY2:Time_Runge_Kutta')                 ; tolerance = 'ABS 2e-11'
    case('SBEND2:Time_Runge_Kutta')                    ; tolerance = 'ABS 1e-11'
    case('SBEND4:Bmad_Standard')                       ; tolerance = 'ABS 1e-11'
    case('SBEND4:Linear')                              ; tolerance = 'ABS 1e-11'
    case('SBEND4:Taylor')                              ; tolerance = 'ABS 1e-12'
    case('SBEND4:Time_Runge_Kutta')                    ; tolerance = 'ABS 4e-13'
    case('SBEND5:Bmad_Standard')                       ; tolerance = 'ABS 5e-13'
    case('SBEND5:Linear')                              ; tolerance = 'ABS 5e-13'
    case('SBEND6:Taylor')                              ; tolerance = 'ABS 2e-11'
    case('SBEND7:Bmad_Standard')                       ; tolerance = 'ABS 2e-13'
    case('SBEND7:Linear')                              ; tolerance = 'ABS 2e-13'
    case('SOLENOID1:Time_Runge_Kutta')                 ; tolerance = 'ABS 1e-11'
    case('SOLENOID2:Symp_Lie_Bmad')                    ; tolerance = 'ABS 2e-14'
    case('SOL_QUAD1:Time_Runge_Kutta')                 ; tolerance = 'ABS 1e-11'
    case('SOL_QUAD2:Time_Runge_Kutta')                 ; tolerance = 'ABS 2e-10'
    case('LCAVITY1:Bmad_Standard')                     ; tolerance = 'ABS 2e-12'
    case('LCAVITY1:Time_Runge_Kutta')                  ; tolerance = 'ABS 3e-11'
    case('LCAVITY1:Runge_Kutta')                       ; tolerance = 'ABS 1e-13'
    case('LCAVITY2:Time_Runge_Kutta')                  ; tolerance = 'ABS 2e-13'
    case('LCAVITY2:Runge_Kutta')                       ; tolerance = 'ABS 1e-13'
    case('LCAVITY3:Runge_Kutta')                       ; tolerance = 'ABS 1e-13'
    case('LCAVITY3:Time_Runge_Kutta')                  ; tolerance = 'ABS 2e-11'
    case('WIGGLER_MAP1:Time_Runge_Kutta')              ; tolerance = 'ABS 2e-13'
    case('WIGGLER_MAP1:Runge_Kutta')                   ; tolerance = 'ABS 1e-13'
    case('WIGGLER_PLANAR1:Bmad_Standard')              ; tolerance = 'ABS 5e-14'
    case('WIGGLER_PLANAR1:Runge_Kutta')                ; tolerance = 'ABS 1e-12'
    case('WIGGLER_PLANAR1:Time_Runge_Kutta')           ; tolerance = 'ABS 2e-12'
    case('WIGGLER_HELICAL1:Runge_Kutta')               ; tolerance = 'ABS 1e-12'
    case('WIGGLER_HELICAL1:Time_Runge_Kutta')          ; tolerance = 'ABS 2e-12'

    case("OCTUPOLE1-Anti:Runge_Kutta")                 ; tolerance = 'ABS 1e-13'
    case("LCAVITY3-Anti:Runge_Kutta")                  ; tolerance = 'ABS 2E-13'
    case("RFCAVITY1-Anti:Runge_Kutta")                 ; tolerance = 'ABS 4E-10'
    case("RFCAVITY1-Anti:Time_Runge_Kutta")            ; tolerance = 'ABS 2e-10'
    case("RFCAVITY2-Anti:Runge_Kutta")                 ; tolerance = 'ABS 4E-10'
    case("RFCAVITY2-Anti:Time_Runge_Kutta")            ; tolerance = 'ABS 2E-10'
    case("SOL_QUAD1-Anti:Symp_Lie_Bmad")               ; tolerance = 'ABS 1E-13'
    case("SOL_QUAD1-Anti:Time_Runge_Kutta")            ; tolerance = 'ABS 2e-12'
    case("SBEND4-Anti:Bmad_Standard")                  ; tolerance = 'ABS 2e-13'
    case("WIGGLER_MAP1-Anti:Runge_Kutta")              ; tolerance = 'ABS 1e-13'
    case("WIGGLER_PLANAR1-Anti:Runge_Kutta")           ; tolerance = 'ABS 5e-13'
    case("WIGGLER_PLANAR1-Anti:Bmad_Standard")         ; tolerance = 'ABS 5e-13'
    case("WIGGLER_PLANAR1-Anti:Time_Runge_Kutta")      ; tolerance = 'ABS 2e-13'                  
    case("WIGGLER_HELICAL1-Anti:Runge_Kutta")          ; tolerance = 'ABS 5e-13'
    case("WIGGLER_HELICAL1-Anti:Time_Runge_Kutta")     ; tolerance = 'ABS 2e-13'                  

    case default 
      if (index(instr, 'Runge_Kutta') /= 0) then
        tolerance = 'ABS 1e-13'
      else
        tolerance = 'ABS 1e-14'
      endif
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
