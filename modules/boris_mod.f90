!+
! Module boris_mod
!
! Module to do Boris integration tracking.
!
! Reference: 
!  "Efficiency of a Boris-like Integration Scheme with Spatial Stepping", 
!   P. H. Stolz et al., Physical Review Special Topics.
!   5, 094001 (2002).
!
! When comparing the paper to this module remember that the paper uses a
! coordinate system:
!   (x, P_x, y, P_y, ct, U)
! While the Bmad coordinate system is:
!   (x, P_x/P_0, y, P_y/P_0, s-beta*c*t, dP/P_0)
!-

#include "CESR_platform.inc"      

module boris_mod

  use em_field_mod

contains

!-----------------------------------------------------------
!-----------------------------------------------------------
!-----------------------------------------------------------
!+
! Subroutine track1_adaptive_boris (start, ele, param, end, track, s_start, s_end)
! 
! Subroutine to do Boris tracking with adaptive step size control.
! This routine is adapted from odeint in Numerical Recipes. 
! See the NR book for more details.
!
! For more information on Boris tracking see the boris_mod documentation.
!
! Note:
!   For each step the error in the orbit must be:
!     error < (|orbit|*rel_tol + abs_tol) / sqrt(N)
! Where N is the number of steps that would be needed at the 
! present step size.
!
! Modules needed:
!   use bmad
!
! Input: 
!   start    -- Coord_struct: Starting coords.
!   ele      -- Ele_struct: Element to track through.
!     %tracking_method -- Determines which subroutine to use to calculate the 
!                         field. Note: BMAD does no supply em_field_custom.
!                           == custom$ then use em_field_custom
!                           /= custom$ then use em_field
!     %value(rel_tol$) -- Real: Relative error tollerance.
!                           Default if zero: 1e-6.
!     %value(abs_tol$) -- Real: Absolute error tollerance.
!                           Default if zero: 1e-7.
!   param    -- Param_struct: Beam parameters.
!     %particle    -- Particle type [positron$, or electron$]
!   track    -- Track_struct: Structure holding the track information.
!     %save_track -- Logical: Set True if track is to be saved.
!   s_start  -- Real, optional: Starting point.
!   s_end    -- Real, optional: Ending point.
!
! Output:
!   end        -- Coord_struct: Ending coords.
!   track      -- Track_struct: Structure holding the track information.
!-

subroutine track1_adaptive_boris (start, ele, param, end, track, s_start, s_end)

  implicit none

  type (coord_struct), intent(in) :: start
  type (coord_struct), intent(out) :: end
  type (ele_struct) ele, loc_ele
  type (param_struct) param
  type (coord_struct) here, orb1, orb2
  type (track_struct) :: track

  real(rp), optional, intent(in) :: s_start, s_end
  real(rp) :: ds, s, s_sav, rel_tol, abs_tol, sqrt_N
  real(rp), parameter :: err_5 = 0.0324, safety = 0.9
  real(rp) :: s1, s2, scale_orb, err_max, ds_temp, rel_tol_N, abs_tol_N

  integer :: n_step

! init

  if (present(s_start)) then
    s1 = s_start
  else
    s1 = 0
  endif

  if (present(s_end)) then
    s2 = s_end
  else 
    s2 = ele%value(l$)
  endif

  if (ele%value(rel_tol$) == 0) then
    rel_tol = 1e-6 
  else
    rel_tol = ele%value(rel_tol$) 
  endif

  if (ele%value(abs_tol$) == 0) then
    abs_tol = 1e-7 
  else
    abs_tol = ele%value(abs_tol$) 
  endif

  s = s1
  ds = sign(track%step0, s2-s1)

  call transfer_ele (ele, loc_ele)
  loc_ele%value(x_offset$) = 0
  loc_ele%value(y_offset$) = 0
  loc_ele%value(s_offset$) = 0
  loc_ele%value(x_pitch$)  = 0
  loc_ele%value(y_pitch$)  = 0
  loc_ele%value(tilt$)     = 0

  here = start
  call boris_energy_correction (ele, param, here)
  call offset_particle (ele, param, here, set$, set_canonical = .false.)
  call track_solenoid_edge (loc_ele, param, set$, here)

! if we are saving the trajectory then allocate enough space in the arrays

  if (track%save_track) then
    s_sav = s - 2 * track%ds_save
    call allocate_saved_orbit (track, int(abs((s2-s1)/track%ds_save))+1)
  endif

! now track

  bmad_status%ok = .true.

  do n_step = 1, track%max_step

    sqrt_N = sqrt(abs((s2-s1)/ds))  ! N = estimated number of steps
    rel_tol_N = rel_tol / sqrt_N
    abs_tol_N = abs_tol / sqrt_N

! record a track if we went far enough.

    if (track%save_track .and. (abs(s-s_sav) > track%ds_save)) &
                    call save_a_step (track, loc_ele, param, s, here%vec, s_sav)

    if ((s+ds-s2)*(s+ds-s1) > 0.0) ds = s2-s

! Make A step. Keep shrinking the step until the error is within bounds.
! The error in a step is estimated by the difference in making one whole step
! or two half steps.

    do

      call track1_boris_partial (here, loc_ele, param, s, ds/2, orb2) 
      call track1_boris_partial (orb2, loc_ele, param, s+ds/2, ds/2, orb2)
      call track1_boris_partial (here, loc_ele, param, s, ds, orb1) 
      scale_orb = maxval((abs(orb1%vec) + abs(orb2%vec))) / 2

      err_max = maxval(abs(orb2%vec - orb1%vec) / &
                                       (scale_orb*rel_tol_N + abs_tol_N))
      if (err_max <= 1) exit

      ds_temp = safety * ds / sqrt(err_max)
      ds = sign(max(abs(ds_temp), 0.1_rp*abs(ds)), ds)

      if (abs(ds) < track%step_min) then
        bmad_status%ok = .false.
        if (bmad_status%type_out) print *, &
            'ERROR IN TRACK1_ADAPTIVE_BORIS: STEPSIZE SMALLER THAN MINIMUM.' 
        if (bmad_status%exit_on_error) call err_exit
        return
      endif

    enddo

! now that we have a good step record the present position and 
! calculate a new step size.

    here = orb2
    s = s + ds

    if (err_max > err_5) then   ! limit increase to no more than a factor of 5
      ds = safety * ds / sqrt(err_max)
    else
      ds = 5 * ds   
    endif

! check if we are done

    if ((s-s2)*(s2-s1) >= 0.0) then
      if (track%save_track) call save_a_step (track, loc_ele, param, s, here%vec, s_sav)
      call track_solenoid_edge (loc_ele, param, unset$, here)
      call offset_particle (ele, param, here, unset$, set_canonical = .false.)
      end = here
      return
    end if

  end do

  bmad_status%ok = .false.
  if (bmad_status%type_out) &
               print *, 'ERROR IN TRACK1_ADAPTIVE_BORIS: TOO MANY STEPS'
  if (bmad_status%exit_on_error) call err_exit

end subroutine track1_adaptive_boris

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine track1_boris (start, ele, param, end, track, s_start, s_end)
!
! Subroutine to do Boris tracking.
! For more information on Boris tracking see the boris_mod documentation.
! 
! Modules needed:
!   use bmad
!
! Input: 
!   start    -- Coord_struct: Starting coords.
!   ele      -- Ele_struct: Element to track through.
!     %tracking_method -- Determines which subroutine to use to calculate the 
!                         field. Note: BMAD does no supply em_field_custom.
!                           == custom$ then use em_field_custom
!                           /= custom$ then use em_field
!     %num_steps      -- number of steps to take
!   param    -- Param_struct: Beam parameters.
!     %particle    -- Particle type [positron$, or electron$]
!   track      -- Track_struct: Structure holding the track information.
!     %save_track -- Logical: Set True if track is to be saved.
!   s_start  -- Real, optional: Starting point.
!   s_end    -- Real, optional: Ending point.
!
! Output:
!   end        -- Coord_struct: Ending coords.
!   track      -- Track_struct: Structure holding the track information.
!     %save_track -- Logical: Set True if track is to be saved.
!-

subroutine track1_boris (start, ele, param, end, track, s_start, s_end)

  implicit none

  type (coord_struct), intent(in) :: start
  type (coord_struct), intent(out) :: end
  type (ele_struct) ele, loc_ele
  type (param_struct) param
  type (track_struct) track
  type (coord_struct) here

  real(rp), optional, intent(in) :: s_start, s_end
  real(rp) s1, s2, s_sav, ds, s

  integer i

! init

  if (present(s_start)) then
    s1 = s_start
  else
    s1 = 0
  endif

  if (present(s_end)) then
    s2 = s_end
  else 
    s2 = ele%value(l$)
  endif

  if (ele%num_steps == 0) then
    print *, 'ERROR IN TRACK1_BORIS: ELEMENT HAS ZERO "NUM_STEPS": ', ele%name
    call err_exit
  endif

  s = s1
  ds = (s2-s1) / ele%num_steps

! go to local coords

  call transfer_ele (ele, loc_ele)
  loc_ele%value(x_offset$) = 0
  loc_ele%value(y_offset$) = 0
  loc_ele%value(s_offset$) = 0
  loc_ele%value(x_pitch$)  = 0
  loc_ele%value(y_pitch$)  = 0
  loc_ele%value(tilt$)     = 0

  here = start
  call boris_energy_correction (ele, param, here)
  call offset_particle (ele, param, here, set$, set_canonical = .false.)
  call track_solenoid_edge (loc_ele, param, set$, here)

! if we are saving the trajectory then allocate enough space in the arrays

  if (track%save_track) then
    s_sav = s - 2.0_rp * track%ds_save
    call allocate_saved_orbit (track, loc_ele%num_steps+1)
    call save_a_step (track, loc_ele, param, s, here%vec, s_sav)
  endif

! track through the body

  do i = 1, loc_ele%num_steps
    call track1_boris_partial (here, loc_ele, param, s, ds, here)
    s = s + ds
    if (track%save_track) call save_a_step (track, loc_ele, param, s, here%vec, s_sav)
  enddo

! back to lab coords

  call track_solenoid_edge (loc_ele, param, unset$, here)
  call offset_particle (ele, param, here, unset$, set_canonical = .false.)

  end = here

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine track1_boris_partial (start, ele, param, s, ds, end)
!
! Subroutine to track 1 step using boris tracking.
! This subroutine is used by track1_boris and track1_adaptive_boris.
!
! Note: Coordinates are with respect to the element coordinate frame.
!
! For more information on Boris tracking see the boris_mod documentation.
!
! Modules needed:
!   use bmad
!
! Input:
!   start -- Coord_struct: Starting coordinates.
!   ele   -- Ele_struct: Element that we are tracking through.
!   param -- Param_struct: 
!     %particle    -- Particle type [positron$, electron$, etc.]
!   s     -- Real(rp): Starting point relative to element beginning.
!   ds    -- Real(rp): step size
!
! Output:
!   end   -- Coord_struct: Ending coordinates.
!-

subroutine track1_boris_partial (start, ele, param, s, ds, end)

  implicit none

  type (ele_struct), intent(in) :: ele
  type (param_struct) param
  type (coord_struct) :: start, end  
  type (em_field_struct) :: field

  real(rp), intent(in) :: s, ds
  real(rp) :: f, p_z, d2, alpha, dxv, dyv, ds2_f, charge, U_rel, p_rel
  real(rp) :: r(3,3), w(3), ex, ey, ex2, ey2, exy, bz, bz2, m_rel, beta

!

  charge = charge_of(param%particle)
  m_rel = mass_of(param%particle) / ele%value(p0c$)

  end = start

! 1) Push the position 1/2 step

  p_rel = 1 + end%vec(6)
  p_z = sqrt(p_rel**2 - end%vec(2)**2 - end%vec(4)**2)
  ds2_f = ds / (2 * p_z)
  U_rel = sqrt (p_rel**2 + m_rel**2)
  beta = p_rel / U_rel  ! particle velocity: v/c

  end%vec(1) = end%vec(1) + ds2_f * end%vec(2) 
  end%vec(3) = end%vec(3) + ds2_f * end%vec(4)
  end%vec(5) = end%vec(5) + ds2_f * (p_z - p_rel) 

! 2) Evaluate the fields.

  call em_field (ele, param, s+ds/2, end, field)
  if (field%type /= em_field$) then
    print *, 'ERROR IN TRACK1_BORIS_PARTIAL: BORIS CAN ONLY TRACK WITH EM FIELDS.'
    print *, '      FOR ELEMENT: ', ele%name
    call err_exit
  endif

! 3) Push the momenta a 1/2 step using only the "b" term.

  f = ds * charge * c_light / (2 * ele%value(p0c$))

  end%vec(2) = end%vec(2) - field%b(2) * f
  end%vec(4) = end%vec(4) + field%b(1) * f
  U_rel = U_rel + field%e(3) * f / c_light
  p_rel = sqrt (U_rel**2 - m_rel**2)
  p_z = sqrt(p_rel**2 - end%vec(2)**2 - end%vec(4)**2)

! 4) Push the momenta a full step using the "R" matrix.

  d2 = ds * charge * c_light / (2 * p_z * ele%value(p0c$)) 

  if (field%e(1) == 0 .and. field%e(2) == 0) then
    if (field%b(3) /= 0) then
      d2 = d2 * field%b(3)
      alpha = 2 * d2 / (1 + d2**2)
      dxv = -d2 * end%vec(2) + end%vec(4)
      dyv = -end%vec(2) - d2 * end%vec(4)
      end%vec(2) = end%vec(2) + alpha * dxv
      end%vec(4) = end%vec(4) + alpha * dyv
    endif
  else
    ex = field%e(1) / c_light;     ex2 = ex**2
    ey = field%e(2) / c_light;     ey2 = ey**2
    bz = field%b(3);               bz2 = bz**2
    exy = ex * ey
    alpha = 2 * d2 / (1 + d2**2 * (bz2 - ex2 - ey2))
    r(1,1:3) = (/ d2 * (ex2 - bz2), bz + d2*exy,      ex + d2*bz*ey    /)
    r(2,1:3) = (/ -bz + d2*exy,     d2 * (ey2 - bz2), ey - d2*bz*ex    /)
    r(3,1:3) = (/ ex - d2*bz*ey,    ey + d2*bz*ex,    d2 * (ex2 + ey2) /)

    w = (/ end%vec(2), end%vec(4), U_rel /)
    w = w + alpha * matmul(r, w)
    end%vec(2:4:2) = w(1:2); U_rel = w(3)
  endif

! 5) Push the momenta a 1/2 step using only the "b" term.

  end%vec(2) = end%vec(2) - field%b(2) * f
  end%vec(4) = end%vec(4) + field%b(1) * f
  U_rel = U_rel + field%e(3) * f / c_light
  p_rel = sqrt (U_rel**2 - m_rel**2)

! 6) Push the position 1/2 step.

  p_z = sqrt(p_rel**2 - end%vec(2)**2 - end%vec(4)**2)
  ds2_f = ds / (2 * p_z)

  end%vec(1) = end%vec(1) + ds2_f * end%vec(2) 
  end%vec(3) = end%vec(3) + ds2_f * end%vec(4)
  end%vec(5) = end%vec(5) + ds2_f * (p_z - p_rel) 
  end%vec(6) = p_rel - 1

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine boris_energy_correction (ele, param, here)
!
! Subroutine to correct the orbit due to a change in the reference energy
! from the start of the element to the end.
!
! Input:
!   ele   -- Ele_struct: Element being tracked through.
!   param -- Param_struct:
!
! Output:
!   here  -- Coord_struct: Coordinates to correct.
!-

subroutine boris_energy_correction (ele, param, here)

  implicit none

  type (ele_struct), intent(in) :: ele
  type (param_struct), intent(in) :: param
  type (coord_struct) :: here

  real(rp) p0, p1
!  character(24) :: r_name = 'boris_energy_correction'

!

  select case (ele%key)
  case (lcavity$) 
    call energy_to_kinetic (ele%value(energy_start$), param%particle, p0c = p0)
    call energy_to_kinetic (ele%value(beam_energy$), param%particle, p0c = p1)
    here%vec(2) = here%vec(2) * p0 / p1
    here%vec(4) = here%vec(4) * p0 / p1
    here%vec(6) = ((1 + here%vec(6)) * p0 - p1) / p1

  case (custom$)
    call energy_to_kinetic (ele%value(beam_energy$)-ele%value(delta_e$), param%particle, p0c = p0)
    call energy_to_kinetic (ele%value(beam_energy$), param%particle, p0c = p1)
    here%vec(2) = here%vec(2) * p0 / p1
    here%vec(4) = here%vec(4) * p0 / p1
    here%vec(6) = ((1 + here%vec(6)) * p0 - p1) / p1

  end select

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine track_solenoid_edge (ele, param, set, orb)

  implicit none

  type (ele_struct), intent(in) :: ele
  type (param_struct), intent(in) :: param
  type (coord_struct) :: orb

  logical, intent(in) :: set

!

  if (ele%key /= solenoid$ .and. ele%key /= sol_quad$) return

  if (set) then
    orb%vec(2) = orb%vec(2) + orb%vec(3) * ele%value(ks$) / 2
    orb%vec(4) = orb%vec(4) - orb%vec(1) * ele%value(ks$) / 2
  else
    orb%vec(2) = orb%vec(2) - orb%vec(3) * ele%value(ks$) / 2
    orb%vec(4) = orb%vec(4) + orb%vec(1) * ele%value(ks$) / 2
  endif

end subroutine

end module
