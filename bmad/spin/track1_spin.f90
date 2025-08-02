!+
! subroutine track1_spin (start_orb, ele, param, end_orb, make_quaternion)
!
! Particle spin tracking through a single element.
!
! Typically this routine should not be directly called. 
! Instead, use track1 which calls this routine.
!
! Input :
!   start_orb        -- Coord_struct: Starting coords.
!   ele              -- Ele_struct: Element to track through.
!   param            -- lat_param_struct: Beam parameters.
!   end_orb          -- Coord_struct: Ending coords.
!     %vec                -- Ending particle position needed for spin_integration spin tracking.
!   make_quaternion  -- logical, optional: If present and true then calculate the 1st
!                          order spin map which is represented as a quaternion.
!
! Output:
!   ele              -- ele_struct: Element to track through
!      %spin_q            -- 1st order spin map made if make_quaternion = True.
!   end_orb          -- Coord_struct: Ending coords.
!      %spin(2)           -- complex(rp): Ending spin
!-

subroutine track1_spin (start_orb, ele, param, end_orb, make_quaternion)

use bmad_routine_interface, dummy => track1_spin
  
implicit none

type (coord_struct) :: start_orb, end_orb, temp_orb
type (ele_struct) :: ele
type (em_field_struct) field
type (lat_param_struct) :: param

real(rp) fc, dxp, dyp, om(3), quat(0:3)
integer method
character(*), parameter :: r_name = 'track1_spin'
logical, optional :: make_quaternion
logical err

! Use spin_integration if spin_tracking_method = tracking$ and particle tracking is not using a RK integration method.

if (start_orb%species == photon$) return

if ((ele%key == drift$ .or. ele%key == marker$ .or. ele%key == fixer$) .and. &
                                                       ele%spin_tracking_method /= custom$) then
  end_orb%spin = start_orb%spin
  if (logic_option(.false., make_quaternion)) ele%spin_q(:,0) = [1, 0, 0, 0]
  return
endif

method = ele%spin_tracking_method
if (method == tracking$) then
  select case (ele%tracking_method)
  case (runge_kutta$, time_runge_kutta$, symp_lie_ptc$, custom$, fixed_step_runge_kutta$, fixed_step_time_runge_kutta$)
    return ! Spin tracking is done at the same time orbital tracking is done
  case (taylor$)
    method = taylor$
  case default
    method = spin_integration$
    if (ele%key == taylor$) method = taylor$
  end select
endif

!

select case (method)
case (off$)
  end_orb%spin = start_orb%spin

case (transverse_kick$)
  fc = start_orb%p0c / charge_of(start_orb%species)
  dxp = (end_orb%vec(2) - start_orb%vec(2)) * fc
  dyp = (end_orb%vec(4) - start_orb%vec(4)) * fc
  field%E = 0
  field%B = [dyp, -dxp, 0.0_rp] / c_light
  om = spin_omega (field, end_orb, +1)
  quat = omega_to_quat(om)
  end_orb%spin = quat_rotate(quat, start_orb%spin)


case (spin_integration$, magnus$)
  call track1_spin_integration (start_orb, ele, param, end_orb)

case (custom$)
  if (.not. associated(track1_spin_custom_ptr)) then
    call out_io (s_error$, r_name, 'TRACK1_SPIN_CUSTOM_PTR HAS NOT BEEN SET IN THIS PROGRAM!', &
                                   'NEEDED FOR CUSTOM ELEMENT OR TRACKING_METHOD = CUSTOM: ' // ele%name)
    end_orb%state = lost$
    return
  endif

  call track1_spin_custom_ptr (start_orb, ele, param, end_orb, err, make_quaternion)

! Notice that PTC spin tracking is only done here only when the (orbital) tracking_method is *not* symp_lie_ptc
case (symp_lie_ptc$)
  temp_orb = start_orb
  call track1_symp_lie_ptc (temp_orb, ele, param)
  end_orb%spin = temp_orb%spin

case (taylor$, sprint$)
  call track1_spin_taylor (start_orb, ele, param, end_orb)

case default
  call out_io (s_fatal$, r_name, 'BAD SPIN_TRACKING_METHOD: ' // spin_tracking_method_name(ele%spin_tracking_method), &
                                 'FOR ELEMENT: ', ele%name)
  if (global_com%exit_on_error) call err_exit
end select

end subroutine track1_spin


