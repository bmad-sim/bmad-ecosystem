!-----------------------------------------------------------
!-----------------------------------------------------------
!-----------------------------------------------------------
!+ 
! Subroutine track1_time_runge_kutta(start_orb, ele, param, end_orb, err_flag, track)
!
! Routine to track a particle through an element using 
! Runge-Kutta time-based tracking. Converts to and from element
! coordinates before and after tracking.
!
! Note that the argument tref in is arbitrary in the function
! convert_particle_coordinates_t_to_s, so the value ele%ref_time is
! used. Consistency is what matters.
!
! Modules Needed:
!   use time_tracker_mod
!
! Input:
!   start_orb                  -- coord_struct: starting position, t-based global
!   ele                    -- ele_struct: element
!    %value                -- real(rp): attribute values
!    %ref_time             -- real(rp): time ref particle passes exit end
!    %s                    -- real(rp): longitudinal ref position at exit end
!   param                  -- lat_param_struct: lattice parameters
!
! Output:
!   end_orb     -- coord_struct: end position, t-based global
!   err_flag    -- Logical: Set True if there is an error. False otherwise
!   track       -- track_struct (optional): particle path
!-


subroutine track1_time_runge_kutta (start_orb, ele, param, end_orb, err_flag, track)

use time_tracker_mod, dummy => track1_time_runge_kutta
use track1_mod, dummy2 => track1_time_runge_kutta

implicit none

type hard_edge_struct
   real(rp) :: s                          ! S position of next hard edge in track_ele frame.
   real(rp) :: s_hard                     ! S-position of next hard edge in hard_ele frame
   type (ele_struct), pointer :: hard_ele ! Points to element with the hard edge.               
   integer  :: hard_end                   ! Integer: Describes hard edge.
end type 


type (coord_struct) :: start_orb, end_orb
type (coord_struct) :: ele_origin
type (lat_param_struct), target, intent(inout) :: param
type (ele_struct), target, intent(inout) :: ele
type (track_struct), optional :: track
type (hard_edge_struct), allocatable :: edge(:)
type (em_field_struct) :: saved_field

real(rp)  dt_step, ref_time, vec6
real(rp)   s_rel, time, s1, s2, del_s, p0c_save, s_save

integer :: i, n_edge

character(30), parameter :: r_name = 'track1_time_runge_kutta'

logical :: local_ref_frame = .true.
logical :: abs_time, err_flag, err

!---------------------------------
! Reset particle lost status

err_flag = .true.

end_orb = start_orb

!If element has zero length, skip tracking
if (ele%value(l$) .eq. 0) then
  
  !If saving tracks, allocate track array and save one point
  if ( present(track) ) then
   !Convert to global-s to local-t coordinates
    call convert_particle_coordinates_s_to_t(end_orb)
    !Tracks use vec(5) = s_rel
    end_orb%vec(5) = 0.0_rp
  if (.not. allocated(track%orb) ) call init_saved_orbit (track, 0)
  call save_a_step (track, ele, param, .false., end_orb%vec(5), end_orb, s_save)
  !Query the local field to save
  call em_field_calc (ele, param, end_orb%vec(5), time, end_orb, local_ref_frame, saved_field, .false., err_flag)
  if (err_flag) return
  track%field(track%n_pt) = saved_field    
    
    
    
  endif
  ! Reset particle to s-coordinates
  end_orb = start_orb
  if (end_orb%direction == +1) then
    end_orb%location = downstream_end$
  else
    end_orb%location = upstream_end$
  endif
  err_flag = .false.
  return
end if

! Specify initial time step.

dt_step = bmad_com%init_ds_adaptive_tracking / c_light

! Get edge array of hard edges
allocate( edge ( 2*max(ele%n_lord, 1) ) )
nullify (edge(1)%hard_ele)
n_edge = 0
do i = 1, size(edge)
  if (i > 1) edge(i)%hard_ele => edge(i-1)%hard_ele
  call calc_next_fringe_edge (ele, +1, edge(i)%s, edge(i)%hard_ele, edge(i)%s_hard, edge(i)%hard_end) 
  if (.not.  associated(edge(i)%hard_ele)) exit
  n_edge = n_edge + 1
enddo 

!------
!Convert particle to element coordinates
! Kicks and multipoles should be turned off in offset_particle

! Particle is moving forward towards the entrance
if (end_orb%direction == +1 .and. end_orb%location /= inside$) then
  call offset_particle (ele, param, set$, end_orb, set_hvkicks = .false., set_multipoles = .false., ds_pos = 0.0_rp ) 

! Interior start, reference momentum is at the end. No edge kicks are given
elseif (end_orb%location == inside$) then
  call offset_particle (ele, param, set$, end_orb, set_hvkicks = .false., set_multipoles = .false., &
                        ds_pos = end_orb%s - (ele%s - ele%value(l$)) )

! Particle is at the exit surface, should be moving backwards
elseif (end_orb%direction == -1 .and. end_orb%location /= inside$) then
  call offset_particle (ele, param, set$, end_orb, set_hvkicks = .false., set_multipoles = .false., &
                        ds_pos = end_orb%s - (ele%s - ele%value(l$)) )

else
  call out_io (s_fatal$, r_name, 'CONFUSED PARTICE ENTERING ELEMENT: ' // ele%name)
  if (global_com%exit_on_error) call err_exit

endif

! ele(s-based) -> ele(t-based)

call convert_particle_coordinates_s_to_t(end_orb)
s_rel =  end_orb%s - (ele%s - ele%value(l$) )
if ( abs(s_rel)  < bmad_com%significant_length ) then
  s_rel = 0
else if ( abs(s_rel - ele%value(l$))  < bmad_com%significant_length ) then
  s_rel = ele%value(l$)
endif

end_orb%vec(5) = s_rel

! Particle time 

time = particle_time(start_orb, ele)

!------
!Check wall

if ( wall3d_d_radius(end_orb%vec, ele) > 0 ) then

  call out_io (s_info$, r_name, "PARTICLE STARTED IN REGION OUTSIDE OF WALL: "//trim(ele%name)//", SETTING TO LOST")
  !Particle won't be tracked, so set end = start with the saved state
  end_orb = start_orb
  end_orb%state = lost$
  return

endif

if ( present(track) ) then
  if (.not. allocated(track%orb)) call init_saved_orbit (track, 1000)   !TODO: pass this from elsewhere
  call save_a_step (track, ele, param, .false., end_orb%vec(5), end_orb, s_save)
  !Query the local field to save
  call em_field_calc (ele, param, end_orb%vec(5), time, end_orb, local_ref_frame, saved_field, .false., err_flag)
  !call em_field_calc (ele, param, orb%vec(5), t_rel, orb, local_ref_frame, saved_field, .false., err_flag)
  if (err_flag) return
  track%field(track%n_pt) = saved_field     
endif

! Track through element

if (n_edge == 0) then
  ! Track whole element with no hard edges
  call odeint_bmad_time(end_orb, ele, param, 0.0_rp, ele%value(l$), time, &
                                                  dt_step, local_ref_frame, err, track)
  if (err) return

! There are hard edges. Track between edges until final exit
else 
  do    
    ! Bracket edges and kick if at edge
    s1 = 0.0_rp
    s2 = ele%value(l$)
    do i = 1, n_edge
      del_s = end_orb%vec(5) - edge(i)%s

      if (abs(del_s) < bmad_com%significant_length) then
        ! At an edge. Kick.         
        ! Get s and ref_time
        if (end_orb%direction == +1) then 
           s1 = edge(i)%s
           ref_time = edge(i)%hard_ele%ref_time - edge(i)%hard_ele%value(delta_ref_time$)
        else 
          s2 = edge(i)%s
          ref_time = edge(i)%hard_ele%ref_time
        end if
        !Convert to s-coordinates for hard edge kick, kick, and convert back
        call convert_particle_coordinates_t_to_s(end_orb, ref_time) 
        call apply_hard_edge_kick (end_orb, edge(i)%s_hard, time, edge(i)%hard_ele, ele, param, edge(i)%hard_end)
        call convert_particle_coordinates_s_to_t(end_orb)
        end_orb%vec(5) = edge(i)%s

      elseif (del_s > 0) then
        if (del_s < end_orb%vec(5) - s1) s1 = edge(i)%s  !new nearest left edge

      else  ! Must be del_s < 0
        if (del_s > end_orb%vec(5) - s2) s2 = edge(i)%s  !new nearest right edge
      endif
    enddo

    if (abs(end_orb%vec(5) - ele%value(L$)) < bmad_com%significant_length .and. end_orb%vec(6) > 0) exit
    if (abs(end_orb%vec(5) - 0.0_rp)        < bmad_com%significant_length .and. end_orb%vec(6) < 0) exit
    if (end_orb%state /= alive$) exit
    
    ! Track
    call odeint_bmad_time(end_orb, ele, param, s1, s2, time, dt_step, local_ref_frame, err, track)
    if (err) return
  
  enddo
endif   


! Cleanup  
deallocate (edge)

!------
!Convert back to s-based coordinates

if (end_orb%location /= inside$ .and. end_orb%vec(6) < 0) then
  !Particle left entrance end going backwards
  !set reference time and momentum
  end_orb%p0c = ele%value(p0c_start$)
  end_orb%direction = -1
  ref_time = ele%ref_time - ele%value(delta_ref_time$)

  !ele(t-based) -> ele(s-based)
  call convert_particle_coordinates_t_to_s(end_orb, ref_time)
  ! call apply_hard_edge_kick (end_orb, ele, param, upstream_end$)
  !unset
  call offset_particle (ele, param, unset$, end_orb, set_hvkicks = .false., set_multipoles = .false.)

elseif (end_orb%state /= alive$) then
  !Particle is lost in the interior of the element.
  !  The reference is a the end of the element
  if (end_orb%vec(6) < 0) then
    end_orb%p0c = ele%value(p0c$)
    end_orb%direction = -1
  else 
    end_orb%p0c = ele%value(p0c$)
    end_orb%direction = 1
  end if

  ref_time = ele%ref_time
  !ele(t-based) -> ele(s-based)
  call convert_particle_coordinates_t_to_s(end_orb, ref_time)
  !unset
  call offset_particle (ele, param, unset$, end_orb, set_hvkicks = .false., set_multipoles = .false., &
                          ds_pos = end_orb%s - (ele%s - ele%value(l$)) )

elseif (end_orb%location /= inside$ .and. end_orb%vec(6) >= 0) then
  !Particle left exit end going forward
  end_orb%p0c = ele%value(p0c$)
  end_orb%direction = 1
  ref_time = ele%ref_time
  !ele(t-based) -> ele(s-based)
  call convert_particle_coordinates_t_to_s(end_orb, ref_time)
  !call apply_hard_edge_kick (end_orb, ele, param, downstream_end$)
  !unset
  call offset_particle (ele, param, unset$, end_orb, set_hvkicks = .false., set_multipoles = .false.)

else
  call out_io (s_fatal$, r_name, 'CONFUSED PARTICE LEAVING ELEMENT: ' // ele%name)
  if (global_com%exit_on_error) call err_exit
endif

!Set relativistic beta
call convert_pc_to (end_orb%p0c * (1 + end_orb%vec(6)), end_orb%species, beta = end_orb%beta)
err_flag = .false.

end subroutine

