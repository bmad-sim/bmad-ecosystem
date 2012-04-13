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
!    %particle             -- integer: positron$, electron$, etc.  
!
! Output:
!   end_orb     -- coord_struct: end position, t-based global
!   err_flag    -- Logical: Set True if there is an error. False otherwise
!   track       -- track_struct (optional): particle path
!   param       -- lat_param_struct: lattice parameters
!    %particle             -- integer: positron$, electron$, etc.  
!    %lost                 -- logical: False only if particle leaves exit end
!-


subroutine track1_time_runge_kutta (start_orb, ele, param, end_orb, err_flag, track)

use time_tracker_mod
use em_field_mod
use bmad_struct
use bmad_interface
use track1_mod

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

real(rp)  dt_step, ref_time, vec6, pc2
real(rp)   s_rel, t_rel, s1, s2, del_s, p0c_save

integer :: i, n_edge

character(30), parameter :: r_name = 'track1_time_runge_kutta'

logical :: local_ref_frame = .true.
logical :: abs_time, err_flag, err

!---------------------------------
! Reset particle lost status

err_flag = .true.
param%lost = .False.

end_orb = start_orb

!If element has zero length, skip tracking
if (ele%value(l$) .eq. 0) then
  
  !If saving tracks, allocate track array and save one point
  if ( present(track) ) then
   !Convert to global-s to local-t coordinates
    call convert_particle_coordinates_s_to_t(end_orb)
    !convert to element coordinates for track_struct
    call init_saved_orbit (track, 0)
    track%n_pt = 0
    track%orb(0) = end_orb
  endif
  ! Reset particle to s-coordinates
  end_orb = start_orb
  end_orb%status = outside$
  err_flag = .false.
  return
end if

!Specify time step; assumes ele%value(ds_step$) has been set
dt_step = ele%value(ds_step$)/c_light


! Get edge array of hard edges
allocate( edge ( 2*max(ele%n_lord, 1) ) )
nullify (edge(1)%hard_ele)
n_edge = 0
do i = 1, size(edge)
  if (i > 1) edge(i)%hard_ele => edge(i-1)%hard_ele
  call calc_next_hard_edge (ele, edge(i)%s, edge(i)%hard_ele, edge(i)%s_hard, edge(i)%hard_end) 
  if (.not.  associated(edge(i)%hard_ele)) exit
  n_edge = n_edge + 1
enddo 


!------
!Convert particle to element coordinates

! Particle is moving forward towards the entrance
if (end_orb%p0c > 0 .and. end_orb%status == outside$) then
  call offset_particle(ele, param, end_orb, set$, set_canonical = .false., ds_pos = 0.0_rp ) 

! Interior start, reference momentum is at the end. No edge kicks are given
elseif (end_orb%status == inside$) then
  call offset_particle(ele, param, end_orb, set$, set_canonical = .false., &
                       ds_pos = end_orb%s - (ele%s - ele%value(l$)) )

! Particle is at the exit surface, should be moving backwards
elseif (end_orb%p0c < 0 .and. end_orb%status == outside$) then
  call offset_particle(ele, param, end_orb, set$, set_canonical = .false., &
                       ds_pos = end_orb%s - (ele%s - ele%value(l$)) )

else
  call out_io (s_fatal$, r_name, 'CONFUSED PARTICE ENTERING ELEMENT: ' // ele%name)
  if (bmad_status%exit_on_error) call err_exit

endif

! ele(s-based) -> ele(t-based)
call convert_particle_coordinates_s_to_t(end_orb)
s_rel =  end_orb%s - (ele%s - ele%value(l$) )
end_orb%vec(5) = s_rel

! Particle time

t_rel = particle_time(start_orb, ele)

!------
!Check wall
ele_origin%vec = [0.0_rp, 0.0_rp,  0.0_rp,  0.0_rp, s_rel,  0.0_rp ]
call  particle_hit_wall_check_time(ele_origin, end_orb, param, ele)

if (param%lost) then

  end_orb = start_orb
  end_orb%status = dead$
   
  !Allocate track array and set value
  if ( present(track) ) then
    call init_saved_orbit (track, 0)
    track%n_pt = 0
    track%orb(0) = end_orb
  endif

  call out_io (s_info$, r_name, "PARTICLE STARTED IN REGION OUTSIDE OF WALL, SKIPPING TRACKING")

else

  if ( present(track) ) then
    call init_saved_orbit (track, 10000)   !TODO: pass this from elsewhere
    track%n_pt = 0
    track%orb(0) = end_orb
  endif


  !If particle passed wall check, track through element
  
  if (n_edge == 0) then
    ! Track whole element with no hard edges
    call odeint_bmad_time(end_orb, ele, param, 0.0_rp, ele%value(l$), t_rel, &
                                      dt_step, local_ref_frame, err, track)
    if (err) return
  else 
    ! There are hard edges. Track between edges until final exit
    do    
      ! Bracket edges and kick if at edge
      s1 = 0.0_rp
      s2 = ele%value(l$)
      do i = 1, n_edge
        del_s = end_orb%vec(5) - edge(i)%s

        if (abs(del_s) < bmad_com%significant_length) then
          ! At an edge. Kick.         
          p0c_save = end_orb%p0c ! Fudge to kick orb in time coordinates by setting p0c = +/- 1
          end_orb%p0c = sign(1.0_rp, p0c_save)
          pc2 = end_orb%vec(2)**2 + end_orb%vec(4)**2 + end_orb%vec(6)**2
          call apply_hard_edge_kick (end_orb, edge(i)%s_hard, t_rel, edge(i)%hard_ele, ele, param, edge(i)%hard_end)
          end_orb%p0c = p0c_save
          end_orb%vec(6) = sign(sqrt(pc2 - end_orb%vec(2)**2 - end_orb%vec(4)**2), end_orb%vec(6))
          if (p0c_save > 0 ) then 
            s1 = edge(i)%s
          else 
            s2 = edge(i)%s
          end if

        elseif (del_s > 0) then
      	  if (del_s < end_orb%vec(5) - s1) s1 = edge(i)%s  !new nearest left edge

        else  ! Must be del_s < 0
          if (del_s > end_orb%vec(5) - s2) s2 = edge(i)%s  !new nearest right edge


        endif
      enddo
      if (abs(end_orb%vec(5) - ele%value(L$)) < bmad_com%significant_length .and. end_orb%vec(6) > 0) exit
      if (abs(end_orb%vec(5) - 0.0_rp)        < bmad_com%significant_length .and. end_orb%vec(6) < 0) exit
      if ( end_orb%status == dead$) exit
      
      ! Track
      call odeint_bmad_time(end_orb, ele, param, s1, s2, t_rel, dt_step, local_ref_frame, err, track)
      if (err) return
    enddo
  endif   
  
end if

! Cleanup  
deallocate (edge)

!------
!Convert back to s-based coordinates

if (end_orb%status == outside$ .and. end_orb%vec(6) < 0) then
  !Particle left entrance end going backwards
  !set reference time and momentum
  end_orb%p0c = -1*ele%value(p0c_start$)
  ref_time = ele%ref_time - ele%value(delta_ref_time$)

  !ele(t-based) -> ele(s-based)
  call convert_particle_coordinates_t_to_s(end_orb, mass_of(param%particle), ref_time)
 ! call apply_hard_edge_kick (end_orb, ele, param, entrance_end$)
  !unset
  call offset_particle(ele, param, end_orb, unset$, set_canonical = .false.)

elseif (end_orb%status == dead$) then
    !Particle is lost in the interior of the element.
    !  The reference is a the end of the element
    if (end_orb%vec(6) < 0) then
      end_orb%p0c = -1*ele%value(p0c$)
    else 
      end_orb%p0c = ele%value(p0c$)
    end if
    ref_time = ele%ref_time
    !ele(t-based) -> ele(s-based)
    call convert_particle_coordinates_t_to_s(end_orb, mass_of(param%particle), ref_time)
    !unset
    call offset_particle(ele, param, end_orb, unset$, set_canonical = .false., &
                                        ds_pos = end_orb%s - (ele%s - ele%value(l$)) )

elseif (end_orb%status == outside$ .and. end_orb%vec(6) .ge. 0) then
  !Particle left exit end going forward
  end_orb%p0c = ele%value(p0c$)
  ref_time = ele%ref_time
  !ele(t-based) -> ele(s-based)
  call convert_particle_coordinates_t_to_s(end_orb, mass_of(param%particle), ref_time)
  !call apply_hard_edge_kick (end_orb, ele, param, exit_end$)
  !unset
  call offset_particle(ele, param, end_orb, unset$, set_canonical = .false.)

else
  call out_io (s_fatal$, r_name, 'CONFUSED PARTICE LEAVING ELEMENT: ' // ele%name)
  if (bmad_status%exit_on_error) call err_exit
endif

!Set relativistic beta
call convert_pc_to (abs(end_orb%p0c) * (1 + end_orb%vec(6)), param%particle, beta = end_orb%beta)


end subroutine

