
! Fudge module to pass arguments into delta_s_target function
! Consider input arguments of: rkck_bmad_T (ele, param, orb, dvec_dt, orb%s, orb%t, new_dt, orb_new, vec_err, local_ref_frame)
module track1_time_runge_kutta_mod

use bmad_struct

type (ele_struct), save, pointer :: ele_com
type (lat_param_struct), save, pointer :: param_com
type (coord_struct), save, pointer ::  orb_com
real(rp), save, pointer :: dvec_dt_com(:)
type (coord_struct), save, pointer :: orb_new_com
real(rp), save, pointer :: vec_err_com(:)
logical, save, pointer :: local_ref_frame_com
real(rp), save, pointer :: s_target_com

contains

! function for zbrent to calculate timestep to exit face surface
function delta_s_target (new_dt)
  real(rp), intent(in)  :: new_dt
  real(rp) :: delta_s_target
  call rkck_bmad_time (ele_com, param_com, orb_com, dvec_dt_com, orb_com%s, orb_com%t, new_dt, orb_new_com, vec_err_com, local_ref_frame_com)
  delta_s_target = orb_new_com%s - s_target_com
end function delta_s_target
  
end module track1_time_runge_kutta_mod

!-----------------------------------------------------------
!-----------------------------------------------------------
!-----------------------------------------------------------
!+ 
! Subroutine track1_time_runge_kutta(start, ele, param, end, track)
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
!   start                  -- coord_struct: starting position, t-based global
!   ele                    -- ele_struct: element
!    %value                -- real(rp): attribute values
!    %ref_time             -- real(rp): time ref particle passes exit end
!    %s                    -- real(rp): longitudinal ref position at exit end
!   param                  -- lat_param_struct: lattice parameters
!    %particle             -- integer: positron$, electron$, etc.  
!
! Output:
!   end     -- coord_struct: end position, t-based global
!   track   -- track_struct (optional): particle path
!   param   -- lat_param_struct: lattice parameters
!    %particle             -- integer: positron$, electron$, etc.  
!    %lost                 -- logical: False only if particle leaves exit end
!-


subroutine track1_time_runge_kutta (start, ele, param, end, track)

use time_tracker_mod
use em_field_mod
use bmad_struct
use bmad_interface
use track1_mod

implicit none

type (coord_struct) :: start, start2, end
type (coord_struct) :: ele_origin
type (lat_param_struct), target, intent(inout) :: param
type (ele_struct), target, intent(inout) :: ele
type (track_struct), optional :: track

real(rp)  dt_step, ref_time, vec6

character(30), parameter :: r_name = 'track1_time_runge_kutta'

logical :: local_ref_frame = .true.

!---------------------------------
!Reset particle lost status
param%lost = .False.

!If element has zero length, skip tracking
if (ele%value(l$) .eq. 0) then
  end = start
  
  !If saving tracks, allocate track array and save one point
  if ( present(track) ) then
   !Convert to global-s to local-t coordinates
    call convert_particle_coordinates_s_to_t(end)
    !convert to element coordinates for track_struct
    end%t = end%t - ele%ref_time
    end%s = end%s - ele%s
    end%vec(5) = end%s
    call init_saved_orbit (track, 0)
    track%n_pt = 0
    track%orb(0) = end
    !Restore s and t to continue tracking
    end = start
  endif

  end%status = outside$

  return
end if

!copy start, because following routines will maniupulate it
start2 = start

!Specify time step; assumes ele%value(ds_step$) has been set
dt_step = ele%value(ds_step$)/c_light


!------
!Convert particle to element coordinates

if (start2%p0c > 0 .and.start2%status == outside$) then
  !Particle is moving forward towards the entrance
  call offset_particle(ele, param, start2, set$, set_canonical = .false., ds_pos = 0.0_rp ) 
  call apply_element_edge_kick (start2, ele, param, entrance_end$)

elseif (start2%status == inside$) then
  !Interior start, reference momentum is at the end. No edge kicks are given
  call offset_particle(ele, param, start2, set$, set_canonical = .false., &
                       ds_pos = start2%s - (ele%s - ele%value(l$)) )

elseif (start2%p0c < 0 .and. start2%status == outside$) then
  !Particle is at the exit surface, should be moving backwards
  call offset_particle(ele, param, start2, set$, set_canonical = .false., &
                       ds_pos = start2%s - (ele%s - ele%value(l$)) )
  call apply_element_edge_kick (start2, ele, param, exit_end$)

else
  call out_io (s_fatal$, r_name, 'CONFUSED PARTICE ENTERING ELEMENT: ' // ele%name)
  if (bmad_status%exit_on_error) call err_exit

endif

! ele(s-based) -> ele(t-based)
call convert_particle_coordinates_s_to_t(start2)
!Shift s and t to ele coordinates
start2%t = start2%t - (ele%ref_time - ele%value(delta_ref_time$))
start2%s = start2%s - (ele%s - ele%value(l$))
start2%vec(5) = start2%s


!------
!Check wall
ele_origin%vec = (/ 0.0_rp, 0.0_rp,  0.0_rp,  0.0_rp, start2%s,  0.0_rp /)
ele_origin%s = start2%s
call  particle_hit_wall_check_time(ele_origin, start2, param, ele)

if (param%lost) then

  end = start2
  end%status = dead$
   
  !Allocate track array and set value
  if ( present(track) ) then
    call init_saved_orbit (track, 0)
    track%n_pt = 0
    track%orb(0) = end
  endif

  call out_io (s_info$, r_name, "PARTICLE STARTED IN REGION OUTSIDE OF WALL, SKIPPING TRACKING")

else

  !If particle passed wall check, track through element
  call odeint_bmad_time(start2, ele, param, end, 0.0_rp, ele%value(l$), &
    dt_step, local_ref_frame, track )

end if

!------
!Convert particle to global curvilinear coordinates

!Shift s and t back to global values
start2%t = start2%t + (ele%ref_time - ele%value(delta_ref_time$))
start2%s = start2%s + (ele%s - ele%value(l$))
end%t = end%t + (ele%ref_time - ele%value(delta_ref_time$))
end%s = end%s + (ele%s - ele%value(l$))


!Convert back to s-based coordinates

if (end%status == outside$ .and. end%vec(6) < 0) then
  !Particle left entrance end going backwards
  !set reference time and momentum
  if (ele_has_constant_reference_energy(ele)) then
    end%p0c = -1*ele%value(p0c$)
  else  ! lcavity, etc.
    end%p0c = -1*ele%value(p0c_start$)
  end if
  ref_time = ele%ref_time - ele%value(delta_ref_time$)

  !ele(t-based) -> ele(s-based)
  call convert_particle_coordinates_t_to_s(end, mass_of(param%particle), ref_time)
  call apply_element_edge_kick (start2, ele, param, entrance_end$)
  !unset
  call offset_particle(ele, param, end, unset$, set_canonical = .false.)

elseif (end%status == dead$) then
    !Particle is lost in the interior of the element.
    !  The reference is a the end of the element
    if (end%vec(6) < 0) then
      end%p0c = -1*ele%value(p0c$)
    else 
      end%p0c = ele%value(p0c$)
    end if
    ref_time = ele%ref_time
    !ele(t-based) -> ele(s-based)
    call convert_particle_coordinates_t_to_s(end, mass_of(param%particle), ref_time)
    !unset
    call offset_particle(ele, param, end, unset$, set_canonical = .false., &
                                        ds_pos = end%s - (ele%s - ele%value(l$)) )

elseif (end%status == outside$ .and. end%vec(6) .ge. 0) then
  !Particle left exit end going forward
  end%p0c = ele%value(p0c$)
  ref_time = ele%ref_time
  !ele(t-based) -> ele(s-based)
  call convert_particle_coordinates_t_to_s(end, mass_of(param%particle), ref_time)
  call apply_element_edge_kick (start2, ele, param, exit_end$)
  !unset
  call offset_particle(ele, param, end, unset$, set_canonical = .false.)

else
  call out_io (s_fatal$, r_name, 'CONFUSED PARTICE LEAVING ELEMENT: ' // ele%name)
  if (bmad_status%exit_on_error) call err_exit
endif

!Set relativistic beta
call convert_pc_to (abs(end%p0c) * (1 + end%vec(6)), param%particle, beta = end%beta)


end subroutine

!-----------------------------------------------------------
!-----------------------------------------------------------
!-----------------------------------------------------------
!+
! Subroutine odeint_bmad_time (start, ele, param, end, s1, s2, &
!                             dt1, local_ref_frame, track)
! 
! Subroutine to do Runge Kutta tracking in time. This routine is adapted from Numerical
! Recipes.  See the NR book for more details.
!
!
! Modules needed:
!   use bmad
!
! Input: 
!   start   -- Coord_struct: Starting coords: (x, px, y, py, s, ps) [t-based]
!   ele     -- Ele_struct: Element to track through.
!     %tracking_method -- Determines which subroutine to use to calculate the 
!                         field. Note: BMAD does no supply em_field_custom.
!                           == custom$ then use em_field_custom
!                           /= custom$ then use em_field_standard
!   param   -- lat_param_struct: Beam parameters.
!     %enegy       -- Energy in GeV
!     %particle    -- Particle type [positron$, or electron$]
!   s1      -- Real: Starting point.
!   s2      -- Real: Ending point.
!   dt1      -- Real: Initial guess for a time step size.
!   local_ref_frame 
!           -- Logical: If True then take the 
!                input and output coordinates as being with 
!                respect to the frame of referene of the element. 
!
!   track   -- Track_struct: Structure holding the track information.
!     %save_track -- Logical: Set True if track is to be saved.
!
! Output:
!   end     -- Coord_struct: Ending coords: (x, px, y, py, s, ps).
!   track   -- Track_struct: Structure holding the track information.
!
!-

subroutine odeint_bmad_time (start, ele, param, end, s1, s2, &
                    dt1, local_ref_frame, track)
use track1_time_runge_kutta_mod
use time_tracker_mod
use em_field_mod

use nr, only: zbrent

implicit none

type (coord_struct), intent(in) :: start
type (coord_struct), intent(out) :: end
type (ele_struct) , target :: ele
type (lat_param_struct), target ::  param
type (track_struct), optional :: track

real(rp), intent(in) :: s1, s2, dt1
real(rp), parameter :: tiny = 1.0e-30_rp, edge_tol = 1e-10
real(rp) :: dt, dt_did, dt_next, s
type (coord_struct), target :: orb
type (coord_struct), target :: orb_new
real(rp), target  :: dvec_dt(6)
real(rp), target  :: vec_err(6)
real(rp), target :: s_target

integer, parameter :: max_step = 100000
integer :: n_step, n_pt, dn_save, n_save_count

logical, target :: local_ref_frame
logical :: exit_flag

character(30), parameter :: r_name = 'odeint_bmad_time'

! init
dt = dt1
orb = start

!Allocate track arrays
n_pt = max_step
if ( present(track) ) then
   call init_saved_orbit (track, n_pt)
   track%n_pt = 0
   track%orb(0) = orb
   !number of save points
   n_save_count = 0
   dn_save = max(floor(track%ds_save/c_light/dt), 1)
endif 


!Now Track
bmad_status%ok = .true.
exit_flag = .false.

do n_step = 1, max_step

  !Get initial kick vector
  !Note that orb%s and orb%t are the relative s and t for element frame
  call em_field_kick_vector_time (ele, param, orb%s, orb%t, orb, local_ref_frame, dvec_dt) 
 
  !Single Runge-Kutta step. Updates orb% vec(6), s, and t to orb_new
  call rkck_bmad_time (ele, param, orb, dvec_dt, orb%s, orb%t, dt, orb_new, vec_err, local_ref_frame)

  !Check entrance and exit faces
  if ( orb_new%s > s2 ) then
   exit_flag = .true.
   s_target = s2
   orb_new%status = outside$
   
   !Set common structures for zbrent's internal functions 
   ele_com => ele
   param_com => param
   orb_com => orb
   dvec_dt_com => dvec_dt
   orb_new_com => orb_new
   vec_err_com => vec_err
   local_ref_frame_com => local_ref_frame
   s_target_com => s_target
   !---
   dt = zbrent (delta_s_target, 0.0_rp, dt, 1d-18)

   !ensure that particle has actually exited after zbrent
   if (orb_new%s < s2) then
      orb_new%s = 2*s2 - orb_new%s
      orb_new%vec(5) = orb_new%s
   end if
   if (abs(s2 - orb_new%s) < edge_tol) then
     orb_new%s = s2 + edge_tol
     orb_new%vec(5) = orb_new%s
   end if
 
  
  else if ( orb_new%s < s1 ) then
    exit_flag = .true. 
    s_target = s1
    orb_new%status = outside$ 
    
    !Set common structures for zbrent's internal functions 
    ele_com => ele
    param_com => param
    orb_com => orb
    dvec_dt_com => dvec_dt
    orb_new_com => orb_new
    vec_err_com => vec_err
    local_ref_frame_com => local_ref_frame
    s_target_com => s_target
    !---
    dt = zbrent (delta_s_target, dt, 0.0_rp, 1d-18)

    ! ensure that particle has actually exited after zbrent
    if (orb_new%s > s1) then
      orb_new%s = 2*s1 - orb_new%s
      orb_new%vec(5) = orb_new%s
    end if
    if (abs(s1 - orb_new%s) < edge_tol) then
      orb_new%s = s1 - edge_tol
      orb_new%vec(5) = orb_new%s
    endif

  else
  !Check wall or aperture at every intermediate step
    call  particle_hit_wall_check_time(orb, orb_new, param, ele)
    !Flag for exit if wall is hit
    if (param%lost) then
      orb_new%status = dead$
      exit_flag = .true.
    end if
  endif
  
  !Update orb
  orb = orb_new

  !Save track
  if ( present(track) ) then
    n_save_count = n_save_count +1
    if (n_save_count == dn_save ) then
      track%n_pt = track%n_pt + 1
      n_pt = track%n_pt
       track%orb(n_pt) = orb
       n_save_count = 0
    end if
  endif

  !Exit when the particle hits surface s1 or s2, or hits wall
  if (exit_flag) then
     end = orb_new   !Return last orb_new that zbrent calculated
     !exit routine
     return
  endif

end do

bmad_status%ok = .false.
if (bmad_status%type_out) then
  call out_io (s_warn$, r_name, 'STEPS EXCEEDED MAX_STEP FOR ELE: '//ele%name )
  !print *, '  Skipping particle; coordinates will not be saved'
  orb_new%status = inside$
  end = orb_new     !Return last coordinate
  return
end if
if (bmad_status%exit_on_error) call err_exit

end subroutine odeint_bmad_time

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
subroutine rkck_bmad_time (ele, param, orb, dr_ds, s, t, h, orb_new, r_err, local_ref_frame)
!Very similar to rkck_bmad, except that em_field_kick_vector_time is called
!  and orb_new%s and %t are updated

use bmad
use time_tracker_mod

implicit none

type (ele_struct) ele
type (lat_param_struct) param
type (coord_struct) orb, orb_new, orb_temp

real(rp), intent(in) :: dr_ds(6)
real(rp), intent(in) :: s, t, h
real(rp), intent(out) :: r_err(6)
real(rp) :: ak2(6), ak3(6), ak4(6), ak5(6), ak6(6), r_temp(6)
real(rp), parameter :: a2=0.2_rp, a3=0.3_rp, a4=0.6_rp, &
    a5=1.0_rp, a6=0.875_rp, b21=0.2_rp, b31=3.0_rp/40.0_rp, &
    b32=9.0_rp/40.0_rp, b41=0.3_rp, b42=-0.9_rp, b43=1.2_rp, &
    b51=-11.0_rp/54.0_rp, b52=2.5_rp, b53=-70.0_rp/27.0_rp, &
    b54=35.0_rp/27.0_rp, &
    b61=1631.0_rp/55296.0_rp, b62=175.0_rp/512.0_rp, &
    b63=575.0_rp/13824.0_rp, b64=44275.0_rp/110592.0_rp, &
    b65=253.0_rp/4096.0_rp, c1=37.0_rp/378.0_rp, &
    c3=250.0_rp/621.0_rp, c4=125.0_rp/594.0_rp, &
    c6=512.0_rp/1771.0_rp, dc1=c1-2825.0_rp/27648.0_rp, &
    dc3=c3-18575.0_rp/48384.0_rp, dc4=c4-13525.0_rp/55296.0_rp, &
    dc5=-277.0_rp/14336.0_rp, dc6=c6-0.25_rp

logical local_ref_frame

!

orb_temp%vec = orb%vec +b21*h*dr_ds
call em_field_kick_vector_time(ele, param, s+a2*h, t, orb_temp, local_ref_frame, ak2)
orb_temp%vec = orb%vec +h*(b31*dr_ds+b32*ak2)
call em_field_kick_vector_time(ele, param, s+a3*h, t, orb_temp, local_ref_frame, ak3) 
orb_temp%vec = orb%vec +h*(b41*dr_ds+b42*ak2+b43*ak3)
call em_field_kick_vector_time(ele, param, s+a4*h, t, orb_temp, local_ref_frame, ak4)
orb_temp%vec = orb%vec +h*(b51*dr_ds+b52*ak2+b53*ak3+b54*ak4)
call em_field_kick_vector_time(ele, param, s+a5*h, t, orb_temp, local_ref_frame, ak5)
orb_temp%vec = orb%vec +h*(b61*dr_ds+b62*ak2+b63*ak3+b64*ak4+b65*ak5)
call em_field_kick_vector_time(ele, param, s+a6*h, t, orb_temp, local_ref_frame, ak6)
!Output new orb and error vector
orb_new%vec = orb%vec +h*(c1*dr_ds+c3*ak3+c4*ak4+c6*ak6)
orb_new%t = orb%t + h
orb_new%s = orb_new%vec(5)
r_err=h*(dc1*dr_ds+dc3*ak3+dc4*ak4+dc5*ak5+dc6*ak6)

end subroutine rkck_bmad_time


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine particle_hit_wall_check_time
!
! Subroutine to check whether particle has collided with element walls,
! and to calculate location of impact if it has
!
! Modules needed:
!   use bmad
!   use capillary_mod
!
! Input
!   orb     -- coord_struct: Previous particle coordinates  [t-based]
!   orb_new -- coord_struct: Current particle coordinates [t-based]
!   param   -- lat_param_struct: Lattice parameters
!    %particle -- integer: Type of particle
!   ele     -- ele_struct: Lattice element
!
! Output
!   param   -- lat_param_struct: Lattice parameters
!    %lost -- logical: True if orbit hit wall
!   orb_new -- coord_struct: Location of hit
!    %phase_x -- real(rp): Used to store hit angle
!-

subroutine  particle_hit_wall_check_time(orb, orb_new, param, ele)
 
use bmad_struct
use capillary_mod

implicit none

type (coord_struct) :: orb, orb_new
type (coord_struct), pointer :: old_orb, now_orb
type (lat_param_struct) :: param
type (ele_struct) :: ele
type (photon_track_struct), target :: particle
integer :: section_ix
real(rp) :: e_tot, perp(3), dummy_real
real(rp) :: edge_tol = 1e-8

!-----------------------------------------------
!Do nothing if there is no wall
if ( .not. associated(ele%wall3d%section) ) then
   return
end if

   
!Prepare coordinate structures for capillary_photon_d_radius
old_orb => particle%old%orb
now_orb => particle%now%orb
old_orb = orb
now_orb = orb_new


!If now_orb is before element, change it
!We can do this because it can't hit the element wall outside of the
!element, and these do not affect tracks
if (now_orb%s < 0) then
   now_orb%s = 0
   now_orb%vec(5) = now_orb%s
end if


!If now_orb is too close to the wall, move it edge_tol away
if (abs(capillary_photon_d_radius(particle%now, ele)) < edge_tol) then
   now_orb%vec(1) = now_orb%vec(1) + sign(edge_tol, now_orb%vec(1))
   now_orb%vec(3) = now_orb%vec(3) + sign(edge_tol, now_orb%vec(3))
end if


!Change from particle coordinates to photon coordinates
! (coord_struct to photon_coord_struct)

!Get e_tot from momentum, calculate beta_i = c*p_i / e_tot
e_tot = sqrt(orb%vec(2)**2 + orb%vec(4)**2 + orb%vec(6)**2 + mass_of(param%particle)**2)
old_orb%vec(2) = orb%vec(2) / e_tot
old_orb%vec(4) = orb%vec(4) / e_tot
old_orb%vec(6) = orb%vec(6) / e_tot

e_tot = sqrt(orb_new%vec(2)**2 + orb_new%vec(4)**2 + orb_new%vec(6)**2 + mass_of(param%particle)**2)
now_orb%vec(2) = orb_new%vec(2) / e_tot
now_orb%vec(4) = orb_new%vec(4) / e_tot
now_orb%vec(6) = orb_new%vec(6) / e_tot


!More coordinate changes
!Equations taken from track_a_capillary in capillary_mod
particle%old%energy = ele%value(e_tot$) * (1 + orb%vec(6))
particle%old%track_len = 0
!particle%old%ix_section = 1

particle%now%energy = ele%value(e_tot$) * (1 + orb_new%vec(6))
particle%now%track_len = sqrt( &
     (now_orb%vec(1) - old_orb%vec(1))**2 + &
     (now_orb%vec(3) - old_orb%vec(3))**2 + &
     (now_orb%vec(5) - old_orb%vec(5))**2)
!particle%now%ix_section = 1


!If particle hit wall, find out where
if (capillary_photon_d_radius(particle%now, ele) > 0) then

   call capillary_photon_hit_spot_calc (particle, ele)

   orb_new = now_orb

   !Calculate perpendicular to get angle of impact
   dummy_real = capillary_photon_d_radius(particle%now, ele, perp)

   !Calculate angle of impact; cos(hit_angle) = norm_photon_vec \dot perp
   !****
   !Notice we store this in orb_new%phase_x; this is so that we don't have
   !to add another argument- yeah, it's a hack.
   !****
   orb_new%phase_x = acos( &
        ((now_orb%vec(1) - old_orb%vec(1)) * perp(1) + &
        (now_orb%vec(3) - old_orb%vec(3)) * perp(2) + &
        (now_orb%vec(5) - old_orb%vec(5)) * perp(3)) / particle%now%track_len)

   !Change back from photon coords to particle coords 
   orb_new%vec(2) = now_orb%vec(2) * e_tot
   orb_new%vec(4) = now_orb%vec(4) * e_tot
   orb_new%vec(6) = now_orb%vec(6) * e_tot

   param%lost = .True.
   orb_new%status = dead$
endif

end subroutine  particle_hit_wall_check_time
