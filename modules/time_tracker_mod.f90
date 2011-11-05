module time_tracker_mod

use bmad_struct
use beam_def_struct
!use bmad_interface
!use write_lat_file_mod


contains

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!+
! Subroutine em_field_kick_vector_time (ele, param, s_rel, t_rel, orbit, local_ref_frame, dvec_dt)
!
! Subroutine to convert particle coordinates from t-based to s-based system. 
!
! Modules needed:
!   use bmad
!
! Input:
!   ele   -- coord_struct: input particle
!   param       -- real: Reference momentum. The sign indicates direction of p_s. 
!   s_rel       -- real: element coordinate system: s
!   t_rel       -- real: element coordinate system: t
!   orbit       -- coord_struct:
!                    %vec(1:6)  in t-based system
!				   				   
!   local_ref_frame
! Output:
!    dvec_dt  -- coord_struct: output particle 
!-

subroutine em_field_kick_vector_time (ele, param, s_rel, t_rel, orbit, local_ref_frame, dvec_dt)

use bmad_struct
use bmad_interface

use em_field_mod

implicit none

type (ele_struct) ele
type (lat_param_struct) param
type (em_field_struct) field

real(rp), intent(in) :: s_rel, t_rel    
type (coord_struct), intent(in) :: orbit
real(rp), intent(out) :: dvec_dt(6)

real(rp) f_bend, kappa_x, kappa_y, ds_dt, h
real(rp) vel(3), force(3)
real(rp), save :: pc, e_tot, mc2, gamma, charge, beta, ds_dt_ref, p0
real(rp), save :: pc_old = -1, particle_old = 0

logical :: local_ref_frame

character(24), parameter :: r_name = 'em_field_kick_vector_time'


! calculate the field. 
! Note that only orbit%vec(1) = x and orbit%vec(3) = y are used in em_field_calc,
!	and they coincide in both coordinate systems,
!	so we can use the 'normal' routine:
call em_field_calc (ele, param, s_rel, t_rel, orbit, local_ref_frame, field, .false.)

mc2 = mass_of(param%particle) ! Note: mc2 is in eV

charge = charge_of(param%particle) ! Note: charge is in units of |e_charge|

!Get e_tot from momentum
e_tot = sqrt( orbit%vec(2)**2 +  orbit%vec(4)**2 +  orbit%vec(6)**2 + mc2**2) 

vel(1:3) = c_light*[  orbit%vec(2),  orbit%vec(4),  orbit%vec(6) ]/ e_tot ! velocities v_x, v_y, v_s:  c*[c*p_x, c*p_y, c*p_s]/e_tot

!Set curvatures kappa_x and kappa_y
if (ele%key == sbend$) then
  if (ele%value(tilt_tot$) /= 0 .and. .not. local_ref_frame) then
    kappa_x = ele%value(g$) * cos(ele%value(tilt_tot$))
    kappa_y = ele%value(g$) * sin(ele%value(tilt_tot$))
  else
    kappa_x = ele%value(g$)
    kappa_y = 0
  endif
  h = 1 + kappa_x *  orbit%vec(1) + kappa_y *  orbit%vec(3) ! h = 1 + \kappa_x * x + \kappa_y * y
endif

! Computation for dr/dt where r(t) = [x, c*p_x, y, c*p_y, s, c*p_s]
! p_x = m c \beta_x \gamma
! p_y = m c \beta_y \gamma
! p_s = m c \beta_s \gamma
!
! h = 1 + \kappa_x * x + \kappa_y * y
!
! dx/dt   = v_x 
!
! dcp_x/dt = cp_s * v_s * h * \kappa_x + c*charge * ( Ex + v_y * Bs - h * v_s * By )
!
! dy/dt   = v_y
!
! dcp_y/dt = cp_s * v_s * h * \kappa_y + c*charge * ( Ey + h * v_s * Bx - v_x * Bs )
!
! ds/dt = v_s
!
! dcp_s/dt = -(2/h) * cp_s * ( v_x * \kappa_x + v_y * \kappa_y ) + c*(charge/h) * ( Es + v_x By - v_y Bx )
!

!Kick vector
if ( ele%key .ne. sbend$  ) then   !Straight coordinate systems have a simple Lorentz force
 
force = charge * (field%E + cross_product(vel, field%B))

dvec_dt(1) = vel(1)

dvec_dt(2) = c_light*force(1)

dvec_dt(3) = vel(2)

dvec_dt(4) = c_light*force(2)

dvec_dt(5) = vel(3)

dvec_dt(6) = c_light*force(3)

else    !Curvilinear coordinates are more complicated

dvec_dt(1) = vel(1)

dvec_dt(2) =  orbit%vec(6) * vel(3) * h * kappa_x + c_light*charge*( field%E(1) + vel(2)* field%B(3) - h*vel(3)*field%B(2) )

dvec_dt(3) = vel(2)

dvec_dt(4) =  orbit%vec(6) * vel(3) * h * kappa_y + c_light*charge*( field%E(2) + h*vel(3)* field%B(1) - vel(1)*field%B(3) )

dvec_dt(5) = vel(3)

dvec_dt(6) = -(2/h)* orbit%vec(6)*(vel(1)*kappa_x + vel(2)*kappa_y) + (c_light*charge/h)*( field%E(3) + vel(1)* field%B(2) - vel(2)*field%B(1) )

endif


end subroutine em_field_kick_vector_time


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine convert_particle_coordinates_t_to_s (particle, p0c, mc2, tref)
!
! Subroutine to convert particle coordinates from t-based to s-based system. 
!
! Modules needed:
!   use bmad
!
! Input:
!   particle   -- coord_struct: input particle
!   p0c        -- real: Reference momentum. The sign indicates direction of p_s. 
!   mc2        -- real: particle rest mass in eV
!   tref       -- real: reference time for z coordinate
! Output:
!    particle   -- coord_struct: output particle 
!-

subroutine convert_particle_coordinates_t_to_s (particle, p0c, mc2, tref)

!use bmad_struct

implicit none

type (coord_struct), intent(inout), target ::particle
real(rp), intent(in) :: p0c
real(rp), intent(in) :: mc2
real(rp), intent(in) :: tref

real(rp) :: pctot

real(rp), pointer :: vec(:)
vec => particle%vec

      !Convert t to s
      pctot = sqrt (vec(2)**2 + vec(4)**2 + vec(6)**2)
      !vec(1) = vec(1)   !this is unchanged
      vec(2) = vec(2)/abs(p0c)
      !vec(3) = vec(3)   !this is unchanged
      vec(4) = vec(4)/abs(p0c)
      vec(5) = -c_light * (pctot/sqrt(pctot**2 +mc2**2)) *  (particle%t - tref) !z \equiv -c \beta(s)  (t(s) -t_0(s)) 
      vec(6) = pctot/abs(p0c) - 1.0_rp

end subroutine convert_particle_coordinates_t_to_s

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine convert_particle_coordinates_s_to_t (particle, p0c)
!
! Subroutine to convert particle coordinates from s-based to t-based system. 
!     The sign of p0c indicates the direction of p_s
!
! Note: t coordinates are:            
!     vec(1) = x                              [m]
!     vec(2) = c*p_x = m c^2 \gamma \beta_x   [eV]
!     vec(3) = y                              [m]
!     vec(4) = c*p_y = m c^2 \gamma beta_y    [eV]
!     vec(5) = s                              [m]
!     vec(6) = c*p_s = m c^2 \gamma \beta_s   [eV]
!
! Modules needed:
!   use bmad
!
! Input:
!   particle   -- coord_struct: input particle
!   p0c        -- real: Reference momentum. The sign indicates direction of p_s 
! Output:
!    particle   -- coord_struct: output particle 
!-

subroutine convert_particle_coordinates_s_to_t (particle, p0c)

!use bmad_struct

implicit none

type (coord_struct), intent(inout), target :: particle
real(rp), intent(in) :: p0c
real(rp), pointer :: vec(:)

vec => particle%vec

      !Convert s to t
      vec(6) = p0c * sqrt( ((1+vec(6)))**2 - vec(2)**2 -vec(4)**2 )
      !vec(1) = vec(1) !this is unchanged
      vec(2) = vec(2)*abs(p0c)
      !vec(3) = vec(3) !this is unchanged
      vec(4) = vec(4)*abs(p0c)
      vec(5) = particle%s
      

end subroutine convert_particle_coordinates_s_to_t



end module
