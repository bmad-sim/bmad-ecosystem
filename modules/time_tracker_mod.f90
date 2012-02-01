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
! Subroutine convert_particle_coordinates_t_to_s (particle, mc2, tref)
!
! Subroutine to convert particle coordinates from t-based to s-based system. 
!
! Modules needed:
!   use bmad
!
! Input:
!   particle   -- coord_struct: input particle coordinates
!                    %vec(:)
!                    %t 
!                    %p0c
!   mc2        -- real: particle rest mass in eV
!   tref       -- real: reference time for z coordinate
! Output:
!    particle   -- coord_struct: output particle 
!-

subroutine convert_particle_coordinates_t_to_s (particle, mc2, tref)

!use bmad_struct

implicit none

type (coord_struct), intent(inout), target ::particle
real(rp) :: p0c
real(rp), intent(in) :: mc2
real(rp), intent(in) :: tref

real(rp) :: pctot

real(rp), pointer :: vec(:)
vec => particle%vec
p0c=abs(particle%p0c)

      !Convert t to s
      pctot = sqrt (vec(2)**2 + vec(4)**2 + vec(6)**2)
      !vec(1) = vec(1)   !this is unchanged
      vec(2) = vec(2)/p0c
      !vec(3) = vec(3)   !this is unchanged
      vec(4) = vec(4)/p0c
      vec(5) = -c_light * (pctot/sqrt(pctot**2 +mc2**2)) *  (particle%t - tref) !z \equiv -c \beta(s)  (t(s) -t_0(s)) 
      vec(6) = pctot/p0c - 1.0_rp

end subroutine convert_particle_coordinates_t_to_s

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine convert_particle_coordinates_s_to_t (particle)
!
! Subroutine to convert particle coordinates from s-based to t-based system. 
!     The sign of particle%p0c indicates the direction of p_s
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
!                       %vec(2), %vec(4), %vec(6)
!                       %s, %p0c
!   p0c        -- real: Reference momentum. The sign indicates direction of p_s 
! Output:
!    particle   -- coord_struct: output particle 
!-

subroutine convert_particle_coordinates_s_to_t (particle)

!use bmad_struct

implicit none

type (coord_struct), intent(inout), target :: particle
real(rp), pointer :: vec(:)

vec => particle%vec

      !Convert s to t
      vec(6) = particle%p0c * sqrt( ((1+vec(6)))**2 - vec(2)**2 -vec(4)**2 )
      !vec(1) = vec(1) !this is unchanged
      vec(2) = vec(2)*abs(particle%p0c)
      !vec(3) = vec(3) !this is unchanged
      vec(4) = vec(4)*abs(particle%p0c)
      vec(5) = particle%s
      

end subroutine convert_particle_coordinates_s_to_t


!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine drift_orbit_time(orbit, mc2, delta_s)
!
! Simple routine to drift a particle orbit in time-based coordinates by a distance delta_s
!
! Modules Needed:
!   use bmad_struct
!
! Input:
!   orbit      -- coord_struct: particle orbit in time-based coordinates
!   mc2        -- real(rp): particle mass in eV
!   delta_s    -- real(rp): s-coordinate distance to drift particle
!                  .
!
! Output:
!   orbit      -- coord_struct: particle orbit in time-based coordinates
!                                     
!-
subroutine drift_orbit_time(orbit, mc2, delta_s)
  use bmad_struct
  
  implicit none
  
  type (coord_struct) :: orbit
  real(rp) :: mc2, delta_s, delta_t, v_s, e_tot, vel(3)
  
  
  e_tot = sqrt( orbit%vec(2)**2 + orbit%vec(4)**2 +  orbit%vec(6)**2 + mc2**2) !Get e_tot from momentum

  vel(1:3) = c_light*[  orbit%vec(2), orbit%vec(4), orbit%vec(6) ]/ e_tot ! velocities v_x, v_y, v_s:  c*[c*p_x, c*p_y, c*p_s]/e_tot

  delta_t = delta_s / vel(3)
  
  !Drift x, y, s
  orbit%vec(1) = orbit%vec(1) + vel(1)*delta_t  !x
  orbit%vec(3) = orbit%vec(3) + vel(2)*delta_t  !y
  orbit%vec(5) = orbit%vec(5) + vel(3)*delta_t  !s
  orbit%s =  orbit%s + delta_s
  orbit%t =  orbit%t + delta_t 

end subroutine drift_orbit_time

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Subroutine write_time_particle_distribution  (time_file_unit, bunch, mc2, err)
!
! Subroutine to write an time-based bunch from a standard Bmad bunch
! 
! Note: The time-based file format is
!       n_particles
!       x/m  m*c^2 \beta_x*\gamma/eV  y/m m*c^2\beta_y*\gamma/eV s/m m*c^2\beta_z*\gamma/eV time/s charge/C
!       . . .
!       all at the same time. 
!       This is very similar to subroutine write_opal_particle_distribution
!
! Input:
!   time_file_unit -- Integer: unit number to write to, if > 0
!   bunch          -- bunch_struct: bunch to be written.
!                            Particles are drifted to bmad_bunch%t_center for output
!   mc2            -- real(rp): particle mass in eV
!
! Output:          
!   err            -- Logical, optional: Set True if, say a file could not be opened.
!-



subroutine write_time_particle_distribution (time_file_unit, bunch, mc2, err)

implicit none

integer			    :: time_file_unit
type (bunch_struct) :: bunch
real(rp)            :: mc2
logical, optional   :: err

type (coord_struct) :: orb
real(rp)        :: dt, pc, gmc
character(40)	:: r_name = 'write_time_particle_distribution'
character(10)   ::  rfmt 
integer n_particle, i


!
if (present(err)) err = .true.

n_particle = size(bunch%particle)

!Format for numbers
  rfmt = 'es13.5'

!Write number of particles to first line
write(time_file_unit, '(i8)') n_particle

!\gamma m c

!Write out all particles to file
do i = 1, n_particle
  orb = bunch%particle(i)
  
  !Get time to track backwards by
  dt = orb%t - bunch%t_center
  
  !Get pc before conversion
  pc = (1+orb%vec(6))*orb%p0c 
  
  !convert to time coordinates
  call convert_particle_coordinates_s_to_t (orb)
  
  !get \gamma m c
  gmc = sqrt(pc**2 + mc2**2) / c_light
  
  !'track' particles backwards in time and write to file
  write(time_file_unit, '(8'//rfmt//')')  orb%vec(1) - dt*orb%vec(2)/gmc, &   !x - dt mc2 \beta_x \gamma / \gamma m c
                                          orb%vec(2), &
                                          orb%vec(3) - dt*orb%vec(4)/gmc, &   !y - dt mc2 \beta_y \gamma / \gamma m c
                                          orb%vec(4), &
                                          orb%vec(5) - dt*orb%vec(6)/gmc, &   !s - dt mc2 \beta_s \gamma / \gamma m c
                                          orb%vec(6), &
										  bunch%t_center,  &!time
										  bunch%particle(i)%charge 
end do 

end subroutine  write_time_particle_distribution


end module
