!+
! Subroutine da_driver (RING, TRACK_INPUT, n_xy_pts, &
!                         point_range, energy, n_energy_pts, in_file,Qx,Qy,Qz, particle, Qp_x, Qp_y)
!
! Subroutine to determine starting point for dynamic aperture tracking.
! The subroutine works by determining where on a radial line y = const * x
! the aperture is.  Here x and y are deviations from the closed orbit.
!
! Modules Needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!   RING    -- lat_struct: Ring containing the lattice.
!   TRACK_INPUT   -- Track_input_struct: Structure holding the input data:
!     %N_TURN         -- Number of turns tracked.
!     %X_INIT         -- Suggested initial x coordinate to start with.
!     %Y_INIT         -- Suggested initial y coordinate to start with.
!     %ACCURACY       -- Accuracy needed of aperture results.
!
!    N_XY_PTS       -- Number of radial lines in a quarter plane.
!    ENERGY         -- Array of energy deviations of initial coordinates.
!    N_ENERGY_PTS   -- Number of off energy runs to do.
!    Qp_x, Qp_y     -- Real: Chromaticity
!
!    IN_FILE        -- name of input file 
!
! Note: The radial lines are spaced equally in angle using coordinates
!       normalized by %X_INIT and %Y_INIT
!-
!........................................................................
! Mod/Commons:
!
! Calls      :
!
! Author     :
!
! Modified   :
!
!........................................................................
!
! $Id$
!
! $Log$
! Revision 1.2  2007/01/30 16:14:31  dcs
! merged with branch_bmad_1.
!
! Revision 1.1.1.1.2.1  2006/12/22 20:30:42  dcs
! conversion compiles.
!
! Revision 1.1.1.1  2005/06/14 14:59:02  cesrulib
! Beam Simulation Code
!
!
!........................................................................
!
#include "CESR_platform.h"


subroutine da_driver (ring, track_input, n_xy_pts, point_range, &
                                           energy, n_energy_pts, in_file, Qx,Qy,Qz, particle, Qp_x, Qp_y)

  use dynamic_aperture_mod                                      
  use bmadz_interface

  implicit none

  type (lat_struct)  ring
  type (lat_param_struct)  ring_param_state
  type (coord_struct) orb0
  type (coord_struct), allocatable, save :: co(:)
  type (aperture_struct)  aperture
  type (track_input_struct)  track_input

  integer i_e, i_xy, it, i, turn_lost, ixr, i_e_max, ix
  integer n_xy_pts, point_range(2), n_energy_pts
  integer int_Q_y, int_Q_x
  integer particle
  integer len 

  real(rp) eps_rel(4), eps_abs(4)
  real(rp) e_init, theta                                   
  real(rp) x0, x1, x2, y0, y1, y2
  real(rp) energy(10)
  real(rp), allocatable :: dk1(:)
  real(rp) Qx, Qy, Qz
  real(rp) Qp_x, Qp_y
  real(rp) phy_x_set, phy_y_set
  real(rp) delta_e/1.e-4/, chrom_x, chrom_y

  logical aperture_bracketed, track_on
  logical ok

  character*60 da_file, in_file
  character date_str*20
  character * 9 particle_type(-1:1)/'electrons',' ','positrons'/
  character*200 file_name

!

  if (track_input%x_init == 0 .or. track_input%y_init == 0) then
    type *, 'ERROR IN DYNAMIC_APERTURE: TRACK_INPUT.X_INIT OR',  &
                                             ' TRACK_INPUT.Y_INIT = 0'
    call err_exit
  endif

  ring_param_state = ring%param

  do i=1,ring%n_ele_track
   if(ring%ele(i)%key /= wiggler$)cycle
   if(ring%ele(i)%key == rfcavity$)cycle
   if(ring%ele(i)%key == quadrupole$)cycle
!   if(ring%ele(i)%tracking_method /= custom$)then
!     ring%ele(i)%tracking_method = taylor$
!     ring%ele(i)%mat6_calc_method = taylor$
!   endif
  end do

  call reallocate_coord (co, ring%n_ele_max)
  allocate(dk1(ring%n_ele_max))

  call twiss_at_start (ring)
  call twiss_propagate_all(ring)

  if(Qz /= 0.)then
   ring%z%tune = Qz * twopi
   call set_z_tune(ring)
  endif
  type *,' Q_z = ',ring%z%tune/twopi

  type *,' Before Qtune: Qx = ',ring%a%tune/twopi,'    Qy = ',ring%b%tune/twopi
  if(Qx /= 0. .and. Qy /= 0.)then 
    int_Q_x = int(ring%ele(ring%n_ele_track)%a%phi / twopi)
    int_Q_y = int(ring%ele(ring%n_ele_track)%b%phi / twopi)
    phy_x_set = (int_Q_x + Qx)*twopi
    phy_y_set = (int_Q_y + Qy)*twopi
    call choose_quads(ring, dk1)
    call custom_set_tune (phy_x_set, phy_y_set, dk1, ring, co, ok) 
    if(.not. ok) type *,' Qtune failed'
  endif


  call twiss_at_start (ring)
  type '(a19,f12.4,a9,f12.4)',' After Qtune: Qx = ',ring%a%tune/twopi,'    Qy = ',ring%b%tune/twopi



  if(Qp_x /= 0. .and. Qp_y /= 0.)then
    call qp_tune(ring, qp_x, qp_y, ok)
    if( .not. ok)then
     type *,' Qp_tune failed '
    endif
  endif
  call chrom_calc(ring, delta_e, chrom_x, chrom_y) 
  type '(a23,f12.4,a11,f12.4)',' After Qp_tune: Qp_x = ',chrom_x,'    Qp_y = ',chrom_y
  call twiss_at_start(ring)
  type '(a21,f12.4,a9,f12.4)',' After Qp_tune: Qx = ',ring%a%tune/twopi,'    Qy = ',ring%b%tune/twopi

    call custom_set_tune (phy_x_set, phy_y_set, dk1, ring, co, ok) 
    if(.not. ok) type *,' Second Qtune failed'

  call twiss_at_start(ring)
  type '(a26,f12.4,a9,f12.4)',' After second qtune: Qx = ',ring%a%tune/twopi,'    Qy = ',ring%b%tune/twopi


  call closed_orbit_calc (ring, co, 4)

!  eps_rel(:) = 0.000001
!  eps_abs(:) = 0.000001
!  call closed_orbit_from_tracking (ring, co, 4, eps_rel, eps_abs, co(0))

  call lat_make_mat6(ring, -1, co)
  call twiss_at_start (ring) 
   ring%param%aperture_limit_on = .true.

! write output
  call file_suffixer (in_file, in_file, '.dat', .true.)
  open (unit = 2, file = in_file,   &
                                        carriagecontrol = 'list')
  write (2, *) 'Lattice  = ', ring.lattice
  write (2, '(a11,i5,3(a10,f7.5))') ' N_turn   =', track_input.n_turn
  write(2,'(a8,a1,f8.4,a1)') '  Q_x = ',"`",ring%a%tune/twopi,"'"
  write(2,'(a8,a1,f8.4,a1)') '  Q_y = ',"`",ring%b%tune/twopi,"'"
  write(2,'(a8,a1,f8.4,a1)') '  Q_z = ',"`",ring%z%tune/twopi,"'"
  write (2, *) 'n_xy_pts =', n_xy_pts
  write (2, *) 'point_range =',point_range
  write (2, *) 'n_energy_pts =', n_energy_pts
  write (2, *) 'x_init   =', track_input.x_init
  write (2, *) 'y_init   =', track_input.y_init
  write (2, *) 'energy    =', &
                       (energy(i),i=1, n_energy_pts)
  write (2, *) 'accuracy =', track_input%accuracy
  write (2, '(a, 6f10.5)') ' Closed_orbit:',  &
                                (co(0)%vec(i), i = 1, 6)
  write(2,'(1x,a16,a1,a9,a1)')' particle_type =',"'", particle_type(particle),"'"  
  write(2,*) 'return'
  write(2,*)

  i_e_max = max(1, n_energy_pts)
  do i_e = 1, i_e_max


     e_init = energy(i_e)

     aperture%closed_orbit = co(0)
     orb0 = co(0)
     orb0%vec(6) = e_init

    call string_trim (in_file, in_file, ix)
    write (da_file, '(a, i1.1, a)') in_file(1:ix-4), i_e, '.dat'
    open (unit = 3, file = da_file,   &
                                        carriagecontrol = 'list')

    file_name = ring%input_file_name
    ix=0
    do while (index(file_name(ix+1:),'/') /= 0)
     ix = index(file_name(ix+1:),'/') +ix
    end do
    len = len_trim(file_name(ix+1:))
    file_name = "`"//file_name(ix+1:ix+len)//"'"
    write (3, '(a11,a)') 'Lattice  = ', file_name
    write (3, *) 'dE_E     = ', e_init
    call date_and_time_stamp (date_str)
    write (3, '(4a)') 'data_date = `', date_str, "'"
    write(3,*)'Q_x = ',"`",'Q_x = ',ring%a%tune/twopi,"'"
    write(3,*)'Q_y = ',"`",'Q_y = ',ring%b%tune/twopi,"'"
    write(3,*)'Q_z = ',"`",'Q_z = ',ring%z%tune/twopi,"'"
    write(2,'(1x,a16,a1,a9,a1)')' particle_type =',"'", particle_type(particle),"'"  
    write (3, *) 'return'
    write (3, *)
    write (3, *)              
    write (3, *) '!        X          Y    Turn  Plane   Lost_At'


    do i_xy = point_range(1), point_range(2)

     theta = (i_xy - 1) * pi / max(1, n_xy_pts - 1)

      call dynamic_aperture (ring, orb0, theta, track_input, aperture)


      write (2, '(2f11.6, i7, 6x, a1, 3x, a)') aperture%x, &
              aperture%y,  &
              aperture%i_turn, &
              plane_name(aperture%plane),  &
              ring%ele(aperture%ix_lat)%name


      write (3, '(2f11.6, i7, 6x, a1, 3x, a)') aperture%x, &
               aperture%y,  &
              aperture%i_turn, &
              plane_name(aperture%plane),  &
              ring%ele(aperture%ix_lat)%name


    enddo ! i_xy

    close (unit = 3)

  enddo ! i_e

  close (unit = 2)

!

  ring%param = ring_param_state
  deallocate(dk1)

  return
  end







































