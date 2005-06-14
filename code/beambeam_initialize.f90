!........................................................................
!+
! Subroutine beambeam_initialize(ring, scan_params, phi_x, phi_y)
!
! Description:
!
! Arguments  :
!
! Mod/Commons:
!
! Calls      :
!
! Author     :
!
! Modified   :
!-
!........................................................................
!
! $Id$
!
! $Log$
! Revision 1.1  2005/06/14 14:59:02  cesrulib
! Initial revision
!
!
!........................................................................
!
#include "CESR_platform.h"

  subroutine beambeam_initialize(ring, scan_params, phi_x, phi_y)
    

   use bmad_struct
   use bmad_interface
   use bmadz_mod
   use bmadz_interface
   use bsim_interface
   use scan_parameters

  implicit none

  type (ring_struct) ring
  type (ring_struct), save :: ring_in, ring_out
  type (coord_struct), allocatable, save ::  start_coord_(:), end_coord_(:) 
  type (coord_struct), allocatable, save :: co_(:), orb_(:)
  type (coord_struct), allocatable, save :: start(:), end(:)
  type (coord_struct) orb, delta_ip
  type (modes_struct) mode
  type (ele_struct) beambeam_ele
  type (scan_params_struct) scan_params

  real(rdef) current
  real(rdef) Qx, Qy
  real(rdef) rgamma, den, xi_v, xi_h
  real(rdef) r_0/2.81794092e-15/
  real(rdef) charge
  real(rdef) phi_x, phi_y
  real(rdef), allocatable :: dk1(:) 
  real(rdef) a_out(1:3)  
  
  integer i, j
  integer ix
  integer ios
  integer particle, i_train, j_car, n_trains_tot, n_cars
  integer slices
  integer n_ok
  integer i_dim
  integer ix_ip/0/
  integer n_typeout/100/

  character date*10, time*10
  character*80 in_file, turns_file

  logical ok
  logical rec_taylor

  call reallocate_coord(co_,ring%n_ele_max)
  call reallocate_coord(orb_, ring%n_ele_max)
  call reallocate_coord(start_coord_, scan_params%n_part)
  call reallocate_coord(end_coord_, scan_params%n_part)
  call reallocate_coord(start, scan_params%n_part)
  call reallocate_coord(end, scan_params%n_part)

  call setup_radiation_tracking(ring, co_, .false., .false.)
  call set_on (rfcavity$, ring, .false.)

  ring%param%particle = scan_params%particle

  call twiss_at_start(ring)
  co_(0)%vec = 0.
  call closed_orbit_at_start(ring, co_(0), 4, .true.)
  call track_all (ring, co_)
  call ring_make_mat6(ring,-1,co_)
  call twiss_at_start(ring)
  call twiss_propagate_all(ring)  

  type *
  type *,' BEAMBEAM_INITIALIZE: Initially '
  type *,'    Qx = ',ring%x%tune/twopi,'    Qy = ',ring%y%tune/twopi
  type '(a15,4e12.4)','  Closed orbit ', co_(0)%vec(1:4)



  particle = scan_params%particle
  i_train = scan_params%i_train
  j_car = scan_params%j_car
  n_trains_tot = scan_params%n_trains_tot
  n_cars = scan_params%n_cars
  current = scan_params%current
  rec_taylor = scan_params%rec_taylor
  slices = scan_params%slices


  if(scan_params%lrbbi)then
   ring_in = ring
   call lrbbi_setup ( ring_in, ring_out, particle, i_train, j_car, n_trains_tot, n_cars, current, rec_taylor)
   ring = ring_out

   call twiss_at_start(ring)
   co_(0)%vec = 0.
   call closed_orbit_at_start(ring, co_(0), 4, .true.)
   call track_all(ring, co_)
   call ring_make_mat6(ring,-1, co_)
   call twiss_at_start(ring)

   type *
   type '(a51,f4.1,a8)', &
        ' BEAMBEAM_INITIALIZE: After parasitic interactions added ', current,'mA/bunch' 
   type *,'    Qx = ',ring%x%tune/twopi,'    Qy = ',ring%y%tune/twopi
   type '(a15,4e12.4)','  Closed orbit ', co_(0)%vec(1:4)

  endif
  
  Qx = ring%x%tune/twopi
  Qy = ring%y%tune/twopi
  ring%z%tune = scan_params%Q_z * twopi
  if(scan_params%beambeam_ip)then
    call beambeam_setup(ring, particle, current, scan_params, slices)

    call twiss_at_start(ring)
    co_(0)%vec = 0.
    call closed_orbit_at_start(ring, co_(0), 4, .true.)
   type '(a15,4e12.4)','  Closed orbit ', co_(0)%vec(1:4)
    call track_all(ring, co_)
    call ring_make_mat6(ring, -1, co_)
    call twiss_at_start(ring)

    type *
    type *,' BEAMBEAM_INITIALIZE: After beambeam added '

    type '(a37,f4.1,a8)', &
         ' BEAMBEAM_INITIALIZE: After beambeam added ', current,'mA/bunch' 
    type *,'    Qx = ',ring%x%tune/twopi,'    Qy = ',ring%y%tune/twopi
    type '(a15,4e12.4)','  Closed orbit ', co_(0)%vec(1:4)
    write(23, *)
    write(23, '(a37,f4.1,a8)') &
         ' BEAMBEAM_INITIALIZE: After beambeam added ', current,'mA/bunch' 

    beambeam_ele = ring%ele_(1)
    write(23,  '(1x,a14)') ' Strong beam: '
    write(23,  '(1x,a12,e12.4)') '  sigma_x = ',beambeam_ele%value(sig_x$)
    write(23,  '(1x,a12,e12.4)') '  sigma_y = ',beambeam_ele%value(sig_y$)
    write(23,  '(1x,a12,e12.4)') '  sigma_z = ',beambeam_ele%value(sig_z$)
    write(23,  '(1x,a14,e12.4,a4,e12.4)') '  Pitch  : x= ',beambeam_ele%value(x_pitch$), &
                                ' y= ',beambeam_ele%value(y_pitch$)
    write(23, '(1x,a14,e12.4,a4,e12.4)') '  Offset : x= ',beambeam_ele%value(x_offset$), &
                                ' y= ',beambeam_ele%value(y_offset$)
    write(23, '(1x,a9,e12.4)') '  Tilt = ', beambeam_ele%value(tilt$)



 endif

  call set_on (rfcavity$, ring, .true.)
  call set_z_tune(ring)


  i_dim = 6


  if(scan_params%close_pretz)then
    call close_pretzel (ring, i_dim) 
    call twiss_at_start(ring)
    call closed_orbit_at_start(ring, co_(0), i_dim, .true.)
    call track_all(ring, co_)
    call ring_make_mat6(ring, -1, co_)
    call twiss_at_start(ring)
    type*,' after close pretzel but before close vertical: '
    type '(1x,3(a9,f12.4))','    Qx = ',ring%x%tune/twopi,'    Qy = ',ring%y%tune/twopi,'   Qz = ',ring%z%tune/twopi
 endif
 if(scan_params%close_vert)call close_vertical(ring,i_dim)

    call twiss_at_start(ring)

    forall( i=0:ring%n_ele_use) co_(i)%vec = 0.
    call closed_orbit_at_start(ring, co_(0), i_dim, .true.)
    call track_all(ring, co_)
    call ring_make_mat6(ring, -1, co_)
    call twiss_at_start(ring)
    call twiss_propagate_all(ring)

    type *
    type *,' After CLOSE VERTICAL ' 
    type '(1x,3(a9,f12.4))','    Qx = ',ring%x%tune/twopi,'    Qy = ',ring%y%tune/twopi,'   Qz = ', ring%z%tune/twopi
    type '(a15,4e12.4)','  Closed orbit ', co_(0)%vec(1:4)

   write(23,'(a36,f6.4,a12,f6.4)')' Beam beam tune shifts:  Delta Qx = ', ring%x%tune/twopi - Qx, &
                                           '  Delta Qy =',ring%y%tune/twopi-Qy

   if(scan_params%beambeam_ip)then
    rgamma = ring%ele_(0)%value(beam_energy$)/mass_of(-1)
    den = twopi*rgamma*beambeam_ele%value(sig_x$)+beambeam_ele%value(sig_y$)
    xi_v = r_0*ring%param%n_part*ring%ele_(0)%y%beta/beambeam_ele%value(sig_y$)/den
    xi_h = r_0*ring%param%n_part*ring%ele_(0)%x%beta/beambeam_ele%value(sig_x$)/den
   write(23,'(2(a10,f7.4))')'   xi_v = ', xi_v,'   xi_h = ',xi_h
   endif

! find beambeam at IP and turn it off
  do i=0,ring%n_ele_max
   if(ring%ele_(i)%name == 'IP_COLLISION')then
     ix_ip = i
     charge = ring%ele_(i)%value(charge$)
     ring%ele_(i)%value(charge$) = 0
     exit
   endif
  end do

  type *,' BEAMBEAM: turn off beam beam at IP'

! Qtune
      allocate(dk1(ring%n_ele_max))
       call choose_quads(ring, dk1)
       call custom_set_tune (phi_x, phi_y, dk1, ring, co_, ok)
      deallocate(dk1)

      if(scan_params%close_pretz)then
        call close_pretzel (ring, i_dim)
!        call twiss_at_start(ring)
        call closed_orbit_at_start(ring, co_(0), i_dim, .true.)
        call track_all(ring, co_)
        call ring_make_mat6(ring, -1, co_)
        call twiss_at_start(ring)
        type *,' beam beam at IP is off'
        type*,' after qtune and after close pretzel but before close vertical: '
    type '(1x,3(a9,f12.4))','    Qx = ',ring%x%tune/twopi,'    Qy = ',ring%y%tune/twopi,'   Qz = ',ring%z%tune/twopi
      endif
      if(scan_params%close_vert)then
        call close_vertical(ring,i_dim)
!        call twiss_at_start(ring)
        call closed_orbit_at_start(ring, co_(0), i_dim, .true.)
        call track_all(ring, co_)
        call ring_make_mat6(ring, -1, co_)
        call twiss_at_start(ring)
        type*,' after qtune with pretzel and vert closed but beam beam at IP off: '
    type '(1x,3(a9,f12.4))','    Qx = ',ring%x%tune/twopi,'    Qy = ',ring%y%tune/twopi,'   Qz = ', ring%z%tune/twopi
      endif
! Turn beambeam back on
  if(ix_ip /= 0)ring%ele_(ix_ip)%value(charge$) = charge
  call twiss_at_start(ring)
  call closed_orbit_at_start(ring, co_(0), i_dim, .true.)
  call track_all(ring, co_)
  call ring_make_mat6(ring, -1, co_)
  call twiss_at_start(ring)
  call beambeam_separation(ring, delta_ip, i_dim)
  type *,' Turn Beambeam on'
  type *,'    Qx = ',ring%x%tune/twopi,'    Qy = ',ring%y%tune/twopi
  type '(a22,4f8.4)', ' dx,dxp,dy,dyp (mm) = ', delta_ip%vec(1:4)*1000.

    beambeam_ele = ring%ele_(1)
    write(23,*)
    write(23,  '(1x,a14)') ' Strong beam: after closing pretzel '
    write(23,  '(1x,a12,e12.4)') '  sigma_x = ',beambeam_ele%value(sig_x$)
    write(23,  '(1x,a12,e12.4)') '  sigma_y = ',beambeam_ele%value(sig_y$)
    write(23,  '(1x,a12,e12.4)') '  sigma_z = ',beambeam_ele%value(sig_z$)
    write(23,  '(1x,a14,e12.4,a4,e12.4)') '  Pitch  : x= ',beambeam_ele%value(x_pitch$), &
                                ' y= ',beambeam_ele%value(y_pitch$)
    write(23, '(1x,a14,e12.4,a4,e12.4)') '  Offset : x= ',beambeam_ele%value(x_offset$), &
                                ' y= ',beambeam_ele%value(y_offset$)
    write(23, '(1x,a9,e12.4)') '  Tilt = ', beambeam_ele%value(tilt$)


  call set_on (rfcavity$, ring, .true.)
  call set_z_tune(ring)
  if(scan_params%radiation)then
    call setup_radiation_tracking(ring, co_, .true., .true.)
    type *,' radiation fluctuations and damping are on'
  endif
  forall(i=0:ring%n_ele_use) orb_(i)%vec = co_(i)%vec
  call radiation_integrals (ring, orb_, mode)

 do i = 1, scan_params%n_part
   start(i)%vec = orb_(0)%vec + scan_params%init(i)%vec
   open(unit=21, file=scan_params%file_name, status='old', access='append')
   call track_a_particle (orb_(0), start(i), scan_params%n_turn, ring, mod(phi_x/twopi,1.), &
                                mod(phi_y/twopi,1.), ring%z%tune/twopi)
   close(unit=21)
 end do
 return

 end
   





















