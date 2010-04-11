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
! Revision 1.3  2007/01/30 16:14:31  dcs
! merged with branch_bmad_1.
!
! Revision 1.2.2.1  2006/12/22 20:30:42  dcs
! conversion compiles.
!
! Revision 1.2  2005/10/17 14:59:02  dlr
! missing optional arguments to close_pretzel and close_vert give error. Include
! arguments
!
! Revision 1.1.1.1  2005/06/14 14:59:02  cesrulib
! Beam Simulation Code
!
!
!........................................................................
!
#include "CESR_platform.inc"

  subroutine beambeam_initialize(ring, scan_params, phi_x, phi_y)
    

   use bmad_struct
   use bmad_interface
   use bmadz_mod
   use bmadz_interface
   use bsim_interface
   use scan_parameters

  implicit none

  type (lat_struct) ring
  type (lat_struct), save :: ring_in, ring_out, ring_save
  type (coord_struct), allocatable, save ::  start_coord(:), end_coord(:) 
  type (coord_struct), allocatable, save :: co(:), orbit(:)
  type (coord_struct), allocatable, save :: start(:), end(:)
  type (coord_struct) orb, delta_ip
  type (coord_struct) final_pos_in, final_pos_out
  type (normal_modes_struct) mode
  type (ele_struct) beambeam_ele
  type (scan_params_struct) scan_params

  real(rp) current
  real(rp) Qx, Qy
  real(rp) rgamma, den, xi_v, xi_h
  real(rp) r_0/2.81794092e-15/
  real(rp) charge
  real(rp) phi_x, phi_y
  real(rp), allocatable :: dk1(:) 
  real(rp) a_out(1:3)  
  
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

  call reallocate_coord(co,ring%n_ele_max)
  call reallocate_coord(orbit, ring%n_ele_max)
  call reallocate_coord(start_coord, scan_params%n_part)
  call reallocate_coord(end_coord, scan_params%n_part)
  call reallocate_coord(start, scan_params%n_part)
  call reallocate_coord(end, scan_params%n_part)

  call set_on_off (rfcavity$, ring, off$)

  ring%param%particle = scan_params%particle

  call twiss_at_start(ring)
  co(0)%vec = 0.
  call closed_orbit_calc(ring, co, 4)
  call lat_make_mat6(ring,-1,co)
  call twiss_at_start(ring)
  call twiss_propagate_all(ring)  

  type *
  type *,' BEAMBEAM_INITIALIZE: Initially '
  type *,'    Qx = ',ring%a%tune/twopi,'    Qy = ',ring%b%tune/twopi
  type '(a15,4e12.4)','  Closed orbit ', co(0)%vec(1:4)



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
   co(0)%vec = 0.
   call closed_orbit_calc(ring, co, 4)
   call lat_make_mat6(ring,-1, co)
   call twiss_at_start(ring)

   type *
   type '(a51,f4.1,a8)', &
        ' BEAMBEAM_INITIALIZE: After parasitic interactions added ', current,'mA/bunch' 
   type *,'    Qx = ',ring%a%tune/twopi,'    Qy = ',ring%b%tune/twopi
   type '(a15,4e12.4)','  Closed orbit ', co(0)%vec(1:4)

  endif
  
  Qx = ring%a%tune/twopi
  Qy = ring%b%tune/twopi
  ring%z%tune = scan_params%Q_z * twopi
  if(scan_params%beambeam_ip)then
    call beambeam_setup(ring, particle, current, scan_params, slices)

    call twiss_at_start(ring)
    co(0)%vec = 0.
    call closed_orbit_calc(ring, co, 4)
   type '(a15,4e12.4)','  Closed orbit ', co(0)%vec(1:4)
    call lat_make_mat6(ring, -1, co)
    call twiss_at_start(ring)

    type *
    type *,' BEAMBEAM_INITIALIZE: After beambeam added '

    type '(a37,f4.1,a8)', &
         ' BEAMBEAM_INITIALIZE: After beambeam added ', current,'mA/bunch' 
    type *,'    Qx = ',ring%a%tune/twopi,'    Qy = ',ring%b%tune/twopi
    type '(a15,4e12.4)','  Closed orbit ', co(0)%vec(1:4)
    write(23, *)
    write(23, '(a37,f4.1,a8)') &
         ' BEAMBEAM_INITIALIZE: After beambeam added ', current,'mA/bunch' 

    beambeam_ele = ring%ele(1)
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

  call set_on_off (rfcavity$, ring, on$)
  call set_z_tune(ring)


  i_dim = 6


  if(scan_params%close_pretz)then
    final_pos_in%vec(:) = 0.
    final_pos_out%vec(:) = 0.
    call close_pretzel (ring, i_dim, final_pos_in, final_pos_out) 
    call twiss_at_start(ring)
    call closed_orbit_calc(ring, co, i_dim)
    call lat_make_mat6(ring, -1, co)
    call twiss_at_start(ring)
    type*,' after close pretzel but before close vertical: '
    type '(1x,3(a9,f12.4))','    Qx = ',ring%a%tune/twopi,'    Qy = ',ring%b%tune/twopi,'   Qz = ',ring%z%tune/twopi
 endif
 if(scan_params%close_vert)call close_vertical(ring,i_dim, final_pos_in, final_pos_out)

    call twiss_at_start(ring)

    forall( i=0:ring%n_ele_track) co(i)%vec = 0.
    call closed_orbit_calc(ring, co, i_dim)
    call lat_make_mat6(ring, -1, co)
    call twiss_at_start(ring)
    call twiss_propagate_all(ring)

    type *
    type *,' After CLOSE VERTICAL ' 
    type '(1x,3(a9,f12.4))','    Qx = ',ring%a%tune/twopi,'    Qy = ',ring%b%tune/twopi,'   Qz = ', ring%z%tune/twopi
    type '(a15,4e12.4)','  Closed orbit ', co(0)%vec(1:4)

   write(23,'(a36,f6.4,a12,f6.4)')' Beam beam tune shifts:  Delta Qx = ', ring%a%tune/twopi - Qx, &
                                           '  Delta Qy =',ring%b%tune/twopi-Qy

   if(scan_params%beambeam_ip)then
    rgamma = ring%ele(0)%value(E_TOT$)/mass_of(-1)
    den = twopi*rgamma*beambeam_ele%value(sig_x$)+beambeam_ele%value(sig_y$)
    xi_v = r_0*ring%param%n_part*ring%ele(0)%b%beta/beambeam_ele%value(sig_y$)/den
    xi_h = r_0*ring%param%n_part*ring%ele(0)%a%beta/beambeam_ele%value(sig_x$)/den
   write(23,'(2(a10,f7.4))')'   xi_v = ', xi_v,'   xi_h = ',xi_h
   endif

! find beambeam at IP and turn it off
  do i=0,ring%n_ele_max
   if(ring%ele(i)%name == 'IP_COLLISION')then
     ix_ip = i
     charge = ring%ele(i)%value(charge$)
     ring%ele(i)%value(charge$) = 0
     exit
   endif
  end do

  type *,' BEAMBEAM: turn off beam beam at IP'

! Qtune
      allocate(dk1(ring%n_ele_max))
       ring_save = ring
       call choose_quads(ring, dk1)
       call custom_set_tune (phi_x, phi_y, dk1, ring, co, ok)
      deallocate(dk1)
      if(.not. ok)then
        ring = ring_save
        print '(a,2f10.3,a)',' Tunes ',phi_x/twopi,phi_y/twopi, ' are unstable. Not tracking ' 
        return
      endif

      if(scan_params%close_pretz)then
        call close_pretzel (ring, i_dim, final_pos_in, final_pos_out)
!        call twiss_at_start(ring)
        call closed_orbit_calc(ring, co, i_dim)
        call lat_make_mat6(ring, -1, co)
        call twiss_at_start(ring)
        type *,' beam beam at IP is off'
        type*,' after qtune and after close pretzel but before close vertical: '
    type '(1x,3(a9,f12.4))','    Qx = ',ring%a%tune/twopi,'    Qy = ',ring%b%tune/twopi,'   Qz = ',ring%z%tune/twopi
      endif
      if(scan_params%close_vert)then
        call close_vertical(ring,i_dim, final_pos_in, final_pos_out)
!        call twiss_at_start(ring)
        call closed_orbit_calc(ring, co, i_dim)
        call lat_make_mat6(ring, -1, co)
        call twiss_at_start(ring)
        type*,' after qtune with pretzel and vert closed but beam beam at IP off: '
    type '(1x,3(a9,f12.4))','    Qx = ',ring%a%tune/twopi,'    Qy = ',ring%b%tune/twopi,'   Qz = ', ring%z%tune/twopi
      endif
! Turn beambeam back on
  if(ix_ip /= 0)ring%ele(ix_ip)%value(charge$) = charge
  call twiss_at_start(ring)
  call closed_orbit_calc(ring, co, i_dim)
  call lat_make_mat6(ring, -1, co)
  call twiss_at_start(ring)
  call beambeam_separation(ring, delta_ip, i_dim)
  type *,' Turn Beambeam on'
  type *,'    Qx = ',ring%a%tune/twopi,'    Qy = ',ring%b%tune/twopi
  type '(a22,4f8.4)', ' dx,dxp,dy,dyp (mm) = ', delta_ip%vec(1:4)*1000.

    beambeam_ele = ring%ele(1)
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


  call set_on_off (rfcavity$, ring, on$)
  call set_z_tune(ring)
  if(scan_params%radiation)then
    bmad_com%radiation_damping_on = .true.
    bmad_com%radiation_fluctuations_on = .true.
    type *,' radiation fluctuations and damping are on'
  endif
  forall(i=0:ring%n_ele_track) orbit(i)%vec = co(i)%vec
  call radiation_integrals (ring, orbit, mode)

 do i = 1, scan_params%n_part
   start(i)%vec = orbit(0)%vec + scan_params%init(i)%vec
   open(unit=21, file=scan_params%file_name, status='old', access='append')
   call track_a_particle (orbit(0), start(i), scan_params%n_turn, ring, mod(phi_x/twopi,1.), &
                                mod(phi_y/twopi,1.), ring%z%tune/twopi)
   close(unit=21)
 end do
 return

 end
   





















