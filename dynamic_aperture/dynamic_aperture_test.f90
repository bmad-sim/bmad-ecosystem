!........................................................................
!+
! program    :  dynamic_aperture_test
!
! Description:  test of dynamic aperture subroutine
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
! Revision 1.2  2007/01/30 16:14:31  dcs
! merged with branch_bmad_1.
!
! Revision 1.1.1.1.2.1  2006/12/22 20:30:43  dcs
! conversion compiles.
!
! Revision 1.1.1.1  2005/06/14 14:59:02  cesrulib
! Beam Simulation Code
!
!
!........................................................................
!

program dynamic_aperture_test

  use dynamic_aperture_mod                                      
  use bmadz_interface
  use cesr_crossings_mod
  use bookkeeper_mod
  use bsim_interface


  implicit none

interface
 subroutine da_driver (ring, track_input, n_xy_pts, point_range, &
                               energy, n_energy_pts, in_file, Qx,Qy,Qz, particle, Qp_x, Qp_y, &
                               delta_fRF, fRF)

!  use bmad_struct
!  use bmad_interface
  use bmad
  use bmadz_interface
!  use bsim_interface
  use dynamic_aperture_mod                                      
  implicit none
  type (lat_struct)  ring
  type (track_input_struct)  track_input

  integer n_xy_pts, point_range(2), n_energy_pts
  integer particle

  real(rdef) energy(10)

  real(rdef) Qx, Qy, Qz
  real(rdef) Qp_x, Qp_y

  real(rdef) delta_fRF, fRF


  character*60 in_file

 end subroutine da_driver
end interface



  type (lat_struct) ring
  type (lat_struct), save :: ring_in, ring_out
  type (track_input_struct) track_input
  type (ele_struct) ele
  type (coord_struct), allocatable, save :: co(:)

  integer n_turn, n_xy_pts, n_energy_pts, i, j, ix, ios, point_range(2), particle
  integer i_train, j_car, n_trains_tot, n_cars
  integer i_bunch
  integer ierr
  integer crossings_tot, bunch_tot, total
  integer, dimension(:,:), allocatable :: ix_lrbbi, ix_master
  integer k
  integer n
  integer int_Q_x, int_Q_y
  integer i_dim/4/  

  real(rp) x_init, y_init, e_max, accuracy, energy(10)
  real(rp) ap_mult, aperture_multiplier, Qx, Qy, Qz, Qx_ini, Qy_ini, Qp_x, Qp_y
  real(rp), dimension(:,:), allocatable :: crossings_1
  real(rp), dimension(:), allocatable :: cross_positions
  real(rp) current
  real(rp) phy_x_set, phy_y_set
  real(rp), allocatable :: dk1(:)
  real(rp) delta_fRF/0./, fRF

  character*60 da_file, in_file
  character*100 lat_file
  character date_str*20

  logical ok
  logical rec_taylor
  logical path_length_patch/.false./

  namelist / input / lat_file,n_turn, n_xy_pts, point_range, n_energy_pts,  &
                     x_init, y_init, e_max, energy, accuracy, &
                     aperture_multiplier, Qx, Qy, Qz, particle, &
                     i_train, j_car, n_trains_tot, n_cars, current, lat_file, &
                     Qx_ini, Qy_ini, Qp_x, Qp_y, rec_taylor, delta_fRF, fRF

! init

  do 

    print '(a, $)', ' Input command file <CR = DA_TEST.IN>: '
    read(*, '(a)') in_file
    call string_trim (in_file, in_file, ix)
    if (ix == 0)in_file = 'da_test.in'
    call file_suffixer (in_file, in_file, '.in', .false.)

    open (unit = 1, file = in_file, status = 'old', iostat = ios)

    if (ios == 0) then
      exit
    else
      print *
      print *, 'ERROR: CANNOT OPEN FILE: ', trim(in_file)
    endif

  enddo

  Qx=0
  Qy=0
  Qx_ini=0
  Qy_ini=0
  aperture_multiplier = 10.
  particle = positron$
  rec_taylor = .true.
  fRF = 5.e8
  read (1, nml = input)
  track_input%n_turn = n_turn
  track_input%x_init = x_init
  track_input%y_init = y_init
  track_input%accuracy = accuracy
  ap_mult = aperture_multiplier
  close (unit = 1)

  call bmad_parser (lat_file, ring)
  if(delta_fRF /= 0.)then
     path_length_patch = .false.
     call implement_pathlength_patch(path_length_patch,ring, delta_fRF, fRF)
     i_dim = 6
  endif
  call reallocate_coord (co, ring%n_ele_track)
  allocate(dk1(ring%n_ele_max))

  call twiss_at_start(ring)
  call twiss_propagate_all(ring)

  type '(a28,f12.6,a9,f12.6)', &
       ' Before initial Qtune: Qx = ',ring%a%tune/twopi,'    Qy = ',ring%b%tune/twopi
  if (Qx_ini /= 0. .and. Qy_ini /= 0.)then 
    int_Q_x = int(ring%ele(ring%n_ele_track)%a%phi / twopi)
    int_Q_y = int(ring%ele(ring%n_ele_track)%b%phi / twopi)
    phy_x_set = (int_Q_x + Qx_ini)*twopi
    phy_y_set = (int_Q_y + Qy_ini)*twopi
    call choose_quads(ring, dk1)
    do i=0,ring%n_ele_track
     co(i)%vec = 0
    end do
   call set_on_off(elseparator$, ring, off$)
   call custom_set_tune (phy_x_set, phy_y_set, dk1, ring, co, ok) 
    if (.not. ok) print *,' Qtune failed'
   call set_on_off(elseparator$, ring, on$)
  endif

  type '(a28,f12.6,a9,f12.6)', &
       '  After initial Qtune: Qx = ',ring%a%tune/twopi,'    Qy = ',ring%b%tune/twopi



  call twiss_at_start(ring)
  call closed_orbit_calc(ring, co, i_dim)
  call lat_make_mat6(ring,-1,co)
  call twiss_at_start(ring)
  call twiss_propagate_all (ring)


  ring%param%particle = particle

! set up for lrbbi

  if (i_train*j_car * n_cars * j_car * current /= 0)then
   ring_in = ring
   call lrbbi_setup (ring_in, ring_out, particle, i_train, j_car, n_trains_tot, n_cars, current, rec_taylor)
   ring = ring_out

   call twiss_at_start(ring)
   co(0)%vec = 0
   call closed_orbit_calc(ring, co, i_dim)
   call lat_make_mat6(ring,-1, co)
   call twiss_at_start(ring)

   print *
   type '(a51,f4.1,a8)', &
        ' BEAMBEAM_SCAN: After parasitic interactions added ', current,'mA/bunch' 
   print *,'    Qx = ',ring%a%tune/twopi,'    Qy = ',ring%b%tune/twopi
   type '(a15,4e12.4)','  Closed orbit ', co(0)%vec(1:4)

  endif

  print *, ' DA_TEST: 1 before '
  call closed_orbit_calc(ring, co, i_dim)
  print *, ' DA_TEST: 1 after '


  call calc_z_tune (ring)
  call twiss_at_start(ring)

  do i = 1, ring%n_ele_track
   if(ring%ele(i)%key == wiggler$)cycle
    ring%ele(i)%value(x1_limit$) = ap_mult * ring%ele(i)%value(x1_limit$)
    ring%ele(i)%value(x2_limit$) = ap_mult * ring%ele(i)%value(x2_limit$)
    ring%ele(i)%value(y1_limit$) = ap_mult * ring%ele(i)%value(y1_limit$)
    ring%ele(i)%value(y2_limit$) = ap_mult * ring%ele(i)%value(y2_limit$)
    if (ring%ele(i)%value(x1_limit$) == 0) ring%ele(i)%value(x1_limit$) = 0.3
    if (ring%ele(i)%value(x2_limit$) == 0) ring%ele(i)%value(x2_limit$) = 0.3
    if (ring%ele(i)%value(y1_limit$) == 0) ring%ele(i)%value(y1_limit$) = 0.3
    if (ring%ele(i)%value(y2_limit$) == 0) ring%ele(i)%value(y2_limit$) = 0.3
  enddo

  ring%param%aperture_limit_on = .true.

! track                                                

  call da_driver (ring, track_input,n_xy_pts, &
                        point_range, energy, n_energy_pts, in_file,Qx,Qy,Qz,particle, Qp_x, Qp_y, &
                        delta_fRF, fRF)

  deallocate(dk1)

end program
                                                            










