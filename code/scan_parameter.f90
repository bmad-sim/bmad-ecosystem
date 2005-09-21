!........................................................................
!+
! module scan_parameters
!
! Description:
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
! Revision 1.3  2005/09/21 20:35:58  dcs
! Put beambeam_setup in scan_parameters module.
!
! Revision 1.2  2005/07/11 14:56:56  sni2
! added damping control to infile
!
! Revision 1.1.1.1  2005/06/14 14:59:02  cesrulib
! Beam Simulation Code
!
!
!........................................................................
!
#include "CESR_platform.h"

module scan_parameters

  use bmad

  type scan_params_struct
    character*80 lat_file
    character*80 file_name
    real(rdef) Q_z
    real(rdef) current  !mA/bunch
    real(rdef) sig_in(3)   !sigx, sigy, sigz, initial distribution
    real(rdef) sig_out(3)  !sigx, sigy, sigz, final distribution
    real(rdef) min_sig     !x/sigx + y/sigy + z/sigz > min_sig
    real(rdef) lum         !luminosity
    real(rdef) coupling_wb    !coupling of horizontal emittance into vertical weak beam
    real(rdef) coupling_sb    !coupling of horizontal emittance into vertical strong beam
    integer n_turn, particle, i_train, j_car, n_trains_tot, n_cars 
    integer n_part !initial distribution
    integer slices 
    integer n_part_out !final distribution
    logical lrbbi, beambeam_ip, close_pretz, close_vert, rec_taylor
    logical radiation
    character*5 fit
    type(coord_struct) final_pos_in ! gives final postion to which closed orbit will converge
    type(coord_struct) init(100)  !initial coordinate for single particle scan
    logical parallel !defines whether code is run in parallel
    integer damping  ! number of turns after which beam resizing begins
 end type scan_params_struct

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine: BEAMBEAM_SETUP (ring, particle,  current, scan_params,  slices)
!
! Description: Subroutine to add beambeam element at IP 
!
! Arguments  :
!   RING -- Ring_struct: Ring containing the lattice to be modified
!   particle  -- Integer : (+-1) for positrons or electrons
!   current   -- real : bunch current in ma
!   slices    -- integer, optional : number of slices, default is 1 
!   SCAN_PARAMS -- scan_params_struct, for emittance coupling or beam sizes
! Output:
!    RING -- Ring_struct: Ring includes Beam-beam element at IP
! 
! Mod/Commons:
!
! Calls      :
!
! Author     :
!
! Modified   :
!-

subroutine beambeam_setup(ring, particle,  current,scan_params, slices)

  use bmadz_struct
  use bmadz_interface
  use bsim_interface

  implicit none

  type ( ring_struct ) ring
  type (ring_struct), save :: ring_oppos
  type (coord_struct), allocatable :: co_(:), co_oppos_(:)
  type (ele_struct) ele, beambeam_ele
  type (modes_struct) mode
  type (scan_params_struct) scan_params

  integer particle
  integer ierr
  integer i,n,ix_ip
  integer, optional, intent(in) :: slices


  real(rdef) current
  real(rdef) major, minor, theta, vert_size_coupled, vert_size_emit
  real(rdef) coupling

  character(20) :: r_name='BEAMBEAM_SETUP'

  allocate( co_(0:ring%n_ele_max))
  allocate( co_oppos_(0:ring%n_ele_max))
! init


     
  ring%param%particle = particle
  call twiss_at_start(ring)
!  type *,' beambeam_setup:1 beta ',ring%ele_(0)%x%beta, ring%ele_(0)%y%beta
  co_(0)%vec = 0.
  call closed_orbit_at_start(ring, co_(0), 4, .true.)
!  type *, ' closed_orbit: co_',co_(0)%vec(1:4)
  call track_all (ring, co_)
  call ring_make_mat6(ring,-1,co_)
  call twiss_at_start(ring)
!  type *,' beambeam_setup:2 beta ',ring%ele_(0)%x%beta, ring%ele_(0)%y%beta
  call twiss_propagate_all (ring)

! set up for lrbbi
  if( current == 0)then
     call out_io(s_abort$,r_name,' bunch specifications incomplete ')
!    type *,' BEAMBEAM_SETUP: bunch specifications incomplete '
    stop
  endif


  ring%param%n_part   = current*0.001 *(ring%param%total_length/c_light)/e_charge
  ring_oppos = ring
  ring_oppos%param%n_part = 0.
  ring_oppos%param%particle = -particle

  call ring_make_mat6(ring_oppos,-1)
  co_oppos_(0)%vec = 0.
  call closed_orbit_at_start(ring_oppos, co_oppos_(0), 4, .true.)
  call track_all (ring_oppos, co_oppos_)
  call ring_make_mat6(ring_oppos,-1,co_oppos_)
  call twiss_at_start(ring_oppos)
!  type *,' beambeam_setup:3 beta ',ring%ele_(0)%x%beta, ring%ele_(0)%y%beta
  call twiss_propagate_all (ring_oppos)
  call set_on (rfcavity$, ring_oppos, .true.)
  call set_z_tune(ring_oppos)
  call radiation_integrals (ring_oppos, co_oppos_, mode)

  call ellipse (ring_oppos%ele_(0), major, minor, theta)

  call init_ele (beambeam_ele)

  beambeam_ele%name = 'IP_COLLISION'
  beambeam_ele%key = beambeam$
  coupling = scan_params%coupling_sb

  if(scan_params%sig_in(1) == 0)then
!    beambeam_ele%value(sig_x$) = sqrt(ring_oppos%ele_(0)%x%beta * mode%a%emittance &
!                                    + (ring_oppos%ele_(0)%x%eta * mode%sige_e)**2)

     beambeam_ele%value(sig_x$) = sqrt( mode%a%emittance * major**2 + &
                                      (ring_oppos%ele_(0)%x%eta * mode%sige_e)**2)

!    beambeam_ele%value(sig_y$) = sqrt(ring_oppos%ele_(0)%y%beta * mode%a%emittance*0.01 &
!                                    + (ring_oppos%ele_(0)%y%eta * mode%sige_e)**2)
     vert_size_coupled = mode%a%emittance * max( minor**2, coupling * ring_oppos%ele_(0)%y%beta)
     vert_size_emit    = mode%b%emittance * ring_oppos%ele_(0)%y%beta
     beambeam_ele%value(sig_y$) = sqrt(vert_size_coupled + vert_size_emit  &
                                    + (ring_oppos%ele_(0)%y%eta * mode%sige_e)**2)
     beambeam_ele%value(sig_z$) = mode%sig_z
  else
     beambeam_ele%value(sig_x$) = scan_params%sig_in(1)
     beambeam_ele%value(sig_y$) = scan_params%sig_in(2)
     beambeam_ele%value(sig_z$) = scan_params%sig_in(3)
  endif


  if(present(slices))then
     beambeam_ele%value(n_slice$) = slices
    else
     beambeam_ele%value(n_slice$) = 1
  endif
  beambeam_ele%value(charge$) = -1 
  beambeam_ele%value(x_pitch$) =  co_oppos_(0)%vec(2)
  beambeam_ele%value(y_pitch$) =  co_oppos_(0)%vec(4)
  beambeam_ele%value(x_offset$) = co_oppos_(0)%vec(1)
  beambeam_ele%value(y_offset$) = co_oppos_(0)%vec(3)
  beambeam_ele%value(tilt$) = theta !tilt of beam ellipse due to coupling

  ix_ip = 1

  call add_superimpose (ring, beambeam_ele, ix_ip)
  call ring_make_mat6(ring, -1)

  call out_io(s_info$,r_name,' Strong beam: ')
  call out_io(s_blank$,r_name,'  sigma_x = \e12.4\ ',beambeam_ele%value(sig_x$))
  call out_io(s_blank$,r_name,'  sigma_y = \e12.4\ ',beambeam_ele%value(sig_y$))
  call out_io(s_blank$,r_name,'  sigma_z = \e12.4\ ',beambeam_ele%value(sig_z$))
  call out_io(s_blank$,r_name,'  Pitch:')
  call out_io(s_blank$,r_name,'    x = \e12.4\ ',beambeam_ele%value(x_pitch$))
  call out_io(s_blank$,r_name,'    y = \e12.4\ ',beambeam_ele%value(y_pitch$))
  call out_io(s_blank$,r_name,'  Offset:')
  call out_io(s_blank$,r_name,'    x = \e12.4\ ',beambeam_ele%value(x_offset$))
  call out_io(s_blank$,r_name,'    y = \e12.4\ ',beambeam_ele%value(y_offset$))
  call out_io(s_blank$,r_name,'  Tilt = \e12.4\ ', beambeam_ele%value(tilt$))

  
  
!  type '(1x,a14)', ' Strong beam: '
!  type '(1x,a12,e12.4)', '  sigma_x = ',beambeam_ele%value(sig_x$)
!  type '(1x,a12,e12.4)', '  sigma_y = ',beambeam_ele%value(sig_y$)
!  type '(1x,a12,e12.4)', '  sigma_z = ',beambeam_ele%value(sig_z$)
!  type '(1x,a14,e12.4,a4,e12.4)', '  Pitch  : x= ',beambeam_ele%value(x_pitch$), &
!                                ' y= ',beambeam_ele%value(y_pitch$)
!  type '(1x,a14,e12.4,a4,e12.4)', '  Offset : x= ',beambeam_ele%value(x_offset$), &
!                                ' y= ',beambeam_ele%value(y_offset$)
!  type '(1x,a9,e12.4)', '  Tilt = ', beambeam_ele%value(tilt$)

!   type '(1x,a1,5x,2x,a12,7x,a1,4x,4a9)','n','Element name','s',' x_pitch ',' y_pitch ', &
!             'x_offset ', &
!             'y_offset'

!  do i=1, ring%n_ele_use
!   ele = ring%ele_(i)
!   if(ele%key == beambeam$)then
!    n=n+1
!   type '(i3,4x,a16,4x,f7.3,2x,4f9.5)',n,ring%ele_(i)%name, ring%ele_(i)%s, &
!          ele%value(x_pitch$), ele%value(y_pitch$), ele%value(x_offset$), ele%value(y_offset$)
!   endif
!   end do


  deallocate( co_)
  deallocate( co_oppos_)

end subroutine                                                            

end module
