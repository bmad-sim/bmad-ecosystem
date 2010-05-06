!........................................................................
!+
!  module    : beambeam_interface 
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
! Revision 1.3  2007/01/30 16:14:31  dcs
! merged with branch_bmad_1.
!
! Revision 1.2.2.1  2006/12/22 20:30:42  dcs
! conversion compiles.
!
! Revision 1.2  2005/09/21 20:59:07  dcs
! more changes to get around compiler bug.
!
! Revision 1.1.1.1  2005/06/14 14:59:02  cesrulib
! Beam Simulation Code
!
!
!........................................................................
!
#include "CESR_platform.inc"

module beambeam_interface


!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------

contains

subroutine size_beam(ring, end, scan_params, transmit, sib_j, n_typeout, orb, phi_x, phi_y, past_params, past_lums, parameters)

  use bmad_struct
  use bmad_interface
  use bmadz_mod
  use bmadz_interface
  use bsim_interface
  use scan_parameters
  
  implicit none

  type(lat_struct) ring
  type(coord_struct) end(:)
  type(scan_params_struct) scan_params
  type(coord_struct) orb(0:)
  
  logical, dimension(2) :: transmit
  
  integer sib_j, n_typeout

  real(rp) :: phi_x, phi_y
  real(rp) :: past_params(:,:), past_lums(:), parameters(:,:)

  type(coord_struct) co_
  integer i,l,n_ok,j
  real(rp) :: param_sum, param_mean(1:3),A(1:3), variances(1:3)
  real(rp) :: t1, t2, t3, t4, t5, t6, t_wait, t_send
  real(rp) :: mean_lum,dev_lum,xsum,x2sum
  real(rp) :: Qx_in , Qy_in
  real(rp), parameter :: j_print_dist=10000
  integer :: t_count, k, p_unit, e_unit, ti_unit, t_unit, g_unit
  integer :: damping
  integer :: mem_mean
  logical :: fit_sig, fit_var
  character*80 params_file, times_file, turns_file, ends_file, dist_file, final_pos_file, in_file
  character*20 wordx, wordy
  type(coord_struct) :: final_pos_out
  type (coord_struct), allocatable, save ::  end_coord(:)

  real(rp) x_offset, y_offset

  call reallocate_coord(end_coord, scan_params%n_part)
  param_mean(:) = 0.0

  damping = scan_params%damping

  fit_sig = .false.
  fit_var = .false.
  if(trim(scan_params%fit)=='sig') fit_sig = .true.
  if(trim(scan_params%fit)=='var') fit_var = .true.

     open(unit=17,file="dist.track")
     do i=1,scan_params%n_part
        if(end(i)%vec(1) == 999. ) cycle  
        write(17,'(6e12.4)')end(i)%vec(:)
     end do
     close(17)

     ! move distribution of particles back to the origin for fitting
     n_ok = 0
     do i = 1, scan_params%n_part
        if(end(i)%vec(1) == 999.)cycle
        n_ok = n_ok +1
        x_offset = orb(0)%vec(2) * end(i)%vec(5)
        y_offset = orb(0)%vec(4) * end(i)%vec(5)
        co_%vec(:) = end(i)%vec(:) - orb(0)%vec(:)
        co_%vec(1) = co_%vec(1) - x_offset
        co_%vec(3) = co_%vec(3) - y_offset
        end_coord(n_ok) = co_
     enddo
     print *,' Finished tracking particle ',n_ok,'  of', scan_params%n_part

     scan_params%n_part_out = n_ok

     ! write distribution after tilting
!     if(mod(sib_j*n_typeout,j_print_dist) == 0.0) then
!        write(wordx,'(f0.2)') scan_params%current
!        write(wordy,'(i6.6)') sib_j*n_typeout
!        dist_file = './dists/dist.calc_' // trim(wordx) // '_' // trim(wordy)
!        d_unit = lunget()
!        open(unit=d_unit,file=dist_file)
!        do i=1,scan_params%n_part_out
!           write(d_unit,'(6e12.4)')end_coord(i)%vec(:)
!        end do
!        close(d_unit)
!     end if
!     
     open(unit=18,file="dist.calc")
     do i=1,n_ok
        write(18,'(6e12.4)')end_coord(i)%vec(:)
     end do
     close(18)

     if(sib_j*n_typeout == scan_params%n_turn)then
        call file_suffixer (scan_params%file_name, in_file, '.end', .true.)
        call histogram(ring%ele(0),end_coord(1:n_ok), in_file, scan_params%sig_out,A)
     else
        call histogram(ring%ele(0),end_coord(1:n_ok), 'junk', scan_params%sig_out,A)
     endif
     variances(:)=scan_params%sig_out(:)

!>>>>>>>>>>>>>>>>>>>>>>>>>>
     ! ANDREW'S CODE--evaluate weak_beam_size & update strong_beam_size
     !run coords though gfit3D, overwrite scan_params%sig_out if fit_sig=.true.
     g_unit=lunget()

     open(unit=g_unit,file="gfit3d.chisq",action="write",access="append")
     write(g_unit,'(i8,2x)',advance='no') sib_j*n_typeout
     close(g_unit)

     call gfit3D(end_coord(1:n_ok),parameters)

!     write(g_unit,'(1x)')
     close(g_unit)

     if(fit_sig) scan_params%sig_out(:) = parameters(3,:)
        
     past_params(sib_j,1:3) = scan_params%sig_out(1:3)
     
     ! determine how many sigmas to use in calculating mean
     ! example: if we have 200000 turns with n_typeout=500, we'll let the mean depend on 
     !   the 10 most recent values of sigma.  But, if n_typeout=100, we need the mean line 
     !   to depend on more values of sigma (50) so that it won't jump around unnecessarily.
     !   We will let the mean depend on the past 1/40th of the total number of turns. 
     mem_mean = 10
!     mem_mean = INT((scan_params%n_turn/40)/n_typeout)
!     if (mem_mean.lt.5) mem_mean = 5
     
     ! use the mean of the sigmas for weak beam to compare beam sizes and adjust accordingly
     param_mean(1:3) = 0.0
     do i=1,3
        param_sum  = 0.0
        if(sib_j - mem_mean >  0) then
           do j=INT(sib_j-(mem_mean-1)),sib_j
              param_sum = param_sum + past_params(j,i)
           end do
           param_mean(i) = param_sum/mem_mean
        else
           do j=1,sib_j
              param_sum = param_sum + past_params(j,i)
           end do
           param_mean(i) = param_sum/sib_j
        end if
     end do
     
     ! write .params file
     p_unit=lunget()
     call file_suffixer (scan_params%file_name, params_file, '.params', .true.)
     open(unit=p_unit, file=params_file,access="append")
     write(p_unit,'(i6,e12.4,e12.4,e12.4,3e12.4,3e12.4,3e12.4)') sib_j*n_typeout,ring%ele(1)%value(sig_x$), &
          ring%ele(1)%value(sig_y$), ring%ele(1)%value(sig_z$), variances(1:3), parameters(3,1:3), param_mean(1:3)
     close(p_unit)
     
     ! we adjust the size of the strong beam to match the weak beam 
     ! when either fit_sig or fit_var == .true.
     if(fit_sig .or. fit_var) then
        ! only start adjusting strong beam after a few damping periods
        if(sib_j*n_typeout > damping) then
           
           if( 1.2*ring%ele(1)%value(sig_x$) .lt. param_mean(1) ) then
              ring%ele(1)%value(sig_x$) = ring%ele(1)%value(sig_x$) + 0.5*ABS( param_mean(1)-ring%ele(1)%value(sig_x$) )
              print *,"Strong Beam x Adjustment"
              transmit(2)=.true.
           end if
           if( 0.8*ring%ele(1)%value(sig_x$) .gt. param_mean(1) ) then
              ring%ele(1)%value(sig_x$) = ring%ele(1)%value(sig_x$) - 0.5*ABS( param_mean(1)-ring%ele(1)%value(sig_x$) )
              print *,"Strong Beam x Adjustment"
              transmit(2)=.true.
           end if
           
           if( 1.2*ring%ele(1)%value(sig_y$) .lt. param_mean(2) ) then
              ring%ele(1)%value(sig_y$) = ring%ele(1)%value(sig_y$) + 0.5*ABS( param_mean(2)-ring%ele(1)%value(sig_y$) )
              print *,"Strong Beam y Adjustment"
              transmit(2)=.true.
           end if
           if( 0.8*ring%ele(1)%value(sig_y$) .gt. param_mean(2) ) then
              ring%ele(1)%value(sig_y$) = ring%ele(1)%value(sig_y$) - 0.5*ABS( param_mean(2)-ring%ele(1)%value(sig_y$) )
              print *,"Strong Beam y Adjustment"
              transmit(2)=.true.
           end if
           
           if( 1.2*ring%ele(1)%value(sig_z$) .lt. param_mean(3) ) then
              ring%ele(1)%value(sig_z$) = ring%ele(1)%value(sig_z$) + 0.5*ABS( param_mean(3)-ring%ele(1)%value(sig_z$) )
              print *,"Strong Beam z Adjustment"
              transmit(2)=.true.
           end if
           if( 0.8*ring%ele(1)%value(sig_z$) .gt. param_mean(3) ) then
              ring%ele(1)%value(sig_z$) = ring%ele(1)%value(sig_z$) - 0.5*ABS( param_mean(3)-ring%ele(1)%value(sig_z$) )
              print *,"Strong Beam z Adjustment"
              transmit(2)=.true.
           end if
        end if
        
        ! overwrite .end file written by histogram with new parameters
        if(( .not.transmit(1)) .and. (fit_sig) )then
           call file_suffixer (scan_params%file_name, in_file, '.end', .true.)
           call writefile(in_file, parameters)
        else
           call writefile('junk', parameters)
        endif
     end if  ! matches if(fit) then
!<<<<<<<<<<<<<<<<<<<<<<<<<<         


!if there are no more good particles or we have reached the designated number of turns, end the simulation
     if(n_ok.eq.0) then 
        transmit(1) = .false.
        print *,"All particles lost, halting simulation"
     endif
     
     if(sib_j*n_typeout.eq.scan_params%n_turn) transmit(1)=.false.
     
     if(ring%ele(1)%key == beambeam$)call luminosity_calc (ring%ele(1), &
          end_coord, ring%param, n_ok, &
          scan_params%lum)
     t_unit=lunget()
     call file_suffixer (scan_params%file_name, turns_file, '.turns', .true.)
     open(unit=t_unit, file= turns_file, status='unknown',access='append')
     write(t_unit,'(2f10.5,2x,2f10.5,3e12.4,2x,3e12.4,2x,i6,2x,e12.4,2x,i5,2x,e12.4,2x,i5)')phi_x/twopi-int(phi_x/twopi), phi_y/twopi-int(phi_y/twopi), &
          ring%a%tune/twopi, ring%b%tune/twopi, &
          ring%ele(1)%value(sig_x$), &
          ring%ele(1)%value(sig_y$), ring%ele(1)%value(sig_z$), &
          scan_params%sig_out(1:3), n_typeout*sib_j, scan_params%lum, scan_params%n_part_out, scan_params%current,scan_params%n_cars
     close(unit=t_unit)
     
     print *, n_typeout*sib_j, scan_params%lum,scan_params%n_part_out, scan_params%current


!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ANDREW CODE--write .ends file
     ! save past luminosities
     past_lums(sib_j)=scan_params%lum
     
     !if this is the last turn: find max & min of mem_mean last lumninocities, write .ends file
     ! recall mem_mean is the integer number of sets of turns that we use to calculate the mean value
     ! of the fitted sigmas.
     if(sib_j*n_typeout.eq.scan_params%n_turn) then
!        max_lum = past_lums(sib_j)
!        min_lum = past_lums(sib_j)
!        do k=sib_j-(mem_mean-1),sib_j
!           if(past_lums(k).gt.max_lum) max_lum = past_lums(k)
!           if(past_lums(k).lt.min_lum) min_lum = past_lums(k)
!        end do

! calculate mean and standard dev of last mem_mean luminosities
        xsum=0
        x2sum=0
        mem_mean = 10
        mem_mean = INT((scan_params%n_turn/10)/n_typeout)
        if (mem_mean.lt.5) mem_mean = 5
        do i=1,mem_mean
           xsum=xsum+past_lums(sib_j-mem_mean+i)
           x2sum=x2sum+past_lums(sib_j-mem_mean+i)**2
        end do
        mean_lum = xsum/mem_mean
        dev_lum = SQRT((x2sum/mem_mean) - (mean_lum)**2)
                

        e_unit=lunget()
        call file_suffixer (scan_params%file_name, ends_file, '.ends', .true.)
        open(e_unit, file=ends_file,access='append')
        write(e_unit,'(I7,2x,4f10.5,2e12.4,2x,2e12.4,2x,i5,2x,e10.3,3e12.4,2x,i3,3e12.4)') scan_params%n_turn, &
             phi_x/twopi-10, phi_y/twopi-9, ring%a%tune/twopi, ring%b%tune/twopi, &
             ring%ele(1)%value(sig_x$), ring%ele(1)%value(sig_y$), &
             scan_params%sig_out(1:2), scan_params%n_part_out, scan_params%current, mean_lum, mean_lum+dev_lum, mean_lum-dev_lum, &
             scan_params%n_cars, ring%ele(0)%b%alpha
        close(e_unit)
     end if
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   end subroutine size_beam

end module beambeam_interface
