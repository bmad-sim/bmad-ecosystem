!........................................................................
!+
! program    :  closed_orbit
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
! Revision 1.2.2.2  2007/01/12 16:59:16  dcs
! beta_x -> beta_a, etc.
!
! Revision 1.2.2.1  2006/12/22 20:30:42  dcs
! conversion compiles.
!
! Revision 1.2  2006/11/17 21:13:14  cesrulib
! changed ele_name size to reflect change in bmad
!
! Revision 1.1.1.1  2005/06/14 14:59:02  cesrulib
! Beam Simulation Code
!
!
!........................................................................
!
#include "CESR_platform.inc"

  program closed_orbit

  use bmadz_interface
  use sim_utils
  use cbar_mod
  use bookkeeper_mod

  implicit none

 type (lat_struct) ring_1, ring_2
 type (lat_struct), save :: ring, ring_two(-1:1)
 type (coord_struct), allocatable, save :: co(:), cot(:), co_high(:), co_low(:)
 type (coord_struct), allocatable, save :: co_off(:)
 type (coord_struct), allocatable, save :: co_electron(:)
 type (coord_struct) traj
 type (coord_struct) dorb
 type (normal_modes_struct) mode
 type (coord_struct) track_start

      integer pgopen, istat1, istat2
      integer i, j, k(6)/1,3,6,2,4,5/
      integer ix
      integer, allocatable :: track_meth(:)
      integer ke
      integer n/0/
      integer nd
      integer plot_flag, last
      integer, parameter :: orbit$=1,beta$=2,cbar$=3,diff$=4, de_beta$=5
      integer, parameter :: eta$=6, de_cbar$=7
      integer ix_cache
      integer, allocatable :: n_ele(:)
      integer i_dim/4/
      integer ix_ele
      integer n_turns
      integer nargs, iargc
      integer n_all
      integer l

      real*4, allocatable :: z(:), x(:), y(:), zz(:,:), xx(:,:), yy(:,:)
      real*4, allocatable :: zz_diff(:), xx_diff(:), yy_diff(:)
      real*4 width/7./, aspect/1./
      real*4 xmax/0./, ymax/0./, xmax0, ymax0
      real*4 xscale, yscale, x_low, y_low
      real*4 xdet(1000), ydet(1000), zdet(1000)
      real(rp) cbar_mat(2,2), cbar_mat1(2,2), cbar_mat2(2,2)
      real(rp) de/8e-4/
      real(rp) rms_x, rms_y
      real(rp) frev
      real*4 length, start, end     
      real(rp)rate_x, rate_y, rate_xq, rate_yq, rate_x_tot, rate_y_tot
      real(rp) d_amp_x, d_amp_y
      real (rp) p1, p2
      real (rp) f,p
      real(rp) axx, axy, ayy
      real(rp) n_part_save
      real(rp) slopes(4)
     
      character*40 lattice
      character*60 lat_file
      character*120 line, last_line, vec_start
      character*16 x_or_y, answer, save_answer
      character*72 comment
      character*20 device_type, last_device_type/' '/
      character*40 ele_names(4)

      logical keep_trying/.true./
      logical write_orbit/.false./                                        
      logical diff
      logical radiation/.true./
      logical transfer_line/.false./
      logical track/.false./
      logical cbarve/.true./


!
      nargs = cesr_iargc()
      if(nargs == 1)then
         call cesr_getarg(1, lat_file)
         print *, 'Using ', trim(lat_file)
       else

      lat_file = 'bmad.'
      print '(a,$)',' Lattice file name ? (default= bmad.) '
      read(5,'(a)') line
       call string_trim(line, line, ix)
       lat_file = line
       if(ix == 0) lat_file = 'bmad.'
       print *, ' lat_file = ', lat_file
     endif

      call bmad_parser (lat_file, ring_1)
      ring  = ring_1
  call reallocate_coord (co, ring%n_ele_max)
  call reallocate_coord (co_electron, ring%n_ele_max)
  call reallocate_coord (cot, ring%n_ele_max)
  call reallocate_coord (co_off, ring%n_ele_max)
  call reallocate_coord (co_high, ring%n_ele_max)
  call reallocate_coord (co_low, ring%n_ele_max)

  co(0)%vec = 0
  co_electron(0)%vec=0
  cot(0)%vec=0
  co_off(0)%vec=0
  co_high(0)%vec=0
  co_low(0)%vec=0
  track_start%vec = 0.

  allocate(track_meth(0:ring%n_ele_max))
  allocate(x(0:ring%n_ele_max))
  allocate(y(0:ring%n_ele_max))
  allocate(z(0:ring%n_ele_max))
  allocate(xx_diff(0:ring%n_ele_max))
  allocate(yy_diff(0:ring%n_ele_max))
  allocate(zz_diff(0:ring%n_ele_max))
  allocate(yy(0:ring%n_ele_max,1:5))
  allocate(xx(0:ring%n_ele_max,1:5))
  allocate(zz(0:ring%n_ele_max,1:5))
  allocate(n_ele(1:5))

   length = ring%ele(ring%n_ele_track)%s
  
    last = 0

 do while (.not. write_orbit)
  do while (keep_trying)

10  print '(a, $)', ' CLOSED_ORBIT: element change or GO> '
      read(5, '(a)',err=10) line
     
  ix = index(line, '!')
  if (ix /= 0) line = line(:ix-1)        ! strip off comments

  call str_upcase(line, line)
  call string_trim(line, line, ix)

  if (ix == 0) then       ! nothing typed. do the same thing
      line = last_line
  endif

   last_line = line

   call str_upcase(line,line)
   if(line(1:1) .eq. 'G')exit
   if(index(line, 'RADIATION') /= 0)then
    if(index(line, 'ON') /= 0) radiation = .true.
    if(index(line, 'OFF') /= 0) radiation = .false.
    exit 
  endif

   if(index(line, 'CBAR_V_E') /= 0)then
    if(index(line, 'ON') /= 0) cbarve = .true.
    if(index(line, 'OFF') /= 0) cbarve = .false.
    exit 
  endif

   if(line(1:2) == 'EX' .or. line(1:2) == 'QU')then
!     if(istat1 > 0)then
!      call pgslct(istat1)
!      call pgclos
!     endif
     if(istat2 > 0)then
      call pgslct(istat2)
      call pgclos
     endif
     stop
   endif
   if(line(1:4) == 'TRAN')then
     if(index(line,'OFF') == 0)then
       transfer_line = .true.
       print *,' Transfer line mode on'
      else
       transfer_line = .false.
       print *,' Transfer line mode off'
     endif
     call set_on_off (rfcavity$, ring, off$)
     cycle
   endif

   if(line(1:4) == 'READ')then
      print '(a,$)',' Lattice file name ? '
      read(5, '(a)') line
       call string_trim(line, line, ix)
       lat_file = line
       print *, ' lat_file = ', lat_file
     call bmad_parser(lat_file, ring_2)
     ring = ring_2
     exit
   endif
   if(line(1:4) == 'RING')then
     if(index(line(6:),'1') /= 0)ring=ring_1
     if(index(line(6:),'2') /= 0)ring=ring_2
     exit
   endif
!  lat_file = line

   if(line(1:2) == '6D')then
     i_dim=6
     print *,' Compute 6-d closed orbit'
     exit
   endif
   if(line(1:2) == '4D')then
     i_dim=4
     co(0)%vec(6)=0.
     co_electron(0)%vec(6)=0.
     print *,' Compute 4-d closed orbit'
     exit
   endif

   if(line(1:4) == 'TRAC')then
     call string_trim(line, line, ix)
     track_start%vec = 0. 
     track = .true.
     call string_trim(line(ix+1:), line, ix)
     if(ix == 0)then
       print *,' TRACK <n_turns> <x> <y> <de/e> <xp> <yp> <dl> '
       cycle 
    else
       read(line(1:ix),*)n_turns
       vec_start = line(ix+1:)
      call num_words(vec_start, j)
       if(ix /= 0) read(line(ix+1:),*)(track_start%vec(k(i)),i=1,j)
     endif
    exit
   endif

   if(line(1:2) == 'PS' .or.  line(1:3) == 'GIF')exit


   call find_change( line, ring)

  if(line(1:2) == 'HE' .or. index(line, '?') /= 0)call list_commands

 end do

   if(line(1:2) /= 'PS' .and.  line(1:3) /= 'GIF')then

!  call lat_make_mat6 (ring, -1)

!  forall( i=0:ring%n_ele_use) co(i)%vec = 0.
     ring%param%particle = positron$

    if(.not. transfer_line)then
     call twiss_at_start(ring)
     call closed_orbit_calc(ring, co, i_dim)
    endif

      print *, ' '
      print *,ring%input_file_name
      print '(a  ,6f10.5)',' positron closed_orbit (mm, mr)         i=0 ', &
        (co(0)%vec(i)*1000.,i=1,6)
      call element_locator('IP_L0_END', ring, ix_ele)
      if(ix_ele >=0)print '(a  ,6f10.5)',' positron closed_orbit (mm, mr) i=ip_l0_end ', &
        (co(ix_ele)%vec(i)*1000.,i=1,6)

      if(track)then
        traj%vec(1:6) = co(0)%vec(1:6) + track_start%vec(1:6)       
        print *
        print '(1x,a6,1x,i5,1x,a6)',  'TRACK ', n_turns, ' turns'
        print '(1x,a16,1x,6e12.4)','TRACK: start    ', track_start%vec
        print '(1x,a16,1x,6e12.4)','TRACK: start+co ', co(0)%vec

        call psp(ring, co(0), traj, n_turns, istat2)

        print '(1x,a16,1x,6e12.4,/)','TRACK: end      ', traj%vec
        track=.false.
      endif

     ring%param%particle = electron$
     n_part_save = ring%param%n_part
     ring%param%n_part = 0.
    if(.not. transfer_line) &
     call closed_orbit_calc(ring, co_electron, i_dim)
      print '(a  ,6f10.5)',' electron closed_orbit (mm, mr) ', &
        (co_electron(0)%vec(i)*1000.,i=1,6)

     ring%param%n_part = n_part_save
     ring%param%particle = positron$

     call lat_make_mat6(ring,-1,co)
    if(.not. transfer_line) &
     call twiss_at_start(ring)
     call calc_z_tune (ring)

!      print *,' Recompute tunes with new matrices'
      print '(3(a11,f7.4))', ' Q_x =     ',ring%a%tune/twopI,    ' Q_y =     ',ring%b%tune/twopi, &
                             ' Q_z =     ', ring%z%tune/twopi 
      print '(2(a11,f7.4))',' Beta_a* = ', ring%ele(0)%a%beta, ' Beta_b* = ', ring%ele(0)%b%beta
      print '(2(a11,f7.4))',' Alpha_a*= ', ring%ele(0)%a%alpha, ' Alpha_y*= ', ring%ele(0)%b%alpha
      call c_to_cbar(ring%ele(0),cbar_mat)
      print '(4(a10,f7.4))',' Cbar11 = ',cbar_mat(1,1),' Cbar12 = ',cbar_mat(1,2), &
              ' Cbar22 = ',cbar_mat(2,2),' Cbar21 = ', cbar_mat(2,1)

      call twiss_propagate_all(ring)
!  calculate off energy beta

      ring_two(1) = ring
      ring_two(-1) = ring
      do i = -1,1,2
       co_off(0)%vec(6) = de *i
      if(.not. transfer_line) &
       call closed_orbit_calc(ring_two(i), co_off,i_dim)
       call lat_make_mat6(ring_two(i), -1, co_off)
      if(.not. transfer_line) &
       call twiss_at_start(ring_two(i))
       call twiss_propagate_all(ring_two(i))
       if(i == -1)forall(j=0:ring%n_ele_track)co_low(j)%vec = co_off(j)%vec
       if(i ==  1)forall(j=0:ring%n_ele_track)co_high(j)%vec = co_off(j)%vec
      end do
      call de_dbeta(ring_two(1), ring_two(-1), de, rms_x, rms_y)

      print '(1x,a12,f12.4,a12,f12.4)',' dB_x*/dE = ',(ring_two(1)%ele(0)%a%beta - ring_two(-1)%ele(0)%a%beta)/2/de, &
                             ' dB_y*/dE = ',(ring_two(1)%ele(0)%b%beta - ring_two(-1)%ele(0)%b%beta)/2/de 

      call c_to_cbar(ring_two(1)%ele(0),cbar_mat1)
      call c_to_cbar(ring_two(-1)%ele(0),cbar_mat2)
      print '(1x,a13,2f12.4)',' dcbar*/dE = ',(cbar_mat1(1,1)-cbar_mat2(1,1))/de, &
                                             (cbar_mat1(1,2)-cbar_mat2(1,2))/de 
      print '(1x,a13,2f12.4)','             ',(cbar_mat1(2,1)-cbar_mat2(2,1))/de, &
                                             (cbar_mat1(2,2)-cbar_mat2(2,2))/de 


      print *,' Chromaticity'
      print '(2(a11,f8.4))', ' dQ_x/dE = ',(ring_two(1)%a%tune - ring_two(-1)%a%tune)/twopi/2/de, &
                             ' dQ_y/dE = ',(ring_two(1)%b%tune - ring_two(-1)%b%tune)/twopi/2/de 
      print '(2(a23,f10.4))', ' sqrt(<(dB_x/dE)^2>) = ',rms_x, ' sqrt(<(dB_y/dE)^2>) = ',rms_y

! sextupole detuning rate
     if( .not. transfer_line)then
      call sext_detune(ring, axx, axy, ayy)
      print *
      print *,' Sextupole detuning rates '
      print '(a12,e12.4)',' Alpha_ax = ',axx 
      print '(a12,e12.4)',' Alpha_ay = ',axy 
      print '(a12,e12.4)',' Alpha_yy = ',ayy
      print * 
     endif

      if(.not. transfer_line .and. cbarve)then   
       ele_names(1) = 'SK_Q03E'
       ele_names(2) = 'SK_Q03W'
       ele_names(3:4) = ' '
       call cbar_v_e(ring, ele_names, slopes)
       print '(a38,a10,a4,a10)',' Energy derivative of cbar for insert ',ele_names(1),' to ',ele_names(2)
       print '(4(a10,f7.3))',' d_Cb11 = ',slopes(1),' d_Cb12 = ',slopes(2), &
               ' d_Cb22 = ',slopes(3),' d_Cr21 = ', slopes(4)
       print '(4(a10,f7.4))',' Etax   = ',ring%ele(0)%a%eta,' Etax_p = ',ring%ele(0)%a%etap, &
               ' Etay   = ',ring%ele(0)%b%eta,' Etay_p = ', ring%ele(0)%b%etap
      endif

      ix_cache = 0
      if(radiation)then
       call radiation_integrals (ring, co, mode, ix_cache)
       print '(a24,e12.4,a25,e12.4)',' horizontal emittance = ', mode%a%emittance, &
                                    '    vertical emittance = ',mode%b%emittance
       print '(a17,e12.4,a18,e12.4)',' Energy spread = ',mode%sige_e,'   Bunch length = ',mode%sig_z
       frev=c_light/ring%ele(ring%n_ele_track)%s
       print '(a11,e12.4)',' Revolution freq    = ', frev
       if(mode%a%alpha_damp /= 0.)then
         print '(a22,e12.4)',' Horiz damping time = ',1/mode%a%alpha_damp/frev
         print '(a22,e12.4)',' Vert damping time =  ',1/mode%b%alpha_damp/frev
         print '(a22,e12.4)',' Long damping time =  ',1/mode%z%alpha_damp/frev
         print '(a7,e12.4)',' i1  =',mode%synch_int(1)
         print '(a7,e12.4)',' i2  =',mode%synch_int(2)
         print '(a7,e12.4)',' i3  =',mode%synch_int(3)
         print '(a7,e12.4)',' i4  =',mode%a%synch_int(4)
         print '(a7,e12.4)',' i5a =',mode%a%synch_int(5)
       endif
      endif

      if(.not. transfer_line)then
       call sextupole_resonance(ring,rate_x, rate_y, rate_xq, rate_yq, rate_x_tot, rate_y_tot, mode%sige_e)
       print '(a45,2e12.4)', &
         ' Synchro-betatron sext growth rate at 0.5+Qs/2, H/V ' , rate_x, rate_y
       print '(a45,2e12.4)', &
         ' Synchro-betatron quad growth rate at 0.5+Qs/2, H/V ' , rate_xq, rate_yq
       print '(a45,2e12.4)', &
         ' Synchro-betatron tot growth rate at 0.5+Qs/2, H/V ' , rate_x_tot, rate_y_tot

       if(ring%z%tune /= 0.)then
         call sync_beta_path(ring, d_amp_x, d_amp_y)
         print '(a31,4e12.4)',' (6d orbit amp- 4d orbit amp)/6d amp/ total volts  x   y', d_amp_x, d_amp_y 
         call sync_beta_volt(ring, d_amp_x, d_amp_y)
         print '(a24,2e12.4)',' (d(amp)/d(volt))  x   y', d_amp_x, d_amp_y
       endif
      endif


     answer='                '

     device_type ='/XSERVE'
     diff=.false.

     print '(a,$)',' Plot ? ([ORBIT,BETA,CBAR, DBETA/DE, ETA, DCBAR/DE, DIFF])', answer
     read(5, '(a)')answer
     save_answer = answer

    else !    if(line(1:2) /= 'PS'.and. line(1:2) /= 'GIF')then
      if(line(1:2) == 'PS')device_type = 'plot.ps/VCPS'
      if(line(1:2) == 'GI')device_type = 'plot.gif/VGIF'
      answer = save_answer
    endif

     call str_upcase(answer,answer)
     call string_trim (answer, answer, ix)
     if(index(answer(1:ix),'OR') /= 0)then
       plot_flag = orbit$
      elseif(index(answer(1:ix),'BE') /= 0)then
       plot_flag=beta$
      elseif(index(answer(1:ix),'DB') /= 0)then
       plot_flag=de_beta$
      elseif(index(answer(1:ix),'CB') /= 0)then
       plot_flag=cbar$
      elseif(index(answer(1:ix),'ET') /=0)then
       plot_flag=eta$
      elseif(index(answer(1:ix),'DC') /=0)then
       plot_flag=de_cbar$
      elseif(index(answer,'DI') /= 0 .or. diff)then
       diff = .true.
      else
       cycle
     endif
    
    if(line(1:2) /= 'PS'.and. line(1:3) /= 'GIF')then
     x_or_y = ' '
     xmax0=0.
     ymax0=0.
     start=0.
     end=length
     print *,answer
     call string_trim(answer(ix+1:),answer,ix)
     type*,answer(1:ix)
     if(ix /= 0.)read(answer(1:ix),*)x_or_y
     call string_trim(answer(ix+1:),answer,ix)
     if(ix /= 0.)read(answer(1:ix),*)p1
     call string_trim(answer(ix+1:),answer,ix)
     type*,answer(1:ix)
     if(ix /= 0.)read(answer(1:ix),*)p2
     
     xmax=0.
     ymax=0.

     xmax0=0.
     ymax0=0.


     if(x_or_y(1:1) == 'Y')then
       xmax0 = p1
       ymax0 = p2
     endif
     if(x_or_y(1:1) == 'X')then
       start = p1
       end = p2
     endif

     print *,' start =', start,'   end = ',end
     print *,'  xmax0 =',xmax0,  '  ymax0 = ',ymax0

     if(plot_flag == last)then
       n=n+1
       if(n>5)n=1
     else
       n=1
     endif
     last=plot_flag

    endif ! if(line(1:2) /= 'PS'.and. line(1:2) /= 'GIF')then

!      device_type = '?'
      if(device_type /= last_device_type)then
       last_device_type = device_type
       istat1 = pgopen(device_type)
       if(istat1 .lt. 1) stop
       call pgpap (width, aspect)
       call pgsubp(1,2)
       call pgask(.false.)
       call pgscr(0, 1., 1., 1.)
       call pgscr(1,0.,0.,0.)
       call pgscr(2, 1., 0., 0.)
       call pgscr(3,0.,0.,0.,1.0)
       call pgsch(2.)
       print '(a,$)', ' Comment ?'
       read(5,'(a)') comment
      endif
      call pgslct(istat1)
      if(device_type(1:7)  == '/XSERVE')then


     if(diff)then
      n=1
      call pgeras

     else


     nd=0
     do i=0,ring%n_ele_track
      z(i) = ring%ele(i)%s

      if(plot_flag == orbit$)then
       x(i)= co(i)%vec(1)*1000.
       y(i)= co(i)%vec(3)*1000.
       if(index(ring_two(1)%ele(i)%name, 'WIG') /= 0 )then
!        print '(1x,a16,3f12.4)',ring_two(1)%ele(i)%name,z(i),x(i),y(i)
       endif
       if(index(ring_two(1)%ele(i)%name, 'IP_L0') /= 0 )then
!        print '(1x,a16,3f12.4)',ring_two(1)%ele(i)%name,z(i),x(i),y(i)
       endif
      endif

      if(plot_flag == beta$)then
       x(i) = ring%ele(i)%a%beta
       y(i) = ring%ele(i)%b%beta
       if(i == 0)print '(1x,a16,3a12)','     Element    ','     z      ','   Beta_a   ','   Beta_b   '
       if(index(ring_two(1)%ele(i)%name, 'WIG') /= 0 )then
!        print '(1x,a16,3f12.4)',ring_two(1)%ele(i)%name,z(i),max(x(i-1),x(i)),max(y(i-1),y(i))
       endif
       if(index(ring_two(1)%ele(i)%name, 'IP_L0') /= 0 )then
!        print '(1x,a16,3f12.4)',ring_two(1)%ele(i)%name,z(i),x(i),y(i)
       endif
      endif

      if(plot_flag == de_beta$)then
       x(i) = (ring_two(1)%ele(i)%a%beta - ring_two(-1)%ele(i)%a%beta)/2/de/ &
                   ring%ele(i)%a%beta
       y(i) = (ring_two(1)%ele(i)%b%beta - ring_two(-1)%ele(i)%b%beta)/2/de/ &
                    ring%ele(i)%b%beta
       if(index(ring_two(1)%ele(i)%name, 'WIG') /= 0 )then
!        print *,ring_two(1)%ele(i)%name,z(i),x(i),y(i)
       endif
       if(index(ring_two(1)%ele(i)%name, 'IP_L0') /= 0 )then
!        print *,ring_two(1)%ele(i)%name,z(i),x(i),y(i)
       endif
      endif

      if(plot_flag == eta$)then
       x(i) = (co_high(i)%vec(1) - co_low(i)%vec(1))/2/de
       y(i) = (co_high(i)%vec(3) - co_low(i)%vec(3))/2/de
      endif

      if(plot_flag == cbar$)then
       call c_to_cbar(ring%ele(i),cbar_mat)
       x(i) = cbar_mat(1,2)
       y(i) = cbar_mat(2,2)
!       if(index(ring_two(1)%ele(i)%name, 'WIG') /= 0 )then
!        print '(1x,a16,3f12.4)',ring_two(1)%ele(ring%ni)%name,z(i),x(i),y(i)
!       endif
!       if(ring%ele(i)%s < 15. .or. ring%ele(ring%n_ele_track)%s - ring%ele(i)%s <15. )then
        print '(1x,a16,3f12.4)',ring_two(1)%ele(i)%name,z(i),x(i),y(i)
!       endif
      endif

      if(plot_flag == de_cbar$)then
       call c_to_cbar(ring_two(1)%ele(i),cbar_mat1)
       call c_to_cbar(ring_two(-1)%ele(i),cbar_mat2)
       x(i) = (cbar_mat1(1,2)-cbar_mat2(1,2))/2/de
       y(i) = (cbar_mat1(2,2)-cbar_mat2(2,2))/2/de
      endif

      if(index(ring%ele(i)%name, 'DET') /= 0)then
        nd = nd+1
        zdet(nd) = z(i)
        xdet(nd) = x(i)
        ydet(nd) = y(i)
      endif

      xmax = max(abs(x(i)),xmax)
      ymax = max(abs(y(i)),ymax)
      if(xmax0 /= 0.)xmax=xmax0
      if(ymax0 /= 0.)ymax=ymax0
     end do
     n_all = ring%n_ele_track
     l = n_all
     endif


     if(diff)then
     nd = 0
     l=0
     do i=0,ring%n_ele_track
      do j = 0, n_all
       if(abs(z(j) -  ring%ele(i)%s) > 0.0001) cycle

       l = l+1

       zz_diff(l) = ring%ele(i)%s 

      if(plot_flag == orbit$)then
       xx_diff(l)= co(i)%vec(1)*1000. -x(j)
       yy_diff(l)= co(i)%vec(3)*1000. -y(j)
      endif
      if(plot_flag == beta$)then
       xx_diff(l) = ring%ele(i)%a%beta -x(j)
       yy_diff(l) = ring%ele(i)%b%beta -y(j)
      endif

      if(plot_flag == de_beta$)then
       xx_diff(l) = (ring_two(1)%ele(i)%a%beta - ring_two(-1)%ele(i)%a%beta)/2/de/ &
                  ring%ele(i)%a%beta - x(j)
       yy_diff(l) = (ring_two(1)%ele(i)%b%beta - ring_two(-1)%ele(i)%b%beta)/2/de/ &
                  ring%ele(i)%b%beta - y(j)
      endif

      if(plot_flag == eta$)then
       xx_diff(l) = (co_high(i)%vec(1) - co_low(i)%vec(1))/2/de - x(j)
       yy_diff(l) = (co_high(i)%vec(3) - co_low(i)%vec(3))/2/de - y(j)
      endif


      if(plot_flag == cbar$)then
       call c_to_cbar(ring%ele(i),cbar_mat)
       xx_diff(l) = cbar_mat(1,2) -x(j)
       yy_diff(l) = cbar_mat(2,2) -y(j)
      endif

      if(plot_flag == de_cbar$)then
       call c_to_cbar(ring_two(1)%ele(i),cbar_mat1)
       call c_to_cbar(ring_two(-1)%ele(i),cbar_mat2)
       xx_diff(l) = (cbar_mat1(1,2)-cbar_mat2(1,2))/2/de - x(j)
       yy_diff(l) = (cbar_mat1(2,2)-cbar_mat2(2,2))/2/de - y(j)
      endif




      xmax = max(abs(xx_diff(l)),xmax)
      ymax = max(abs(yy_diff(l)),ymax)
      if(xmax0 /= 0.)xmax=xmax0
      if(ymax0 /= 0.)ymax=ymax0

      if(index(ring%ele(i)%name, 'DET') /= 0)then
        nd = nd+1
        zdet(nd) = zz_diff(l)
        xdet(nd) = xx_diff(l)
        ydet(nd) = yy_diff(l)
      endif
       print *,ring%ele(i)%name, l, zz_diff(l), xx_diff(l), yy_diff(l)       
      end do

     end do
    
     x(0:l) = xx_diff(0:l)
     y(0:l) = yy_diff(0:l)
     z(0:l) = zz_diff(0:l)
     endif


  endif !      if(device_type(1:7)  == '/XSERVE')then

!      call pgpage

!       if(n==1)then
       call pgsci(1)
       if(plot_flag==orbit$)then
         p = int(log10(xmax))
         if(p<=0)p=p-1
         f=xmax/10**p
         xscale=(int(f*2+1)/2.)*10**p
         print *,' p,f, xmax, xscale ',p,f, xmax, xscale
         call pgenv(start, end,-xscale,xscale,0,1)
         call pglab('z (m)','x(mm)',' Closed orbit')
       endif
       if(plot_flag == beta$)then
         xscale=(int(xmax/10.)+1)*10
         x_low = 0.
         if(diff)x_low = -xscale
         call pgenv(start, end,x_low,xscale,0,1)
         call pglab('z (m)','Bx(m)',' Beta')
       endif
       if(plot_flag == de_beta$)then
         xscale=(int(xmax/5.)+1)*5

         call pgenv(start, end,-xscale,xscale,0,1)
         call pglab('z (m)','dBx/dE(m)',' dBeta/dE')
       endif
       if(plot_flag == eta$)then
         xscale=(int(xmax/5.)+1)*5

         call pgenv(start, end,-xscale,xscale,0,1)
         call pglab('z (m)','etax',' eta')
       endif
       if(plot_flag == cbar$)then
         xscale=(int(xmax/0.1)+1)*0.1

         call pgenv(start, end,-xscale,xscale,0,1)
         call pglab('z (m)','cbar12',' cbar')
       endif

       if(plot_flag == de_cbar$)then
         xscale=(int(xmax/0.1)+1)*0.1

         call pgenv(start, end,-xscale,xscale,0,1)
         call pglab('z (m)','d(cbar12)/dE',' cbar')
       endif
!       endif

!       do i=1,ring%n_ele_track
       do i=1,l
         zz(i,n)=z(i)
         xx(i,n)=x(i)
       end do
!         n_ele(n) = ring%n_ele_track
         n_ele(n) = l

       do j =1,n
         call pgsci(j)
         forall(i=1:n_ele(j))z(i)=zz(i,j)
         forall(i=1:n_ele(j))x(i)=xx(i,j)
         call pgline(n_ele(j), z, x)
       end do
       call pgpt(nd, zdet, xdet, 18)
       call pgmtxt('T',3.,0.,0.,comment)

       call pgsci(1)
       if(plot_flag == orbit$)then
         p = int(log10(ymax))
         if(p<=0)p=p-1
         f=ymax/10**p
         yscale=(int(f*2+1)/2.)*10**p
         print *,' p,f, ymax, yscale ',p,f, ymax, yscale
         call pgenv(start, end,-yscale,yscale,0,1)
         call pglab('z (m)','y(mm)',' Closed orbit')
       endif
       if(plot_flag == beta$)then
         yscale=(int(ymax/10.)+1)*10
         y_low = 0.
         if(diff)y_low=-yscale
         call pgenv(start, end, y_low, yscale,0,1)
         call pglab('z (m)','By(m)',' Beta')
       endif
       if(plot_flag == de_beta$)then
         yscale=(int(ymax/5.)+1)*5

         call pgenv(start, end,-yscale,yscale,0,1)
         call pglab('z (m)','dBy/dE(m)',' dBeta/dE')
       endif
       if(plot_flag == eta$)then
         yscale=(int(ymax/0.5)+1)*0.5

         call pgenv(start, end,-yscale,yscale,0,1)
         call pglab('z (m)','etay',' eta')
       endif
       if(plot_flag == cbar$)then 
         yscale=(int(ymax/0.1)+1)*0.1

         call pgenv(start, end,-yscale,yscale,0,1)
         call pglab('z (m)','cbar22',' cbar')
       endif
       if(plot_flag == de_cbar$)then 
         yscale=(int(ymax/0.1)+1)*0.1

         call pgenv(start, end,-yscale,yscale,0,1)
         call pglab('z (m)','d(cbar22)/dE',' cbar')
       endif

!       do i=1,ring%n_ele_track
       do i=1,l
         zz(i,n)=z(i)
         yy(i,n)=y(i)
       end do

       do j =1,n
         call pgsci(j)
         forall(i=1:n_ele(j))z(i)=zz(i,j)
         forall(i=1:n_ele(j))y(i)=yy(i,j)
         call pgline(n_ele(j), z, y)
       end do

         call pgpt(nd, zdet, ydet, 18)
!     endif
    

!     answer = ' '
!     print '(a,$)',' Write orbit ?', answer
!     accept *,answer
!     if(answer(1:1) == 'y' .or. answer(1:1) == 'Y')exit

      if(device_type(1:7) /= '/XSERVE') then
       print *,' write ',device_type(1:index(device_type,'/')-1)
       if(istat1 > 0)call pgslct(istat1)
       call pgclos
      endif

 end do
 


      print *,' write orbit data to fort.33'
      print *,' plot orbit with orbit.pcm and orbit_ir.pcm'
      open(unit=33)
      write (33,*) ' z(meters), x,xp,y,yp (mm,mrad)'
      write(33,2)
2     format(1x,'ele',14x,'z',9x,'x',5x,'xp',6x,'y',7x,'yp',10x,'l',7x,'energy')

       do i=1,ring%n_ele_track
  write(33,1)ring%ele(i)%name,ring%ele(i)%s,(co(i)%vec(j)*1000.,j=1,4),   &
                       co(i)%vec(5), co(i)%vec(6)
       end do
1      format(1x,a13,f8.3,4f8.2,2e12.4) 
       close(unit=33)



      end

  subroutine de_dbeta(ring_high, ring_low, de, rms_x, rms_y)
  use bmad_struct
  use bmad_interface

 implicit none

 type (lat_struct) ring_high, ring_low

 real(rp) de, rms_x, rms_y, avg_x, avg_y,sum_x,sum_y

  integer i

   sum_x=0.
   sum_y=0.
   do i=1,ring_high%n_ele_track
      sum_x = sum_x + (ring_high%ele(i)%a%beta - ring_low%ele(i)%a%beta)/2/de
      sum_y = sum_y + (ring_high%ele(i)%b%beta - ring_low%ele(i)%b%beta)/2/de
   end do
    avg_x = sum_x/ring_high%n_ele_track
    avg_y = sum_y/ring_high%n_ele_track

    sum_x=0.
    sum_y=0.

   do i=1,ring_high%n_ele_track
      sum_x = sum_x + ( (ring_high%ele(i)%a%beta - ring_low%ele(i)%a%beta)/2/de -  avg_x)**2
      sum_y = sum_y + ( (ring_high%ele(i)%b%beta - ring_low%ele(i)%b%beta)/2/de -  avg_y)**2
   end do

   rms_x = sqrt(sum_x/ring_high%n_ele_track)
   rms_y = sqrt(sum_y/ring_high%n_ele_track)

   return
   end

  subroutine list_commands
    implicit none

    print *
    print *,' "READ" : to read another lattice into ring_2'
    print *,' "RING_1(2)" : switch to ring_1(2). '
    print *,'        Ring_1(2) is the ring structure for the lattice read at startup'
    print *,' "RADIATION ON(OFF)" : turn radiation damping and fluctuations on(off)'
    print *,' "CBAR_V_E ON(OFF)" : turn cbar vs energy calc on(off)'
    print *,' "TRANSFER ON(OFF)" :transfer line mode (ON) or closed ring (OFF)'
    print *,' "6D" :compute 6-dimensional closed orbit'
    print *,' "4D" :compute 4-dimensional closed orbit'
    print *,' "TRACK" :track and plot phase space'
    print *,' At <Plot> prompt:'
    print *,'                  type "PS" or "GIF" for hardcopy of last plot'
    print *,'                  type "<data_type>  X  <x_min>  <x_max> " to'
    print *,'                   set xrange' 
    print *,'                  type "<data_type>  Y  <y_up>  <y_low> " to'
    print *,'                   set absolute yrange for upper and lower plots' 

    print *
    return
  end
 
  subroutine sextupole_resonance(ring, rate_x, rate_y, rate_xq, rate_yq, rate_x_tot, rate_y_tot, delta_e)
  use bmad_struct
  use bmad_interface

 implicit none

 type (lat_struct) ring

 real(rp) rate_x, rate_y, delta_e, rate_xq, rate_yq
 real(rp) xf, yf,sum_x_real, sum_x_imagine, sum_y_real, sum_y_imagine
 real(rp) xfq, yfq,sum_x_realq, sum_x_imagineq, sum_y_realq, sum_y_imagineq
 real(rp) rate_x_tot, rate_y_tot
 real(rp) frev

  integer i

   frev=c_light/ring%ele(ring%n_ele_track)%s
   sum_x_real =0.
   sum_x_imagine =0.
   sum_y_real =0.
   sum_y_imagine =0.
   sum_x_realq =0.
   sum_x_imagineq =0.
   sum_y_realq =0.
   sum_y_imagineq =0.

   do i=1,ring%n_ele_track
     if(ring%ele(i)%key == sextupole$)then
     xf = ring%ele(i)%value(k2$) * ring%ele(i)%a%eta *ring%ele(i)%a%beta
     yf = ring%ele(i)%value(k2$) * ring%ele(i)%a%eta *ring%ele(i)%b%beta

     sum_x_real = sum_x_real+ xf*cos(2*ring%ele(i)%a%phi)
     sum_x_imagine = sum_x_imagine +xf*sin(2*ring%ele(i)%a%phi)
     sum_y_real = sum_y_real + yf*cos(2*ring%ele(i)%b%phi)
     sum_y_imagine = sum_y_imagine + yf*sin(2*ring%ele(i)%b%phi)
!!     print '(a16,a16,i,2e12.4)',' name, i, xf, yf', ring%ele(i)%name, i, xf, yf
   
    elseif (ring%ele(i)%key == quadrupole$ .or. index(ring%ele(i)%name,'Q01')/= 0) then
     xfq = ring%ele(i)%value(k1$) *ring%ele(i)%a%beta
     yfq = ring%ele(i)%value(k1$) *ring%ele(i)%b%beta
     sum_x_realq = sum_x_realq+ xfq*cos(2*ring%ele(i)%a%phi)
     sum_x_imagineq = sum_x_imagineq +xfq*sin(2*ring%ele(i)%a%phi)
     sum_y_realq = sum_y_realq + yfq*cos(2*ring%ele(i)%b%phi)
     sum_y_imagineq = sum_y_imagineq + yfq*sin(2*ring%ele(i)%b%phi)
 !    print '( a16,2e12.4,i)', ring%ele(i)%name, xfq, sum_x_realq, ring%ele(i)%key
   
   endif   
  
   end do

   rate_x = sqrt(sum_x_real**2 + sum_x_imagine**2)* frev * delta_e
   rate_y = sqrt(sum_y_real**2 + sum_y_imagine**2)* frev * delta_e
   rate_xq = sqrt(sum_x_realq**2 + sum_x_imagineq**2)* frev * delta_e
   rate_yq = sqrt(sum_y_realq**2 + sum_y_imagineq**2)* frev * delta_e

   rate_x_tot = sqrt((sum_x_real+sum_x_realq)**2 + (sum_x_imagine+sum_x_imagineq)**2)* frev * delta_e
   rate_y_tot = sqrt((sum_y_real+sum_y_realq)**2 + (sum_y_imagine+sum_y_imagineq)**2)* frev * delta_e

   return
   end

   subroutine  num_words(line, ix)
   implicit none
   character*(*) line
   integer i, ix

   ix=0
   call string_trim(line, line, i)
   do while (i /= 0 )
     ix = ix+1
     call string_trim(line(i+1:), line, i)
   end do
   return
   end

  subroutine psp(ring, co, traj, n_turns, istat2)
  use bmad
  implicit none
  type (lat_struct) ring
  type (coord_struct) traj, co
  type (coord_struct), allocatable, save :: co_lat(:)

  real*4, allocatable, save :: psx(:), psxp(:), psy(:), psyp(:)
  real*4 xscale/0./, xscalep/0./, yscale/0./,yscalep/0./

  integer n_turns, i, pgopen, istat2
  integer number_turns
  integer icall
  logical first1/.true./,first2/.true./
  logical already_open/.false./

  call reallocate_coord(co_lat, ring%n_ele_max)

  allocate(psx(n_turns), psxp(n_turns), psy(n_turns), psyp(n_turns))


  co_lat(0)%vec = traj%vec
  do i=1,n_turns
    call track_all(ring, co_lat)
    if(ring%param%lost)exit
    co_lat(0)%vec = co_lat(ring%n_ele_track)%vec
    psx(i) = (co_lat(0)%vec(1) - co%vec(1))*1000.
    psxp(i) = (co_lat(0)%vec(2) - co%vec(2))*1000.
    psy(i) = (co_lat(0)%vec(3) - co%vec(3))*1000.
    psyp(i) = (co_lat(0)%vec(4) - co%vec(4))*1000.
  end do
  number_turns = i-1
  print *,' Number of turns = ', number_turns
  traj%vec = co_lat(0)%vec

  xscale=0.
  xscalep=0.
  do i=1,number_turns
    xscale = max(abs(psx(i)),xscale)
    xscalep = max(abs(psxp(i)),xscalep)
  end do

  yscale=0.
  yscalep=0.
  do i=1,number_turns
    yscale = max(abs(psy(i)),yscale)
    yscalep = max(abs(psyp(i)),yscalep)
  end do

  if(.not. already_open)then
     istat2 =  pgopen('/XSERVE')
     if(istat2 >= 0)already_open = .true.
     icall = 1
  endif
!  icall = icall +1
  call pgslct(istat2)
  call pgpap(4.5,1.6)
  call pgsubp(1,2)
  call pgsch(2.)     
  call pgask(.false.)
  call pgsci(icall)

  call pgenv(-yscale, yscale, -yscalep, yscalep,0,1)
  call pglab('y(mm)','yp',' vertical phase space')
  call pgpt(number_turns, psy, psyp, 18)

  call pgenv(-xscale, xscale, -xscalep, xscalep,0,1)
  call pglab('x(mm)','xp',' horizontal phase space')
  call pgpt(number_turns, psx, psxp, 18)

  deallocate(psx, psxp, psy, psyp)

  return
  end







