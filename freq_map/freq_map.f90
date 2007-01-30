!........................................................................
!+
! program    :  freq_map
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
!
! $Log$
! Revision 1.3  2007/01/30 16:14:31  dcs
! merged with branch_bmad_1.
!
!
!
!........................................................................
!
#include "CESR_platform.h"

  program freq_map

  use bmad
!  use bmad_interface
  use bmadz_interface
  use bsim_interface
  use nr
  use nrutil, only: dpc
  use scan_parameters

  implicit none

  type (lat_struct) ring, ring_in
  type (coord_struct), allocatable :: co(:), orbit(:), data(:)
  type (coord_struct) start_coord, orb
  type (normal_modes_struct) mode
  type (scan_params_struct) scan_params

  integer i, j, k, isign, r, s, g
  integer ix
  integer n_turn, particle, i_train, j_car, n_trains_tot, n_cars, n_part, slices
  integer ios
  integer i_x, i_y, i_z
  integer n_x, n_y, n_z
  integer ix_pn, ix_dot

  real(rdef) eps_rel(4), eps_abs(4)
  real(rdef), allocatable :: tune(:,:,:)
  real(rdef) Qx, Qy, Qz, current, Q, d, b, c, A
  real(rdef) sig_in(3) !sigx, sigy, sigz, initial distribution
  real(rdef) coupling_sb, coupling_wb
  real(rdef) x0, y0, e0, x1, y1, e1, dx, dy, de
  real(rp) fft1(4096), ft1a(4096), ft1b(4096), fft2(4096), ft2a(4096),  ft2b(4096)
  real(rp) fft3(4096), ft3a(4096), ft3b(4096), fft4(4096), ft4a(4096),  ft4b(4096)
  real(rp) fft5(4096), ft5a(4096), ft5b(4096), fft6(4096), ft6a(4096), ft6b(4096)
  complex(dpc) n1(1, 4096), m1(1, 4096),  n2(1, 4096), m2(1, 4096), n3(1, 4096), m3(1, 4096), n4(1, 4096)
  complex(dpc) cfftx(1, 4096), cffty(1, 4096), cffte(1, 4096)
  real(rp) cftxa(4096), cftxb(4096), cftya(4096), cftyb(4096), cftea(4096), cfteb(4096)

  character*40 lattice
  character*60 lat_file
  character*80 line, last_line, file_name, temp_name, prefix_name
  character*1 answer
  character*60 in_file
  character*2 wordx, wordy, wordz

  logical keep_trying/.true./
  logical write_orbit/.false./                                        
  logical beambeam_ip, close_pretz, close_vert, lrbbi, rec_taylor, radiation,  go
  logical ok

  namelist / parameters /lat_file, Qx, Qy, Qz, &
                         x0, y0, e0, x1, y1, e1, dx, dy, de, &
                n_turn, particle, &
                i_train, j_car, n_trains_tot, n_cars, current, &
                lrbbi, rec_taylor, beambeam_ip, close_pretz, close_vert, &
                slices, radiation, sig_in, coupling_sb, go

!
  do   !read file with input parameters

    type '(a, $)', ' Input command file <CR=freq_map.in>: '
    read  '(a)', file_name
    call string_trim (file_name, file_name, ix)
    if (ix .eq. 0) file_name = 'freq_map.in'

    open (unit= 1, file = file_name, status = 'old', iostat = ios)

    if (ios == 0) then
      exit
    else
      type *
      type *, 'ERROR: CANNOT OPEN FILE: ', trim(file_name)
    endif
 
  enddo

  particle = positron$
  lrbbi = .false.
  beambeam_ip = .false.
  close_pretz = .false.
  close_vert = .false.
  read(1, nml = parameters)
  
  scan_params%Q_z    =        Qz 
  scan_params%lat_file =      lat_file 
  scan_params%n_turn =        n_turn 
  scan_params%particle =      particle
  scan_params%i_train =       i_train 
  scan_params%j_car =         j_car 
  scan_params%n_trains_tot =  n_trains_tot 
  scan_params%n_cars =        n_cars 
  scan_params%current =       current
  scan_params%lrbbi =         lrbbi
  scan_params%beambeam_ip =   beambeam_ip 
  scan_params%close_pretz =   close_pretz
  scan_params%close_vert =    close_vert 
  scan_params%slices =        slices 
  scan_params%rec_taylor =    rec_taylor
 ! scan_params%radiation  =    radiation
  scan_params%sig_in     =    sig_in
  scan_params%coupling_sb = coupling_sb

   type *, ' lat_file = ', lat_file

   call bmad_parser (lat_file, ring)

   call reallocate_coord(orbit, ring%n_ele_max)
   call reallocate_coord(co, ring%n_ele_max)

! in the next do loop, tune lattice parameters

   if(go)keep_trying=.false.
   do while (keep_trying)

    type '(a, $)', ' FREQ_MAP: element change or GO> '
     read  '(a)', line
    
     ix = index(line, '!')
     if (ix /= 0) line = line(:ix-1)        ! strip off comments

     call str_upcase(line, line)
     call string_trim(line, line, ix)

     if (ix == 0) then       ! nothing typed. do the same thing
       line = last_line
     endif

     last_line = line

     if(line(1:1) .eq. 'G')exit

     call find_change( line, ring)

  end do

! after tuning

  call set_on (rfcavity$, ring, .false.)

  ring%param%particle = particle

  call twiss_at_start(ring)
  co(0)%vec = 0.
  call closed_orbit_at_start(ring, co(0), 4, .true.)
  call track_all (ring, co)
  call lat_make_mat6(ring,-1,co)


  call twiss_at_start(ring)
  call twiss_propagate_all (ring)
  type *
  type *,' FREQ_MAP: Before parasitic added '
  type *,'    Qx = ',ring%a%tune/twopi,'    Qy = ',ring%b%tune/twopi
  type '(a15,4e12.4)','  Closed orbit ', co(0)%vec(1:4)


  if(lrbbi)then
   ring_in = ring
   call lrbbi_setup ( ring_in, ring, particle, i_train, j_car, n_trains_tot, n_cars, current, rec_taylor)
  endif
   
  call twiss_at_start(ring)
  co(0)%vec = 0.
  call closed_orbit_at_start(ring, co(0), 4, .true.)

  type *
  type *,' FREQ_MAP: After parasitic added '
  type *,'    Qx = ',ring%a%tune/twopi,'    Qy = ',ring%b%tune/twopi
  type '(a15,4e12.4)','  Closed orbit ', co(0)%vec(1:4)
  
  ring%z%tune = Qz * twopi
  if(beambeam_ip)call beambeam_setup(ring, particle, current, scan_params, slices)

  call twiss_at_start(ring)
  co(0)%vec = 0.
  call closed_orbit_at_start(ring, co(0), 4, .true.)

  type *
  type *,' FREQ_MAP: After beambeam added '
  type *,'    Qx = ',ring%a%tune/twopi,'    Qy = ',ring%b%tune/twopi
  type '(a15,4e12.4)','  Closed orbit ', co(0)%vec(1:4)


  if(close_pretz)call close_pretzel (ring, 6)
!  if(close_vert)call close_vertical(ring)

  forall( i=1:ring%n_ele_track) orbit(i)%vec = 0.

  call twiss_at_start(ring)

   type *
   type *,' After CLOSE PRETZEL ' 
  type *,'    Qx = ',ring%a%tune/twopi,'    Qy = ',ring%b%tune/twopi
  call closed_orbit_at_start(ring, orbit(0), 4, .true.)
  call track_all (ring, orbit)
  call lat_make_mat6(ring, -1, orbit)
  call twiss_at_start(ring)
  call twiss_propagate_all(ring)

  call q_tune(ring, Qx, Qy, ok)
  if(ok)print '(1x,a7,a6,f10.4,a6,f10.4)','Qtune: ',' Qx = ',ring%a%tune/twopi,' Qy = ', ring%b%tune/twopi

  call set_on (rfcavity$, ring, .true.)
  call set_z_tune(ring)

  do i = 1, ring.n_ele_max
    if(ring.ele(i).value(x_limit$) == 0.)ring.ele(i).value(x_limit$) = 0.05
    if(ring.ele(i).value(y_limit$) == 0.)ring.ele(i).value(y_limit$) = 0.05
  enddo

 
  n_x = nint(abs((x1 - x0)/dx))+1
  n_y = nint(abs((y1 - y0)/dy))+1
  if(de /= 0.)n_z = nint(abs((e1 - e0)/de))+1
  if(de == 0.)n_z=1

  allocate(tune(n_x*n_y*n_z, 3, 2))

  dx = (x1-x0)/(n_x-1)
  dy = (y1-y0)/(n_y-1)
  if(n_z /= 1)de = (e1-e0)/(n_z-1)

  ix_dot = index(file_name,'.')
  prefix_name = file_name(1:ix_dot-1)
  call string_trim(prefix_name, prefix_name,ix_pn)
  call file_suffixer (file_name, in_file, '.out', .true.)
!  open(unit=14, file = in_file)

!  write(14,'(a14,a)')'  Input file :',file_name
!  write(14,'(a10,f10.3,a11,f11.4)') &
!          ' beta_a = ',ring%ele(0)%a%beta, ' alpha_a = ', ring%ele(0)%a%alpha 
!  write(14,'(a10,f10.3,a11,f11.4)') &
!          ' beta_b = ',ring%ele(0)%b%beta, ' alpha_y = ', ring%ele(0)%b%alpha 
!  write(14,'(a10,f10.3,a11,f11.4)') &
!          '  eta_a = ',ring%ele(0)%a%eta, '  eta_ap = ', ring%ele(0)%a%etap
!  write(14,'(a10,f10.3,a11,f11.4)') &
!          '  eta_b = ',ring%ele(0)%b%eta, '  eta_bp = ', ring%ele(0)%b%etap
  open(unit=13, file= 'uber.out')


  print '(/,1x,3(a6,e12.4,4x))',' x0 = ',x0, ' y0 = ',y0, ' e0 = ',e0
  print '(1x,3(a6,e12.4,4x))',' dx = ',dx, ' dy = ',dy, ' de = ',de
  print '(1x,3(a6,i4,6x,4x))','n_x = ',n_x,'n_y = ',n_y,'n_z = ',n_z
  print '(1x,(a9,i),/)','n_turn = ',n_turn
  g= 0
  do i_x = 1,n_x
   start_coord%vec(1:6) = 0.
   start_coord%vec(1) = x0 + dx * i_x
   do i_y = 1,n_y
     start_coord%vec(3) = y0 + dy * i_y
     do i_z = 1,n_z
       start_coord%vec(6) = e0 + de * i_z
    g= g+1
    write(wordx,'(i2.2)')i_x
    write(wordy,'(i2.2)')i_y
    write(wordz,'(i2.2)')i_z
    temp_name = prefix_name(1:ix_pn)//'_'//wordx//'_'//wordy//'_'//wordz//'.out'

!     open(unit=16, file=temp_name)

       co(0)%vec = orbit(0)%vec + start_coord%vec
!       write(14,*)' !'
!       write(14,'(23x,a26,23x,10x,23x,a27)')'phase space at IP at start', &
!                                'phase space at IP at turn n'
!       write(14,'(1x,6a12,6x,a4,6a12)')'    x       ','     xp     ','     y      ', &
!                                       '   yp       ','     dl     ','    dE/E    ', &
!                             'turn','    x       ','     xp     ','     y      ', &
!                                       '   yp       ','     dl     ','    dE/E    '

       do j =1, n_turn
         call track_all(ring, co)
         if(ring%param%lost)then
            type *,' Particle lost in turn ',j
            exit
         endif
         co(0)%vec = co(ring%n_ele_track)%vec

   fft1(j)= co(0)%vec(1) -orbit(0)%vec(1)
   fft2(j)= co(0)%vec(2) -orbit(0)%vec(2)
   fft3(j)= co(0)%vec(3) -orbit(0)%vec(3)
   fft4(j)= co(0)%vec(4) -orbit(0)%vec(4)
   fft5(j)= co(0)%vec(5) -orbit(0)%vec(5)
   fft6(j)= co(0)%vec(6) -orbit(0)%vec(6)
   
   cfftx(1, j)= cmplx(fft1(j), fft2(j))*(2*(sin((pi*j)/n_turn))*(sin((pi*j)/n_turn)))
   cffty(1, j)= cmplx(fft3(j), fft4(j))*(2*(sin((pi*j)/n_turn))*(sin((pi*j)/n_turn)))
   cffte(1, j)= cmplx(fft5(j), -fft6(j))*(2*(sin((pi*j)/n_turn))*(sin((pi*j)/n_turn)))

end do

  if(.not. ring%param%lost)then

   n1(1, 1:n_turn/2)= cfftx(1, 1:n_turn/2)
   m1(1, 1:n_turn/2)= cfftx(1, (n_turn/2+1):n_turn)
   n2(1, 1:n_turn/2)= cffty(1, 1:n_turn/2)
   m2(1, 1:n_turn/2)= cffty(1, (n_turn/2+1):n_turn)
   n3(1, 1:n_turn/2)= cffte(1, 1:n_turn/2)
   m3(1, 1:n_turn/2)= cffte(1, (n_turn/2+1):n_turn)



isign= 1

 call fourrow(n1(:,1:n_turn/2), isign)
     forall(i=1:n_turn/2)cftxa(i)=sqrt(n1(1,i)*conjg(n1(1,i)))

 call fourrow(m1(:,1:n_turn/2), isign)
     forall(i=1:n_turn/2)cftxb(i)=sqrt(m1(1,i)*conjg(m1(1,i)))

 call fourrow(n2(:,1:n_turn/2), isign)
     forall(i=1:n_turn/2)cftya(i)=sqrt(n2(1,i)*conjg(n2(1,i)))

 call fourrow(m2(:,1:n_turn/2), isign)
     forall(i=1:n_turn/2)cftyb(i)=sqrt(m2(1,i)*conjg(m2(1,i)))

 call fourrow(n3(:,1:n_turn/2), isign)
     forall(i=1:n_turn/2)cftea(i)=sqrt(n3(1,i)*conjg(n3(1,i)))

 call fourrow(m3(:,1:n_turn/2), isign)
     forall(i=1:n_turn/2)cfteb(i)=sqrt(m3(1,i)*conjg(m3(1,i)))

  Q= 1

  do k= 2, n_turn/2
   if (cftxa(Q)<cftxa(k)) Q= k
  end do

  d= cftxa(Q)
  b= cftxa(Q+1)
  c= cos(twopi/(n_turn/2))
  A= (-(d+b*c)*(d-b)+b*sqrt(c*c*(d+b)*(d+b)-2*d*b*(2*c*c-c-1)))/(d*d+b*b+2*d*b*c)

  tune(g, 1, 1)= Q/(n_turn/2)+(1/twopi)*asin(A*sin(twopi/(n_turn/2)))

!tune of x for 2nd half etc.
Q= 1
do k= 2, n_turn/2

if (cftxb(Q)<cftxb(k)) Q= k
end do
d= cftxb(Q)
b= cftxb(Q+1)
c= cos(twopi/(n_turn/2))
A= (-(d+b*c)*(d-b)+b*sqrt(c*c*(d+b)*(d+b)-2*d*b*(2*c*c-c-1)))/(d*d+b*b+2*d*b*c)

tune(g, 1, 2)= Q/(n_turn/2)+(1/twopi)*asin(A*sin(twopi/(n_turn/2)))

!find tune for y first half etc.
Q= 1
do k= 2, n_turn/2

if (cftya(Q)<cftya(k)) Q= k
end do
d= cftya(Q)
b= cftya(Q+1)
c= cos(twopi/(n_turn/2))
A= (-(d+b*c)*(d-b)+b*sqrt(c*c*(d+b)*(d+b)-2*d*b*(2*c*c-c-1)))/(d*d+b*b+2*d*b*c)

tune(g, 2, 1)= Q/(n_turn/2)+(1/twopi)*asin(A*sin(twopi/(n_turn/2)))

!find tune for y second half etc.
Q= 1
do k= 2, n_turn/2

if (cftyb(Q)<cftyb(k)) Q= k
end do
d= cftyb(Q)
b= cftyb(Q+1)
c= cos(twopi/(n_turn/2))
A= (-(d+b*c)*(d-b)+b*sqrt(c*c*(d+b)*(d+b)-2*d*b*(2*c*c-c-1)))/(d*d+b*b+2*d*b*c)

tune(g, 2, 2)= Q/(n_turn/2)+(1/twopi)*asin(A*sin(twopi/(n_turn/2)))

!find tune for z first half etc.
Q= 1
do k= 2, n_turn/4

if (cftea(Q)<cftea(k)) Q= k
end do
d= cftea(Q)
b= cftea(Q+1)
c= cos(twopi/(n_turn/2))
A= (-(d+b*c)*(d-b)+b*sqrt(c*c*(d+b)*(d+b)-2*d*b*(2*c*c-c-1)))/(d*d+b*b+2*d*b*c)

tune(g, 3, 1)= Q/(n_turn/2)+(1/twopi)*asin(A*sin(twopi/(n_turn/2)))

!find tune for z second half etc.
Q= 1
do k= 2, n_turn/2

if (cfteb(Q)<cfteb(k)) Q= k
end do
d= cfteb(Q)
b= cfteb(Q+1)
c= cos(twopi/(n_turn/2))
A= (-(d+b*c)*(d-b)+b*sqrt(c*c*(d+b)*(d+b)-2*d*b*(2*c*c-c-1)))/(d*d+b*b+2*d*b*c)

tune(g, 3, 2)= Q/(n_turn/2)+(1/twopi)*asin(A*sin(twopi/(n_turn/2)))

write(13, '(12e14.5)') start_coord%vec(1), start_coord%vec(3), start_coord%vec(6), &
                               tune(g, :, 1), tune(g, :, 2), tune(g,:,2)-tune(g,:,1)


!  close(unit=16)

 endif !.not.ring%param%lost
   ring%param%lost = .false.

     end do !i_z
   end do !i_y
  end do !i_x
close(unit=13)

 end













