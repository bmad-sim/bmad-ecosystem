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

  program freq_map

  use bmad_struct
  use bmad_interface
!  use bmadz_interface
!  use cesr_utils
  use nr
!  use bmad
  use nrutil, only: dpc

  implicit none

  type (ring_struct) ring, ring_in
  type (coord_struct), allocatable :: co_(:), orb_(:), data_(:)
  type (coord_struct) start_coord, orb
  type (modes_struct) mode

  integer i, j, k, isign, r, s, g
  integer ix
  integer n_turn, particle, i_train, j_car, n_trains_tot, n_cars, n_part
  integer ios
  integer i_x, i_y, i_z
  integer n_x, n_y, n_z
  integer ix_pn, ix_dot

  real(rdef) eps_rel(4), eps_abs(4)
  real(rdef), allocatable :: tune(:,:,:)
  real(rdef) Qx, Qy, Qz, current, Q, d, b, c, A
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
  logical beambeam_ip, close_pretz, close_vert, lrbbi, go
  logical ok

  namelist / parameters /lat_file, Qx, Qy, Qz, &
                         x0, y0, e0, x1, y1, e1, dx, dy, de, &
                n_turn, particle, &
                i_train, j_car, n_trains_tot, n_cars, current, &
                lrbbi, beambeam_ip, close_pretz, close_vert, go

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


   type *, ' lat_file = ', lat_file

   call bmad_parser (lat_file, ring)

   call reallocate_coord(orb_, ring%n_ele_max)
   call reallocate_coord(co_, ring%n_ele_max)

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
  co_(0)%vec = 0.
  call closed_orbit_at_start(ring, co_(0), 4, .true.)
  call track_all (ring, co_)
  call ring_make_mat6(ring,-1,co_)


  call twiss_at_start(ring)
  call twiss_propagate_all (ring)
  type *
  type *,' FREQ_MAP: Before parasitic added '
  type *,'    Qx = ',ring%x%tune/twopi,'    Qy = ',ring%y%tune/twopi
  type '(a15,4e12.4)','  Closed orbit ', co_(0)%vec(1:4)


  if(lrbbi)then
   ring_in = ring
   call lrbbi_setup ( ring_in, ring, particle, i_train, j_car, n_trains_tot, n_cars, current)
  endif
   
  call twiss_at_start(ring)
  co_(0)%vec = 0.
  call closed_orbit_at_start(ring, co_(0), 4, .true.)

  type *
  type *,' FREQ_MAP: After parasitic added '
  type *,'    Qx = ',ring%x%tune/twopi,'    Qy = ',ring%y%tune/twopi
  type '(a15,4e12.4)','  Closed orbit ', co_(0)%vec(1:4)
  
  ring%z%tune = Qz * twopi
  if(beambeam_ip)call beambeam_setup(ring, particle, current)

  call twiss_at_start(ring)
  co_(0)%vec = 0.
  call closed_orbit_at_start(ring, co_(0), 4, .true.)

  type *
  type *,' FREQ_MAP: After beambeam added '
  type *,'    Qx = ',ring%x%tune/twopi,'    Qy = ',ring%y%tune/twopi
  type '(a15,4e12.4)','  Closed orbit ', co_(0)%vec(1:4)


  if(close_pretz)call close_pretzel (ring, 6)
!  if(close_vert)call close_vertical(ring)

  forall( i=1:ring%n_ele_use) orb_(i)%vec = 0.

  call twiss_at_start(ring)

   type *
   type *,' After CLOSE PRETZEL ' 
  type *,'    Qx = ',ring%x%tune/twopi,'    Qy = ',ring%y%tune/twopi
  call closed_orbit_at_start(ring, orb_(0), 4, .true.)
  call track_all (ring, orb_)
  call ring_make_mat6(ring, -1, orb_)
  call twiss_at_start(ring)
  call twiss_propagate_all(ring)

  call q_tune(ring, Qx, Qy, ok)
  if(ok)print '(1x,a7,a6,f10.4,a6,f10.4)','Qtune: ',' Qx = ',ring%x%tune/twopi,' Qy = ', ring%y%tune/twopi

  call set_on (rfcavity$, ring, .true.)
  call set_z_tune(ring)

  do i = 1, ring.n_ele_max
    if(ring.ele_(i).value(x_limit$) == 0.)ring.ele_(i).value(x_limit$) = 0.05
    if(ring.ele_(i).value(y_limit$) == 0.)ring.ele_(i).value(y_limit$) = 0.05
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
!          ' beta_x = ',ring%ele_(0)%x%beta, ' alpha_x = ', ring%ele_(0)%x%alpha 
!  write(14,'(a10,f10.3,a11,f11.4)') &
!          ' beta_y = ',ring%ele_(0)%y%beta, ' alpha_y = ', ring%ele_(0)%y%alpha 
!  write(14,'(a10,f10.3,a11,f11.4)') &
!          '  eta_x = ',ring%ele_(0)%x%eta, '  eta_xp = ', ring%ele_(0)%x%etap
!  write(14,'(a10,f10.3,a11,f11.4)') &
!          '  eta_y = ',ring%ele_(0)%y%eta, '  eta_yp = ', ring%ele_(0)%y%etap
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

       co_(0)%vec = orb_(0)%vec + start_coord%vec
!       write(14,*)' !'
!       write(14,'(23x,a26,23x,10x,23x,a27)')'phase space at IP at start', &
!                                'phase space at IP at turn n'
!       write(14,'(1x,6a12,6x,a4,6a12)')'    x       ','     xp     ','     y      ', &
!                                       '   yp       ','     dl     ','    dE/E    ', &
!                             'turn','    x       ','     xp     ','     y      ', &
!                                       '   yp       ','     dl     ','    dE/E    '

       do j =1, n_turn
         call track_all(ring, co_)
         if(ring%param%lost)then
            type *,' Particle lost in turn ',j
            exit
         endif
         co_(0)%vec = co_(ring%n_ele_use)%vec

   fft1(j)= co_(0)%vec(1) -orb_(0)%vec(1)
   fft2(j)= co_(0)%vec(2) -orb_(0)%vec(2)
   fft3(j)= co_(0)%vec(3) -orb_(0)%vec(3)
   fft4(j)= co_(0)%vec(4) -orb_(0)%vec(4)
   fft5(j)= co_(0)%vec(5) -orb_(0)%vec(5)
   fft6(j)= co_(0)%vec(6) -orb_(0)%vec(6)
   
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













