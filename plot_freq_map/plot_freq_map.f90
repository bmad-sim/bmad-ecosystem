   program plot_freq_map

   use column_mod
   implicit none
   
   type(column_struct), allocatable :: column(:)

   interface
    subroutine read_scan (file_name, column,size) 
     use column_mod
     implicit none
     type(column_struct), allocatable :: column(:)
     character*200 file_name
     integer size
    end  subroutine
   end interface

   character*200 file_name, title
   character*22 wedge_label
   character*20 device_type

   integer pgopen
   integer ix, iy, iz, ix_dim, iy_dim
   integer n_x, n_y, n_eng
   integer ix_j, iy_j
   integer j,n, k
   integer j_left, j_right, k_up, k_down
   integer ieng, ix_energy

   character*20 data_type
   character*6 ans/'N'/

      INTEGER I, L, C1, C2, NC
      REAL CONTRA, BRIGHT, ANGLE, C, S, ALEV(1)
      CHARACTER*16 VAL
 
!   parameter (n_x=31, n_y=31)
!    integer n_x, n_y
   real xmin, xmax, ymin, ymax, energy_min, energy_max
   real tr(6), fmin,fmax
   real*8 fmin_dp
   real*8 fmax_dp
   real*8 xmin_dp, xmax_dp, ymin_dp, ymax_dp
   real, allocatable :: z(:,:)
   real delta_x, delta_y, delta_eng
   real factor
   real, allocatable:: saved_x(:), saved_y(:), saved_eng(:)
   real*8 sigma_x, sigma_y
   real aspect
   real z_off_scale

   logical write_ps/.false./, write_gif/.false./
   logical fill_gaps/.false./   

    print '(a, $)', ' Freq map data file <CR=uber.out>: '
    read(*, '(a)') file_name
    call string_trim (file_name, file_name, ix)
    if(ix == 0)file_name = 'uber.out'

   call read_scan(file_name, column, n)

   allocate(saved_x(n), saved_y(n), saved_eng(n))

!   n = size(column)

   print *,' Column      parameter'
   print *,'   1           x[m]   '
   print *,'   2           y[m]   '
   print *,'   3       Delta E/E  '
   print *,'   4           Q_x    '
   print *,'   5           Q_y    '
   ix=1
   call query_int(' column for horizontal axis ', ix, 'i3')
   iy=2
   call query_int(' column for vertical axis ', iy, 'i3')
   ieng = 3  !energies

   sigma_x = 1.

!   call query_real(' sigma_x ? ', sigma_x, 'e12.4')
   column(1:n)%row(ix) = column(1:n)%row(ix)/sigma_x

   sigma_y = 1.
!   call query_real(' sigma_y ? ', sigma_y, 'e12.4')
   column(1:n)%row(iy) = column(1:n)%row(iy)/sigma_y

10 print '(a)', ' Plot: '
   print '(a)', ' Delta Qx (DQX) '
   print '(a)', ' Delta Qy (DQY) '
   print '(a)', ' Delta Qz (DQZ) '
   print '(a)', ' Qx (QX) '
   print '(a)', ' Qy (QY) '
   print '(a,$)', ' Plot ? '
   read (*,'(a)') data_type
   call str_upcase(data_type, data_type)
   iz=10
   if(data_type(1:3) == 'DQX')then
      iz=10
      wedge_label = 'Delta Qx'
      title = 'Delta Qx'
      factor = 1.  
   elseif(data_type(1:3) == 'DQY')then
      iz=11
      wedge_label = 'Delta Qy'
      title = 'Delta Qy'
      factor =1. 
   elseif(data_type(1:3) == 'DQZ')then
      iz=9
      wedge_label = 'Delta Qs'
      title = 'Delta Qs'
      factor = 1.
   elseif(data_type(1:2) == 'QX')then
      iz=4
      wedge_label = 'Qx'
      title = 'Qx'
      factor = 1.
   elseif(data_type(1:3) == 'QY')then
      iz=5
      wedge_label = 'Qy'
      title = 'Qy'
      factor = 1.
   else
      print *,' try again'
      goto 10
   endif


   if(iz >= 9) then
      z_off_scale = log(1.e-19)
      column(1:n)%row(iz)= log(abs(column(1:n)%row(iz)))
     else
      z_off_scale = 0.
   endif

!   call query_int(' column for z ', iz, 'i3')

! find max and min in x
   xmin = minval(column(1:n)%row(ix))
   xmax = maxval(column(1:n)%row(ix))
   ymin = minval(column(1:n)%row(iy))
   ymax = maxval(column(1:n)%row(iy))
   energy_min = minval(column(1:n)%row(ieng))
   energy_max = maxval(column(1:n)%row(ieng))

! find number of steps

   j=1
   k=1
   l=1

   saved_x(1)=column(1)%row(ix)
   saved_y(1)=column(1)%row(iy)
   saved_eng(1)=column(1)%row(ieng)
   do i = 1,n

    if(all(saved_x(1:j)-column(i)%row(ix) /= 0.))then
      j=j+1
      saved_x(j) = column(i)%row(ix)
!      type *,' j, saved_x(j-1), column(i)%row(ix)', j, saved_x(j-1), column(i)%row(ix)
     endif

    if(all(saved_y(1:k)-column(i)%row(iy) /= 0.))then
     k=k+1
     saved_y(k) = column(i)%row(iy)
!     type *,' k, saved_y(k) ', k, saved_y(k)
    endif

    if(all(saved_eng(1:l)-column(i)%row(ieng) /= 0.))then
     l=l+1
     saved_eng(l) = column(i)%row(ieng)
!     type *,' l, saved_eng(l) ', l, saved_eng(l)
    endif

!    if(j > 5000 .or. k > 5000 .or. l > 5000)then
!      print *,' More than 5000 different x,y,or energy values.'
!      print *,' Increase size of saved_ array'
!    endif

   end do
   n_x = j
   n_y = k
   n_eng = l

   n_x = min(200,n_x)
   n_y = min(200,n_y)

   print *,' n_x, n_y,_eng', n_x, n_y, n_eng
   print *,' Energies '
   print '(a16)',' n_eng, energy '
   print '(i5,3x,e12.4)',(i,saved_eng(i), i=1,n_eng)
   print '(a,$)', ' Choose an energy (integer) ?'
   read (*,'(i)') ,ix_energy

!   delta_x = (xmax-xmin)/(n_x-1)
!   delta_y = (ymax-ymin)/(n_y-1)
!   delta_eng = (energy_max-energy_min)/(n_eng-1)
  print *, ' xmin,xmax,ymin,ymax ,energy_min,energy_max', xmin,xmax,ymin,ymax,energy_min, energy_max

   xmin_dp = xmin
   xmax_dp = xmax
   ymin_dp = ymin
   ymax_dp = ymax
   call query_real(' xmin ? ', xmin_dp, 'f10.3')
   call query_real(' xmax ? ', xmax_dp, 'f10.3')
   call query_real(' ymin ? ', ymin_dp, 'f10.3')
   call query_real(' ymax ? ', ymax_dp, 'f10.3')
   xmin = xmin_dp
   xmax = xmax_dp
   ymin = ymin_dp
   ymax = ymax_dp
   print *,' xmin,xmax,ymin,ymax ', xmin,xmax,ymin,ymax

   delta_x = (xmax-xmin)/(n_x-1)
   delta_y = (ymax-ymin)/(n_y-1)
   delta_eng = (energy_max-energy_min)/(n_eng-1)

   print *,' delta_x, delta_y ', delta_x, delta_y

  allocate(z(n_x, n_y))
   do j=1,n_x
    do k=1,n_y
      z(j,k)=z_off_scale
    end do
   end do
   do j=1,n
    if(column(j)%row(ieng) /= saved_eng(ix_energy))cycle
     ix_j = (column(j)%row(ix)-xmin)/delta_x +1+ .5
      if(ix_j <= 0 .or. ix_j > n_x)cycle
     iy_j = (column(j)%row(iy)-ymin)/delta_y  +1+ .5
      if(iy_j <= 0 .or. iy_j > n_y)cycle
     z(ix_j,iy_j) = column(j)%row(iz)
      write(11, '(a18,3i8,1x,e12.4)')' j, ix_j, iy_j, z ',j, ix_j,iy_j,z(ix_j,iy_j)
     if( ix_j >n_x .or. iy_j > n_y)then
      print *,' ix_j, iy_j, z ', ix_j,iy_j,z(ix_j,iy_j)
      print *, column(j)%row(ix), xmin,delta_x, column(j)%row(iy), ymin, delta_y
    endif 

   end do 

    print '(a,$)', ' Fill gaps (Y/N) ?'
    read(5,*) data_type
    if(data_type(1:1) == 'y' .or. data_type(1:1) == 'Y')fill_gaps = .true.

 if(fill_gaps)then
! fill in gaps
   do j=1,n_x
    do k=1,n_y
     if(z(j,k) == z_off_scale)then
       do l = j,1,-1
         if(z(l,k) /= z_off_scale)then
            j_left = l
            exit
         endif
         j_left = 1
       end do

       do l=j,n_x
         if(z(l,k) /= z_off_scale)then
            j_right = l
            exit
         endif
         j_right = n_x
       end do

       if(j_right == n_x .and. j_left == 1)then  !look in y
        do l = k,1,-1
          if(z(j,l) /= z_off_scale)then
            k_down = l
            exit
          endif
          k_down = 1
        end do

        do l=k,n_y
          if(z(j,l) /= z_off_scale)then
            k_up = l
            exit
          endif
          k_up = n_y
        end do

        z(j,k) = (z(j,k_up) - z(j,k_down))/(k_up-k_down)*(k-k_down) +z(j,k_down)

       else
        z(j,k) = (z(j_right,k) - z(j_left,k))/(j_right-j_left)*(j-j_left) +z(j_left,k)
       endif
      endif
     end do
    end do
  endif !if fill_gaps

   fmin = minval(column(1:n)%row(iz))
   fmax = maxval(column(1:n)%row(iz))
   print *, ' fmin, fmax ', fmin, fmax
   fmin_dp = fmin
   fmax_dp = fmax
   call query_real(' fmin ? ', fmin_dp, 'f10.3')
   call query_real(' fmax ? ', fmax_dp, 'f10.3')
   fmin=fmin_dp
   fmax=fmax_dp
!
! Open device for graphics.
!
!      IF (PGOPEN('?') .LT. 1) STOP

      device_type ='/XSERVE'
!      device_type = '?'
 20  continue

      IF (PGOPEN(device_type) .LT. 1) STOP
      CALL PGQINF('TYPE', VAL, L)
      WRITE (*,*) 'PGPLOT device type: ', VAL(1:L)

      call pgpap(6.,1.)

      C1=0
      CALL PGQCIR(c1, C2)
      NC = MAX(0, C2-C1+1)
      WRITE (*,*) 'Number of color indices used for image: ', NC
      IF (NC .LT.8) THEN 
         WRITE (*,*) 'Not enough colors available on this device'
         STOP
      ELSE
         WRITE (*,*)
      END IF
!
!
! Set the coordinate transformation matrix: 
! world coordinate = pixel number.
!

      aspect = xmax/ymax
!      aspect = 1
!      if(ix == 1 .and. iy == 2) aspect = 10.
      print *,' aspect = ',aspect
      TR(1) = xmin -delta_x
      TR(2) = delta_x
      TR(3) = 0.
      TR(4) = (ymin-delta_y) * aspect
      TR(5) = 0.0
      TR(6) = delta_y * aspect
!
! Clear the screen. Set up window and viewport.
!
      CALL PGPAGE
      CALL SETVP
!      CALL PGWNAD(0.0, 1.0+n_x, 0.0, 1.0+n_y)
      CALL PGWNAD(xmin-delta_x/2,xmax+delta_x/2, &
               (ymin-delta_y/2)* aspect,(ymax+delta_x/2)* aspect)
!
! Set up the color map.
!
      BRIGHT = 0.5
      CONTRA  = 1.0
      CALL PALETT(2, CONTRA, BRIGHT)
!
! Draw the map with PGIMAG.  
!
!      call pgenv(xmin, xmax, ymin,ymax,1,0)
      CALL PGIMAG(z,n_x,n_y,1,n_x,1,n_y,FMIN,FMAX,TR)
!
! Annotate the plot.
!
      print '(a,$)',' Scan input file name ?'
      read(*,'(a)')file_name

      CALL PGMTXT('t',1.0,0.5,0.0,file_name)
      CALL PGMTXT('t', 1.0, 0.0, 0.0, title)
      CALL PGSCH(1.)
      CALL PGBOX('bcntsi1',0.0,0,'bcntsiv1',0.0,0)
      if(ix == 1)CALL PGMTXT('b',3.0,1.0,1.0,'horizontal displacement (m)')
      if(iy == 2)CALL PGMTXT('l',5.5,1.0,1.0,'vertical displacement X 10 (m)')
      if(ix == 4)CALL PGMTXT('b',3.0,1.0,1.0,'horizontal tune')
      if(iy == 5)CALL PGMTXT('l',5.5,1.0,1.0,'vertical tune')

! Draw a wedge.
!

      CALL PGWEDG('BI', 4.0, 5.0, FMIN, FMAX, wedge_label)
      CALL PGSCH(1.0)
!      CALL PGIDEN    
    
    if(index(ans,'N') /= 0)then
     print '(a, $)', ' post script[p] and/or gif[g] output ? : '
!     print '(a, $)', ' post script[p] ? : '
     read(*, '(a)')ans
     call str_upcase(ans, ans)
    endif
     if(index(ans,'P') /=0 .and. .not. write_ps)then
      device_type = 'scan.ps/VCPS'
      print *,' write scan.ps'
      write_ps = .true.
      goto 20
     endif
     if(index(ans,'G') /=0 .and. .not. write_gif)then
      device_type = 'scan.gif/VGIF'
      print *,' write scan.gif'
      write_gif = .true.
      goto 20
     endif
    call pgclos
end

!---------------------------------------------------
   subroutine read_scan (file_name, column, size) 

   use column_mod
   implicit none


   type(column_struct), allocatable, intent(out) :: column(:)
   character*200 file_name
   character*200 line
   character*20 word
   integer i, n, ix, j, size
   integer k

   open(unit=2, file=file_name, type='OLD')
   n=0
   i=0
   do while(.true.)
    read(2,'(a200)', end=99)line
!    type '(a6,a200)', ' line ',line
    i=i+1
    if(i > 0)n=n+1
   end do
99 continue
   close(unit=2)
   size = n 

   allocate(column(n))

   open(unit=2, file=file_name, type='OLD')
   n=0
   type *,' i ',i
!   i=0
   do k= 1,i !while(.true.)   
    read(2,'(a200)', end=100, err=100)line
!    i=i+1
    if(k > 0)then
     n=n+1
!    if(n > 180) type *, ' line n ',n,line
     ix=0
!     print *,' line',line
      do j=1,12
        call string_trim(line(ix+1:), line, ix)
        word=line(1:ix)
!        type *,'n, j, word', n, j, word
        if(index(word,'NaN') /= 0)then
          column(n)%row(j) = 0.
         else
          read(word,*)column(n)%row(j)
        endif
      end do
    endif
   end do
100 continue
  close(unit=2)
  return
 end
   
