   program plot_scan

   use column_mod
   implicit none
   
   type(column_struct), allocatable :: column(:)

   interface
    subroutine read_scan (file_name, column,size) 
     use column_mod
     implicit none
     type(column_struct), allocatable :: column(:)
     character*80 file_name
     integer size
    end  subroutine
   end interface

   character*80 file_name, title
   character*22 wedge_label
   character*20 device_type

   integer pgopen
   integer ix, iy, iz, ix_dim, iy_dim
   integer n_x, n_y
   integer ix_j, iy_j
   integer j,n, k
   integer j_left, j_right, k_up, k_down
   integer turns, max_turns
   integer word_len
   integer itest
   integer nx,ny

   character*20 data_type,hv,plane
   character*6 ans/'N'/

      INTEGER I, L, C1, C2, NC
      REAL CONTRA, BRIGHT, ANGLE, C, S, ALEV(1)
      CHARACTER*16 VAL
 
   parameter (n_x=26, n_y=26)

   real xmin, xmax, ymin, ymax
   real tr(6), fmin,fmax
   real*8 fmin_dp
   real*8 fmax_dp
   real, allocatable :: z(:,:)
   real delta_x, delta_y
   real qx, qy, qs

   logical write_ps/.false./, write_gif/.false./
   

    print '(a, $)', ' Scan data file <CR=sb_combined.out>: '
    read(*, '(a)') file_name
    call string_trim (file_name, file_name, ix)
    if(ix == 0)file_name = 'sb_combined.out'

   call read_scan(file_name, column, n)
  print '(a,$)',' Plot in terms of frequency (F) or fractional tune (T) ?'
  read (*,'(a)') data_type
  call str_upcase(data_type, data_type)

  if(data_type(1:1) == 'F')then  
   column(1:n)%row(1) = column(1:n)%row(1)*390.1
   column(1:n)%row(2) = column(1:n)%row(2)*390.1
   column(1:n)%row(3) = column(1:n)%row(3)*390.1
  endif

!   n = size(column)

   ix=1
   iy=2
   
10 print '(a,$)', ' Plot amplitude  (A), or log of amplitude (L) ?'
   read (*,'(a)') data_type
   call str_upcase(data_type, data_type)
   iz=12
   if(index(data_type, 'L') /= 0)then
      wedge_label = ''
      title = 'Log(Amplitude)'
   elseif(index(data_type,'A') /= 0)then
      wedge_label = ' '
      title = 'Amplitude'
   else
      print *,' try again'
      goto 10
   endif
   do while(.true.)
    print '(a,$)', ' Horizontal (H) or vertical (V) amplitude ?'
    read (*,'(a)') hv
    call str_upcase(hv,hv)
    call string_trim(title, title, word_len)
    if(index(hv,'H') /= 0)then
      iz  = 6
      title = title(1:word_len)//' Horizontal'
      exit
     elseif(index(hv,'V') /= 0)then
      iz = 7
      title = title(1:word_len)//' Vertical'
      exit
     else
      print *,' try again '
     endif
    end do     


   do while(.true.)

    print '(a,$)', ' Choose plane HV, LV, or LH ?'
    read (*,'(a)') plane
    call str_upcase(plane, plane)
    if(index(plane,'HV') /= 0)then
      ix =1
      iy=2
      qs=column(1)%row(3)
      call query_real (' synchrotron tune ', qs,'f5.3')
      exit
     elseif(index(plane,'LV') /= 0)then
      ix=2
      iy=3
      qx=column(1)%row(1)
      call query_real (' horizontal tune ', qx,'f5.3')
      exit
     elseif(index(plane,'LH') /= 0)then
      ix=1
      iy=3
      qy=column(1)%row(2)
      call query_real (' vertical tune ', qy,'f5.3')
      exit
     else
      print *,' try again '
     endif
    end do

    turns = 8     
  

!   call query_int(' column for z ', iz, 'i3')

  nx=0
  ny=0
! find max and min in x
   do i =1,n
     if(any(column(i)%row(ix) == column(1:i-1)%row(ix)))cycle
     nx = nx+1
   end do
   do i =1,n
     if(any(column(i)%row(iy) == column(1:i-1)%row(iy)))cycle
     ny = ny+1
   end do
   allocate(z(1:nx,1:ny))
   xmin = minval(column(1:n)%row(ix))
   xmax = maxval(column(1:n)%row(ix))
   ymin = minval(column(1:n)%row(iy))
   ymax = maxval(column(1:n)%row(iy))
   max_turns = maxval(column(1:n)%row(turns))
   delta_x = (xmax-xmin)/(nx-1)
   delta_y = (ymax-ymin)/(ny-1)
  print *, ' xmin,xmax,ymin,ymax,nx,ny', xmin,xmax,ymin,ymax,nx,ny

   fmin = minval(column(1:n)%row(iz))
   fmax = maxval(column(1:n)%row(iz))
   if(index(data_type,'L') /= 0)then
     fmin = minval(log(column(1:n)%row(iz)))
     fmax = maxval(log(column(1:n)%row(iz)))
   endif
   print *, ' fmin, fmax ', fmin, fmax
   fmin_dp = fmin
   fmax_dp = fmax
   call query_real(' fmin ? ', fmin_dp, 'f10.3')
   call query_real(' fmax ? ', fmax_dp, 'f10.3')
   fmin=fmin_dp
   fmax=fmax_dp

   do j=1,nx
    do k=1,ny
      z(j,k)=0.
    end do
   end do
   do j=1,n
     ix_j = (column(j)%row(ix)-xmin)/delta_x  + .5 + 1
     iy_j = (column(j)%row(iy)-ymin)/delta_y  + .5 + 1
     z(ix_j,iy_j) = column(j)%row(iz)
     if(index(data_type,'L') /= 0)z(ix_j,iy_j) = log(column(j)%row(iz))
     if(column(j)%row(turns) > 0 .and. column(j)%row(turns) < max_turns) &
       z(ix_j,iy_j) = fmax   
   if( ix_j >nx .or. iy_j > ny)then
      print '(6f10.5)', column(j)%row(ix), xmin,delta_x, column(j)%row(iy), ymin, delta_y
      print *,' ix_j, iy_j, z ', ix_j,iy_j,z(ix_j,iy_j)
    endif 

   end do 

 print '(a,$)',' Fill in the gaps (Y/N) ?'
 read(*,'(a)')ans
 call str_upcase(ans,ans) 
 if(index(ans,'Y') /= 0)then 
! fill in gaps
   do j=1,nx
    do k=1,ny
     if(z(j,k) == 0.)then
       do l = j,1,-1
         if(z(l,k) /= 0.)then
            j_left = l
            exit
         endif
         j_left = 1
       end do

       do l=j,nx
         if(z(l,k) /= 0.)then
            j_right = l
            exit
         endif
         j_right = nx
       end do

       if(j_right == nx .and. j_left == 1)then  !look in y
        do l = k,1,-1
          if(z(j,l) /= 0.)then
            k_down = l
            exit
          endif
          k_down = 1
        end do

        do l=k,ny
          if(z(j,l) /= 0.)then
            k_up = l
            exit
          endif
          k_up = ny
        end do

        z(j,k) = (z(j,k_up) - z(j,k_down))/(k_up-k_down)*(k-k_down) +z(j,k_down)

       else
        z(j,k) = (z(j_right,k) - z(j_left,k))/(j_right-j_left)*(j-j_left) +z(j_left,k)
       endif
      endif
     end do
    end do
 endif
  ans='N'
  print *, ' fmin, fmax ', fmin, fmax
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
      CALL PGQCIR(C1, C2)
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
      TR(1) = xmin
      TR(2) = delta_x
      TR(3) = 0.
      TR(4) = ymin
      TR(5) = 0.0
      TR(6) = delta_y
!
! Clear the screen. Set up window and viewport.
!
      CALL PGPAGE


       call pgscr(0, 1., 1., 1.)
       call pgscr(1,0.,0.,0.)
       call pgscr(2, 1., 0., 0.)
       call pgscr(3,0.,0.,0.,1.0)
       call pgsch(2.)

      CALL SETVP
!      CALL PGWNAD(0.0, 1.0+nx, 0.0, 1.0+ny)
      CALL PGWNAD(xmin+delta_x/2,xmax+delta_x/2,ymin+delta_y/2,ymax+delta_x/2)



!
! Annotate the plot.
!
      print '(a,$)',' Scan input file name ?'
      read(*,'(a)')file_name
      print '(2a)', ' file_name = ', file_name
      if(device_type == '/XSERVE')then
        CALL PGSCI(0)
       else
        CALL PGSCI(1)
      endif
      ix = index(title,'Vertical')
      title = title(1:ix+7)//'   '//file_name(1:len(file_name)) 
      print *, title
      CALL PGSCH(1.5)
!      CALL PGMTXT('t',1.0,0.5,0.0,file_name)
      CALL PGMTXT('t', 1.0, 0.0, 0.0, title)
      CALL PGSCH(1.)
      CALL PGBOX('bcntsia1',0.0,0,'bcntsiva1',0.0,0)
      if(ix == 1)CALL PGMTXT('b',3.0,1.0,1.0,'horizontal tune (single beam)')
      if(ix == 3)CALL PGMTXT('b',3.0,1.0,1.0,'horizontal tune including lrbbi')
      if(iy == 2)CALL PGMTXT('l',3.5,1.0,1.0,'vertical tune (single beam)')
      if(iy == 4)CALL PGMTXT('l',3.5,1.0,1.0,'vertical tune including lrbbi')


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
      CALL PGIMAG(z,nx,ny,1,nx,1,ny,FMIN,FMAX,TR)

! Draw a wedge.
!

      CALL PGWEDG('BI', 4.0, 5.0, FMIN, FMAX, wedge_label)
      CALL PGSCH(1.0)
!      CALL PGIDEN    
      call PGUPDT

    
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


   type(column_struct), allocatable :: column(:)
   character*80 file_name
   character*140 line
   character*20 word
   integer i, n, ix, j, size
   integer n_columns/8/, n_skip/9/
   integer io

   open(unit=2, file=file_name, status='OLD')
   n=0
   i=0
   do while(.true.)
    read(2,'(a140)', end=99)line
    i=i+1
    if(i > n_skip)n=n+1
   end do
99 continue
   close(unit=2)
   size = n 
   allocate(column(n))

   open(unit=2, file=file_name, status='OLD')
   n=0
   i=0
   do while(.true.)   
    read(2,'(a140)', end=100)line
    i=i+1
    if(i > n_skip .and. line(1:15) /= '               ' .and. &
         n_columns >=1 )then
     n=n+1
     ix=0
      do j=1,n_columns
        if(j == n_columns .and. ix == 0)then
          word = '-1'
         else
          call string_trim(line(ix+1:), line, ix)
          word=line(1:ix)
        endif
        if(index(word,'NaN') /= 0)then
          column(n)%row(j) = 0.
         else
          read(word,*, iostat = io)column(n)%row(j)
       endif
      end do
    endif
   end do
100 continue
  size = n
!  print *,' at return from read scan, n=',n
  close(unit=2)
  return
 end
   


