
subroutine plot_pretz(lat,ncross, cross)

 use bmad

  implicit none

 type(lat_struct) lat
 type(coord_struct), allocatable :: co_positron(:), co_electron(:)
 type(ele_struct) ele

  integer trains, bunches, spacing
  integer ncross
  integer n_turns, i, pgopen, istat2
  integer i_dim/4/
  integer j,k, l

  real(rp)  cross(1000)
  real(rp) scale/1000./
  real(rp) theta
  real(rp) r,xp,yp,xe,ye, x, y
  real(rp) xlast, ylast
  real(rp) slast,s, z, xc, yc,xlow,xhigh,ylow,yhigh,bar/1.5/,mag/3.0/
  real(rp) magg/1.5/
  real(rp) spast, xpast, ypast
  real(rp) inj/0.02/
  real(rp) zero/0.0/
  real(rp) xx(4), yy(4)
  real(rp) sfeed(3)/55.091335,56.443835,700./,xfeed(3),yfeed(3)

  character*5 sources(10)/'B04W ','B05W ','B06W ','B07W ','B04E ','B05E ','B06E ','B07E ','WIG_W','WIG_E'/

  call reallocate_coord (co_positron, lat%n_ele_max)
  call reallocate_coord (co_electron, lat%n_ele_max)
  call closed_orbit_calc(lat, co_positron, i_dim)
  lat%param%particle = electron$
  call closed_orbit_calc(lat, co_electron, i_dim)

! IP_L0
     write(37,'(2f10.3)')zero, zero
     write(37,*)

   do i = 2, lat%n_ele_track

   xlast = lat%ele(i-1)%floor%z
   ylast = lat%ele(i-1)%floor%x

   ele = lat%ele(i)
   theta =  ele%floor%theta

   
   x = ele%floor%z
   y = ele%floor%x
   s = lat%ele(i)%s
   if(s >0. .and. x==0. .and. y == 0.)cycle
   j=1
   do while((x-xlast == 0. .and. y-ylast == 0.) .or. (slast >0. .and. xlast == 0. .and. ylast == 0.))
    j=j+1
    ylast = lat%ele(i-j)%floor%x
    xlast = lat%ele(i-j)%floor%z
    slast = lat%ele(i-j)%s
   end do

   theta = atan2(ylast - y, xlast-x)

   r = co_positron(i)%vec(1)
   xp = r * sin(-theta)*scale
   yp = r * cos(theta)*scale

   r = co_electron(i)%vec(1)
   xe = r * sin(-theta)*scale
   ye = r * cos(theta)*scale

   if(ele%s /= 0. .and. (x == 0. .and. y == 0.))cycle
   write(35,'(9f10.3)')ele%s,-x,-y,-xp,-yp, -xe,-ye,r, theta
   if(lat%ele(i-1)%key == patch$ .or. lat%ele(i)%key == patch$)cycle
   if(ele%name(1:4) == 'Q34W') then
     write(37,'(2f10.3)')-(x+inj*sin(-theta)*scale), -(y+inj*cos(theta)*scale)
     write(37,*)
   endif
   if(ele%name(1:4) == 'Q34E') then
     write(37,'(2f10.3)')-(x+inj*sin(-theta)*scale),-( y+inj*cos(theta)*scale)
     write(37,*)
   endif
   if(ele%name(1:5) == 'H_SEP')then
     write(40,'(2f10.3)')-( x+sin(-theta)*mag),-(y+cos(theta)*mag)
     write(40,'(2f10.3)')-( x-sin(-theta)*mag),-(y-cos(theta)*mag)
     write(40,*)
   endif

   do k=1,ncross
    spast = lat%ele(i-1)%s
    if(s-spast == 0.)cycle
    if(cross(k) >= spast .and. cross(k) < s)then
    ypast = lat%ele(i-1)%floor%x
    xpast = lat%ele(i-1)%floor%z
     
      z = (cross(k)-spast)/(s-spast)
      xc = (x-xpast)* z + xpast
      yc = (y-ypast) * z + ypast
      xlow = xc - bar *( sin(-theta))
      xhigh = xc + bar * ( sin(-theta))
      ylow = yc - bar * ( cos(theta))
      yhigh = yc + bar * ( cos(theta))
     write(36,'(2f10.3,i,6f10.3)') -xlow, -ylow, k, z, theta, x, xpast, y, ypast
     write(36,'(2f10.3,i,2f10.3)') -xhigh, -yhigh, k, z, theta
     write(36,*) 
     write(36,*) 
     write(36,*) 
   endif
   end do

   do k=1,2
    spast = lat%ele(i-1)%s
    if(s-spast == 0.)cycle
    if(sfeed(k) >= spast .and. sfeed(k) < s)then
      ypast = lat%ele(i-1)%floor%x
      xpast = lat%ele(i-1)%floor%z

      z=(sfeed(k)-spast)/(s-spast)
      xfeed(k) = (x-xpast) * z + xpast
      yfeed(k) = (y-ypast) * z + ypast
      write(39,'(2f10.3,/)') -xfeed(k), -yfeed(k)
    endif
   end do

    do k=1,10
     if(ele%name(1:5) == sources(k)(1:5))then
       xpast = lat%ele(i-1)%floor%z
       ypast = lat%ele(i-1)%floor%x
       xx(1) = xpast -mag*sin(-theta)
       yy(1) = ypast -mag*cos(theta)

       xx(2) = xpast +mag*sin(-theta)
       yy(2) = ypast +mag*cos(theta)

       xx(4) = x -mag*sin(-theta)
       yy(4) = y -mag*cos(theta)

       xx(3) = x +mag*sin(-theta)
       yy(3) = y +mag*cos(theta)

      do l=1,4
        write(38,'(2f10.3)')-xx(l),-yy(l)
      end do
        write(38,'(2f10.3)')-xx(1),-yy(1)
      write(38,*)
     endif
    end do 

     if(ele%key == quadrupole$ .or. ele%key == sbend$ .or.ele%key == rbend$ &
                             .or. ele%key == wiggler$)then
       xpast = lat%ele(i-1)%floor%z
       ypast = lat%ele(i-1)%floor%x
       xx(1) = xpast -magg*sin(-theta)
       yy(1) = ypast -magg*cos(theta)

       xx(2) = xpast +magg*sin(-theta)
       yy(2) = ypast +magg*cos(theta)

       xx(4) = x -magg*sin(-theta)
       yy(4) = y -magg*cos(theta)

       xx(3) = x +magg*sin(-theta)
       yy(3) = y +magg*cos(theta)

      do l=1,4
        write(41,'(2f10.3)')-xx(l),-yy(l)
      end do
        write(41,'(2f10.3)')-xx(1),-yy(1)
      write(41,*)
     endif
  end do

  return
  end



  






