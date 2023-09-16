module adts_mod

contains 

subroutine tune_ele_by_ele(ring,orb0,n_turns,track_state,nu_x,nu_y,std_nu_x,std_nu_y,dump_coords)
  use bmad
  use four_mod
  !use sls_lib

  implicit none

  type(lat_struct) ring
  type(coord_struct) :: orb0
  integer n_turns
  integer track_state
  real(rp) nu_x, nu_y
  real(rp), optional :: std_nu_x, std_nu_y
  integer, optional :: dump_coords

  type(coord_struct), allocatable :: orb(:)
  real(rp) d_nu_x, d_nu_y
  integer j
  real(rp), allocatable :: nu_x_stash(:), nu_y_stash(:) 

  allocate(orb(0:ring%n_ele_track))
  allocate(nu_x_stash(1:n_turns))
  allocate(nu_y_stash(1:n_turns))

  orb(0) = orb0
  do j=1,n_turns
    call track_all(ring,orb,track_state=track_state)
    if(track_state == moving_forward$) then 
      call accumulate_phase_advance(ring, orb, d_nu_x, d_nu_y, dump_coords)
      nu_x_stash(j) = d_nu_x
      nu_y_stash(j) = d_nu_y
      orb(0)=orb(ring%n_ele_track)
    else
      ! write(*,'(a,i4,a,6es11.2)') "ADTS: Particle lost at turn ", j,  ". Initial coordiantes: ", orb0%vec
      ! write(*,*) "ADTS: Exiting and returning nu_x = nu_y = -1"
      exit
    endif
  enddo

  if(track_state == moving_forward$) then
    nu_x = sum(nu_x_stash(:))/twopi/n_turns
    nu_y = sum(nu_y_stash(:))/twopi/n_turns
    if(present(std_nu_x)) std_nu_x = sqrt(sum((nu_x_stash/twopi-nu_x)**2)/(n_turns-1))
    if(present(std_nu_y)) std_nu_y = sqrt(sum((nu_y_stash/twopi-nu_y)**2)/(n_turns-1))
  else
    nu_x = -1.0
    nu_y = -1.0
    if(present(std_nu_x)) std_nu_x = -1.0
    if(present(std_nu_y)) std_nu_y = -1.0
  endif

  deallocate(nu_x_stash)
  deallocate(nu_y_stash)
  deallocate(orb)
end subroutine

subroutine accumulate_phase_advance(ring, orb, d_nu_x, d_nu_y, dump_coords)
  use bmad
  use four_mod
  !use sls_lib

  implicit none

  type(lat_struct) ring
  type(coord_struct) orb(0:)
  real(rp) d_nu_x, d_nu_y
  integer, optional :: dump_coords

  type(coord_struct) norm_coords
  integer i
  real(rp) x1prev, xp1prev, y1prev, yp1prev
  real(rp) x1, xp1, y1, yp1
  real(rp) arg_x, arg_y, angle_x, angle_y
  logical ok

  d_nu_x = 0.0d0
  d_nu_y = 0.0d0

  do i=0,ring%n_ele_track
    !stash normal coordiantes
    if(i .gt. 0) then
      x1prev = x1
      xp1prev = xp1
      y1prev = y1
      yp1prev = yp1
    endif
    call xy_to_action(ring, i, orb(i)%vec, norm_coords%vec, ok)
    if(.not. ok) then
      write(*,*) "xy_to_action is not ok."
      error stop
    endif
    x1 = norm_coords%vec(1)
    xp1 = norm_coords%vec(2)
    y1 = norm_coords%vec(3)
    yp1 = norm_coords%vec(4)

    if( is_linear_ele(ring%ele(i)) ) then
      if( i .gt. 0 ) then
        d_nu_x = d_nu_x + ( ring%ele(i)%a%phi - ring%ele(i-1)%a%phi )
        d_nu_y = d_nu_y + ( ring%ele(i)%b%phi - ring%ele(i-1)%b%phi )
      endif
    else
      if(i .gt. 0) then
        arg_x = (x1*x1prev + xp1*xp1prev)/sqrt(x1*x1 + xp1*xp1)/sqrt(x1prev*x1prev + xp1prev*xp1prev)
        if( abs(arg_x) .lt. 1.0d0 ) then  !protect against rounding errors sending arg_x above 1.0
          angle_x = acos( arg_x )
        else
          angle_x = 0.0d0
        endif
        if( x1*xp1prev - xp1*x1prev < 0 ) angle_x = -angle_x
        d_nu_x = d_nu_x + angle_x

        arg_y = (y1*y1prev + yp1*yp1prev)/sqrt(y1*y1+yp1*yp1)/sqrt(y1prev*y1prev+yp1prev*yp1prev)
        if( abs(arg_y) .lt. 1.0d0 ) then  !protect against rounding errors sending arg_y above 1.0
          angle_y = acos( arg_y )
        else
          angle_y = 0.0d0
        endif
        if( y1*yp1prev - yp1*y1prev < 0 ) angle_y = -angle_y
        d_nu_y = d_nu_y + angle_y
      endif
    endif
    if( present(dump_coords) ) write(dump_coords-1,'(i6,18es17.9)') i, orb(i)%vec
    if( present(dump_coords) ) write(dump_coords,  '(i6,18es17.9)') i, norm_coords%vec, angle_x, d_nu_x, arg_x, angle_y, d_nu_y, arg_y
  enddo

  contains

    function is_linear_ele(ele)
    use bmad
    implicit none
    logical is_linear_ele
    type(ele_struct) ele

    if( (ele%key == sbend$) .or. &
        (ele%key == sextupole$) .or. &
        (ele%key == multipole$) ) then
      is_linear_ele = .false.
    else
      is_linear_ele = .true.
    endif
    end function

end subroutine

end module
