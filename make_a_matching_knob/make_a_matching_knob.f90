program make_a_matching_knob

  use bmad

  implicit none

  type(lat_struct) ring, ring0
  type(coord_struct), allocatable :: co(:)
  character*10 mags(20)
  type (ele_pointer_struct), allocatable :: eles(:)

  integer, parameter :: max_it = 200
  integer i,j,k,l
  integer n_vars, n_cons
  integer n_loc
  integer nQx, nQy
  integer match_point
  integer xix, yix
  integer status
  logical converged
  logical set_x, set_y
  character*100 in_file
  character*100 lat_file
  character*20 use_line
  character*16 var_str
  character*25 set_str
  real(rp), allocatable :: jac_delta(:)
  real(rp), allocatable :: var(:), var0(:), delta_var(:)
  real(rp), allocatable :: con0(:), con(:)
  real(rp), allocatable :: Jac(:,:)
  real(rp), allocatable :: Jacp(:,:)
  real(rp) dk
  real(rp) eps
  real(rp) alpha_n
  real(rp) Qxmin, Qxmax, Qymin, Qymax, dQx, dQy

  namelist /knob/ lat_file,use_line,match_point,eps,alpha_n,dk,set_x,set_y,mags,nQx,nQy,Qxmin,Qxmax,Qymin,Qymax 

  mags = ''

  call getarg(1,in_file)
  open (unit = 10, file = in_file, action='read')
  read (10, nml = knob)
  close (10)

  i = 0
  do while(.true.)
    if(mags(i+1) .ne. '') then
      i = i + 1
    else
      exit
    endif
  enddo
  n_vars = i

  n_cons = 5
  if(set_x .and. .not. set_y) then
    n_cons = 6
    xix = 6
    nQy = 1
  elseif(set_y .and. .not. set_x) then
    n_cons = 6
    yix = 6
    nQx = 1
  elseif(set_x .and. set_y) then
    n_cons = 7
    xix = 6
    yix = 7
  else
    write(*,*) "Neither set_x nor set_y true.  Aborting."
    stop
  endif

  call bmad_parser(lat_file, ring, use_line=use_line)  !parameter
  call set_on_off(rfcavity$, ring, off$)
  bmad_com%radiation_damping_on = .false.
  bmad_com%radiation_fluctuations_on = .false.
  call twiss_and_track(ring,co,status)

  allocate(jac_delta(n_vars))
  allocate(var0(n_vars))
  allocate(var(n_vars))
  allocate(delta_var(n_vars))
  allocate(con0(n_cons))
  allocate(con(n_cons))
  allocate(Jac(n_cons,n_vars))
  allocate(Jacp(n_vars,n_cons))
  var0 = 0.0d0
  var = 0.0d0

  ! Constraints
  con = 0.0d0
  con0 = 0.0d0
  con(1) = ring%ele(match_point)%a%beta
  con(2) = ring%ele(match_point)%a%alpha
  con(3) = ring%ele(match_point)%b%beta
  con(4) = ring%ele(match_point)%b%alpha
  con(5) = 0.0d0  !horizontal dispersion at element 0

  call get_mag_str(mags,ring,var0)

  write(*,*) "Starting simulation..."

  ! Loop to find solution var
  ring0 = ring
  open(11,file='knob.grid')
  write(11,'(10a13)') '# dQx', 'dQy', (trim(mags(i)),i=1,n_vars)
  do i = 1,nQx
    do j = 1,nQy
      converged = .false.
      ring = ring0
      call get_mag_str(mags,ring,var)
      ! Desired tune adjustment
      if(set_x) then
        dQx = (Qxmax-Qxmin)/(nQx-1.0)*(i-1.0) + Qxmin
        con(xix) = ring%ele(ring%n_ele_track)%a%phi/twopi + dQx
      endif
      if(set_y) then
        dQy = (Qymax-Qymin)/(nQy-1.0)*(j-1.0) + Qymin
        con(yix) = ring%ele(ring%n_ele_track)%b%phi/twopi + dQy
      endif

      ! Intermediate values of constraints.
      con0(1) = ring%ele(match_point)%a%beta
      con0(2) = ring%ele(match_point)%a%alpha
      con0(3) = ring%ele(match_point)%b%beta
      con0(4) = ring%ele(match_point)%b%alpha
      con0(5) = ring%ele(0)%a%eta
      if(set_x) con0(xix) = ring%ele(ring%n_ele_track)%a%phi/twopi
      if(set_y) con0(yix) = ring%ele(ring%n_ele_track)%b%phi/twopi
      do k=1,max_it
        !make jacobian
        do l=1,n_vars  !calculate response for family l
          jac_delta(:) = 0.0d0
          jac_delta(l) = dk
          call set_mag_str(mags,ring,var+jac_delta)
          call twiss_and_track(ring,co,status)
          if(status /= ok$) then
            write(*,'(a,a)') "Could not calculate ring during Jacobian calculation for ", mags(l)
            stop
          endif

          Jac(1,l) = (con0(1) - ring%ele(match_point)%a%beta) / dk
          Jac(2,l) = (con0(2) - ring%ele(match_point)%a%alpha) / dk
          Jac(3,l) = (con0(3) - ring%ele(match_point)%b%beta) / dk
          Jac(4,l) = (con0(4) - ring%ele(match_point)%b%alpha) / dk
          Jac(5,l) = (con0(5) - ring%ele(0)%a%eta) / dk
          if(set_x) Jac(xix,l) = (con0(xix) - ring%ele(ring%n_ele_track)%a%phi/twopi) / dk
          if(set_y) Jac(yix,l) = (con0(yix) - ring%ele(ring%n_ele_track)%b%phi/twopi) / dk
        enddo

        call mat_pseudoinverse(Jac,Jacp)

        delta_var = matmul(Jacp,(con-con0))
        var = var - alpha_n*delta_var
        call set_mag_str(mags,ring,var)
        call twiss_and_track(ring,co,status)
        if(status /= ok$) then
          write(*,'(a)') "Step failed."
          exit
        endif

        con0(1) = ring%ele(match_point)%a%beta
        con0(2) = ring%ele(match_point)%a%alpha
        con0(3) = ring%ele(match_point)%b%beta
        con0(4) = ring%ele(match_point)%b%alpha
        con0(5) = ring%ele(0)%a%eta
        if(set_x) con0(xix) = ring%ele(ring%n_ele_track)%a%phi/twopi
        if(set_y) con0(yix) = ring%ele(ring%n_ele_track)%b%phi/twopi

        if( maxval(abs(con-con0)) .lt. eps ) then
          converged = .true.
          exit
        endif
      enddo
      if(converged) then
        write(11,'(2f13.5,10f13.6)') dQx, dQy, var-var0
      endif

      write(*,'(2i3,a,l1,a,i5,a)') i, j, " converged? ", converged, " in ", k, " iterations."
    enddo
  enddo
  close(11)

  contains

  subroutine get_mag_str(mags,ring,strengths)
    use bmad

    implicit none

    character(*) mags(:)
    type(lat_struct) ring
    real(rp) strengths(:)

    integer n_mags, n_loc
    type (ele_pointer_struct), allocatable :: eles(:)
    logical err
    integer i

    n_mags = size(mags)

    do i=1, n_mags
      if(mags(i) == '') exit
      call lat_ele_locator(mags(i), ring, eles, n_loc, err)
      if(err .or. (n_loc .lt. 1)) then
        write(*,*) "Get ele attribute error.  Terminating."
        call err_exit
      endif
      strengths(i) = eles(1)%ele%value(k1$)
    enddo

    if(allocated(eles)) deallocate(eles)
  end subroutine

  subroutine set_mag_str(mags,ring,strengths)
    use bmad

    implicit none

    character(*) mags(:)
    type(lat_struct) ring
    real(rp) strengths(:)

    integer n_mags, n_loc
    type (ele_pointer_struct), allocatable :: eles(:)
    logical err
    character*18 var_str
    character*30 set_str
    integer i,j

    n_mags = size(mags)

    do i=1, n_mags
      if(mags(i) == '') exit
      call lat_ele_locator(mags(i), ring, eles, n_loc, err)
      write(var_str,'(f18.8)') strengths(i)
      set_str = 'k1='//trim(adjustl(var_str))

      do j=1, n_loc
        call set_ele_attribute (eles(j)%ele, set_str, err)
        if(err) then
          write(*,*) "Set ele attribute error.  Terminating.", set_str
          !error stop
          call err_exit
        endif
      enddo
    enddo

    if(allocated(eles)) deallocate(eles)
    call lattice_bookkeeper(ring)
  end subroutine

end program









