!+
! Subroutine phase_space_fit (x, xp, twiss, tune, emitt, x_0, xp_0, chi, tol)
!
! Subroutine to fit one turn data (x, x') and extract the twiss parameters,
! etc. 
!
! modules needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!   x(:)  -- Real: Array of x data.
!   xp(:) -- Real: Array of x' data.
!   tol   -- Real: Optional. Tolerance for the fit. If chi < tol 
!               then the subroutine declares success and exits 
!               from the field of battle. Default: tol = 1e-3
!
! Output:
!   twiss -- Twiss_struct: 
!     %phi   -- phase of the oscillations. 
!               x = sqrt(emitt*beta) * sin(tune*n + phi)
!               n = turn #, first data point is n = 1.
!     %beta  -- Beta parameter.
!     %alpha -- alpha parameter.
!   tune  -- Fractional tune 0 < tune < twopi.
!   emitt -- Real: Sqrt(emitt*beta) is the amplitude.
!   x_0   -- Real: x Closed orbit offset.
!   xp_0  -- Real: x' Closed orbit offset.
!   chi   -- Real: Figure of merit for the fit. 
!                 0 => Perfect.
!                 1 => Terrible.
!-

subroutine phase_space_fit (x, xp, twiss, tune, emitt, x_0, xp_0, chi, tol)

  use bmad_struct
  use bmad_interface
  use nr, only: mrqmin
  implicit none

  type (twiss_struct) twiss

  real, optional :: tol
  real x(:), xp(:), tune, emitt, target
  real A_x, A_xp, phi_x, phi_xp, x_0, xp_0
  real chisq, alamda, dphi, sin_del, chi, x1, x2

  real, allocatable :: a(:), x_in(:), y_in(:), sig(:)
  real, allocatable :: covar(:,:), alpha(:,:)

  integer i, n_cross, n_pts, n_in, n1, n2

  logical, allocatable :: maska(:)
  logical init_needed

! interface

  interface
    subroutine phase_fit_func (x, a, yfit, dyda)
      use nrtype
      real(sp), intent(in) :: x(:), a(:)
      real(sp), intent(out) :: yfit(:)
      real(sp), intent(out) :: dyda(:,:)
    end subroutine
  end interface

! init
! we use a "Hann" window function (see the Num. Rec. Book) so that it is less
! likely that the fitting algorithm will not get stuck in a local minimum

  if (.not. allocated(a)) then
    allocate (a(7), covar(7,7), alpha(7,7), maska(7))
    maska = .true.
  endif

  n_pts = size(x)         
  n_in = 2 * n_pts

  init_needed = .false.

  if (allocated(x_in)) then
    if (size(x_in) /= n_in) then
      deallocate (x_in, y_in, sig, maska, covar, alpha)
      init_needed = .true.
    endif
  else
    init_needed = .true.
  endif

  if (init_needed) then
    allocate (x_in(n_in), y_in(n_in), sig(n_in))
    do i = 1, n_pts
      sig(i) = 2 / (1 - cos(2*pi*i/(n_pts+1.0))) 
      sig(i+n_pts) = sig(i)
    enddo
  endif

! error check

  if (n_pts /= size(xp)) then
    print *, 'ERROR IN PHASE_SPACE_FIT: X AND XP ARRAY SIZES DO NOT MATCH.'
    print *, '      SIZES:', size(x), size(xp)
    call err_exit
  endif

! we fit to the form:
!   x  = A_x  * sin(tune*n + phi_x)  + x_0
!   xp = A_xp * sin(tune*n + phi_xp) + xp_0

! as a starting point make a rough guess what the parameters are

  A_x = (maxval(x) - minval(x)) / 2
  A_xp = (maxval(xp) - minval(xp)) / 2
  x_0 = sum(x) / n_pts
  xp_0 = sum(xp) / n_pts
                                        
  n_cross = -1
  do i = 1, n_pts-1
    if (x(i) > x_0 .and. x(i+1) <= x_0 .or. &
                        x(i) <= x_0 .and. x(i+1) > x_0) then 
      if (n_cross == -1) n1 = i
      n2 = i
      n_cross = n_cross + 1
    endif
  enddo

  x1 = (n1 * x(n1+1) - (n1+1) * x(n1)) / (x(n1+1) - x(n1))
  x2 = (n2 * x(n2+1) - (n2+1) * x(n2)) / (x(n2+1) - x(n2))
  tune = n_cross * pi / (x2 - x1)
 
  if (x(n1) <= 0) then
    phi_x = -tune * x1 
  else
    phi_x = -tune * x1 + pi
  endif

  do i = 1, n_pts-1
    if (xp(i) > xp_0 .and. xp(i+1) <= xp_0) then
      x1 = (i * x(i+1) - (i+1) * x(i)) / (x(i+1) - x(i))
      phi_xp = -tune * x1 + pi
    elseif (xp(i) <= xp_0 .and. xp(i+1) > xp_0) then
      x1 = (i * x(i+1) - (i+1) * x(i)) / (x(i+1) - x(i))
      phi_xp = -tune * x1 
    endif
  enddo

! now use the Numerical Recipes routine mrqmin (Levenberg-Marquardt) to
! find the best fit to the data

  a(1) = A_x
  a(2) = A_xp     
  a(3) = phi_x
  a(4) = phi_xp
  a(5) = x_0
  a(6) = xp_0
  a(7) = tune

  x_in = 0
  y_in(1:n_pts) = x 
  y_in(n_pts+1:n_in) = xp 

  alamda = -1
  if (present(tol)) then
    target =  tol**2 * (A_x**2 + A_xp**2) * n_pts
  else
    target =  1e-6 * (A_x**2 + A_xp**2) * n_pts
  endif

  do i = 1, 50
    call mrqmin (x_in, y_in, sig, a, maska, covar, alpha, chisq, &
                                                 phase_fit_func, alamda)
    if (chisq < target) exit
    if (alamda > 1e20 .or. alamda < 1e-20) exit
  enddo

  alamda = 0
  call mrqmin (x_in, y_in, sig, a, maska, covar, alpha, chisq, &
                                                   phase_fit_func, alamda)

! calc chi

  A_x    = a(1)
  A_xp   = a(2)
  phi_x  = a(3)
  phi_xp = a(4)
  x_0    = a(5)
  xp_0   = a(6)
  tune   = a(7)
             
  chi = 0            
  do i = 1, n_pts
    chi = chi + (x(i) - A_x*sin(i*tune + phi_x) - x_0)**2/A_x**2 + &
                (xp(i) - A_xp*sin(i*tune + phi_xp) - xp_0)**2/A_xp**2 
  enddo
  chi = sqrt(chi / (2 * n_pts))

! Convert to Twiss parameters.

  dphi = phi_xp - phi_x
  sin_del = sin(dphi)
                               
  emitt = A_x * A_xp * sin_del
  twiss%beta = A_x / (A_xp * sin_del)
  twiss%phi = phi_x           
  twiss%alpha = - cos(dphi) / dphi

  if (emitt < 0) then
    emitt = -emitt
    twiss%beta = -twiss%beta
    twiss%phi = -twiss%phi
    twiss%alpha = -twiss%alpha
    tune = twopi - tune
  endif

  tune = modulo(tune, twopi)

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine phase_fit_func (x, a, yfit, dyda)
!
! Subroutine for mrqmin to evaluate the fitting function. 
! See the Num Rec book for more details.
!-

subroutine phase_fit_func (x, a, yfit, dyda)

  use nrtype

  implicit none

  real(sp), intent(in) :: x(:), a(:)
  real(sp), intent(out) :: yfit(:)
  real(sp), intent(out) :: dyda(:, :)

  real A_x, phi_x, x_0, A_xp, phi_xp, xp_0, tune

  real sin_x, cos_x, sin_xp, cos_xp

  integer i, k, n_pts

  logical init_needed

! init
                  
  A_x    = a(1)
  A_xp   = a(2)
  phi_x  = a(3)
  phi_xp = a(4)
  x_0    = a(5)
  xp_0   = a(6)
  tune   = a(7)

  n_pts = size(x)/2

! calc yfit and dyda

  do i = 1, n_pts

    k = i + n_pts

    sin_x  = sin(tune*i + phi_x)
    cos_x  = cos(tune*i + phi_x)
    sin_xp = sin(tune*i + phi_xp)
    cos_xp = cos(tune*i + phi_xp)

    yfit(i) = (A_x * sin_x + x_0) 
    yfit(k) = (A_xp * sin_xp + xp_0) 

    dyda(i, 1) = sin_x 
    dyda(i, 2) = 0
    dyda(i, 3) = A_x * cos_x 
    dyda(i, 4) = 0
    dyda(i, 5) = 1
    dyda(i, 6) = 0
    dyda(i, 7) = i * A_x * cos_x 

    dyda(k, 1) = 0
    dyda(k, 2) = sin_xp 
    dyda(k, 3) = 0
    dyda(k, 4) = A_xp * cos_xp 
    dyda(k, 5) = 0
    dyda(k, 6) = 1
    dyda(k, 7) = i * A_xp * cos_xp 

  enddo

end subroutine
                              
