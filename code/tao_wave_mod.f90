module tao_wave_mod

use tao_utils

real(rp), private, allocatable, save :: f_com(:,:)
type (tao_d1_data_struct), private, pointer, save :: d1_ptr

contains

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine tao_wave_cmd (curve_name, plot_place)
!
! Routine to setup wave plotting.
!
! Input:
!   curve_name -- Character(*) curve for wave analysis.
!   place      -- Character(*) place on plot page to put the plot.
!-

subroutine tao_wave_cmd (curve_name, plot_place)

use quick_plot

implicit none

type (tao_curve_array_struct), allocatable, target, save :: curve_array(:)
type (tao_curve_struct), pointer :: curve
type (tao_plot_region_struct), pointer :: region
type (tao_plot_struct), save, target :: wave_plot
type (tao_plot_struct), pointer :: plot
type (tao_graph_struct), pointer :: wg, wg1
type (tao_curve_struct), pointer :: wc, wc1
type (tao_d2_data_struct), pointer :: d2_ptr
type (tao_d1_data_array_struct), allocatable, save :: d1_array(:)

integer i, p1, p2

character(*) curve_name, plot_place
character(20) :: r_name = 'tao_wave_cmd'

logical :: init_needed = .true.
logical err

! Find the curve

call tao_find_plots (err, curve_name, 'REGION', curve = curve_array, always_allocate = .true.)
if (err) return

if (size(curve_array) == 0) then
  call out_io (s_error$, r_name, 'No curve found.')
  return
endif

if (size(curve_array) > 1) then
  call out_io (s_error$, r_name, 'Multiple curves match the given curve name. Nothing done.')
  return
endif

curve => curve_array(1)%c
plot => curve%g%p

! Find the plot region

if (plot_place /= '') then
  call tao_find_plot_region (err, plot_place, region)
  if (err) return
else
  region => plot%r
endif

! Initiate the wave plot.
! There are four graphs. Each graph has a single curve.
! The first graph holds the original curve
! The other three graphs are: 
!   2) extended original curve 
!   3) curve - a_region_fit
!   4) curve - b_region_fit

if (init_needed) then
  allocate(wave_plot%graph(4))
  init_needed = .false.
endif

! Transfer curve information to the wave plot

wave_plot%name        = 'wave'
wave_plot%x_axis_type = plot%x_axis_type
wave_plot%x           = plot%x
wave_plot%autoscale_gang_x = .false.

do i = 1, 4
  wg => wave_plot%graph(i)
  wg = curve%g
  wg%x = plot%x
  wg%box = [1, 5-i, 1, 3]

  if (size(wg%curve) /= 1) then
    deallocate (wg%curve)
    allocate (wg%curve(1))
    wg%curve(1) = curve
  endif
  wg%curve(1)%g => wg
enddo

wave_plot%graph(1)%visible = .false.
wave_plot%graph(2)%type = 'wave.0'  ! Extended origninal curve.
wave_plot%graph(3)%type = 'wave.a'  ! Curve - A region fit.
wave_plot%graph(4)%type = 'wave.b'  ! Curve - B region fit.

wave_plot%graph(2)%name = '0'  ! Extended origninal curve.
wave_plot%graph(3)%name = 'a'  ! Curve - A region fit.
wave_plot%graph(4)%name = 'b'  ! Curve - B region fit.

! Set up curves to wrap around the lattice, etc.

wg1 => wave_plot%graph(1)
wc1 => wg1%curve(1)

wc1%y_axis_scale_factor = 1

if (wc1%data_source == 'data_array') then
  call tao_find_data (err, wc1%data_type, d2_ptr, d1_array, ix_uni = wc1%ix_universe)
  d1_ptr => d1_array(1)%d1
else
  call out_io (s_error$, r_name, &
                    'ANALYSIS FOR THIS DATA_SOURCE NOT YET IMPLEMENTED: ' // wc1%data_source)
  return
endif

do i = 2, 4
  wg => wave_plot%graph(i)
  ! wg%x%major_div_nominal = nint(1.5wg1%x%major_div_nominal)
  p1 = nint(0.7 * wg%x%major_div_nominal)
  p2 = nint(1.3 * wg%x%major_div_nominal)
  call qp_calc_and_set_axis ('X', 0.0_rp, 1.5*wg1%x%max, p1, p2, 'GENERAL', wg%x%type)
  call qp_get_axis ('X', wg%x%min, wg%x%max, wg%x%major_div, wg%x%places)
enddo

! Place the wave plot in the region

call tao_plot_struct_transfer (wave_plot, region%plot)
region%visible = .true.
region%plot%r => region

! Init some wave parameters

s%wave%ix_a1 = -1
s%wave%ix_a2 = -1
s%wave%ix_b1 = -1
s%wave%ix_b2 = -1

end subroutine tao_wave_cmd

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine tao_wave_analysis (plot)
!
!-

subroutine tao_wave_analysis (plot)

implicit none

type (tao_plot_struct), target :: plot
type (tao_graph_struct), pointer :: wg, wg1
type (tao_curve_struct), pointer :: wc, wc1

integer i, j, n, n_dat_arr, n_pt_curve, n_tot

character(20) :: r_name = 'tao_wave_analysis'
character(40) dtype

logical err

!

plot%graph(2)%x%draw_label = .false.
plot%graph(2)%x%draw_numbers = .false.
plot%graph(3)%x%draw_label = .false.
plot%graph(3)%x%draw_numbers = .false.

plot%graph(3)%title = 'After A-Region Subtraction'
plot%graph(4)%title = 'After B-Region Subtraction'

plot%graph(3)%title_suffix = ''
plot%graph(4)%title_suffix = ''

! Transfer information from the base curve to the wave curves.
! First count how many points we will have

wg1 => plot%graph(1)
wc1 => wg1%curve(1)

n_dat_arr = size(d1_ptr%d)
n_pt_curve = size(wc1%x_symb)
s%wave%i_wrap_pt = n_pt_curve

do i = 1, n_pt_curve
  if (wc1%ix_symb(i) > nint(0.5 * n_dat_arr)) then
    n_tot = n_pt_curve + i - 1
    exit
  endif
enddo

! Init the region placements if needed.

if (s%wave%ix_a1 == -1) then
  s%wave%ix_a1 = 5
  s%wave%ix_a2 = 15
  s%wave%ix_b1 = n_tot - 15
  s%wave%ix_b2 = n_tot - 5
endif

! Transfer the data to th wave arrays

do i = 2, 4
  wg => plot%graph(i)
  wc => wg%curve(1)
  call re_allocate (wc%x_line, n_tot)
  call re_allocate (wc%y_line, n_tot)
  call re_allocate (wc%x_symb, n_tot)
  call re_allocate (wc%y_symb, n_tot)
  call re_allocate (wc%ix_symb, n_tot)
  call re_allocate (s%wave%ix_data, n_tot)

  wc%x_symb(1:n_pt_curve)  = wc1%x_symb
  wc%y_symb(1:n_pt_curve)  = wc1%y_symb
  wc%ix_symb(1:n_pt_curve) = wc1%ix_symb

  wc%x_symb(n_pt_curve+1:n_tot)  = wc1%x_symb(1:n_tot-n_pt_curve) + n_dat_arr
  wc%y_symb(n_pt_curve+1:n_tot)  = wc1%y_symb(1:n_tot-n_pt_curve)
  wc%ix_symb(n_pt_curve+1:n_tot) = wc1%ix_symb(1:n_tot-n_pt_curve)

  wc%x_line = wc%x_symb
  wc%y_line = wc%y_symb

  s%wave%ix_data(1:n_pt_curve) = wc1%ix_symb
  s%wave%ix_data(n_pt_curve+1:n_tot) = wc1%ix_symb(1:n_tot-n_pt_curve) + n_dat_arr  

enddo

! Wave region setup

s%wave%i_a1 = -1
s%wave%i_b1 = -1
s%wave%n_a = 0
s%wave%n_b = 0

do i = 1, n_tot
  if (s%wave%ix_a1 <= s%wave%ix_data(i) .and. s%wave%ix_data(i) <= s%wave%ix_a2) &
                                                       s%wave%n_a = s%wave%n_a + 1
  if (s%wave%ix_b1 <= s%wave%ix_data(i) .and. s%wave%ix_data(i) <= s%wave%ix_b2) &
                                                       s%wave%n_b = s%wave%n_b + 1
  if (s%wave%ix_data(i) >= s%wave%ix_a1 .and. s%wave%i_a1 == -1) s%wave%i_a1 = i
  if (s%wave%ix_data(i) >= s%wave%ix_b1 .and. s%wave%i_b1 == -1) s%wave%i_b1 = i
enddo

s%wave%i_a2 = s%wave%i_a1 + s%wave%n_a - 1
s%wave%i_b2 = s%wave%i_b1 + s%wave%n_b - 1

if (s%wave%n_a < 5 .or. s%wave%n_b < 5) then
  call out_io (s_error$, r_name, 'Not enough points in one of the regions to do the analysis!')
  return
endif

! Specific analysis

dtype = plot%graph(1)%curve(1)%data_type
s%wave%data_type = dtype

if (dtype(1:5) == 'orbit') then
  call tao_orbit_wave_anal (plot)
elseif (dtype(1:5) == 'phase') then
  call tao_phase_wave_anal (plot)
elseif (dtype(1:4) == 'cbar') then
  call tao_cbar_wave_anal (plot)
else
  call out_io (s_error$, r_name, 'INVALID CURVE DATA_SOURCE FOR THIS OPERATION: ' // dtype)
  return
endif

end subroutine tao_wave_analysis

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine tao_orbit_wave_anal (plot)
!
!-

subroutine tao_orbit_wave_anal (plot)

implicit none

type (tao_plot_struct), target :: plot
type (tao_universe_struct), pointer :: u
type (lat_struct), pointer :: lat
type (tao_curve_struct), pointer :: curve

real(rp) dphi, coef_a(2), coef_b(2), rms_a(3), rms_b(3)
real(rp), allocatable :: phi(:), sin_phi(:), cos_phi(:)
real(rp) coef_ba(2), rms_ba(3), phi0, amp_a, amp_b, amp_ba, beta, tune
real(rp) phi1

integer i, j, ix, n_track, n_curve_pt, ix_ele, m_min, m_max, nc

character(40) :: r_name = 'tao_orbit_wave_anal'

! Init

curve => plot%graph(2)%curve(1)
u => tao_pointer_to_universe (curve%ix_universe)
lat => u%model%lat

n_track = lat%n_ele_track
n_curve_pt = size(curve%x_symb)

allocate (phi(n_curve_pt), sin_phi(n_curve_pt), cos_phi(n_curve_pt))

if (curve%data_type == 'orbit.x') then
  tune = lat%ele(n_track)%a%phi
elseif (curve%data_type == 'orbit.y') then
  tune = lat%ele(n_track)%b%phi
else
  call err_exit
endif

do i = 1, n_curve_pt

  ix = curve%ix_symb(i)

  dphi = 0
  if (i > s%wave%i_wrap_pt) dphi = tune

  ix_ele = d1_ptr%d(ix)%ix_ele

  if (curve%data_type == 'orbit.x') then
    phi(i) = lat%ele(ix_ele)%a%phi + dphi
    beta = lat%ele(ix_ele)%a%beta
  elseif (curve%data_type == 'orbit.y') then
    phi(i) = lat%ele(ix_ele)%b%phi + dphi
    beta = lat%ele(ix_ele)%b%beta
  endif

  sin_phi(i) = sin(phi(i)) * sqrt(beta)
  cos_phi(i) = cos(phi(i)) * sqrt(beta)

enddo

! Fit

s%wave%n_func = 2
call tao_wave_fit (curve, s%wave%i_a1, s%wave%n_a, coef_a, rms_a, sin_phi, cos_phi)
call tao_wave_fit (curve, s%wave%i_b1, s%wave%n_b, coef_b, rms_b, sin_phi, cos_phi)

! Put the results into the plot

plot%graph(3)%curve(1)%y_symb = plot%graph(2)%curve(1)%y_symb - &
                                  coef_a(1) * sin_phi - coef_a(2) * cos_phi
plot%graph(4)%curve(1)%y_symb = plot%graph(2)%curve(1)%y_symb - &
                                  coef_b(1) * sin_phi - coef_b(2) * cos_phi

plot%graph(3)%curve(1)%y_line = plot%graph(3)%curve(1)%y_symb
plot%graph(4)%curve(1)%y_line = plot%graph(4)%curve(1)%y_symb

! Compute the possible places between the regions where there could be a kick

coef_ba = coef_b - coef_a
rms_ba = sqrt(rms_a**2 + rms_b**2)

amp_a  = sqrt(sum(coef_a**2))
amp_b  = sqrt(sum(coef_b**2))
amp_ba = sqrt(sum(coef_ba**2))

if (coef_ba(1) == 0 .and. coef_ba(2) == 0) then
  phi0 = 0
else
  phi0 = -atan2(coef_ba(2), coef_ba(1))
endif

m_min = int((phi(s%wave%i_a2) - phi0)/pi) + 1
m_max = int((phi(s%wave%i_b1) - phi0)/pi)

nc = m_max - m_min + 1
s%wave%n_kick = nc
if (.not. allocated (s%wave%kick) .or. size(s%wave%kick) < nc) then
  if (allocated(s%wave%kick)) deallocate(s%wave%kick)
  allocate (s%wave%kick(nc))
endif

nc = 0
do i = m_min, m_max
  nc = nc + 1
  phi1 = i * pi + phi0
  s%wave%kick(nc)%amp = coef_ba(1) * cos(phi1) - coef_ba(2) * sin(phi1)
  do j = s%wave%i_a2, s%wave%i_b1
    if (phi(j) < phi1) s%wave%kick(nc)%ix_dat = s%wave%ix_data(j)
  enddo
  if (phi1 > tune) phi1 = phi1 - tune
  s%wave%kick(nc)%phi = phi1
enddo

! Compute the fit rms

if (amp_a == 0 .or. amp_b == 0 .or. amp_ba == 0) then
  s%wave%rms_rel_a = 0
  s%wave%rms_rel_b = 0
  s%wave%rms_rel_k = 0
  s%wave%rms_phi   = 0
else
  s%wave%rms_rel_a = sqrt(rms_a(1)**2 + rms_a(2)**2) / amp_a
  s%wave%rms_rel_b = sqrt(rms_b(1)**2 + rms_b(2)**2) / amp_b
  s%wave%rms_rel_k = sqrt(rms_ba(1)**2 * coef_ba(1)**2 + &
                            rms_ba(2)**2 * coef_ba(2)**2) / amp_ba**2
  s%wave%rms_phi   = sqrt(rms_ba(2)**2 * coef_ba(1)**2 +  &
                            rms_ba(1)**2 * coef_ba(2)**2) / amp_ba**2
endif

!s%wave%coef_a(1:2)  = coef_a
!s%wave%coef_b(1:2)  = coef_b
!s%wave%coef_ba(1:2) = coef_ba
!s%wave%rms_a(1:2)  = rms_a
!s%wave%rms_b(1:2)  = rms_b
!s%wave%rms_ba(1:2) = rms_ba
!s%wave%amp_a  = amp_a
!s%wave%amp_b  = amp_b
!s%wave%amp_ba = amp_ba

! Cleanup

deallocate (phi, sin_phi, cos_phi)

end subroutine tao_orbit_wave_anal

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine tao_phase_wave_anal (plot)
!
!-

subroutine tao_phase_wave_anal (plot)

implicit none

type (tao_plot_struct) plot

character(40) :: r_name = 'tao_phase_wave_anal'

!

call err_exit

end subroutine tao_phase_wave_anal

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine tao_cbar_wave_anal (plot)
!
!-

subroutine tao_cbar_wave_anal (plot)

implicit none

type (tao_plot_struct) plot

character(40) :: r_name = 'tao_cbar_wave_anal'

!

call err_exit

end subroutine tao_cbar_wave_anal

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine tao_wave_fit (curve, ix1, n_dat, coef, rms, f1, f2, f3, f4)
!
! Routine for fitting the curve data to up to four functions using a least squares
! SVD fit.
!
! Input:
!   curve -- Tao_curve_struct: Curve containing the data.
!   ix1   -- Integer: Index of first point in the data array.
!   n_dat -- Integer: Number of data points.
!   f1(:) -- Real(rp): First fit function.
!   f2(:) -- Real(rp), optional: Second fit function.
!   f3(:) -- Real(rp), optional: third fit function.
!   f4(:) -- Real(rp), optional: fourth fit function.
!
! Output:
!   coef(1:n_func)  -- Real(rp): Fit coefficients.
!   rms(1:n_func+1) -- Real(rp): Variances with rms(n_func+1) = sqrt(chi^2/n_dat).
!-

subroutine tao_wave_fit (curve, ix1, n_dat, coef, rms, f1, f2, f3, f4)

use nr, only: svdfit, svdvar

implicit none

!interface
!  function tao_wave_funcs(x,n)
!    use precision_def
!    implicit none
!    real(rp), intent(in) :: x
!    integer, intent(in) :: n
!    real(rp), dimension(n) :: tao_wave_funcs
!  end function tao_wave_funcs
!end interface


type (tao_curve_struct) curve

real(rp) f1(:), coef(:), rms(:)
real(rp), optional :: f2(:), f3(:), f4(:)
real(rp), allocatable :: cvm(:,:), v(:,:)
real(rp), allocatable :: x(:), sig(:), w(:)
real(rp) chisq

integer i, ix1, ix2, n_dat, n_tot, n_f

! If the number of data points is 1 or less then cannot do the analysis

if (n_dat < 2) then
  coef = 0
  rms = 0
  return
endif

! Setup

ix2 = ix1 + n_dat - 1

n_f = s%wave%n_func
allocate (x(n_dat), f_com(n_f, n_dat), v(n_f,n_f), w(n_f), cvm(n_f,n_f), sig(n_dat))

f_com(1,:) = f1(ix1:ix2)
if (present(f2)) f_com(2,:) = f2(ix1:ix2)
if (present(f3)) f_com(3,:) = f3(ix1:ix2)
if (present(f4)) f_com(4,:) = f4(ix1:ix2)

sig = 1
x = [(i, i = 1, n_dat)]

! Fit

sig = 1
call svdfit (x, curve%y_symb(ix1:ix2), sig, coef, v, w, chisq, tao_wave_funcs)
call svdvar (v, w, cvm)

rms(n_f+1) = sqrt(chisq/n_dat)
forall (i = 1:n_f) rms(i) = sqrt(cvm(i,i)) * rms(n_f+1)

! Cleanup

deallocate (x, f_com, cvm, v, sig, w)

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function tao_wave_funcs (x, ma) result (value)
!
! Function called by NR routine SVDFIT
!-

function tao_wave_funcs (x, ma) result (value)

implicit none

integer, intent(in) :: ma
real(rp), intent(in) :: x
real(rp) value(ma)

!

value = f_com(1:ma, nint(x))

end function tao_wave_funcs


end module
