module tao_wave_mod

use tao_graph_setup_mod

private ele_at_curve_point

contains

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine tao_wave_cmd (curve_name, plot_place, err_flag)
!
! Routine to do the initial setup for wave plotting.
! The wave analysis is done by the routine tao_wave_analysis.
!
! Input:
!   curve_name  -- Character(*) curve for wave analysis.
!   plot_place  -- Character(*) place on plot page to put the wave plot.
!-

subroutine tao_wave_cmd (curve_name, plot_place, err_flag)

use quick_plot

implicit none

type (tao_curve_array_struct), allocatable, target, save :: curve_array(:)
type (tao_curve_struct), pointer :: curve
type (tao_plot_struct), save, target :: wave_plot
type (tao_plot_struct), pointer :: plot
type (tao_graph_struct), pointer :: wg, wg0
type (tao_curve_struct), pointer :: wc, wc0
type (tao_d1_data_array_struct), allocatable :: d1_array(:)
type (tao_universe_struct), pointer :: u
type (branch_struct), pointer :: branch

integer i, p1, p2, ix_curve

character(*) curve_name, plot_place
character(*), parameter :: r_name = 'tao_wave_cmd'

logical :: init_needed = .true.
logical err_flag, err

! Find the curve

err_flag = .true.

if (any(curve_name == tao_wave_data_name)) then
  call tao_find_plots (err, '*', 'REGION', curve = curve_array, blank_means_all = .true.)
  if (err) return

  ix_curve = 0
  do i = 1, size(curve_array)
    if (curve_array(i)%c%data_type /= curve_name) cycle
    if (ix_curve /= 0) then
      call out_io (s_error$, r_name, 'Multiple curves have the same data_type as: ' // curve_name, 'Nothing done.')
      return
    endif
    ix_curve = i
  enddo

else
  call tao_find_plots (err, curve_name, 'REGION', curve = curve_array, blank_means_all = .true.)
  if (err) return

  ix_curve = 0
  do i = 1, size(curve_array) 
    if (.not. any(curve_array(i)%c%data_type == tao_wave_data_name)) cycle
    if (ix_curve /= 0) then
      call out_io (s_error$, r_name, 'Name does not resolve to a unique curve: ' // curve_name, 'Nothing done.')
      return
    endif
    ix_curve = i
  enddo
endif

if (ix_curve == 0) then
  call out_io (s_error$, r_name, 'No displayed curve found with this name: ' // curve_name)
  return
endif

curve => curve_array(ix_curve)%c
plot => curve%g%p

! Find the plot region

if (plot_place /= '') then
  call tao_find_plot_region (err, plot_place, s%wave%region)
  if (err) return
else
  s%wave%region => plot%r
endif

! Initiate the wave plot.
! There are three graphs. Each graph has a single curve.
!   Graph 1: extended original curve 
!   Graph 2: curve - a_region_fit
!   Graph 3: curve - b_region_fit
! Also the original graph is saved in the s%wave structure

if (init_needed) then
  allocate(wave_plot%graph(3))
  init_needed = .false.
endif

! Transfer curve information to the wave plot

wave_plot%name             = 'wave'
wave_plot%type             = 'wave'
wave_plot%x_axis_type      = plot%x_axis_type
wave_plot%autoscale_gang_x = .false.

do i = 1, 4
  if (i == 4) then
    wg => s%wave%base_graph
  else
    wg => wave_plot%graph(i)
  endif
  wg = curve%g
  wg%x = plot%graph(1)%x
  wg%box = [1, 4-i, 1, 3]

  if (size(wg%curve) /= 1) then
    deallocate (wg%curve)
    allocate (wg%curve(1))
    wg%curve(1) = curve
  endif
  wg%curve(1)%g => wg
enddo

!! wave_plot%graph(1)%visible = .false.
wave_plot%graph(1)%type = 'wave.0'  ! Extended origninal curve.
wave_plot%graph(2)%type = 'wave.a'  ! Curve - A region fit.
wave_plot%graph(3)%type = 'wave.b'  ! Curve - B region fit.

wave_plot%graph(1)%name = '0'  ! Extended origninal curve.
wave_plot%graph(2)%name = 'a'  ! Curve - A region fit.
wave_plot%graph(3)%name = 'b'  ! Curve - B region fit.

! Set up curves to wrap around the lattice, etc.

wg0 => s%wave%base_graph
wc0 => wg0%curve(1)

wc0%y_axis_scale_factor = 1

u => tao_pointer_to_universe (tao_curve_ix_uni(curve))
branch => u%model%lat%branch(0)

if (wc0%data_source /= 'data') then
  call out_io (s_error$, r_name, 'PLOT DATA_SOURCE FOR WAVE ANALYSIS MUST BE "data"')
  return
endif

call tao_find_data (err, wc0%data_type, d1_array = d1_array, ix_uni = tao_curve_ix_uni(wc0))
s%wave%d1_dat => d1_array(1)%d1
wg0%x%min = 0
if (plot%x_axis_type == 's') then
  wg0%x%max = branch%ele(branch%n_ele_track)%s
else
  wg0%x%max = ubound(s%wave%d1_dat%d, 1)
endif

do i = 1, 3
  wg => wave_plot%graph(i)
  p1 = nint(0.7 * wg%x%major_div_nominal)
  p2 = nint(1.3 * wg%x%major_div_nominal)
  if (branch%param%geometry == closed$) then
    call qp_calc_axis_params (0.0_rp, 1.5*wg0%x%max, p1, p2, wg%x)
  else
    call qp_calc_axis_params (0.0_rp, wg0%x%max, p1, p2, wg%x)
  endif
enddo

! Place the wave plot in the region

call tao_plot_struct_transfer (wave_plot, s%wave%region%plot)
s%wave%region%visible = .true.
s%wave%region%plot%r => s%wave%region
err_flag = .false.

end subroutine tao_wave_cmd

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine tao_wave_analysis (plot)
!
! Routine to do a wave anaysis. 
!
! Input:
!   plot      -- tao_plot_struct: Plot region setup by tao_wave_cmd.
!
! Output:
!   plot      -- tao_plot_struct: Plot with wave analysis curves.
!-

subroutine tao_wave_analysis (plot)

implicit none

type (tao_plot_struct), target :: plot
type (tao_graph_struct), pointer :: wg, wg0
type (tao_curve_struct), pointer :: wc, wc0
type (tao_universe_struct), pointer :: u
type (tao_d1_data_struct), pointer :: d1_dat
type (branch_struct), pointer :: branch

real(rp) dx
integer i, j, n, n_dat_arr, n_pt0, n_tot

character(*), parameter :: r_name = 'tao_wave_analysis'
character(40) dtype

logical err

!

call tao_graph_setup (plot, s%wave%base_graph)

plot%graph(1)%x%draw_label = .false.
plot%graph(1)%x%draw_numbers = .false.
plot%graph(2)%x%draw_label = .false.
plot%graph(2)%x%draw_numbers = .false.

plot%graph(2)%title = 'After A-Region Subtraction'
plot%graph(3)%title = 'After B-Region Subtraction'

plot%graph(2)%title_suffix = ''
plot%graph(3)%title_suffix = ''

! Transfer information from the base curve to the wave curves.
! First count how many points we will have

wg0 => s%wave%base_graph
wc0 => wg0%curve(1)

if (.not. allocated(wc0%x_line)) then
  n_pt0 = size(wc0%x_symb)
  allocate (wc0%x_line(n_pt0), wc0%y_line(n_pt0))
  wc0%x_line = wc0%x_symb
  wc0%y_line = wc0%y_symb
  if (allocated(wc0%ix_symb)) then
    call re_allocate(wc0%ix_line, n_pt0)
    wc0%ix_line = wc0%ix_symb
  endif
elseif (.not. allocated(wc0%x_symb)) then
  n_pt0 = size(wc0%x_line)
  allocate (wc0%x_symb(n_pt0), wc0%y_symb(n_pt0))
  wc0%x_symb  = wc0%x_line
  wc0%y_symb  = wc0%y_line
  if (allocated(wc0%ix_line)) then
    call re_allocate(wc0%ix_symb, n_pt0)
    wc0%ix_symb = wc0%ix_line
  endif
else
  n_pt0 = size(wc0%x_symb)
endif

u => tao_pointer_to_universe (tao_curve_ix_uni(wc0))
branch => u%model%lat%branch(0)

s%wave%i_curve_wrap_pt = n_pt0

d1_dat => s%wave%d1_dat
n_dat_arr = size(d1_dat%d)

if (branch%param%geometry == closed$ .and. wc0%data_source == 'data') then
  do i = 1, n_pt0
    if (wc0%ix_symb(i) > nint(0.5 * n_dat_arr)) then
      n_tot = n_pt0 + i - 1
      exit
    endif
  enddo
else
  n_tot = n_pt0
endif

! Init the region placements if needed.

if (s%wave%ix_a1 == -1) then
  s%wave%ix_a1 = 5
  s%wave%ix_a2 = 15
  s%wave%ix_b1 = n_pt0 - 15
  s%wave%ix_b2 = n_pt0 - 5
endif

! Transfer the data to the wave arrays

if (plot%x_axis_type == 's') then
  dx = u%model%lat%param%total_length
else
  dx = n_dat_arr
endif

do i = 1, 3
  wg => plot%graph(i)
  wc => wg%curve(1)
  call re_allocate (wc%x_line, n_tot)
  call re_allocate (wc%y_line, n_tot)
  call re_allocate (wc%x_symb, n_tot)
  call re_allocate (wc%y_symb, n_tot)
  call re_allocate (wc%ix_symb, n_tot)
  call re_allocate (s%wave%ix_data, n_tot)

  wc%x_symb(1:n_pt0)  = wc0%x_symb
  wc%y_symb(1:n_pt0)  = wc0%y_symb
  wc%ix_symb(1:n_pt0) = wc0%ix_symb

  wc%x_symb(n_pt0+1:n_tot)  = wc0%x_symb(1:n_tot-n_pt0) + dx
  wc%y_symb(n_pt0+1:n_tot)  = wc0%y_symb(1:n_tot-n_pt0)
  wc%ix_symb(n_pt0+1:n_tot) = wc0%ix_symb(1:n_tot-n_pt0)

  wc%x_line = wc%x_symb
  wc%y_line = wc%y_symb

  s%wave%ix_data(1:n_pt0) = wc0%ix_symb
  s%wave%ix_data(n_pt0+1:n_tot) = wc0%ix_symb(1:n_tot-n_pt0) + dx 

enddo

! Wave region setup

s%wave%i_a1 = -1
s%wave%i_b1 = -1
s%wave%n_a = 0
s%wave%n_b = 0

do i = 1, n_tot
  if (s%wave%ix_a1 <= s%wave%ix_data(i) .and. s%wave%ix_data(i) <= s%wave%ix_a2) s%wave%n_a = s%wave%n_a + 1
  if (s%wave%ix_b1 <= s%wave%ix_data(i) .and. s%wave%ix_data(i) <= s%wave%ix_b2) s%wave%n_b = s%wave%n_b + 1
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

dtype = s%wave%base_graph%curve(1)%data_type
s%wave%data_type = dtype

select case (dtype)
case ('orbit.x', 'orbit.y', 'beta.a', 'beta.b', 'eta.x', 'eta.y', 'ping_a.amp_x', 'ping_b.amp_y')
  call tao_orbit_beta_wave_anal (plot)

case ('phase.a', 'phase.b', 'ping_a.phase_x', 'ping_b.phase_y')
  call tao_phase_wave_anal (plot)

case ('ping_a.amp_sin_rel_y', 'ping_a.amp_cos_rel_y', 'ping_b.amp_sin_rel_x', 'ping_b.amp_cos_rel_x', &
      'ping_a.amp_sin_y', 'ping_a.amp_cos_y', 'ping_b.amp_sin_x', 'ping_b.amp_cos_x', &
      'ping_a.amp_y', 'ping_a.phase_y', 'ping_b.amp_x', 'ping_b.phase_x', 'cbar.12', 'cbar.11', 'cbar.22')
  call tao_cbar_wave_anal (plot)

case ('cbar.21')
  call out_io (s_error$, r_name, 'WAVE ANALYSIS IS NOT POSSIBLE FOR CBAR.21')

case default
  call out_io (s_error$, r_name, 'INVALID CURVE DATA_SOURCE FOR THIS OPERATION: ' // dtype)
end select

!

do i = 1, 3
  wg => plot%graph(i)
  wc => wg%curve(1)
  wc%x_symb = wc%x_symb * wg%x_axis_scale_factor
  wc%y_symb = wc%y_symb * wc%y_axis_scale_factor
  wc%x_line = wc%x_line * wg%x_axis_scale_factor
  wc%y_line = wc%y_line * wc%y_axis_scale_factor
enddo

end subroutine tao_wave_analysis

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine tao_orbit_beta_wave_anal (plot)
!
!-

subroutine tao_orbit_beta_wave_anal (plot)

implicit none

type (tao_plot_struct), target :: plot
type (tao_universe_struct), pointer :: u
type (lat_struct), pointer :: lat
type (tao_curve_struct), pointer :: curve
type (ele_struct), pointer :: ele, ele0

real(rp) dphi, coef_a(2), coef_b(2), coef_ba(2)
real(rp) rms_a(3), rms_b(3), rms_ba(3)
real(rp), allocatable :: phi(:), sin_phi(:), cos_phi(:)
real(rp) amp_a, amp_b, amp_ba, beta, tune
real(rp) phi0, phi_kick, s0, s1, phi1, phi2

integer i, j, n_track, n_curve_pt, m_min, m_max, nc

character(20) data_type
character(*), parameter :: r_name = 'tao_orbit_beta_wave_anal'

! Init

curve => plot%graph(1)%curve(1)
u => tao_pointer_to_universe (tao_curve_ix_uni(curve))
lat => u%model%lat

n_track = lat%n_ele_track
n_curve_pt = size(curve%x_symb)

! Eta and orbit analysis use the same formulas

data_type = curve%data_type
if (data_type(1:3) == 'eta') data_type = 'orbit' // trim(data_type(4:))

!

allocate (phi(n_curve_pt), sin_phi(n_curve_pt), cos_phi(n_curve_pt))

if (data_type == 'orbit.x' .or. data_type == 'beta.a' .or. data_type(1:6) == 'ping_a') then
  tune = lat%ele(n_track)%a%phi
elseif (data_type == 'orbit.y' .or. data_type == 'beta.b' .or. data_type(1:6) == 'ping_b') then
  tune = lat%ele(n_track)%b%phi
else
  call out_io (s_error$, r_name, 'UNKNOWN DATA_TYPE: ' // data_type)
  return
endif

do i = 1, n_curve_pt

  dphi = 0
  if (i > s%wave%i_curve_wrap_pt) dphi = tune

  ele => ele_at_curve_point(plot, i)

  if (data_type == 'orbit.x' .or. data_type == 'beta.a' .or. data_type(1:6) == 'ping_a') then
    phi(i) = ele%a%phi + dphi
    beta = ele%a%beta
  elseif (data_type == 'orbit.y' .or. data_type == 'beta.b' .or. data_type(1:6) == 'ping_b') then
    phi(i) = ele%b%phi + dphi
    beta = ele%b%beta
  endif

  if (data_type(1:5) == 'orbit') then
    sin_phi(i) = sin(phi(i)) * sqrt(beta)
    cos_phi(i) = cos(phi(i)) * sqrt(beta)
  elseif (data_type(1:4) == 'beta') then
    sin_phi(i) = sin(2*phi(i)) * beta
    cos_phi(i) = cos(2*phi(i)) * beta
    phi(i) = 2 * phi(i)
  else
    sin_phi(i) = sin(2*phi(i)) * sqrt(beta)
    cos_phi(i) = cos(2*phi(i)) * sqrt(beta)
    phi(i) = 2 * phi(i)    
  endif

enddo

if (data_type(1:4) == 'beta' .or. data_type(1:4) == 'ping') tune = 2 * tune

! Fit

s%wave%n_func = 2
call tao_wave_fit (curve, s%wave%i_a1, s%wave%n_a, coef_a, rms_a, sin_phi, cos_phi)
call tao_wave_fit (curve, s%wave%i_b1, s%wave%n_b, coef_b, rms_b, sin_phi, cos_phi)

! Put the results into the plot

plot%graph(2)%curve(1)%y_symb = plot%graph(1)%curve(1)%y_symb - coef_a(1) * sin_phi - coef_a(2) * cos_phi
plot%graph(3)%curve(1)%y_symb = plot%graph(1)%curve(1)%y_symb - coef_b(1) * sin_phi - coef_b(2) * cos_phi

plot%graph(2)%curve(1)%y_line = plot%graph(2)%curve(1)%y_symb
plot%graph(3)%curve(1)%y_line = plot%graph(3)%curve(1)%y_symb

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
  phi_kick = i * pi + phi0
  s%wave%kick(nc)%amp = coef_ba(1) * cos(phi_kick) - coef_ba(2) * sin(phi_kick)
  do j = s%wave%i_b1, s%wave%i_a2, -1
    if (phi(j) < phi_kick) exit
  enddo

  if (s%wave%ix_data(j) <= size(s%wave%d1_dat%d)) then
    s%wave%kick(nc)%ix_dat_before_kick = s%wave%ix_data(j)
  else
    s%wave%kick(nc)%ix_dat_before_kick = s%wave%ix_data(j) - size(s%wave%d1_dat%d)
  endif

  if (phi_kick > tune) phi_kick = phi_kick - tune
  if (data_type(1:4) == 'beta' .or. data_type(1:4) == 'ping') phi_kick = phi_kick / 2
  s%wave%kick(nc)%phi = phi_kick

  ele => ele_at_curve_point(plot, j)
  ele0 => pointer_to_next_ele(ele, -1)

  do    
    if (data_type == 'orbit.x' .or. data_type == 'beta.a' .or. data_type(1:6) == 'ping_a') then
      phi1 = ele0%a%phi
      phi2 = ele%a%phi
    elseif (data_type == 'orbit.y' .or. data_type == 'beta.b' .or. data_type(1:6) == 'ping_b') then
      phi1 = ele0%b%phi
      phi2 = ele%b%phi
    endif

    if (phi1 <= phi_kick .and. phi_kick <= phi2) exit
    ele => pointer_to_next_ele(ele, -1, .true.)
    ele0 => pointer_to_next_ele(ele, -1)
  enddo

  s%wave%kick(nc)%s = (ele0%s * (phi_kick - phi1) + ele%s * (phi2 - phi_kick)) / (phi2 - phi1)
  s%wave%kick(nc)%ele => ele
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
  if (data_type(1:4) == 'beta' .or. data_type(1:4) == 'ping') s%wave%rms_phi = s%wave%rms_phi / 2

endif

! Cleanup

deallocate (phi, sin_phi, cos_phi)

end subroutine tao_orbit_beta_wave_anal

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine tao_phase_wave_anal (plot)
!
!-

subroutine tao_phase_wave_anal (plot)

implicit none

type (tao_plot_struct), target :: plot
type (tao_universe_struct), pointer :: u
type (lat_struct), pointer :: lat
type (tao_curve_struct), pointer :: curve
type (ele_struct), pointer :: ele, ele0

real(rp) dphi, coef_a(3), coef_b(3), coef_ba(3)
real(rp) rms_a(4), rms_b(4), rms_ba(4)
real(rp), allocatable :: phi(:), sin_2phi(:), cos_2phi(:), one(:)
real(rp) amp_a, amp_b, amp_ba, beta, tune
real(rp) phi2_0, phi_kick, kick_amp, phi1, phi2

integer i, j, n_track, n_curve_pt, m_min, m_max, nc

character(*), parameter :: r_name = 'tao_phase_wave_anal'

! Init

curve => plot%graph(1)%curve(1)
u => tao_pointer_to_universe (tao_curve_ix_uni(curve))
lat => u%model%lat

n_track = lat%n_ele_track
n_curve_pt = size(curve%x_symb)

allocate (phi(n_curve_pt), sin_2phi(n_curve_pt), cos_2phi(n_curve_pt), one(n_curve_pt))

select case (curve%data_type)
case ('phase.a', 'ping_a.phase_x')
  tune = lat%ele(n_track)%a%phi
case ('phase.b', 'ping_b.phase_y')
  tune = lat%ele(n_track)%b%phi
case default
  call out_io (s_error$, r_name, 'UNKNOWN CURVE DATA_TYPE: ' // curve%data_type)
  return
end select

do i = 1, n_curve_pt

  dphi = 0
  if (i > s%wave%i_curve_wrap_pt) dphi = tune

  ele => ele_at_curve_point(plot, i)

  select case (curve%data_type)
  case ('phase.a', 'ping_a.phase_x')
    phi(i) = ele%a%phi + dphi
  case ('phase.b', 'ping_b.phase_y')
    phi(i) = ele%b%phi + dphi
  end select

  sin_2phi(i) = sin(2*phi(i))
  cos_2phi(i) = cos(2*phi(i))
  one(i)      = 1

enddo

! Fit

s%wave%n_func = 3
call tao_wave_fit (curve, s%wave%i_a1, s%wave%n_a, coef_a, rms_a, sin_2phi, cos_2phi, one)
call tao_wave_fit (curve, s%wave%i_b1, s%wave%n_b, coef_b, rms_b, sin_2phi, cos_2phi, one)

! Put the results into the plot

plot%graph(2)%curve(1)%y_symb = plot%graph(1)%curve(1)%y_symb - coef_a(1) * sin_2phi - coef_a(2) * cos_2phi - coef_a(3) * one
plot%graph(3)%curve(1)%y_symb = plot%graph(1)%curve(1)%y_symb - coef_b(1) * sin_2phi - coef_b(2) * cos_2phi - coef_b(3) * one

plot%graph(2)%curve(1)%y_line = plot%graph(2)%curve(1)%y_symb
plot%graph(3)%curve(1)%y_line = plot%graph(3)%curve(1)%y_symb

! Compute the possible places between the regions where there could be a kick

coef_ba = coef_b - coef_a
rms_ba = sqrt(rms_a**2 + rms_b**2)

amp_a  = sqrt(sum(coef_a(1:2)**2))
amp_b  = sqrt(sum(coef_b(1:2)**2))
amp_ba = sqrt(sum(coef_ba(1:2)**2))

! kick_amp is k of quad so positive k means horizontal focusing.

kick_amp = -sign(amp_ba, coef_ba(3)) - coef_ba(3)
if (curve%data_type == 'phase.a' .or. curve%data_type == 'ping_a.phase_x') kick_amp = -kick_amp

if (coef_ba(1) == 0 .and. coef_ba(2) == 0) then
  phi2_0 = 0
else
  phi2_0 = atan2(-coef_ba(1), -coef_ba(2))
endif

if (coef_ba(3) < 0) phi2_0 = phi2_0 + pi

!

m_min = int((2*phi(s%wave%i_a2) - phi2_0)/pi) + 1
m_max = int((2*phi(s%wave%i_b1) - phi2_0)/pi)

nc = m_max - m_min + 1
s%wave%n_kick = nc
if (.not. allocated (s%wave%kick) .or. size(s%wave%kick) < nc) then
  if (allocated(s%wave%kick)) deallocate(s%wave%kick)
  allocate (s%wave%kick(nc))
endif

nc = 0
do i = m_min, m_max
  nc = nc + 1
  phi_kick = (i * pi + phi2_0) / 2
  do j = s%wave%i_b1, s%wave%i_a2, -1
    if (phi(j) < phi_kick) exit
  enddo

  if (s%wave%ix_data(j) <= size(s%wave%d1_dat%d)) then
    s%wave%kick(nc)%ix_dat_before_kick = s%wave%ix_data(j)
  else
    s%wave%kick(nc)%ix_dat_before_kick = s%wave%ix_data(j) - size(s%wave%d1_dat%d)
  endif

  if (phi_kick > tune) phi_kick = phi_kick - tune
  s%wave%kick(nc)%phi = phi_kick
  s%wave%kick(nc)%amp = kick_amp * (-1)**i

  ele => ele_at_curve_point(plot, j)
  ele0 => pointer_to_next_ele(ele, -1)

  do
    select case (curve%data_type)
    case ('phase.a', 'ping_a.phase_x')
      phi1 = ele0%a%phi
      phi2 = ele%a%phi
    case ('phase.b', 'ping_b.phase_y')
      phi1 = ele0%b%phi
      phi2 = ele%b%phi
    end select

    if (phi1 <= phi_kick .and. phi_kick <= phi2) exit
    ele => pointer_to_next_ele(ele, -1, .true.)
    ele0 => pointer_to_next_ele(ele, -1)
  enddo

  s%wave%kick(nc)%s = (ele0%s * (phi_kick - phi1) + ele%s * (phi2 - phi_kick)) / (phi2 - phi1)
  s%wave%kick(nc)%ele => ele
enddo

! Compute the fit rms

if (amp_a == 0 .or. amp_b == 0 .or. amp_ba == 0) then
  s%wave%chi_c     = 0
  s%wave%rms_rel_a = 0
  s%wave%rms_rel_b = 0
  s%wave%rms_rel_k = 0
  s%wave%rms_phi   = 0
else
  s%wave%chi_c     = abs(abs(coef_ba(3)) - amp_ba) / abs(kick_amp)
  s%wave%rms_rel_a = sqrt(rms_a(1)**2 + rms_a(2)**2) / amp_a
  s%wave%rms_rel_b = sqrt(rms_b(1)**2 + rms_b(2)**2) / amp_b
  s%wave%rms_rel_k = sqrt(rms_ba(1)**2 * coef_ba(1)**2 + rms_ba(2)**2 * coef_ba(2)**2 + rms_ba(3)**2 * coef_ba(3)**2) / kick_amp**2
  s%wave%rms_phi   = sqrt(rms_ba(2)**2 * coef_ba(1)**2 + rms_ba(1)**2 * coef_ba(2)**2) / (2 * amp_ba**2)
endif

! Cleanup

deallocate (phi, sin_2phi, cos_2phi, one)

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

type (tao_plot_struct), target :: plot
type (tao_universe_struct), pointer :: u
type (lat_struct), pointer :: lat
type (tao_curve_struct), pointer :: curve
type (ele_struct), pointer :: ele, ele0

real(rp) coef_a(4), coef_b(4), coef_ba(4)
real(rp) rms_a(5), rms_b(5), rms_ba(5)
real(rp), allocatable :: phi_s(:), phi_r(:), sin_s(:), sin_r(:), cos_s(:), cos_r(:)
real(rp) amp_as, amp_ar, amp_bs, amp_br, amp_ba_s, amp_ba_r, beta, tune_s, tune_r
real(rp) phi0_s, phi0_r, phi1_s, phi1_r, phi_r_ave, sqrt_beta, phi1, phi2

integer i, j, m, n, n_track, n_curve_pt, m_min, m_max, nc

character(*), parameter :: r_name = 'tao_cbar_wave_anal'

! Init

curve => plot%graph(1)%curve(1)
u => tao_pointer_to_universe (tao_curve_ix_uni(curve))
lat => u%model%lat

n_track = lat%n_ele_track
n_curve_pt = size(curve%x_symb)

n = n_curve_pt
allocate (phi_s(n), phi_r(n), sin_s(n), cos_s(n), sin_r(n), cos_r(n))

tune_s = lat%ele(n_track)%a%phi + lat%ele(n_track)%b%phi
tune_r = lat%ele(n_track)%a%phi - lat%ele(n_track)%b%phi

do i = 1, n_curve_pt
  ele => ele_at_curve_point(plot, i)

  phi_s(i) = ele%a%phi + ele%b%phi
  phi_r(i) = ele%a%phi - ele%b%phi

  if (i > s%wave%i_curve_wrap_pt) then
    phi_s(i) = ele%a%phi + ele%b%phi + tune_s
    phi_r(i) = ele%a%phi - ele%b%phi + tune_r
  endif

  if (substr(curve%data_type,1,4) == 'cbar') then
    sin_s(i) = sin(phi_s(i))
    cos_s(i) = cos(phi_s(i))
    sin_r(i) = sin(phi_r(i))
    cos_r(i) = cos(phi_r(i))
  elseif (substr(curve%data_type,1,6) == 'ping_a') then
    sqrt_beta = sqrt(ele%b%beta)
    sin_s(i) = sqrt_beta * sin(phi_s(i))
    cos_s(i) = sqrt_beta * cos(phi_s(i))
    sin_r(i) = sqrt_beta * sin(phi_r(i))
    cos_r(i) = sqrt_beta * cos(phi_r(i))
  elseif (substr(curve%data_type,1,6) == 'ping_b') then
    sqrt_beta = sqrt(ele%a%beta)
    sin_s(i) = sqrt_beta * sin(phi_s(i))
    cos_s(i) = sqrt_beta * cos(phi_s(i))
    sin_r(i) = sqrt_beta * sin(phi_r(i))
    cos_r(i) = sqrt_beta * cos(phi_r(i))
  endif
enddo

! Fit

s%wave%n_func = 4
call tao_wave_fit (curve, s%wave%i_a1, s%wave%n_a, coef_a, rms_a, sin_s, cos_s, sin_r, cos_r)
call tao_wave_fit (curve, s%wave%i_b1, s%wave%n_b, coef_b, rms_b, sin_s, cos_s, sin_r, cos_r)

! Put the results into the plot

plot%graph(2)%curve(1)%y_symb = plot%graph(1)%curve(1)%y_symb - coef_a(1) * sin_s - coef_a(2) * cos_s - coef_a(3) * sin_r - coef_a(4) * cos_r
plot%graph(3)%curve(1)%y_symb = plot%graph(1)%curve(1)%y_symb - coef_b(1) * sin_s - coef_b(2) * cos_s - coef_b(3) * sin_r - coef_b(4) * cos_r

plot%graph(2)%curve(1)%y_line = plot%graph(2)%curve(1)%y_symb
plot%graph(3)%curve(1)%y_line = plot%graph(3)%curve(1)%y_symb

! Compute the possible places between the regions where there could be a kick

coef_ba = coef_b - coef_a
rms_ba = sqrt(rms_a**2 + rms_b**2)

amp_as  = sqrt(sum(coef_a(1:2)**2))
amp_bs  = sqrt(sum(coef_b(1:2)**2))
amp_ba_s = sqrt(sum(coef_ba(1:2)**2))

amp_ar   = sqrt(sum(coef_a(3:4)**2))
amp_br   = sqrt(sum(coef_b(3:4)**2))
amp_ba_r = sqrt(sum(coef_ba(3:4)**2))

if (coef_ba(1) == 0 .and. coef_ba(2) == 0) then
  phi0_s = 0
else
  select case (curve%data_type)
  case ('cbar.11', 'cbar.22', 'ping_a.amp_cos_rel_y', 'ping_b.amp_cos_rel_x')
    phi0_s = atan2(-coef_ba(2), coef_ba(1))
  case ('cbar.12', 'ping_a.amp_sin_rel_y', 'ping_b.amp_sin_rel_x')
    phi0_s = atan2(coef_ba(1), coef_ba(2))
  end select
endif

if (coef_ba(3) == 0 .and. coef_ba(4) == 0) then
  phi0_r = 0
else
  select case (curve%data_type)
  case ('cbar.11', 'cbar.22', 'ping_a.amp_cos_rel_y', 'ping_b.amp_cos_rel_x')
    phi0_r = atan2(-coef_ba(4), coef_ba(3))
  case ('cbar.12', 'ping_a.amp_sin_rel_y', 'ping_b.amp_sin_rel_x')
    phi0_r = atan2(coef_ba(3), coef_ba(4))
  end select
endif

m_min = int((phi_s(s%wave%i_a2) - phi0_s)/pi) + 1
m_max = int((phi_s(s%wave%i_b1) - phi0_s)/pi)

nc = m_max - m_min + 1
s%wave%n_kick = nc
if (.not. allocated (s%wave%kick) .or. size(s%wave%kick) < nc) then
  if (allocated(s%wave%kick)) deallocate(s%wave%kick)
  allocate (s%wave%kick(nc))
endif

nc = 0
do i = m_min, m_max
  nc = nc + 1
  phi1_s = i * pi + phi0_s
  s%wave%kick(nc)%amp = -2 * (coef_ba(1) * sin(phi1_s) + coef_ba(2) * cos(phi1_s))
  do j = s%wave%i_b1, s%wave%i_a2, -1
    if (phi_s(j) < phi1_s) exit
  enddo

  if (s%wave%ix_data(j) <= size(s%wave%d1_dat%d)) then
    s%wave%kick(nc)%ix_dat_before_kick = s%wave%ix_data(j)
  else
    s%wave%kick(nc)%ix_dat_before_kick = s%wave%ix_data(j) - size(s%wave%d1_dat%d)
  endif

  phi_r_ave = (phi_r(j) + phi_r(j+1)) / 2
  m = nint((phi_r_ave - phi0_r) / pi)
  phi1_r = phi0_r + m * pi

  if (phi1_s > tune_s) then
    phi1_s = phi1_s - tune_s
    phi1_r = phi1_r - tune_r
  endif

  s%wave%kick(nc)%phi_s = phi1_s
  s%wave%kick(nc)%phi_r = phi1_r

  ele => ele_at_curve_point(plot, j)
  ele0 => pointer_to_next_ele(ele, -1)

  do
    phi1 = ele0%a%phi + ele0%b%phi
    phi2 = ele%a%phi + ele%b%phi

    if (phi1 <= phi1_s .and. phi1_s <= phi2) exit
    ele => pointer_to_next_ele(ele, -1, .true.)
    ele0 => pointer_to_next_ele(ele, -1)
  enddo

  s%wave%kick(nc)%s = (ele0%s * (phi1_s - phi1) + ele%s * (phi2 - phi1_s)) / (phi2 - phi1)
  s%wave%kick(nc)%ele => ele
enddo

! Compute the fit rms

if (amp_ar == 0 .or. amp_as == 0 .or. amp_br == 0 .or. amp_bs == 0) then
  s%wave%chi_a      = 0
  s%wave%rms_rel_ar = 0
  s%wave%rms_rel_as = 0
  s%wave%rms_rel_br = 0
  s%wave%rms_rel_bs = 0
  s%wave%rms_rel_ks = 0
  s%wave%rms_rel_kr = 0
  s%wave%rms_phi_s  = 0
  s%wave%rms_phi_r  = 0
else
  s%wave%chi_a      = abs(amp_ba_s - amp_ba_r) / max(amp_ba_s + amp_ba_r, 1e-10)
  s%wave%rms_rel_ar = sqrt(rms_a(3)**2 + rms_a(4)**2) / max(amp_ar, 1e-10)
  s%wave%rms_rel_as = sqrt(rms_a(1)**2 + rms_a(2)**2) / max(amp_as, 1e-10)
  s%wave%rms_rel_bs = sqrt(rms_b(1)**2 + rms_b(2)**2) / max(amp_bs, 1e-10)
  s%wave%rms_rel_br = sqrt(rms_b(3)**2 + rms_b(4)**2) / max(amp_br, 1e-10)
  s%wave%rms_rel_ks = sqrt(rms_ba(1)**2 * coef_ba(1)**2 + &
                        rms_ba(2)**2 * coef_ba(2)**2) / max(amp_ba_s**2, 1e-10)
  s%wave%rms_rel_kr = sqrt(rms_ba(3)**2 * coef_ba(3)**2 + &
                        rms_ba(4)**2 * coef_ba(4)**2) / max(amp_ba_r**2, 1e-10)
  s%wave%rms_phi_s = sqrt(rms_ba(2)**2 * coef_ba(1)**2 + &
                        rms_ba(1)**2 * coef_ba(2)**2) / max(amp_ba_s**2, 1e-10)
  s%wave%rms_phi_r = sqrt(rms_ba(4)**2 * coef_ba(3)**2 + &
                        rms_ba(3)**2 * coef_ba(4)**2) / max(amp_ba_r**2, 1e-10)
endif

! Cleanup

deallocate (phi_s, phi_r, sin_s, sin_r, cos_s, cos_r)

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

implicit none

type (tao_curve_struct) curve
real(rp), allocatable :: A(:,:), A2(:,:)

real(rp) f1(:), coef(:), rms(:)
real(rp), optional :: f2(:), f3(:), f4(:)
real(rp), allocatable :: cvm(:,:), v(:,:), w(:), b(:), wti(:)
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
allocate (A(n_dat, n_f), v(n_f,n_f), w(n_f), wti(n_f), cvm(n_f,n_f), b(n_dat))

A(:,1) = f1(ix1:ix2)
if (present(f2)) A(:,2) = f2(ix1:ix2)
if (present(f3)) A(:,3) = f3(ix1:ix2)
if (present(f4)) A(:,4) = f4(ix1:ix2)
A2 = A

b = curve%y_symb(ix1:ix2)

! Fit

call svd_fit(A2, b, 1.0e-5_rp, coef, chisq, w, v)

wti = 0
where (w /= 0) wti = 1.0_rp / (w*w)

cvm = v * spread(wti, dim=1, ncopies=n_f)
cvm = matmul(cvm, transpose(v))

rms(n_f+1) = sqrt(chisq/n_dat)
forall (i = 1:n_f) rms(i) = sqrt(cvm(i,i)) * rms(n_f+1)

end subroutine tao_wave_fit

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

function ele_at_curve_point (plot, ix_pt) result (ele)

use twiss_and_track_mod, only: twiss_and_track_at_s

implicit none

type (tao_plot_struct), target :: plot
type (ele_struct), pointer :: ele
type (tao_curve_struct), pointer :: curve, curve2
type (tao_d1_data_array_struct), allocatable :: d1_array(:)
type (tao_universe_struct), pointer :: u
type (tao_d1_data_struct), pointer :: d1_dat
type (branch_struct), pointer :: branch

real(rp) s_pos

integer ix_pt
integer ix, ix_ele

logical err
character(*), parameter :: r_name = 'ele_at_curve_point'

!

curve => s%wave%base_graph%curve(1)
curve2 => plot%graph(1)%curve(1)

u => tao_pointer_to_universe (tao_curve_ix_uni(curve))
branch => u%model%lat%branch(0)

if (curve%data_source == 'data') then
  call tao_find_data (err, curve%data_type, d1_array = d1_array, ix_uni = tao_curve_ix_uni(curve))
  d1_dat => d1_array(1)%d1
  ix = curve2%ix_symb(ix_pt)
  ix_ele = d1_dat%d(ix)%ix_ele
  ele => branch%ele(ix_ele)
else
  call out_io (s_error$, r_name, 'ANALYSIS FOR THIS DATA_SOURCE NOT YET IMPLEMENTED: ' // curve%data_source)
  nullify(ele)
  return
endif

end function ele_at_curve_point

end module
