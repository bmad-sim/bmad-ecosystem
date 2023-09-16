!+
! Program ptc_profiler
!
! Program to profile (show computation times) for PTC.
!-

program ptc_profiler

use bmad
use pointer_lattice, only: c_linear_map, SPIN0, operator(+), track_probe, bmadl, probe, &
                      c_damap, assignment(=), alloc, kill

implicit none

type profile_struct
  type (spin_orbit_map1_struct) :: map
  real(rp) t_taylor, t_probe, t_fibre
end type

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (taylor_struct), pointer :: st
type (coord_struct), allocatable :: orbit(:)
type (spin_orbit_map1_struct), target :: map0
type (profile_struct), target :: new(20), old(20)

real(rp) stop_time, spin_diff_old, spin_diff_new
integer ie, ii, i, k, n, p, is, ns, n_time_calc
integer num_steps(20), integrator_order(3)

logical err

character(100) lat_file
character(2), parameter :: q_name(0:3) = ['q0', 'qx', 'qy', 'qz']

namelist / params / num_steps, lat_file, stop_time, n_time_calc, integrator_order

!

lat_file = 'lat.bmad'
num_steps = -1
stop_time = -1
integrator_order = -1

open (1, file = 'ptc_profiler.init', status = 'old')
read (1, nml = params)
close (1)

call bmad_parser (lat_file, lat, .false.)

bmad_com%spin_tracking_on = .true.
call reallocate_coord(orbit, lat)

do ie = 1, lat%n_ele_track
  ele => lat%ele(ie)
  if (ele%name == 'END') cycle

  call init_coord (orbit(ie-1), lat%particle_start, ele, upstream_end$)

  do ii = 1, size(integrator_order)
    if (integrator_order(ii) == -1) exit
    ele%value(integrator_order$) = integrator_order(ii)

    do is = 1, size(num_steps)
      if (num_steps(is) == -1) exit

      ele%value(num_steps$) = num_steps(is)
      call set_flags_for_changed_attribute(ele, ele%value(num_steps$))

      ptc_com%old_integrator = 1   ! True
      call kill_taylor (ele%taylor);  call kill_taylor(ele%spin_taylor)
      call spin_concat_linear_maps(err, old(is)%map, lat%branch(0), ie-1, ie, orbit = orbit)

      old(is)%t_taylor = taylor_timer(ele, orbit(ie-1), stop_time, n_time_calc)
      old(is)%t_probe  = probe_timer(ele, orbit(ie-1), stop_time, n_time_calc)
      old(is)%t_fibre  = fibre_timer(ele, orbit(ie-1), stop_time, n_time_calc)

      ptc_com%old_integrator = -1  ! False
      call kill_taylor (ele%taylor);  call kill_taylor(ele%spin_taylor)
      call spin_concat_linear_maps(err, new(is)%map, lat%branch(0), ie-1, ie, orbit = orbit)

      new(is)%t_taylor = taylor_timer(ele, orbit(ie-1), stop_time, n_time_calc)
      new(is)%t_probe  = probe_timer(ele, orbit(ie-1), stop_time, n_time_calc)
      new(is)%t_fibre  = fibre_timer(ele, orbit(ie-1), stop_time, n_time_calc)
    enddo

    map0 = new(is-1)%map

    print *
    print '(2a)',     'Element: ', trim(ele%name)
    print '(a, 2i3)', 'Taylor_order:     ', bmad_com%taylor_order, ptc_private%taylor_order_ptc
    print '(a, i3)',  'Integrator_order: ', nint(ele%value(integrator_order$))
    print '(a, l1)',  'Exact_model:      ', ptc_com%exact_model
    print '(a)', '      |                    old spin integrator               |                    new spin integrator'
    print '(a)', 'Steps | %spin_diff    orb_diff   t_probe   t_fibre  t_taylor | %spin_diff    orb_diff   t_probe   t_fibre  t_taylor'

    do is = 1, size(num_steps)
      if (num_steps(is) == -1) exit
      spin_diff_old = spin_diff(old(is)%map, map0)
      spin_diff_new = spin_diff(new(is)%map, map0)
      print '(i5, 2(f13.6, f12.8, 3f10.1))', num_steps(is), &
                    spin_diff_old, orbit_diff(old(is)%map, map0), old(is)%t_probe, old(is)%t_fibre, old(is)%t_taylor, &
                    spin_diff_new, orbit_diff(new(is)%map, map0), new(is)%t_probe, new(is)%t_fibre, new(is)%t_taylor
    enddo
  enddo

  print *
  print '(14x, a, 7x, 6(2x, a, 7x))', '0th order', 'dx  ', 'dpx', 'dy ', 'dpy', 'dz ', 'dpz'
  do i = 0, 3
    print '(i4, 2x, a, 2x, f12.6, 4x, 6f12.6)', i, q_name(i), map0%spin_q(i,:)
  enddo
enddo

!-----------------------------------------------------------------
contains

! Returns time to do n_time_calc ele_to_taylor calculations.

function taylor_timer(ele, orb0, stop_time, n_time_calc) result (t_taylor)

type (ele_struct) ele
type (coord_struct) orb0

real(rp) stop_time, t_taylor, this_t
integer n_time_calc, ic

!

call run_timer ('START')
ic = 0

do
  call ele_to_taylor (ele, ele%branch%param, orb0)
  ic = ic + 1
  call run_timer ('READ', this_t)
  if (this_t > stop_time) exit
enddo

t_taylor = (this_t * n_time_calc) / ic

end function taylor_timer

!-----------------------------------------------------------------
! contains

! Returns time to do n_time_calc track_probe calculations.

function probe_timer(ele, orb0, stop_time, n_time_calc) result (t_probe)

type (ele_struct) ele
type (coord_struct) orb0
type (probe) ptc_probe
type (probe_8) ptc_probe8
type (c_damap) ptc_cdamap

real(rp) stop_time, t_probe, this_t
integer n_time_calc, ic

!

call run_timer ('START')
ic = 0

do
  call alloc(ptc_cdamap)
  call alloc(ptc_probe8)
  ptc_probe = 0
  ptc_probe = orb0%vec
  ptc_cdamap = 1
  ptc_probe8 = ptc_cdamap + ptc_probe

  call track_probe (ptc_probe8, ptc_private%base_state+SPIN0, fibre1 = bmadl%start)
  ic = ic + 1

  call kill (ptc_probe8)
  call kill (ptc_cdamap)

  call run_timer ('READ', this_t)
  if (this_t > stop_time) exit
enddo

t_probe = (this_t * n_time_calc) / ic

end function probe_timer

!-----------------------------------------------------------------
! contains

! Returns time to do n_time_calc ele_to_fibre calculations.

function fibre_timer(ele, orb0, stop_time, n_time_calc) result (t_fibre)

type (ele_struct) ele
type (coord_struct) orb0
type (fibre), pointer :: ptc_fibre

real(rp) stop_time, t_fibre, this_t
integer n_time_calc, ic
logical err_flag

!

call run_timer ('START')
ic = 0

do
  call ele_to_fibre (ele, ptc_fibre, ele%branch%param, .false., err_flag, ref_in = orb0)
  ic = ic + 1
  call run_timer ('READ', this_t)
  if (this_t > stop_time) exit
enddo

t_fibre = (this_t * n_time_calc) / ic

end function fibre_timer

!-----------------------------------------------------------------
! contains

function spin_diff (map, map0) result (diff)

type (spin_orbit_map1_struct) map, map0
real(rp) diff, denom, vec(0:3), vec0(0:3), dvec
integer i, j

!

diff = 0
do i = 1, 6
  vec = map%spin_q(:,i)
  vec0 = map0%spin_q(:,i)
  dvec = 0
  do j = 0, 3
    if (vec(j) == 0 .and. vec0(j) == 0) cycle
    dvec = max(dvec, abs(vec(j)-vec0(j))) 
  enddo
  if (dvec == 0) cycle
  if (mod(i,2) == 0) then
    diff = max(diff, 3*dvec/(norm2(map0%spin_q(:,2)) + norm2(map0%spin_q(:,4)) + norm2(map0%spin_q(:,6))))
  else
    diff = max(diff, 3*dvec/(norm2(map0%spin_q(:,1)) + norm2(map0%spin_q(:,3)) + norm2(map0%spin_q(:,5))))
  endif
enddo

end function spin_diff

!-----------------------------------------------------------------
! contains

function orbit_diff (map, map0) result (diff)

type (spin_orbit_map1_struct) map, map0
real(rp) diff, dmat(6,6)
integer i, j

!

dmat = map%orb_mat - map0%orb_mat
diff = maxval(abs(dmat))

end function orbit_diff

end program
