module tao_wave_mod

use tao_utils

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

use tao_struct
use tao_interface

implicit none

type (tao_curve_array_struct), allocatable, target, save :: curve_array(:)
type (tao_curve_struct), pointer :: curve
type (tao_plot_region_struct), pointer :: region
type (tao_plot_struct), save :: wave_plot
type (tao_plot_struct), pointer :: plot

integer i

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

do i = 1, 4
  wave_plot%graph(1) = curve%g
  wave_plot%graph(i)%x = plot%x
  if (size(wave_plot%graph(i)%curve) /= 1) then
    deallocate (wave_plot%graph(i)%curve)
    allocate (wave_plot%graph(i)%curve(1))
    wave_plot%graph(i)%curve(1) = curve
  endif
enddo

wave_plot%graph(1)%visible = .false.

wave_plot%graph(2)%type = 'wave.0'  ! Extended origninal curve.
wave_plot%graph(3)%type = 'wave.a'  ! Curve - A region fit.
wave_plot%graph(4)%type = 'wave.b'  ! Curve - B region fit.

! Place the wave plot in the region

call tao_plot_struct_transfer (curve%g%p, region%plot)
region%visible = .true.
region%plot%r => region

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

type (tao_plot_struct) plot

character(20) :: r_name = 'tao_wave_analysis'
character(40) dtype

!

dtype = plot%graph(1)%curve(1)%data_type

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

type (tao_plot_struct) plot

character(40) :: r_name = 'tao_orbit_wave_anal'


!



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

end module
