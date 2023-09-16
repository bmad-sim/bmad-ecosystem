!+
! Subroutine lat_to_ptc_layout (lat)
!
! Subroutine to create a PTC layout from a Bmad lattice branch.
! Note: If lat_to_ptc_layout has already been setup, you should first do a
!           call kill_ptc_layouts(lat)
! This deallocates the pointers in PTC
!
! Note: If not already done, before you call this routine you need to first call:
!    call set_ptc (...)
! [This is normally done in bmad_parser.]
!
! Note: If a Bmad element is using a hard edge model (EG: RFcavity element), there 
! will be three corresponding PTC fibre elements: (drift, RF. drift) for example.
! In this case, ele%ptc_fibre will be set to point to the last PTC fibre. That is the 
! exit end of ele will correspond to the exit end of ele%ptc_fibre.
!
! Note: Photon branches will not be included in the layout.
!
! Input:
!   lat                   -- lat_struct: Input lattice
!
! Output:
!   lat%branch(:)%ptc              -- Pointers to generated layouts.
!   lat%branch(:)%ele(:)%ptc_fibre -- Pointer to PTC fibres
!-

subroutine lat_to_ptc_layout (lat)

use ptc_layout_mod, dummy => lat_to_ptc_layout
use madx_ptc_module, only: m_u, m_t, fibre, append_empty_layout, survey, make_node_layout, &
                           append_point, set_up, ring_l, survey

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele, ele0
type (fibre), pointer :: save_fib
type (layout), pointer :: lay

real(rp) ptc_orientation(3,3), ang(3)

integer i, j, ix_pass, i_save
logical logic

character(20), parameter :: r_name = 'lat_to_ptc_layout'

! Setup m_u

do i = 0, ubound(lat%branch, 1)
  branch => lat%branch(i)
  if (branch%param%particle == photon$) cycle
  call branch_to_ptc_m_u (branch)
enddo

! setup m_t

do i = 0, ubound(lat%branch, 1)

  branch => lat%branch(i)
  if (branch%param%particle == photon$) cycle

  call append_empty_layout(m_t)
  call set_up (m_t%end)
  write (m_t%end%name, '(a, i4)') 'm_t Bmad branch:', i  ! For diagnostic purposes

  branch%ptc%m_t_layout => m_t%end   ! Save layout

  lay => m_t%end

  do j = 0, branch%n_ele_track
    ele => branch%ele(j)

    if (.not. associated(ele%ptc_fibre)) then
      call out_io (s_fatal$, r_name, 'NO FIBRE ASSOCIATED WITH ELEMENT: ' // ele%name)
      if (global_com%exit_on_error) call err_exit
    endif

    call append_this_fibre(ele%ptc_fibre, .true.)

    ! Must add an energy patch if the reference energy shifts.

    if (j > 0) then
      ele0 => branch%ele(j-1)
      if (ele0%ptc_fibre%mag%p%p0c /= ele%ptc_fibre%mag%p%p0c) then
        ele%ptc_fibre%patch%energy = 1
      endif
    endif

  enddo

  lay%closed = .true.

  logic = .true.
  call ring_l(lay, logic)

  ele => branch%ele(0)
  ang = [ele%floor%phi, ele%floor%theta, ele%floor%psi]
  call bmad_patch_parameters_to_ptc (ang, ptc_orientation)

  !! call survey (m_t%end, ptc_orientation, ele%floor%r)

  ! Create the integration node arrays for the fibres.

  call set_ptc_quiet (12, set$, i_save)
  call make_node_layout(m_t%end)
  call survey (m_t%end%start, a = ele%floor%r, ent = ptc_orientation)
  call set_ptc_quiet (12, unset$, i_save)

  ! Beambeam elements are special. 
  ! Now that the integration node arrays have been created, any beambeam elements may be setup.

  do j = 0, branch%n_ele_track
    ele => branch%ele(j)
    call misalign_ptc_fibre (ele, .true., ele%ptc_fibre, .true.)
    if (ele%key == beambeam$ .and. ele%is_on) call beambeam_fibre_setup(ele, ele%ptc_fibre)
  enddo

enddo

!-----------------------------------------------------------------------------
contains

subroutine append_this_fibre(ele_fib, do_point)

type (fibre), pointer :: ele_fib
type (fibre), pointer :: this_fib
logical, optional :: do_point

!

call append_point(lay, ele_fib)
this_fib => lay%end

if (ele%key == patch$ .or. ele%key == floor_shift$) then
  this_fib%dir = ele%value(upstream_coord_dir$)
else
  this_fib%dir = ele%orientation
endif

this_fib%charge = charge_of(default_tracking_species(branch%param)) / charge_of(branch%param%particle)

if (logic_option(.false., do_point)) ele%ptc_fibre => this_fib

end subroutine append_this_fibre

end subroutine lat_to_ptc_layout 

