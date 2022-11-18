!+
! Subroutine ptc_linear_isf_calc (branch, ele_isf)
!
! Routine to calculate the G-matrix and other stuff relavent to spin matching.
!
! Note: A call to lat_to_ptc_layout must be done before calling this routine.
!
! Input:
!   branch          -- branch_struct: Lattice branch to analyze.
!
! Output:
!   ele_isf(0:)     -- linear_ele_isf_struct, allocatable: ISF at every element.
!-

subroutine ptc_linear_isf_calc (branch, ele_isf)

use bmad_routine_interface, dummy => ptc_linear_isf_calc
use pointer_lattice

implicit none

type (branch_struct), target :: branch
type (linear_ele_isf_struct), allocatable, target :: ele_isf(:)
type (ele_struct), pointer :: ele
type (linear_ele_isf_struct), pointer :: eisf
type (linear_isf1_struct), pointer :: isf1
type (fibre), pointer :: ptc_fibre, fib_now, fib_next
type (integration_node), pointer :: node_now
type (internal_state) ptc_state
type (c_damap) cdamap, u, cdamap_1
type (probe) probe_orb
type (probe_8) p8_a
type (c_normal_form) cc_norm
type (c_linear_map) q0, q2, q_invar, q_y, q_1turn

integer n, ie, ii
logical rf_on

!

rf_on = rf_is_on(branch)
if (rf_on) then
  ptc_state = ptc_private%base_state - NOCAVITY0 + SPIN0
else
  ptc_state = ptc_private%base_state + NOCAVITY0 + SPIN0
endif

! Allocate ele_isf

if (allocated(ele_isf)) then
  if (ubound(ele_isf, 1) /= branch%n_ele_track) deallocate(ele_isf)
endif
if (.not. allocated(ele_isf)) allocate(ele_isf(0:branch%n_ele_track))

do ie = 0, branch%n_ele_track
  ele => branch%ele(ie)
  eisf => ele_isf(ie)
  n = ele%ptc_fibre%t2%pos_in_fibre
  if (allocated(eisf%node)) then
    if (size(eisf%node) /= n) deallocate(eisf%node)
  endif
  if (.not. allocated(eisf%node)) allocate (eisf%node(n))
enddo

!

call init_all(ptc_state, 1, 0)   ! Only need first order map for this analysis

call alloc(cdamap, cdamap_1, u)
call alloc(p8_a)
call alloc(cc_norm)

q_y = 2
cdamap_1 = 1   ! Unit map

ele => branch%ele(0)
eisf => ele_isf(0)

ptc_fibre => pointer_to_fibre(ele)

!

if (branch%param%geometry == closed$) then
  call find_orbit_x (eisf%node(1)%orb0, ptc_state, 1.e-8_rp, fibre1 = ptc_fibre) 

  probe_orb = eisf%node(1)%orb0

  p8_a = probe_orb + cdamap_1

  call track_probe(p8_a, ptc_state, fibre1 = ptc_fibre)

  cdamap = p8_a

  call ptc_set_rf_state_for_c_normal(ptc_state%nocavity)
  call c_normal(cdamap, cc_norm, dospin = .true.)

  cc_norm%n = cc_norm%atot**(-1) * cdamap * cc_norm%atot
  p8_a = probe_orb + cc_norm%atot

  u = cc_norm%atot
  p8_a = probe_orb + u

else
  !! p8_a = probe_orb + u
endif

fib_now => ptc_fibre
node_now => fib_now%t1
q_invar = 1   ! Set %mat = unit matrix

! Track

do ie = 0, branch%n_ele_track
  do ii = 1, fib_now%t2%pos_in_fibre
    if (ie /= 0 .and. ii /= 1) then
      call track_probe(p8_a, ptc_state, node1 = node_now, node2 = node_now%next)
      node_now => node_now%next
    endif

    eisf => ele_isf(ie)
    eisf%node(ii)%orb0 = p8_a%x
    eisf%node(ii)%s = node_now%s(2)  ! Offset from start of element.

    u = p8_a
    q_invar = u

    q2 = q_invar * q_y * q_invar**(-1)
    eisf%node(ii)%isf = q2%q

    ! This can be slow for nonlinear
    !  cdamap = u * cc_norm%n * u**(-1)
    !  q_1turn = cdamap
    !  eisf%node(ii)%m_1turn(1:6,1:6) = q_1turn%mat
  enddo
enddo

!

call kill (cdamap_1)
call kill (cdamap)
call kill (u)
call kill (p8_a)
call kill (cc_norm)

call init_all (ptc_private%base_state, ptc_private%taylor_order_ptc, 0)

end subroutine ptc_linear_isf_calc
