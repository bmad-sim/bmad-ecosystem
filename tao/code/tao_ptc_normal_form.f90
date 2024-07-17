!+
! Subroutine tao_ptc_normal_form (do_calc, tao_lat, ix_branch, rf_on)
!
! Routine to calculate normal form quantities.
!
! Input:
!   do_calc     -- logical: Set True to do the calculation. 
!   tao_lat     -- tao_lattice_struct: Lattice to work on.
!   ix_branch   -- integer: Branch of lattice to work on.
!   rf_on       -- integer, optional: RF state for calculation. yes$, no$, or maybe$ (default)
!                   maybe$ means that RF state in branch is used.
!-

subroutine tao_ptc_normal_form (do_calc, tao_lat, ix_branch, rf_on)

use tao_interface, dummy => tao_ptc_normal_form
use ptc_layout_mod
use pointer_lattice

implicit none

type (tao_lattice_struct), target :: tao_lat
type (tao_lattice_branch_struct), pointer :: tao_branch
type (ptc_normal_form_struct), pointer :: ptc_nf
type (branch_struct), pointer :: branch
type (c_taylor) beta, c_tay(6)
type (c_damap) c_da

real(rp) mm, beta0
integer ix_branch
integer, optional :: rf_on
logical do_calc, this_rf_on

!

branch => tao_lat%lat%branch(ix_branch)
tao_branch => tao_lat%tao_branch(ix_branch)
ptc_nf => tao_branch%ptc_normal_form

if (ptc_nf%valid_map) then
  call kill(ptc_nf%one_turn_map)
  call kill(ptc_nf%normal_form)
  call kill(ptc_nf%phase)
  call kill(ptc_nf%spin_tune)
  call kill(ptc_nf%path_length)
  call kill(ptc_nf%isf)
  ptc_nf%valid_map = .false.
endif

nullify(tao_branch%bmad_normal_form%ele_origin)

if (.not. do_calc .or. branch%param%geometry == open$) return

!

call set_ptc_verbose(.false.)

if (.not. associated(ptc_nf%ele_origin)) ptc_nf%ele_origin => branch%ele(0)
if (.not. associated(branch%ptc%m_t_layout)) call lat_to_ptc_layout (tao_lat%lat)

call ptc_one_turn_map_at_ele (ptc_nf%ele_origin, ptc_nf%orb0, ptc_nf%one_turn_map, ptc_nf%state, pz = 0.0_rp, rf_on = rf_on)

call alloc(ptc_nf%normal_form)
call alloc(ptc_nf%phase)
call alloc(ptc_nf%spin_tune)
call alloc(ptc_nf%path_length)
call alloc(ptc_nf%isf)
call alloc(beta)
call alloc(c_da)
call alloc(c_tay)

ptc_nf%valid_map = .true.

call ptc_map_to_normal_form (ptc_nf%one_turn_map, ptc_nf%state%nocavity, &
                                      ptc_nf%normal_form, ptc_nf%phase, ptc_nf%spin_tune)

mm = mass_of(branch%param%particle) / branch%ele(0)%value(p0c$)
beta0 = branch%ele(0)%value(p0c$) / branch%ele(0)%value(E_tot$)
 
c_da = 1
beta = c_da%v(6)
beta = (1 + beta) / sqrt((1+beta)**2 + mm**2)
ptc_nf%path_length = (branch%param%total_length - ptc_nf%phase(3)) * beta / beta0

c_da%q%x = 0.0d0
c_da%q%x(2) = 1.0d0
c_da = ptc_nf%normal_form%atot*c_da*ptc_nf%normal_form%atot**(-1)
ptc_nf%isf = c_da%q


! this_rf_on = rf_is_on(branch)
! call normal_form_taylors(normal_form%m, this_rf_on, dhdj = normal_form%dhdj, &
!                                  A = normal_form%A, A_inverse = normal_form%A_inv)  ! Get A, A_inv, dhdj
! call normal_form_complex_taylors(normal_form%m, this_rf_on, F = normal_form%F, L = normal_form%L)  ! Get complex L and F

call set_ptc_verbose(.true.)
call kill (beta)
call kill (c_da)
call kill (c_tay)

end subroutine
