!+
! Function momentum_compaction(branch) result (mom_comp)
!
! Routine to compute the momentum compaction factor = dL/L / dp/p.
!
! Note:
!   slip_factor      = mom_comp - 1 / gamma^2
!   gamma_transition = 1 / sqrt(mom_comp)
!
! Input:
!   branch      -- branch_struct, optional: Lattice branch to calculate on.
!
! Output:
!   mom_comp    -- real(rp): Momentum compaction.
!-

function momentum_compaction(branch) result (mom_comp)

use bmad_routine_interface, dummy => momentum_compaction

implicit none

type (branch_struct), target :: branch
type (ele_struct), pointer :: ele
type (ele_struct) ele0
type (coord_struct) orbit

real(rp) mom_comp
real(rp) m(6,6), mat4(4,4), eta_vec(4), length, pc, mc2
integer ie

! Form 1-turn matrix

length = branch%param%total_length

call mat_make_unit(m)
do ie = 1, branch%n_ele_track
  ele => branch%ele(ie)
  if (ele%key == rfcavity$) then
    orbit = ele%map_ref_orb_in
    call track_a_drift(orbit, ele%value(l$), m, .true.)
  else
    m = matmul(ele%mat6, m)
  endif
enddo

! Calc momentum compaction from: L = c*beta*t_ref - z

call mat_make_unit(mat4)
mat4 = mat4 - m(1:4,1:4)
call mat_inverse (mat4, mat4)
eta_vec = matmul(mat4, m(1:4,6))

mom_comp = -(sum(eta_vec * m(5,1:4)) + m(5,6)) / length  ! z contribution

orbit = branch%ele(1)%map_ref_orb_in
pc = (1+orbit%vec(6)) * orbit%p0c 
mc2 = mass_of(orbit%species)
mom_comp = mom_comp + mc2**2 / (pc**2 + mc2**2)  ! c*beta*t_ref contribution

end function
