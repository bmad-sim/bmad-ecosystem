!+
! Subroutine ele_at_s (lat, s, ix_ele, ix_branch, err_flag, s_eff)
!
! Subroutine to return the index of the element at position s.
! That is, ix_ele is choisen such that:
!     lat%ele(ix_ele-1)%s < s <= lat%ele(ix_ele)%s
! Exception: s = 0 --> ix_ele = 1.
!
! Note: For a circular lattice s is evaluated at the effective s which
! is modulo the lattice length:
!     s_eff = s - lat_length * floor(s/lat_length)
!
! Modules needed:
!   use bmad
!
! Input:
!   lat       -- lat_struct: Lattice of elements.
!   s         -- Real(rp): Longitudinal position.
!   ix_branch -- Integer, optional: Branch index. Default is 0.
!
! Output:
!   ix_ele    -- Integer: Index of element at s.
!   err_flag  -- logical, optional: Set True if s is out of bounds. False otherwise.
!   s_eff     -- Real(rp), optional: Effective s. Equal to s with a linear lattice.
!-

subroutine ele_at_s (lat, s, ix_ele, ix_branch, err_flag, s_eff)

use bmad, except_dummy => ele_at_s

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch

real(rp) s, ss, ll, ds_fudge
real(rp), optional :: s_eff

integer ix_ele, n1, n2, n3
integer, optional :: ix_branch

character(16), parameter :: r_name = 'ele_at_s'
logical, optional :: err_flag

!

if (present(err_flag)) err_flag = .true.

ds_fudge = bmad_com%significant_longitudinal_length
branch => lat%branch(integer_option(0, ix_branch))
ll = branch%param%total_length

if (branch%param%lattice_type == circular_lattice$) then
  ss = s - ll * floor((s-branch%ele(0)%s)/ll)
else
  ss = s
  if (s < branch%ele(0)%s - ds_fudge .or. s > branch%ele(branch%n_ele_track)%s + ds_fudge) then
    call out_io (s_fatal$, r_name, 'S POSITION OUT OF BOUNDS \f10.2\ ' , s)
    if (bmad_status%exit_on_error) call err_exit
    return
  endif
endif

n1 = 0
n3 = branch%n_ele_track

do

  if (n3 == n1 + 1) then
    ix_ele = n3
    exit
  endif

  n2 = (n1 + n3) / 2

  if (ss < branch%ele(n2)%s) then
    n3 = n2
  else
    n1 = n2
  endif

enddo

if (present(s_eff)) s_eff = ss
if (present(err_flag)) err_flag = .false.

end subroutine
