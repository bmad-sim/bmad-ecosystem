!+
! Subroutine lat_make_mat6_hook (finished, lat, ix_ele, ref_orb, ix_branch, err_flag)
!
! Prototype routine that can be customized for constructing element transfer matrices.
!
! Subroutine to make the first order transfer map for an element or elements:
!   r_out = M * r_in + vec0
! M is the 6x6 linear transfer matrix (Jacobian) about the 
! reference orbit ref_orb.
!
! If the element lat%ele(ix_ele) is a lord element then the martices of 
! all the slave elements will be recomputed.
!
! Input:
!   lat         -- lat_struct: Lat containing the elements.
!   ix_ele      -- Integer, optional: Index of the element. If not present
!                    or negative, the matrices for all elements will be calculated.
!   ref_orb(0:) -- Coord_struct, optional: Coordinates of the reference orbit
!                   around which the matrix is calculated. If not present 
!                   then the referemce is taken to be the origin.
!   ix_branch   -- Integer, optional: Branch index. Default is 0 (main lattice).
!                   -1 => All branches/all elements (ref_orb & ix_ele will be ignored).
!
! Output:
!   finished    -- logical: When set True, the standard lat_make_mat6 code will not be called.
!   lat         -- lat_struct:
!     ele(:)%mat6  -- Real(rp): 1st order (Jacobian) 6x6 transfer matrix.
!     ele(:)%vec0  -- Real(rp): 0th order transfer vector.
!   err_flag    -- Logical, optional: True if there is an error. False otherwise.
!-

recursive subroutine lat_make_mat6_hook (finished, lat, ix_ele, ref_orb, ix_branch, err_flag)

use da_program_mod

implicit none
                                       
type (lat_struct), target :: lat
type (coord_struct), optional :: ref_orb(0:)
type (coord_struct), allocatable :: orbit(:)
type (branch_struct), pointer :: branch
type (ltt_params_struct), pointer :: lttp
type (ltt_com_struct), pointer :: ltt_com

real(rp), pointer :: mat6(:,:), vec0(:)

integer, optional :: ix_ele, ix_branch
integer state, n

logical, optional :: err_flag
logical finished

character(*), parameter :: r_name = 'lat_make_mat6_hook'

!

branch => lat%branch(integer_option(0, ix_branch))
finished = .false.

branch => lat%branch(ix_branch)
lttp => ltt_params_global
ltt_com => ltt_com_global

! Normal Bmad tracking?

finished = .false.
if (lttp%tracking_method == 'OLD' .or. ltt_com%track_bypass) return

! This is to get around the problem that the other tracking_methods will not calculate ref_orb so
! approximate the ref_orb by normal tracking 

call reallocate_coord (orbit, lat, ix_branch)
orbit(0) = ref_orb(0)
call track_all(lat, orbit, ix_branch, state, err_flag)
n = branch%n_ele_max
ref_orb(1:n) = orbit(1:n)

end subroutine
