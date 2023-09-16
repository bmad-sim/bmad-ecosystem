!+
! Subroutine g_integrals_calc (lat)
!
! Routine to calculate approximate g (bending force) integrals to be used in setting tolerances
! for radiation damping and stochastic matrices.
!
! Input:
!   lat       -- lat_struct: Lattice to integrate through.
!
! Output:
!   lat       -- lat_struct:
!     %branch(:)%g1_integral   -- Integral of |g| branch-by-branch.
!     %branch(:)%g2_integral   -- Integral of g^2 branch-by-branch. Same as I_2 radiation integral.
!     %branch(:)%g3_integral   -- Integral of g^3 branch-by-branch. Same as I_3 radiation integral.
!-

subroutine g_integrals_calc (lat)

use coord_mod, dummy => g_integrals_calc

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele
type (lat_param_struct), pointer :: p
type (em_field_struct) field
type (coord_struct) orbit

real(rp) ds, g_vec(3), g
integer ib, ie, j

!

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  p => branch%param

  p%g1_integral = 0
  p%g2_integral = 0
  p%g3_integral = 0

  do ie = 1, branch%n_ele_track
    ele => branch%ele(ie)
    select case (ele%key)
    case (sbend$, rf_bend$)
      g = ele%value(g$)
      p%g1_integral = p%g1_integral + g    * ele%value(l$) 
      p%g2_integral = p%g2_integral + g**2 * ele%value(l$) 
      p%g3_integral = p%g3_integral + g**3 * ele%value(l$) 

    case (wiggler$, undulator$, em_field$)
      call init_coord(orbit, ele, upstream_end$)
      ds = 1e-3_rp
      do j = 0, int(ele%value(l$) / ds)
        call g_bending_strength_from_em_field (ele, p, j*ds,  orbit, .true., g_vec)
        g = norm2(g_vec)
        p%g1_integral = p%g2_integral + g    * ele%value(l$) 
        p%g2_integral = p%g2_integral + g**2 * ele%value(l$) 
        p%g3_integral = p%g3_integral + g**3 * ele%value(l$) 
      enddo
    end select
  enddo
enddo

end subroutine
