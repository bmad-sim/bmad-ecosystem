!+
! Subroutine i2_i3_rad_int_calc (lat)
!
! Routine to calculate I2 and I3 radiation integrals be used in setting tolerances
! for radiation damping and stochastic matrices.
!
! Input:
!   lat       -- lat_struct: Lattice to ingegrate through.
!
! Output:
!   lat       -- lat_struct:
!     %branch(:)%I2_rad_int   -- I2 radiation integral branch-by-branch.
!     %branch(:)%I3_rad_int   -- I3 radiation integral branch-by-branch.
!-

subroutine i2_i3_rad_int_calc (lat)

use coord_mod, dummy => i2_i3_rad_int_calc

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

  p%i2_rad_int = 0
  p%i3_rad_int = 0

  do ie = 1, branch%n_ele_track
    ele => branch%ele(ie)
    select case (ele%key)
    case (sbend$)
      g = ele%value(g$)
      p%i2_rad_int = p%i2_rad_int + g**2 * ele%value(l$) 
      p%i3_rad_int = p%i3_rad_int + g**3 * ele%value(l$) 

    case (wiggler$, undulator$, em_field$)
      call init_coord(orbit, ele, upstream_end$)
      ds = 1e-3_rp
      do j = 0, int(ele%value(l$) / ds)
        call g_bending_strength_from_em_field (ele, p, j*ds,  orbit, .true., g_vec)
        g = norm2(g_vec)
        p%i2_rad_int = p%i2_rad_int + g**2 * ele%value(l$) 
        p%i3_rad_int = p%i3_rad_int + g**3 * ele%value(l$) 
      enddo
    end select
  enddo
enddo

end subroutine
