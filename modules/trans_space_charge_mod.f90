#include "CESR_platform.inc"

module trans_space_charge_mod

use bmad_struct
use bmad_interface
use make_mat6_mod

contains

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!+
! Subroutine setup_trans_space_charge_calc (calc_on, lattice, mode, closed_orb)
!
! Subroutine to initialize constants needed by the transverse space charge 
! tracking routine track1_space_charge. This routine must be called if 
! the lattice or any of the other input parameters are changed.
!
! Modules needed:
!   use trans_space_charge_mod
!
! Input:
!   calc_on    -- Logical: True turns on the space charge calculation.
!   lattice    -- lat_struct: Lattice for tracking.
!     %param%n_part -- Number of particles in a bunch
!   mode       -- normal_modes_struct: Structure holding the beam info.
!     %a%emittance  -- a-mode emitance.
!     %b%emittance  -- b-mode emittance.
!     %sig_z        -- Real(rp): Bunch length.
!     %sigE_E       -- Real(rp): Sigma_E/E relative energy spread
!   closed_orb(0:) -- Coord_struct, optional: Closed orbit. If not present
!                       the closed orbit is taken to be zero. 
!-

subroutine setup_trans_space_charge_calc (calc_on, lattice, mode, closed_orb)

  implicit none

  type (lat_struct), target :: lattice
  type (coord_struct), optional :: closed_orb(0:)
  type (normal_modes_struct) mode
  type (trans_space_charge_struct), pointer :: sc
  type (ele_struct), pointer :: ele
  type (twiss_struct) a, b

  real(rp) c11, c12, c22, g, g2, xx_ave, xy_ave, yy_ave, phi
  real(rp) xx_rot_ave, yy_rot_ave, a_emit, b_emit, length, g3

  integer i, m
  logical calc_on

! Transfer some data to the common block for later use

  bmad_com%trans_space_charge_on = calc_on

! Allocate space in the common block.

!------------------------------
! Loop over all lattice elements

  do i = 1, lattice%n_ele_track

    ele => lattice%ele(i)
    if (.not. associated(ele%trans_sc)) allocate(ele%trans_sc)
    sc => ele%trans_sc

    sc%sig_z = mode%sig_z

! Save the reference closed orbit.

    if (present(closed_orb)) then
      sc%closed_orb = closed_orb(i)
    else
      sc%closed_orb%vec = 0
    endif

! Due to coupling the beam ellipse may be rotated in the x-y plane.
! phi is this rotation angle.
! In the rotated frame the beam, by construction is decoupled.
! sc%sig_x and sc%sig_y are the x and y sigmas in the rotated frame.

    c11 = ele%c_mat(1,1); c12 = ele%c_mat(1,2); c22 = ele%c_mat(2,2)
    a = ele%a
    b = ele%b
    g = ele%gamma_c
    g2 = ele%gamma_c**2
    a_emit = mode%a%emittance
    b_emit = mode%b%emittance

    xx_ave = g2 * a_emit * a%beta + b_emit * (c11**2 * b%beta - &
                     2 * c11 * c12 * b%alpha + c12**2 * b%gamma) + &
                     (a%eta_lab * mode%sigE_E)**2

    xy_ave = g * (a_emit * (-c22 * a%beta - c12 * a%alpha) + &
                  b_emit * ( c11 * b%beta - c12 * b%alpha)) + &
                   a%eta_lab * b%eta_lab * mode%sigE_E**2

    yy_ave = g2 * b_emit * b%beta + a_emit * (c22**2 * a%beta + &
                  2 * c22 * c12 * a%alpha + c12**2 * a%gamma) + &
                  (b%eta_lab * mode%sigE_E)**2

    phi = atan2(2 * xy_ave, xx_ave - yy_ave) / 2

    sc%phi = phi
    sc%cos_phi = cos(phi)
    sc%sin_phi = sin(phi)

    xx_rot_ave = (xx_ave + yy_ave + cos(2*phi) * (xx_ave - yy_ave) + &
                                               2 * sin(2*phi) * xy_ave) / 2
    yy_rot_ave = xx_ave + yy_ave - xx_rot_ave

    sc%sig_x = sqrt(xx_rot_ave)
    sc%sig_y = sqrt(yy_rot_ave)

! The length over which the space charge acts is taken to be half 
! the length of the element + half the length of the next element.

    length = (ele%value(l$) + lattice%ele(i+1)%value(l$)) / 2
    if (i == 1) length = ele%value(l$) + lattice%ele(i+1)%value(l$) / 2
    if (i == lattice%n_ele_track) length = ele%value(l$) / 2

! Calculate the kick constant.
! Taken from:
!   W. Decking, R. Brinkmann
!   "Space Charge Problems in the TESLA Damping Ring"
!   EPAC 2000, Vienna.
! The extra factor of 4pi comes from the normalization of 
!   the bbi_kick routine used in track1_trans_space_charge.

    g3 = (ele%value(p0c$) / mass_of(lattice%param%particle))**3
    sc%kick_const = length * r_e *  lattice%param%n_part / &
        (sqrt(twopi**3) * g3 * (sc%sig_x + sc%sig_y) * mode%sig_z)

  enddo

end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!+
! Subroutine track1_trans_space_charge (start, ele, param, end)
!
! Routine to apply the space charge kick to a particle at the end of 
! an element. The routine setup_trans_space_charge_calc must be called
! initially before any tracking is done. This routine assumes a Gaussian 
! bunch and is only valid with relativistic particles where the effect
! of the space charge is small.
!
! Modules needed:
!   use trans_space_charge_mod
!
! Input:
!   start  -- Coord_struct: Starting position
!   ele    -- Ele_struct: Element tracked through.
!   param  -- lat_param_struct:
!
! Output:
!   end   -- Coord_struct: End position
!-

subroutine track1_trans_space_charge (start, ele, param, end)

  implicit none

  type (coord_struct), intent(in)  :: start
  type (coord_struct), intent(out) :: end
  type (ele_struct), target, intent(inout)  :: ele
  type (lat_param_struct), intent(inout) :: param
  type (trans_space_charge_struct), pointer :: sc

  real(rp) x, y, x_rel, y_rel, kx_rot, ky_rot
  real(rp) kx, ky, kick_const

! Init

  end = start
  if (.not. associated(ele%trans_sc)) return

  sc => ele%trans_sc

! Rotate into frame where beam is not tilted.

  x = end%vec(1) - sc%closed_orb%vec(1)
  y = end%vec(3) - sc%closed_orb%vec(3)

  x_rel =  (x * sc%cos_phi + y * sc%sin_phi) / sc%sig_x
  y_rel = (-x * sc%sin_phi + y * sc%cos_phi) / sc%sig_y

  call bbi_kick (x_rel, y_rel, sc%sig_y/sc%sig_x, kx_rot, ky_rot)

! Transform the kick back to the lab coords and apply.

  kx = kx_rot * sc%cos_phi - ky_rot * sc%sin_phi
  ky = kx_rot * sc%sin_phi + ky_rot * sc%cos_phi

  kick_const = sc%kick_const * exp(-0.5 * (end%vec(5)/sc%sig_z)**2) 

! The negative sign is due to the bbi kick assuming beams of opposite sign.

  end%vec(2) = end%vec(2) - kick_const * kx
  end%vec(4) = end%vec(4) - kick_const * ky

end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!+
! Subroutine make_mat6_trans_space_charge (ele, param)
!
! Routine to add the space charge kick to the element transfer matrix.
! The routine setup_trans_space_charge_calc must be called
! initially before any tracking is done. This routine assumes a Gaussian 
! bunch and is only valid with relativistic particles where the effect
! of the space charge is small.
!
! Modules needed:
!   use trans_space_charge_mod
!
! Input:
!   start  -- Coord_struct: Starting position
!   ele    -- Ele_struct: Element tracked through.
!   param  -- lat_param_struct:
!
! Output:
!   end   -- Coord_struct: End position
!-

subroutine make_mat6_trans_space_charge (ele, param)

  implicit none

  type (ele_struct), target, intent(inout)  :: ele
  type (lat_param_struct), intent(inout) :: param
  type (trans_space_charge_struct), pointer :: sc

  real(rp) kx_rot, ky_rot, kick_const, sc_kick_mat(6,6)

! Setup the space charge kick matrix and concatenate it with the 
! existing element transfer matrix.

  if (.not. associated(ele%trans_sc)) return
  sc => ele%trans_sc

  call mat_make_unit (sc_kick_mat)

  kx_rot = 4 * pi * sc%kick_const / sc%sig_x
  ky_rot = 4 * pi * sc%kick_const / sc%sig_y
  sc_kick_mat(2,1) = kx_rot 
  sc_kick_mat(4,3) = ky_rot 

  call tilt_mat6(sc_kick_mat, -sc%phi)

  ele%mat6 = matmul (sc_kick_mat, ele%mat6)

end subroutine

end module
