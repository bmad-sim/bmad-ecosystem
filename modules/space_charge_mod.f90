#include "CESR_platform.inc"

module space_charge_mod

use bmad_struct
use bmad_interface

type space_charge_struct
  type (coord_struct) closed_orb
  real(rp) kick_const
  real(rp) sig_x
  real(rp) sig_y
  real(rp) phi      ! Rotation angle to go from lab frame to rotated frame.
  real(rp) sin_phi
  real(rp) cos_phi
  real(rp) sig_z
endtype    

type space_charge_common_struct
  type (space_charge_struct), allocatable :: v(:)
end type

! sc_com%v(i) holds the parameters for the space_charge calculation
! at lattice element i.

type (space_charge_common_struct), save, private, target :: sc_com

contains

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!+
! Subroutine setup_space_charge_calc (calc_on, lattice, mode, closed_orb)
!
! Subroutine to initialize constants needed by the space charge 
! tracking routine track1_space_charge. This routine must be called if 
! the lattice or any of the other input parameters are changed.
!
! Modules needed:
!   use space_charge_mod
!
! Input:
!   calc_on    -- Logical: True turns on the space_charge calculation.
!   lattice    -- Ring_struct: Lattice for tracking.
!     %param%n_part -- Number of particles in a bunch
!   mode       -- modes_struct: Structure holding the beam info.
!     %a%emittance  -- a-mode emitance.
!     %b%emittance  -- b-mode emittance.
!     %sig_z        -- Real(rp): Bunch length.
!     %sigE_E       -- Real(rp): Sigma_E/E relative energy spread
!   closed_orb(0:) -- Coord_struct, optional: Closed orbit. If not present
!                       the closed orbit is taken to be zero. 
!-

subroutine setup_space_charge_calc (calc_on, lattice, mode, closed_orb)

  implicit none

  type (ring_struct), target :: lattice
  type (coord_struct), optional :: closed_orb(0:)
  type (modes_struct) mode
  type (space_charge_struct), pointer :: v
  type (ele_struct), pointer :: ele
  type (twiss_struct) a, b

  real(rp) c11, c12, c22, g, g2, xx_ave, xy_ave, yy_ave, phi
  real(rp) xx_rot_ave, yy_rot_ave, a_emit, b_emit, length, g4

  integer i, m, n_use
  logical calc_on

! Transfer some data to the common block for later use

  bmad_com%space_charge_on = calc_on

! Allocate space in the common block.

  n_use = lattice%n_ele_use
  m = 0
  if (allocated(sc_com%v)) m = size(sc_com%v)
  if (m /= n_use) then
    if (allocated(sc_com%v)) deallocate(sc_com%v)
    allocate (sc_com%v(n_use))
  endif

!------------------------------
! Loop over all lattice elements

  do i = 1, n_use

    ele => lattice%ele_(i)
    v => sc_com%v(i)

    v%sig_z = mode%sig_z

! Save the reference closed orbit.

    if (present(closed_orb)) then
      v%closed_orb = closed_orb(i)
    else
      v%closed_orb%vec = 0
    endif

! Due to coupling the beam ellipse may be rotated in the x-y plane.
! phi is this rotation angle.
! In the rotated frame the beam, by construction is decoupled.
! v%sig_x and v%sig_y are the x and y sigmas in the rotated frame.

    c11 = ele%c_mat(1,1); c12 = ele%c_mat(1,2); c22 = ele%c_mat(2,2)
    a = ele%x
    b = ele%y
    g = ele%gamma_c
    g2 = ele%gamma_c**2
    a_emit = mode%a%emittance
    b_emit = mode%b%emittance

    xx_ave = g2 * a_emit * a%beta + b_emit * (c11**2 * b%beta - &
                     2 * c11 * c12 * b%alpha + c12**2 * b%gamma) + &
                     (a%eta_lab * mode%sigE_E)**2

    xy_ave = g * (a_emit * (c22 * a%beta + c12 * a%alpha) + &
                  b_emit * (c11 * b%beta - c12 * b%alpha)) + &
                   a%eta_lab * b%eta_lab * mode%sigE_E**2

    yy_ave = g2 * a_emit * a%gamma + b_emit * (c22**2 * b%beta + &
                  2 * c22 * c12 * b%alpha + c12**2 * b%gamma) + &
                  (b%eta_lab * mode%sigE_E)**2

    phi = atan2(2 * xy_ave, xx_ave - yy_ave) / 2

    v%phi = phi
    v%cos_phi = cos(phi)
    v%sin_phi = sin(phi)

    xx_rot_ave = (xx_ave + yy_ave + cos(2*phi) * (xx_ave - yy_ave) + &
                                               2 * sin(2*phi) * xy_ave) / 2
    yy_rot_ave = xx_ave + yy_ave - xx_rot_ave

    v%sig_x = sqrt(xx_rot_ave)
    v%sig_y = sqrt(yy_rot_ave)

! The length over which the space charge acts is taken to be half 
! the length of the element + half the length of the next element.

    length = (ele%value(l$) + lattice%ele_(i+1)%value(l$)) / 2
    if (i == 1) length = ele%value(l$) + lattice%ele_(i+1)%value(l$) / 2
    if (i == n_use) length = ele%value(l$) / 2

! Calculate the kick constant.
! Taken from:
!   W. Decking, R. Brinkmann
!   "Space Charge Problems in the TESLA Damping Ring"
!   EPAC 2000, Vienna.
! The extra factor of 4pi comes from the normalization of 
!   the bbi_kick routine used in track1_space_charge.

    g4 = (ele%value(beam_energy$) / mass_of(lattice%param%particle))**4
    v%kick_const = length * r_e *  lattice%param%n_part / &
        (sqrt(twopi**3) * g4 * (v%sig_x + v%sig_y) * mode%sig_z)

  enddo

end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!+
! Subroutine track1_space_charge (start, ele, param, end)
!
! Routine to apply the space charge kick to a particle at the end of 
! an element. The routine setup_space_charge_calc must be called
! initially before any tracking is done. This routine assumes a Gaussian 
! bunch and is only valid with relativistic particles where the effect
! of the space charge is small.
!
! Modules needed:
!   use space_charge_mod
!
! Input:
!   start  -- Coord_struct: Starting position
!   ele    -- Ele_struct: Element tracked through.
!   param  -- Param_struct:
!
! Output:
!   end   -- Coord_struct: End position
!-

subroutine track1_space_charge (start, ele, param, end)

  implicit none

  type (coord_struct), intent(in)  :: start
  type (coord_struct), intent(out) :: end
  type (ele_struct),   intent(inout)  :: ele
  type (param_struct), intent(inout) :: param
  type (space_charge_struct), pointer :: v

  real(rp) x, y, x_rel, y_rel, kx_rot, ky_rot
  real(rp) kx, ky, kick_const

! Init

  end = start
  v => sc_com%v(ele%ix_ele)

! Rotate into frame where beam is not tilted.

  x = end%vec(1) - v%closed_orb%vec(1)
  y = end%vec(3) - v%closed_orb%vec(3)

  x_rel =  (x * v%cos_phi + y * v%sin_phi) / v%sig_x
  y_rel = (-x * v%sin_phi + y * v%cos_phi) / v%sig_y

  call bbi_kick (x_rel, y_rel, v%sig_y/v%sig_x, kx_rot, ky_rot)

! Transform the kick back to the lab coords and apply.

  kx = kx_rot * v%cos_phi - ky_rot * v%sin_phi
  ky = kx_rot * v%sin_phi - ky_rot * v%cos_phi

  kick_const = v%kick_const * exp(-(end%vec(5)/v%sig_z)**2) 

  end%vec(1) = end%vec(1) + kick_const * kx
  end%vec(3) = end%vec(3) + kick_const * ky

end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!+
! Subroutine make_mat6_space_charge (ele, param)
!
! Routine to add the space charge kick to the element transfer matrix.
! The routine setup_space_charge_calc must be called
! initially before any tracking is done. This routine assumes a Gaussian 
! bunch and is only valid with relativistic particles where the effect
! of the space charge is small.
!
! Modules needed:
!   use space_charge_mod
!
! Input:
!   start  -- Coord_struct: Starting position
!   ele    -- Ele_struct: Element tracked through.
!   param  -- Param_struct:
!
! Output:
!   end   -- Coord_struct: End position
!-

subroutine make_mat6_space_charge (ele, param)

  implicit none

  type (ele_struct),   intent(inout)  :: ele
  type (param_struct), intent(inout) :: param
  type (space_charge_struct), pointer :: v

  real(rp) kx_rot, ky_rot, kick_const, sc_kick_mat(6,6)

! Setup the space charge kick matrix and concatenate it with the 
! existing element transfer matrix.

  v => sc_com%v(ele%ix_ele)

  kx_rot = -4 * pi / v%sig_x
  ky_rot = -4 * pi / v%sig_y

  kick_const = v%kick_const 

  call mat_make_unit (sc_kick_mat)
  sc_kick_mat(2,1) = kick_const * kx_rot 
  sc_kick_mat(4,3) = kick_const * ky_rot 
  call tilt_mat6(sc_kick_mat, -v%phi)

  ele%mat6 = matmul (sc_kick_mat, ele%mat6)

end subroutine

end module
