module high_energy_space_charge_mod

use bmad_interface

contains

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!+
! Subroutine setup_high_energy_space_charge_calc (calc_on, branch, n_part, mode, closed_orb)
!
! Routine to initialize constants needed by the ultra relativistic space charge 
! tracking routine track1_high_energy_space_charge. This setup routine must be called if 
! the lattice or any of the other input parameters are changed.
!
! Input:
!   calc_on    -- Logical: Turns on or off the space charge calculation.
!   branch     -- branch_struct: Lattice for tracking.
!   n_part     -- Real(rp): Number of actual particles in a bunch. Used to compute the bunch charge.
!   mode       -- normal_modes_struct: Structure holding the beam info.
!     %a%emittance  -- a-mode unnormalized emittance.
!     %b%emittance  -- b-mode unnormalized emittance.
!     %sig_z        -- Real(rp): Bunch length.
!     %sigE_E       -- Real(rp): Sigma_E/E relative energy spread
!   closed_orb(0:) -- Coord_struct, optional: Closed orbit. If not present
!                       the closed orbit is taken to be zero. 
!-

subroutine setup_high_energy_space_charge_calc (calc_on, branch, n_part, mode, closed_orb)

implicit none

type (branch_struct), target :: branch
type (coord_struct), optional :: closed_orb(0:)
type (normal_modes_struct) mode
type (high_energy_space_charge_struct), pointer :: sc
type (ele_struct), pointer :: ele
type (twiss_struct), pointer :: a, b
type (xy_disp_struct), pointer :: x, y

real(rp) c11, c12, c22, g, g2, xx_ave, xy_ave, yy_ave, phi, n_part
real(rp) xx_rot_ave, yy_rot_ave, a_emit, b_emit, length, g3, mc2, q2

integer i, m
logical calc_on

! Transfer some data to the common block for later use

bmad_com%high_energy_space_charge_on = calc_on

if (present(closed_orb)) then
  mc2 = mass_of(closed_orb(0)%species)
  q2 = charge_of(closed_orb(0)%species)
else
  mc2 = mass_of(branch%param%particle)
  q2 = charge_of(branch%param%particle)
endif

! Loop over all branch elements

do i = 1, branch%n_ele_track

  ele => branch%ele(i)
  if (.not. associated(ele%high_energy_space_charge)) allocate(ele%high_energy_space_charge)
  sc => ele%high_energy_space_charge

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
  a => ele%a
  b => ele%b
  x => ele%x
  y => ele%y
  g = ele%gamma_c
  g2 = ele%gamma_c**2
  a_emit = mode%a%emittance
  b_emit = mode%b%emittance

  xx_ave = g2 * a_emit * a%beta + b_emit * (c11**2 * b%beta - &
                   2 * c11 * c12 * b%alpha + c12**2 * b%gamma) + &
                   (x%eta * mode%sigE_E)**2

  xy_ave = g * (a_emit * (-c22 * a%beta - c12 * a%alpha) + &
                b_emit * ( c11 * b%beta - c12 * b%alpha)) + &
                 x%eta * y%eta * mode%sigE_E**2

  yy_ave = g2 * b_emit * b%beta + a_emit * (c22**2 * a%beta + &
                2 * c22 * c12 * a%alpha + c12**2 * a%gamma) + &
                (y%eta * mode%sigE_E)**2

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

  length = (ele%value(l$) + branch%ele(i+1)%value(l$)) / 2
  if (i == 1) length = ele%value(l$) + branch%ele(i+1)%value(l$) / 2
  if (i == branch%n_ele_track) length = ele%value(l$) / 2

! Calculate the kick constant.
! Taken from:
!   W. Decking, R. Brinkmann
!   "Space Charge Problems in the TESLA Damping Ring"
!   EPAC 2000, Vienna.
! The extra factor of 4pi comes from the normalization of 
!   the bbi_kick routine used in track1_space_charge.

  g3 = (ele%value(p0c$) / mc2)**3
  sc%kick_const = length * classical_radius_factor * n_part * q2 / &
                   (sqrt(twopi**3) * g3 * mc2 * (sc%sig_x + sc%sig_y) * mode%sig_z)

enddo

end subroutine setup_high_energy_space_charge_calc 

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!+
! Subroutine track1_high_energy_space_charge (ele, param, orbit)
!
! Routine to apply the ultra-relative space charge kick to a particle at the end of an element. 
! The routine setup_high_energy_space_charge_calc must be called initially before any tracking is done. 
! This routine assumes a Gaussian bunch and is only valid with relativistic particles where the 
! effect of the space charge is small.
!
! Input:
!   orbit   -- Coord_struct: Starting position
!   ele     -- Ele_struct: Element tracked through.
!   param   -- lat_param_struct:
!
! Output:
!   orbit   -- Coord_struct: End position
!-

subroutine track1_high_energy_space_charge (ele, param, orbit)

implicit none

type (coord_struct) :: orbit
type (ele_struct), target, intent(inout)  :: ele
type (lat_param_struct), intent(inout) :: param
type (high_energy_space_charge_struct), pointer :: sc

real(rp) x, y, x_rel, y_rel, kx, ky
real(rp) nk(2), dnk(2,2), kick_const

! Init

if (.not. associated(ele%high_energy_space_charge)) return

sc => ele%high_energy_space_charge

! Rotate into frame where beam is not tilted.

x = orbit%vec(1) - sc%closed_orb%vec(1)
y = orbit%vec(3) - sc%closed_orb%vec(3)

x_rel =  x * sc%cos_phi + y * sc%sin_phi 
y_rel = -x * sc%sin_phi + y * sc%cos_phi

call bbi_kick (x_rel, y_rel, [sc%sig_y, sc%sig_x], nk, dnk)

! Transform the kick back to the lab coords and apply.

kx = nk(1) * sc%cos_phi - nk(2) * sc%sin_phi
ky = nk(1) * sc%sin_phi + nk(2) * sc%cos_phi

kick_const = sc%kick_const * exp(-0.5 * (orbit%vec(5)/sc%sig_z)**2) / (1 + orbit%vec(6))**3

! The negative sign is due to the bbi kick assuming beams of opposite sign.

orbit%vec(2) = orbit%vec(2) - kick_const * kx
orbit%vec(4) = orbit%vec(4) - kick_const * ky

end subroutine track1_high_energy_space_charge 

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!+
! Subroutine make_mat6_high_energy_space_charge (ele, param)
!
! Routine to add the ultra relativistic space charge kick to the element transfer matrix.
! The routine setup_space_charge_calc must be called
! initially before any tracking is done. This routine assumes a Gaussian 
! bunch and is only valid with relativistic particles where the effect
! of the space charge is small.
!
! Input:
!   ele    -- Ele_struct: Element tracked through.
!   param  -- lat_param_struct:
!
! Output:
!   end   -- Coord_struct: End position
!-

subroutine make_mat6_high_energy_space_charge (ele, param)

implicit none

type (ele_struct), target, intent(inout)  :: ele
type (lat_param_struct), intent(inout) :: param
type (high_energy_space_charge_struct), pointer :: sc

real(rp) kx_rot, ky_rot, kick_const, sc_kick_mat(6,6)

! Setup the space charge kick matrix and concatenate it with the 
! existing element transfer matrix.

if (.not. associated(ele%high_energy_space_charge)) return
sc => ele%high_energy_space_charge

call mat_make_unit (sc_kick_mat)

kx_rot = 4 * pi * sc%kick_const / sc%sig_x
ky_rot = 4 * pi * sc%kick_const / sc%sig_y
sc_kick_mat(2,1) = kx_rot 
sc_kick_mat(4,3) = ky_rot 

call tilt_mat6(sc_kick_mat, -sc%phi)

ele%mat6 = matmul (sc_kick_mat, ele%mat6)

end subroutine make_mat6_high_energy_space_charge

end module
