!+
! Subroutine radiation_integrals (ring, orb_, mode)
!
! Subroutine to calculate the synchrotron radiation integrals along with the
! emittance, and energy spread.
!
! Note: A negative emittance is possible and just means that the beam is
! unstable. That is, you have a negative damping partition number.
!
! Modules needed:
!   use bmad
!
! Input:
!   ring      -- Ring_struct: Ring to use. The calculation assumes that 
!                    the Twiss parameters have been calculated.
!   orb_(0:)  -- Coord_struct: Closed orbit.
!
! Output:
!   mode -- Modes_struct: Parameters for the ("horizontal like") a-mode,
!                              ("vertical like") b-mode, and the z-mode
!     %synch_int(1:3) -- Synchrotron integrals.
!     %sig_e          -- Sigma_E/E energy spread
!     %sig_z          -- Bunch Length
!     %e_loss         -- Energy loss in eV per turn
!     %a, %b, %z      -- Amode_struct: Substructure
!       %emittance      -- Emittance
!       %synch_int(4:5) -- Synchrotron integrals
!       %j_damp         -- Damping partition factor
!       %alpha_damp     -- Exponential damping coefficient per turn
!
! Notes:
!   1) %synch_int(1) = momentum_compaction * ring_length
!
!   2) There is a common block where the integrals for the individual elements
!      are saved. To access this common block a use statement is needed:
!         use rad_int_common
!      In the common block:
!         ric%i1_(:)  -- I1 integral for each element.
!         ric%i2_(:)  -- I2 integral for each element.
!         ric%i3_(:)  -- I3 integral for each element.
!         ric%i4a_(:) -- "A" mode I4 integral for each element.
!         ric%i4b_(:) -- "B" mode I4 integral for each element.
!         ric%i5a_(:) -- "A" mode I5 integral for each element.
!         ric%i5b_(:) -- "B" mode I5 integral for each element.
!-       

!$Id$
!$Log$
!Revision 1.7  2003/03/04 16:03:29  dcs
!VMS port
!
!Revision 1.6  2002/11/04 16:49:26  dcs
!*** empty log message ***
!
!Revision 1.5  2002/07/16 20:44:01  dcs
!*** empty log message ***
!
!Revision 1.4  2002/06/13 14:54:28  dcs
!Interfaced with FPP/PTC
!
!Revision 1.3  2002/02/23 20:32:22  dcs
!Double/Single Real toggle added
!
!Revision 1.2  2001/09/27 18:31:56  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine radiation_integrals (ring, orb_, mode)
                     
  use precision_def
  use nr
  use rad_int_common

  implicit none

  type (ring_struct), target :: ring
  type (coord_struct), target :: orb_(0:), start, end
  type (modes_struct) mode

  real(rdef), parameter :: c_gam = 4.425e-5, c_q = 3.84e-13
  real(rdef), save :: i1, i2, i3, i4a, i4b, i4z, i5a, i5b, m65, G_max, g3_ave
  real(rdef) theta, energy, gamma2_factor, energy_loss, arg, ll
  real(rdef) v(4,4), v_inv(4,4)

  integer i, ix, ir, key

!---------------------------------------------------------------------
! init

  m65 = 0

  ric%i1_ = 0;   ric%i2_ = 0;  ric%i3_ = 0
  ric%i4a_ = 0;  ric%i4b_ = 0
  ric%i5a_ = 0;  ric%i5b_ = 0

!---------------------------------------------------------------------
! loop over all elements

  ric%ring => ring

  do ir = 1, ring%n_ele_use     

    ric%ele => ring%ele_(ir)         
    ric%orb1 => orb_(ir)

    if (.not. ric%ele%is_on) cycle

    ric%ele0 => ring%ele_(ir-1)
    ric%orb0 => orb_(ir-1)

    ll = ric%ele%value(l$) 
    if (ll == 0) cycle

    key = ric%ele%key

    ric%g_x0 = -ric%ele%value(hkick$) / ll
    ric%g_y0 = -ric%ele%value(vkick$) / ll

    if (key == rfcavity$) m65 = m65 + ric%ele%mat6(6,5) 

! custom

    if (key == custom$) then
      call custom_radiation_integrals (ring, ir, orb_)
      cycle
    endif

! exact calculation involves runge_kutta tracking.

    if (ric%ele%exact_rad_int_calc) then

      ric%d_orb%vec = (/ 1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-5 /)

      rk_com%n_pts = 0
      call track1_runge_kutta (ric%orb0, ric%ele, ric%ring%param, ric%orb1)
      call transfer_rk_track (rk_com, ric%rk_track(0))

      do i = 1, 6
        start = ric%orb0
        start%vec(i) = start%vec(i) + ric%d_orb%vec(i)
        call track1_runge_kutta (start, ric%ele, ric%ring%param, end)
        call transfer_rk_track (rk_com, ric%rk_track(i))
      enddo

      ric%i1_(ir)  =   qromb_rad_int (eval_i1,  0.0_rdef, ll, i1, 'I1')
      ric%i2_(ir)  =   qromb_rad_int (eval_i2,  0.0_rdef, ll, i2, 'I2')
      ric%i3_(ir)  =   qromb_rad_int (eval_i3,  0.0_rdef, ll, i3, 'I3')
      ric%i4a_(ir) = ric%i4a_(ir) + &
                         qromb_rad_int (eval_i4a, 0.0_rdef, ll, i4a, 'I4A')
      ric%i4b_(ir) = ric%i4b_(ir) + &
                         qromb_rad_int (eval_i4b, 0.0_rdef, ll, i4b, 'I4B')
      ric%i5a_(ir) =   qromb_rad_int (eval_i5a, 0.0_rdef, ll, i5a, 'I5A')
      ric%i5b_(ir) =   qromb_rad_int (eval_i5b, 0.0_rdef, ll, i5b, 'I5B')

      cycle

    endif

! for an old style wiggler we make the approximation that the variation of G is
! fast compaired to the variation in eta.

    if (key == wiggler$ .and. ric%ele%sub_key == periodic_type$) then
      G_max = sqrt(2*abs(ric%ele%value(k1$)))       ! 1/rho at max B
      ric%i2_(ir) = ll * G_max**2 / 2
      g3_ave = 4 * G_max**3 / (3 * pi)
      ric%i3_(ir) = ll * g3_ave
      ric%g_x0 = g3_ave**(1.0/3)
      ric%g_y0 = 0
      ric%k1 = 0
      ric%s1 = 0
      ric%i5a_(ir) = qromb_rad_int (eval_i5a, 0.0_rdef, ll, i5a, 'I5A')
      ric%i5b_(ir) = qromb_rad_int (eval_i5b, 0.0_rdef, ll, i5b, 'I5B')
      cycle
    endif
   

    if (key == wiggler$ .and. ric%ele%sub_key == map_type$) then

!      call track1_runge_kutta (ric%orb0, ric%ele, ric%ring%param, end)
!      call transfer_rk_track (rk_com, ric%rk_track(0))

      call make_v_mats (ric%ele0, v, v_inv)
      ric%eta_a0 = &      
          matmul(v, (/ ric%ele0%x%eta, ric%ele0%x%etap, 0.0_rdef, 0.0_rdef /))
      ric%eta_b0 = &
          matmul(v, (/ 0.0_rdef, 0.0_rdef, ric%ele0%y%eta, ric%ele0%y%etap /))

      call make_v_mats (ric%ele, v, v_inv)
      ric%eta_a1 = &      
          matmul(v, (/ ric%ele%x%eta, ric%ele%x%etap, 0.0_rdef, 0.0_rdef /))
      ric%eta_b1 = &
          matmul(v, (/ 0.0_rdef, 0.0_rdef, ric%ele%y%eta, ric%ele%y%etap /))

      ric%i1_(ir)  =   qromb_rad_int (eval_i1,  0.0_rdef, ll, i1, 'I1')
      ric%i2_(ir)  =   qromb_rad_int (eval_i2,  0.0_rdef, ll, i2, 'I2')
      ric%i3_(ir)  =   qromb_rad_int (eval_i3,  0.0_rdef, ll, i3, 'I3')
      ric%i4a_(ir) = ric%i4a_(ir) + &
                         qromb_rad_int (eval_i4a, 0.0_rdef, ll, i4a, 'I4A')
      ric%i4b_(ir) = ric%i4b_(ir) + &
                         qromb_rad_int (eval_i4b, 0.0_rdef, ll, i4b, 'I4B')
      ric%i5a_(ir) =   qromb_rad_int (eval_i5a, 0.0_rdef, ll, i5a, 'I5A')
      ric%i5b_(ir) =   qromb_rad_int (eval_i5b, 0.0_rdef, ll, i5b, 'I5B')

      cycle

    endif

!
 
    if (ric%g_x0 == 0 .and. ric%g_y0 == 0 .and. &
          key /= quadrupole$ .and. key /= sol_quad$ .and. key /= sbend$) cycle

    if (key == sbend$) then
      theta = ric%ele%value(tilt$) + ric%ele%value(roll$)
      ric%g_x0 = ric%g_x0 + cos(theta) * ric%ele%value(g$) 
      ric%g_y0 = ric%g_y0 - sin(theta) * ric%ele%value(g$) 
    endif

    ric%g2 = ric%g_x0**2 + ric%g_y0**2
    ric%g = sqrt(ric%g2)

    ric%i2_(ir)  = ric%g2 * ll
    ric%i3_(ir)  = ric%g2 * ric%g * ll

    if (key == quadrupole$ .or. key == sol_quad$) then
      theta = ric%ele%value(tilt$)
      ric%k1 = ric%ele%value(k1$) * cos(2*theta)
      ric%s1 = ric%ele%value(k1$) * sin(2*theta)
    elseif (key == sbend$) then
      theta = ric%ele%value(tilt$) + ric%ele%value(roll$)
      ric%k1 = ric%ele%value(k1$) * cos(2*theta)
      ric%s1 = ric%ele%value(k1$) * sin(2*theta)
    else
      ric%k1 = 0
      ric%s1 = 0
    endif

! edge effects for a bend. In this case we ignore any rolls.

    if (key == sbend$) then
      call propagate_part_way(0.0_rdef)
      ric%i4a_(ir) = -ric%eta_a(1) * ric%g2 * tan(ric%ele%value(e1$))
      ric%i4b_(ir) = -ric%eta_b(1) * ric%g2 * tan(ric%ele%value(e1$))
      call propagate_part_way(ll)
      ric%i4a_(ir) = ric%i4a_(ir) - &
                           ric%eta_a(1) * ric%g2 * tan(ric%ele%value(e2$))
      ric%i4b_(ir) = ric%i4a_(ir) - &
                           ric%eta_b(1) * ric%g2 * tan(ric%ele%value(e2$))
    endif

! integrate 

    ric%i1_(ir)  =   qromb_rad_int (eval_i1,  0.0_rdef, ll, i1, 'I1')
    ric%i4a_(ir) = ric%i4a_(ir) + &
                         qromb_rad_int (eval_i4a, 0.0_rdef, ll, i4a, 'I4A')
    ric%i4b_(ir) = ric%i4b_(ir) + &
                         qromb_rad_int (eval_i4b, 0.0_rdef, ll, i4b, 'I4B')
    ric%i5a_(ir) =   qromb_rad_int (eval_i5a, 0.0_rdef, ll, i5a, 'I5A')
    ric%i5b_(ir) =   qromb_rad_int (eval_i5b, 0.0_rdef, ll, i5b, 'I5B')

  enddo

!---------------------------------------------------------------------
! now put everything together

  i1   = sum(ric%i1_(1:ring%n_ele_use))
  i2   = sum(ric%i2_(1:ring%n_ele_use))
  i3   = sum(ric%i3_(1:ring%n_ele_use))
  i4a  = sum(ric%i4a_(1:ring%n_ele_use))
  i4b  = sum(ric%i4b_(1:ring%n_ele_use))
  i5a  = sum(ric%i5a_(1:ring%n_ele_use))
  i5b  = sum(ric%i5b_(1:ring%n_ele_use))

  i4z = i4a + i4b

  energy = ring%param%beam_energy
  gamma2_factor = (energy * 1956.95e-9)**2
  energy_loss = 1e9 * c_gam * (1e-9 * energy)**4 * i2 / pi

  mode%synch_int(1) = i1
  mode%synch_int(2) = i2
  mode%synch_int(3) = i3

  mode%a%synch_int(4) = i4a
  mode%b%synch_int(4) = i4b
  mode%z%synch_int(4) = i4z

  mode%a%synch_int(5) = i5a
  mode%b%synch_int(5) = i5b

  if (i2 /= 0) then

    mode%a%emittance = c_q * gamma2_factor * i5a / (i2 - i4a)
    mode%b%emittance = c_q * gamma2_factor * i5b / (i2 - i4b)

    mode%a%j_damp = 1 - i4a / i2
    mode%b%j_damp = 1 - i4b / i2
    mode%z%j_damp = 2 + i4z / i2

    arg = (c_q * i3 * gamma2_factor / (2*i2 + i4z))
    mode%sig_e = sqrt(max(0.0_rdef, arg))

  endif

  mode%a%alpha_damp = energy_loss * mode%a%j_damp / energy
  mode%b%alpha_damp = energy_loss * mode%b%j_damp / energy
  mode%z%alpha_damp = energy_loss * mode%z%j_damp / energy

  mode%e_loss = energy_loss

  if(abs(m65) > 0. ) then
    mode%sig_z = sqrt( mode%synch_int(1)/abs(m65) ) * mode%sig_e
  else
    mode%sig_z = 0.
  endif

end subroutine
