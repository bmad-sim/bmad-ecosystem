!+
! Function attribute_name (ele, ix_att)
!
! Function to return the name of an attribute for a particular type of 
! BMAD element. 
!
! Modules Needed:
!   use bmad
!
! Input:
!   ele    -- Ele_struct: Integer: Key name of element type (e.g. SBEND$, etc.)
!   ix_att -- Integer: Index of attribute (e.g. k1$)
!
! Output:
!   attribute_name -- Character(16): Name of attribute.
!
! Example:
!   ele%key = sbend$
!   name = attribute_name (ele, K1$)
! Result:
!   name -> 'k1'
!-

#include "CESR_platform.inc"

function attribute_name (ele, ix_att) result (at_name)

  use bmad_struct
  use bmad_interface, except => attribute_name

  implicit none

  type (ele_struct) ele

  integer i, key, ix_att

  character(16) attrib_array(n_key, n_attrib_special_maxx)
  character(16) at_name

  logical, save :: init_needed = .true.

!---------------------------------------------------------------------
! Init

  if (init_needed) then

    attrib_array = null_name

    do i = 1, n_key
                                    
      if (i == def_beam$)      cycle
      if (i == def_parameter$) cycle

      attrib_array(i, type$)     = 'TYPE'
      attrib_array(i, alias$)    = 'ALIAS'
      attrib_array(i, descrip$)  = 'DESCRIP'

      if (i == init_ele$) cycle
      if (i == group$)    cycle
      if (i == overlay$)  cycle
      if (i == i_beam$)   cycle

      attrib_array(i, x_limit$)     = 'X_LIMIT'
      attrib_array(i, y_limit$)     = 'Y_LIMIT'
      attrib_array(i, aperture$)    = 'APERTURE'
      attrib_array(i, aperture_at$) = 'APERTURE_AT'
      attrib_array(i, beam_energy$) = 'BEAM_ENERGY'
      attrib_array(i, p0c$)         = 'P0C'
      attrib_array(i, mat6_calc_method$)  = 'MAT6_CALC_METHOD'
      attrib_array(i, tracking_method$)   = 'TRACKING_METHOD'

      if (i == marker$)       cycle

      attrib_array(i, is_on$)       = 'IS_ON'

      if (i == match$)        cycle
      if (i == patch$)        cycle
      if (i == beambeam$)     cycle
      if (i == hom$)          cycle
      if (i == multipole$)    cycle 
      if (i == ab_multipole$) cycle
      if (i == taylor$)       cycle

      attrib_array(i, integration_ord$)   = 'INTEGRATION_ORD'
      attrib_array(i, num_steps$)         = 'NUM_STEPS'
      attrib_array(i, symplectify$)       = 'SYMPLECTIFY'
      attrib_array(i, rel_tol$)           = 'REL_TOL'
      attrib_array(i, abs_tol$)           = 'ABS_TOL'

      if (i == hkicker$)      cycle
      if (i == vkicker$)      cycle
      if (i == custom$)       cycle
      if (i == drift$)        cycle
      if (i == monitor$)      cycle
      if (i == instrument$)   cycle
      if (i == kicker$)       cycle

      attrib_array(i, x_offset$) = 'X_OFFSET'
      attrib_array(i, y_offset$) = 'Y_OFFSET'
      attrib_array(i, s_offset$) = 'S_OFFSET'

      if (i == rcollimator$) cycle
      if (i == ecollimator$) cycle

      attrib_array(i, hkick$)  = 'HKICK'
      attrib_array(i, vkick$)  = 'VKICK'

      attrib_array(i, x_pitch$)   = 'X_PITCH'
      attrib_array(i, y_pitch$)   = 'Y_PITCH'

    enddo

!

    do i = 1, n_key
      select case (i)
      case (elseparator$, kicker$, octupole$, quadrupole$, sbend$, rbend$, &
             sextupole$, solenoid$, sol_quad$, ab_multipole$, wiggler$)
        attrib_array(i, a0$:a20$) = (/ 'A0 ', &
                                       'A1 ', 'A2 ', 'A3 ', 'A4 ', 'A5 ', & 
                                       'A6 ', 'A7 ', 'A8 ', 'A9 ', 'A10', &
                                       'A11', 'A12', 'A13', 'A14', 'A15', &
                                       'A16', 'A17', 'A18', 'A19', 'A20' /)
        attrib_array(i, b0$:b20$) = (/ 'B0 ', &
                                       'B1 ', 'B2 ', 'B3 ', 'B4 ', 'B5 ', & 
                                       'B6 ', 'B7 ', 'B8 ', 'B9 ', 'B10', &
                                       'B11', 'B12', 'B13', 'B14', 'B15', &
                                       'B16', 'B17', 'B18', 'B19', 'B20' /)
      end select
    enddo

    attrib_array(def_parameter$, beam_energy$)  = 'BEAM_ENERGY'
    attrib_array(def_parameter$, lattice_type$) = 'LATTICE_TYPE'
    attrib_array(def_parameter$, taylor_order$) = 'TAYLOR_ORDER'
    attrib_array(def_parameter$, ran_seed$)     = 'RAN_SEED'

    attrib_array(taylor$, symplectify$) = 'SYMPLECTIFY'

    attrib_array(match$, l$)             = 'L'
    attrib_array(match$, beta_x0$)       = 'BETA_X0'
    attrib_array(match$, alpha_x0$)      = 'ALPHA_X0'
    attrib_array(match$, beta_y0$)       = 'BETA_Y0'
    attrib_array(match$, alpha_y0$)      = 'ALPHA_Y0'
    attrib_array(match$, beta_x1$)       = 'BETA_X1'
    attrib_array(match$, alpha_x1$)      = 'ALPHA_X1'
    attrib_array(match$, beta_y1$)       = 'BETA_Y1'
    attrib_array(match$, alpha_y1$)      = 'ALPHA_Y1'
    attrib_array(match$, dphi_x$)        = 'DPHI_X'
    attrib_array(match$, dphi_y$)        = 'DPHI_Y'
    attrib_array(match$, eta_x0$)        = 'ETA_X0'
    attrib_array(match$, etap_x0$)       = 'ETAP_X0'
    attrib_array(match$, eta_y0$)        = 'ETA_Y0'
    attrib_array(match$, etap_y0$)       = 'ETAP_Y0'
    attrib_array(match$, eta_x1$)        = 'ETA_X1'
    attrib_array(match$, etap_x1$)       = 'ETAP_X1'
    attrib_array(match$, eta_y1$)        = 'ETA_Y1'
    attrib_array(match$, etap_y1$)       = 'ETAP_Y1'

    attrib_array(i_beam$, x_offset$)     = 'X_OFFSET'
    attrib_array(i_beam$, y_offset$)     = 'Y_OFFSET'
    attrib_array(i_beam$, s_offset$)     = 'S_OFFSET'
    attrib_array(i_beam$, x_pitch$)      = 'X_PITCH'
    attrib_array(i_beam$, y_pitch$)      = 'Y_PITCH'
    attrib_array(i_beam$, s_center$)     = 'S_CENTER'
    attrib_array(i_beam$, tilt$)         = 'TILT'

    attrib_array(lcavity$, l$)             = 'L'
    attrib_array(lcavity$, energy_start$)  = 'ENERGY_START'
    attrib_array(lcavity$, phi0$)          = 'PHI0'
    attrib_array(lcavity$, gradient$)      = 'GRADIENT'
    attrib_array(lcavity$, rf_frequency$)  = 'RF_FREQUENCY'
    attrib_array(lcavity$, e_loss$)        = 'E_LOSS'
    attrib_array(lcavity$, delta_e$)       = 'DELTA_E'
    attrib_array(lcavity$, sr_wake_file$)  = 'SR_WAKE_FILE'
    attrib_array(lcavity$, lr_wake_file$)  = 'LR_WAKE_FILE'
    attrib_array(lcavity$, field_calc$)    = 'FIELD_CALC'
    attrib_array(lcavity$, freq_spread$)   = 'FREQ_SPREAD'

    attrib_array(group$, command$)        = 'COMMAND'
    attrib_array(group$, old_command$)    = 'OLD_COMMAND'
    attrib_array(group$, coef$)           = 'COEF'
    attrib_array(group$, start_edge$)     = 'START_EDGE'
    attrib_array(group$, end_edge$)       = 'END_EDGE'
    attrib_array(group$, accordion_edge$) = 'ACCORDION_EDGE'
    attrib_array(group$, symmetric_edge$) = 'SYMMETRIC_EDGE'

    attrib_array(drift$, l$)          = 'L'
    attrib_array(drift$, is_on$)      =  null_name    
    attrib_array(drift$, field_calc$) = 'FIELD_CALC'

    attrib_array(monitor$,    l$) = 'L'
    attrib_array(instrument$, l$) = 'L'

    attrib_array(rcollimator$, l$)     = 'L'
    attrib_array(ecollimator$, l$)     = 'L'

    attrib_array(hkicker$, l$)          = 'L'
    attrib_array(hkicker$, tilt$)       = 'TILT'
    attrib_array(hkicker$, kick$)       = 'KICK'
    attrib_array(hkicker$, field_calc$) = 'FIELD_CALC'

    attrib_array(vkicker$, l$)          = 'L'
    attrib_array(vkicker$, tilt$)       = 'TILT'
    attrib_array(vkicker$, kick$)       = 'KICK'
    attrib_array(vkicker$, field_calc$) = 'FIELD_CALC'

    attrib_array(kicker$, l$)           = 'L'
    attrib_array(kicker$, hkick$)       = 'HKICK'
    attrib_array(kicker$, vkick$)       = 'VKICK'
    attrib_array(kicker$, tilt$)        = 'TILT'
    attrib_array(kicker$, h_displace$)  = 'H_DISPLACE'
    attrib_array(kicker$, v_displace$)  = 'V_DISPLACE'
    attrib_array(kicker$, radius$)      = 'RADIUS'
    attrib_array(kicker$, field_calc$)  = 'FIELD_CALC'

    attrib_array(sbend$, l$)          = 'L'
    attrib_array(sbend$, angle$)      = 'ANGLE'
    attrib_array(sbend$, e1$)         = 'E1'
    attrib_array(sbend$, e2$)         = 'E2'
    attrib_array(sbend$, h1$)         = 'H1'
    attrib_array(sbend$, h2$)         = 'H2'
    attrib_array(sbend$, k1$)         = 'K1'
    attrib_array(sbend$, g$)          = 'G'
    attrib_array(sbend$, delta_g$)    = 'DELTA_G'
    attrib_array(sbend$, tilt$)       = 'TILT'
    attrib_array(sbend$, roll$)       = 'ROLL'
    attrib_array(sbend$, hgap$)       = 'HGAP'
    attrib_array(sbend$, hgapx$)      = 'HGAPX'
    attrib_array(sbend$, fint$)       = 'FINT'
    attrib_array(sbend$, fintx$)      = 'FINTX'
    attrib_array(sbend$, rho$)        = 'RHO'
    attrib_array(sbend$, l_chord$)    = 'L_CHORD'
    attrib_array(sbend$, b_field$)    = 'B_FIELD'
    attrib_array(sbend$, radius$)     = 'RADIUS'
    attrib_array(sbend$, field_calc$) = 'FIELD_CALC'

    attrib_array(bend_sol_quad$, l$)          = 'L'
    attrib_array(bend_sol_quad$, angle$)      = 'ANGLE'
    attrib_array(bend_sol_quad$, k1$)         = 'K1'
    attrib_array(bend_sol_quad$, g$)          = 'G'
    attrib_array(bend_sol_quad$, ks$)         = 'KS'
    attrib_array(bend_sol_quad$, dks_ds$)     = 'DKS_DS'
    attrib_array(bend_sol_quad$, quad_tilt$)  = 'QUAD_TILT'
    attrib_array(bend_sol_quad$, bend_tilt$)  = 'BEND_TILT'
    attrib_array(bend_sol_quad$, x_quad$)     = 'X_QUAD'
    attrib_array(bend_sol_quad$, y_quad$)     = 'Y_QUAD'
    attrib_array(bend_sol_quad$, tilt$)       = 'TILT'
    attrib_array(bend_sol_quad$, rho$)        = 'RHO'
    attrib_array(bend_sol_quad$, radius$)     = 'RADIUS'
    attrib_array(bend_sol_quad$, field_calc$) = 'FIELD_CALC'

    attrib_array(patch$, x_pitch$)    = 'X_PITCH'
    attrib_array(patch$, y_pitch$)    = 'Y_PITCH'
    attrib_array(patch$, x_offset$)   = 'X_OFFSET'
    attrib_array(patch$, y_offset$)   = 'Y_OFFSET'
    attrib_array(patch$, z_offset$)   = 'Z_OFFSET'
    attrib_array(patch$, dE_offset$)  = 'DE_OFFSET'
    attrib_array(patch$, tilt$)       = 'TILT'


    attrib_array(rbend$, l$)          = 'L'
    attrib_array(rbend$, angle$)      = 'ANGLE'
    attrib_array(rbend$, e1$)         = 'E1'
    attrib_array(rbend$, e2$)         = 'E2'
    attrib_array(rbend$, h1$)         = 'H1'
    attrib_array(rbend$, h2$)         = 'H2'
    attrib_array(rbend$, k1$)         = 'K1'
    attrib_array(rbend$, g$)          = 'G'
    attrib_array(rbend$, delta_g$)    = 'DELTA_G'
    attrib_array(rbend$, tilt$)       = 'TILT'
    attrib_array(rbend$, roll$)       = 'ROLL'
    attrib_array(rbend$, hgap$)       = 'HGAP'
    attrib_array(rbend$, hgapx$)      = 'HGAPX'
    attrib_array(rbend$, fint$)       = 'FINT'
    attrib_array(rbend$, fintx$)      = 'FINTX'
    attrib_array(rbend$, rho$)        = 'RHO'
    attrib_array(rbend$, l_chord$)    = 'L_CHORD'
    attrib_array(rbend$, b_field$)    = 'B_FIELD'
    attrib_array(rbend$, radius$)     = 'RADIUS'
    attrib_array(rbend$, field_calc$) = 'FIELD_CALC'

    attrib_array(quadrupole$, l$)          = 'L'
    attrib_array(quadrupole$, tilt$)       = 'TILT'
    attrib_array(quadrupole$, k1$)         = 'K1'
    attrib_array(quadrupole$, B_gradient$) = 'B_GRADIENT'
    attrib_array(quadrupole$, radius$)     = 'RADIUS'
    attrib_array(quadrupole$, field_calc$) = 'FIELD_CALC'

    attrib_array(sextupole$, l$)          = 'L'
    attrib_array(sextupole$, tilt$)       = 'TILT'
    attrib_array(sextupole$, k2$)         = 'K2'
    attrib_array(sextupole$, B_gradient$) = 'B_GRADIENT'
    attrib_array(sextupole$, radius$)     = 'RADIUS'
    attrib_array(sextupole$, field_calc$) = 'FIELD_CALC'

    attrib_array(octupole$, l$)          = 'L'
    attrib_array(octupole$, tilt$)       = 'TILT'
    attrib_array(octupole$, k3$)         = 'K3'
    attrib_array(octupole$, B_gradient$) = 'B_GRADIENT'
    attrib_array(octupole$, radius$)     = 'RADIUS'
    attrib_array(octupole$, field_calc$) = 'FIELD_CALC'

    attrib_array(solenoid$, l$)          = 'L'
    attrib_array(solenoid$, ks$)         = 'KS'
    attrib_array(solenoid$, b_field$)    = 'B_FIELD'
    attrib_array(solenoid$, radius$)     = 'RADIUS'
    attrib_array(solenoid$, field_calc$) = 'FIELD_CALC'

    attrib_array(rfcavity$, l$)             = 'L'
    attrib_array(rfcavity$, voltage$)       = 'VOLTAGE'
    attrib_array(rfcavity$, rf_frequency$)  = 'RF_FREQUENCY'
    attrib_array(rfcavity$, phi0$)          = 'PHI0'
    attrib_array(rfcavity$, harmon$)        = 'HARMON'
    attrib_array(rfcavity$, field_calc$)    = 'FIELD_CALC'

    attrib_array(elseparator$, l$)          = 'L'
    attrib_array(elseparator$, gap$)        = 'GAP'
    attrib_array(elseparator$, e_field$)    = 'E_FIELD'
    attrib_array(elseparator$, voltage$)    = 'VOLTAGE'
    attrib_array(elseparator$, tilt$)       = 'TILT'
    attrib_array(elseparator$, radius$)     = 'RADIUS'
    attrib_array(elseparator$, field_calc$) = 'FIELD_CALC'

    attrib_array(beambeam$, sig_x$)       = 'SIG_X'
    attrib_array(beambeam$, sig_y$)       = 'SIG_Y'
    attrib_array(beambeam$, sig_z$)       = 'SIG_Z'
    attrib_array(beambeam$, bbi_const$)   = 'BBI_CONSTANT'
    attrib_array(beambeam$, charge$)      = 'CHARGE'
    attrib_array(beambeam$, n_slice$)     = 'N_SLICE'
    attrib_array(beambeam$, symplectify$) = 'N_SLICE'
    attrib_array(beambeam$, x_offset$)    = 'X_OFFSET'
    attrib_array(beambeam$, y_offset$)    = 'Y_OFFSET'
    attrib_array(beambeam$, s_offset$)    = 'S_OFFSET'
    attrib_array(beambeam$, x_pitch$)     = 'X_PITCH'
    attrib_array(beambeam$, y_pitch$)     = 'Y_PITCH'
    attrib_array(beambeam$, tilt$)        = 'TILT'
    attrib_array(beambeam$, field_calc$)  = 'FIELD_CALC'

    attrib_array(wiggler$, l$)          = 'L'
    attrib_array(wiggler$, k1$)         = 'K1'
    attrib_array(wiggler$, b_max$)      = 'B_MAX'
    attrib_array(wiggler$, rho$)        = 'RHO'
    attrib_array(wiggler$, n_pole$)     = 'N_POLE'
    attrib_array(wiggler$, tilt$)       = 'TILT'
    attrib_array(wiggler$, radius$)     = 'RADIUS'
    attrib_array(wiggler$, term$)       = 'TERM'
    attrib_array(wiggler$, polarity$)   = 'POLARITY'
    attrib_array(wiggler$, z_patch$)    = 'Z_PATCH'
    attrib_array(wiggler$, radius$)     = 'RADIUS'
    attrib_array(wiggler$, field_calc$) = 'FIELD_CALC'

    attrib_array(sol_quad$, l$)          = 'L'
    attrib_array(sol_quad$, k1$)         = 'K1'
    attrib_array(sol_quad$, ks$)         = 'KS'
    attrib_array(sol_quad$, tilt$)       = 'TILT'
    attrib_array(sol_quad$, radius$)     = 'RADIUS'
    attrib_array(sol_quad$, field_calc$) = 'FIELD_CALC'

    attrib_array(def_beam$, particle$)   = 'PARTICLE'
    attrib_array(def_beam$, energy_gev$) = 'ENERGY'
    attrib_array(def_beam$, n_part$)     = 'N_PART'

    attrib_array(multipole$, tilt$)     = 'TILT'
    attrib_array(multipole$, k0l$:k20l$) = (/ 'K0L ', &
                                   'K1L ', 'K2L ', 'K3L ', 'K4L ', 'K5L ', & 
                                   'K6L ', 'K7L ', 'K8L ', 'K9L ', 'K10L', &
                                   'K11L', 'K12L', 'K13L', 'K14L', 'K15L', &
                                   'K16L', 'K17L', 'K18L', 'K19L', 'K20L' /)
    attrib_array(multipole$, t0$:t20$) = (/ 'T0 ', &
                                   'T1 ', 'T2 ', 'T3 ', 'T4 ', 'T5 ', & 
                                   'T6 ', 'T7 ', 'T8 ', 'T9 ', 'T10', &
                                   'T11', 'T12', 'T13', 'T14', 'T15', &
                                   'T16', 'T17', 'T18', 'T19', 'T20' /)
    attrib_array(multipole$, x_offset$) = 'X_OFFSET'
    attrib_array(multipole$, y_offset$) = 'Y_OFFSET'
    attrib_array(multipole$, s_offset$) = 'S_OFFSET'
    attrib_array(multipole$, lrad$)     = 'LRAD'

    attrib_array(ab_multipole$, tilt$)     = 'TILT'
    attrib_array(ab_multipole$, x_offset$) = 'X_OFFSET'
    attrib_array(ab_multipole$, y_offset$) = 'Y_OFFSET'
    attrib_array(ab_multipole$, s_offset$) = 'S_OFFSET'
    attrib_array(ab_multipole$, lrad$)     = 'LRAD'

    attrib_array(accel_sol$, l$)             = 'L'
    attrib_array(accel_sol$, voltage$)       = 'VOLTAGE'
    attrib_array(accel_sol$, phi0$)          = 'PHI0'
    attrib_array(accel_sol$, rf_wavelength$) = 'RF_WAVELENGTH'
    attrib_array(accel_sol$, b_z$)           = 'B_Z'
    attrib_array(accel_sol$, b_x1$)          = 'B_X1'
    attrib_array(accel_sol$, b_y1$)          = 'B_Y1'
    attrib_array(accel_sol$, s_st1$)         = 'S_ST1'
    attrib_array(accel_sol$, l_st1$)         = 'L_ST1'
    attrib_array(accel_sol$, b_x2$)          = 'B_X2'
    attrib_array(accel_sol$, b_y2$)          = 'B_Y2'
    attrib_array(accel_sol$, s_st2$)         = 'S_ST2'
    attrib_array(accel_sol$, l_st2$)         = 'L_ST2'
    attrib_array(accel_sol$, x_beg_limit$)   = 'X_BEG_LIMIT'
    attrib_array(accel_sol$, y_beg_limit$)   = 'Y_BEG_LIMIT'
    attrib_array(accel_sol$, field_calc$)    = 'FIELD_CALC'

    attrib_array(custom$, l$)          = 'L'
    attrib_array(custom$, tilt$)       = 'TILT'
    attrib_array(custom$,  val1$)      = 'VAL1'
    attrib_array(custom$,  val2$)      = 'VAL2'
    attrib_array(custom$,  val3$)      = 'VAL3'
    attrib_array(custom$,  val4$)      = 'VAL4'
    attrib_array(custom$,  val5$)      = 'VAL5'
    attrib_array(custom$,  val6$)      = 'VAL6'
    attrib_array(custom$,  val7$)      = 'VAL7'
    attrib_array(custom$,  val8$)      = 'VAL8'
    attrib_array(custom$,  val9$)      = 'VAL9'
    attrib_array(custom$, val10$)      = 'VAL10'
    attrib_array(custom$, val11$)      = 'VAL11'
    attrib_array(custom$, val12$)      = 'VAL12'
    attrib_array(custom$, x_offset$)   = 'X_OFFSET'
    attrib_array(custom$, y_offset$)   = 'Y_OFFSET'
    attrib_array(custom$, s_offset$)   = 'S_OFFSET'
    attrib_array(custom$, x_pitch$)    = 'X_PITCH'
    attrib_array(custom$, y_pitch$)    = 'Y_PITCH'
    attrib_array(custom$, field_calc$) = 'FIELD_CALC'
    attrib_array(custom$, delta_e$)    = 'DELTA_E'

    init_needed = .false.

  endif

!--------------------------------------------------------------------
! Main part of subroutine

  key = ele%key

  if (key <= 0 .or. key > n_key) then
    at_name = 'BAD ELE KEY'
  elseif (ix_att <= 0 .or. ix_att > n_attrib_special_maxx) then
    at_name = 'BAD INDEX'
  elseif (ele%control_type == overlay_lord$) then
    if (ix_att == ele%ix_value) then
      at_name = ele%attribute_name
    else
      at_name = 'INVALID INDEX'
    endif
  else
    at_name = attrib_array(key, ix_att)
  endif

end function




