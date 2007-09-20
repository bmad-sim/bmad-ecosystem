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
!   ele    -- Ele_struct: 
!     %key    -- Integer: Key name of element type (e.g. SBEND$, etc.)
!   ix_att -- Integer: Index of attribute (e.g. k1$)
!
! Output:
!   attribute_name -- Character(40): Name of attribute. 
!              If %key is invalid then                   attribute_name = "!BAD ELE KEY"
!              If ix_att is invalid then                 attribute_name = "!BAD INDEX"
!              If ix_att is invalid for an overlay then  attribute_name = "!INVALID INDEX"
!              If ix_att does not correspond to an attribute for the given key then
!                                                        attribute_name = "!NULL"
!
! Example:
!   ele%key = sbend$
!   name = attribute_name (ele, k1$)
! Result:
!   name -> 'K1'
!-

#include "CESR_platform.inc"

function attribute_name (ele, ix_att) result (at_name)

  use bmad_struct
  use bmad_interface, except_dummy => attribute_name

  implicit none

  type (ele_struct) ele

  integer i, key, ix_att

  character(40) attrib_array(n_key, n_attrib_special_maxx)
  character(40) at_name

  logical, save :: init_needed = .true.

!---------------------------------------------------------------------
! Init

  if (init_needed) then

    attrib_array = null_name

    do i = 1, n_key
                                    
      if (i == def_beam$)       cycle
      if (i == def_parameter$)  cycle
      if (i == def_beam_start$) cycle

      attrib_array(i, type$)     = 'TYPE'
      attrib_array(i, alias$)    = 'ALIAS'
      attrib_array(i, descrip$)  = 'DESCRIP'

      if (i == init_ele$) cycle
      if (i == group$)    cycle
      if (i == overlay$)  cycle
      if (i == girder$)   cycle

      attrib_array(i, x_limit$)               = 'X_LIMIT'
      attrib_array(i, x1_limit$)              = 'X1_LIMIT'
      attrib_array(i, x2_limit$)              = 'X2_LIMIT'
      attrib_array(i, y_limit$)               = 'Y_LIMIT'
      attrib_array(i, y1_limit$)              = 'Y1_LIMIT'
      attrib_array(i, y2_limit$)              = 'Y2_LIMIT'
      attrib_array(i, aperture$)              = 'APERTURE'
      attrib_array(i, aperture_at$)           = 'APERTURE_AT'
      attrib_array(i, offset_moves_aperture$) = 'OFFSET_MOVES_APERTURE'
      attrib_array(i, E_tot$)                 = 'E_TOT'
      attrib_array(i, p0c$)                   = 'P0C'
      attrib_array(i, mat6_calc_method$)      = 'MAT6_CALC_METHOD'
      attrib_array(i, tracking_method$)       = 'TRACKING_METHOD'

      attrib_array(i, is_on$)       = 'IS_ON'

      if (i == marker$)       cycle
      if (i == match$)        cycle
      if (i == patch$)        cycle
      if (i == beambeam$)     cycle
      if (i == hom$)          cycle
      if (i == multipole$)    cycle 
      if (i == ab_multipole$) cycle

      attrib_array(i, symplectify$)       = 'SYMPLECTIFY'
      attrib_array(i, map_with_offsets$)    = 'MAP_WITH_OFFSETS'

      if (i == taylor$)       cycle

      attrib_array(i, integrator_order$)  = 'INTEGRATOR_ORDER'
      attrib_array(i, num_steps$)         = 'NUM_STEPS'
      attrib_array(i, ds_step$)           = 'DS_STEP'
      attrib_array(i, rel_tol$)           = 'REL_TOL'
      attrib_array(i, abs_tol$)           = 'ABS_TOL'
      attrib_array(i, csr_calc_on$)       = 'CSR_CALC_ON'

      if (i == hkicker$)      cycle
      if (i == vkicker$)      cycle
      if (i == custom$)       cycle

      attrib_array(i, hkick$)    = 'HKICK'
      attrib_array(i, vkick$)    = 'VKICK'
      attrib_array(i, bl_hkick$) = 'BL_HKICK'
      attrib_array(i, bl_vkick$) = 'BL_VKICK'

      if (i == drift$)        cycle
      if (i == kicker$)       cycle
      if (i == monitor$)      cycle
      if (i == instrument$)   cycle

      attrib_array(i, x_offset$) = 'X_OFFSET'
      attrib_array(i, y_offset$) = 'Y_OFFSET'
      attrib_array(i, s_offset$) = 'S_OFFSET'

      if (i == rcollimator$) cycle
      if (i == ecollimator$) cycle

      attrib_array(i, x_pitch$)   = 'X_PITCH'
      attrib_array(i, y_pitch$)   = 'Y_PITCH'

    enddo

!

    do i = 1, n_key
      select case (i)
      case (elseparator$, kicker$, octupole$, quadrupole$, sbend$, rbend$, &
             sextupole$, solenoid$, sol_quad$, ab_multipole$, wiggler$, bend_sol_quad$, &
             hkicker$, vkicker$)
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

    attrib_array(def_parameter$, E_TOT$)        = 'E_TOT'
    attrib_array(def_parameter$, lattice_type$) = 'LATTICE_TYPE'
    attrib_array(def_parameter$, lattice$)      = 'LATTICE'
    attrib_array(def_parameter$, taylor_order$) = 'TAYLOR_ORDER'
    attrib_array(def_parameter$, ran_seed$)     = 'RAN_SEED'
    attrib_array(def_parameter$, n_part$)       = 'N_PART'
    attrib_array(def_parameter$, particle$)     = 'PARTICLE'

    attrib_array(def_beam$, particle$)   = 'PARTICLE'
    attrib_array(def_beam$, energy_gev$) = 'ENERGY'
    attrib_array(def_beam$, n_part$)     = 'N_PART'

    attrib_array(def_beam_start$, x$)     = 'X'
    attrib_array(def_beam_start$, p_x$)   = 'P_X'
    attrib_array(def_beam_start$, y$)     = 'Y'
    attrib_array(def_beam_start$, p_y$)   = 'P_Y'
    attrib_array(def_beam_start$, z$)     = 'Z'
    attrib_array(def_beam_start$, p_z$)   = 'P_Z'

    attrib_array(taylor$, l$)           = 'L'
    attrib_array(taylor$, x_offset$)    = 'X_OFFSET'   
    attrib_array(taylor$, y_offset$)    = 'Y_OFFSET'   
    attrib_array(taylor$, s_offset$)    = 'S_OFFSET'   
    attrib_array(taylor$, x_pitch$)     = 'X_PITCH'   
    attrib_array(taylor$, y_pitch$)     = 'Y_PITCH'   
    attrib_array(taylor$, tilt$)        = 'TILT' 

    attrib_array(match$, l$)             = 'L'
    attrib_array(match$, beta_a0$)       = 'BETA_A0'
    attrib_array(match$, alpha_a0$)      = 'ALPHA_A0'
    attrib_array(match$, beta_b0$)       = 'BETA_B0'
    attrib_array(match$, alpha_b0$)      = 'ALPHA_B0'
    attrib_array(match$, beta_a1$)       = 'BETA_A1'
    attrib_array(match$, alpha_a1$)      = 'ALPHA_A1'
    attrib_array(match$, beta_b1$)       = 'BETA_B1'
    attrib_array(match$, alpha_b1$)      = 'ALPHA_B1'
    attrib_array(match$, dphi_a$)        = 'DPHI_A'
    attrib_array(match$, dphi_b$)        = 'DPHI_B'
    attrib_array(match$, eta_a0$)        = 'ETA_A0'
    attrib_array(match$, etap_a0$)       = 'ETAP_A0'
    attrib_array(match$, eta_b0$)        = 'ETA_B0'
    attrib_array(match$, etap_b0$)       = 'ETAP_B0'
    attrib_array(match$, eta_a1$)        = 'ETA_A1'
    attrib_array(match$, etap_a1$)       = 'ETAP_A1'
    attrib_array(match$, eta_b1$)        = 'ETA_B1'
    attrib_array(match$, etap_b1$)       = 'ETAP_B1'

    attrib_array(girder$, x_offset$)     = 'X_OFFSET'
    attrib_array(girder$, y_offset$)     = 'Y_OFFSET'
    attrib_array(girder$, s_offset$)     = 'S_OFFSET'
    attrib_array(girder$, x_pitch$)      = 'X_PITCH'
    attrib_array(girder$, y_pitch$)      = 'Y_PITCH'
    attrib_array(girder$, s_center$)     = 'S_CENTER'
    attrib_array(girder$, tilt$)         = 'TILT'

    attrib_array(lcavity$, l$)             = 'L'
    attrib_array(lcavity$, tilt$)          = 'TILT'
    attrib_array(lcavity$, lrad$)          = 'LRAD'   ! This is for felv testing.
    attrib_array(lcavity$, p0c_start$)     = 'P0C_START'
    attrib_array(lcavity$, E_TOT_START$)  = 'E_TOT_START'
    attrib_array(lcavity$, dphi0$)         = 'DPHI0'
    attrib_array(lcavity$, phi0$)          = 'PHI0'
    attrib_array(lcavity$, gradient$)      = 'GRADIENT'
    attrib_array(lcavity$, rf_frequency$)  = 'RF_FREQUENCY'
    attrib_array(lcavity$, e_loss$)        = 'E_LOSS'
    attrib_array(lcavity$, delta_e$)       = 'DELTA_E'
    attrib_array(lcavity$, sr_wake_file$)  = 'SR_WAKE_FILE'
    attrib_array(lcavity$, lr_wake_file$)  = 'LR_WAKE_FILE'
    attrib_array(lcavity$, field_calc$)    = 'FIELD_CALC'
    attrib_array(lcavity$, freq_spread$)   = 'FREQ_SPREAD'
    attrib_array(lcavity$, coupler_at$)       = 'COUPLER_AT'
    attrib_array(lcavity$, coupler_strength$) = 'COUPLER_STRENGTH'
    attrib_array(lcavity$, coupler_angle$)    = 'COUPLER_ANGLE'
    attrib_array(lcavity$, coupler_phase$)    = 'COUPLER_PHASE'
    attrib_array(lcavity$, gradient_err$)     = 'GRADIENT_ERR'
    attrib_array(lcavity$, phi0_err$)         = 'PHI0_ERR'


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

    attrib_array(monitor$, l$)        = 'L'
    attrib_array(monitor$, x_offset$) = 'X_OFFSET'
    attrib_array(monitor$, y_offset$) = 'Y_OFFSET'
    attrib_array(monitor$, x_pitch$)  = 'X_PITCH'
    attrib_array(monitor$, y_pitch$)  = 'Y_PITCH'
    attrib_array(monitor$, tilt$)     = 'TILT'

    attrib_array(instrument$, l$) = 'L'
    attrib_array(instrument$, x_offset$) = 'X_OFFSET'
    attrib_array(instrument$, y_offset$) = 'Y_OFFSET'
    attrib_array(instrument$, x_pitch$)  = 'X_PITCH'
    attrib_array(instrument$, y_pitch$)  = 'Y_PITCH'
    attrib_array(instrument$, tilt$)     = 'TILT'

    attrib_array(marker$, x_offset$) = 'X_OFFSET'
    attrib_array(marker$, y_offset$) = 'Y_OFFSET'
    attrib_array(marker$, tilt$)     = 'TILT'

    attrib_array(rcollimator$, l$)     = 'L'
    attrib_array(rcollimator$, is_on$) = null_name

    attrib_array(ecollimator$, l$)     = 'L'
    attrib_array(ecollimator$, is_on$) = null_name

    attrib_array(hkicker$, l$)          = 'L'
    attrib_array(hkicker$, tilt$)       = 'TILT'
    attrib_array(hkicker$, kick$)       = 'KICK'
    attrib_array(hkicker$, field_calc$) = 'FIELD_CALC'
    attrib_array(hkicker$, bl_kick$)    = 'BL_KICK'
    attrib_array(Hkicker$, s_offset$)   = 'S_OFFSET'

    attrib_array(vkicker$, l$)          = 'L'
    attrib_array(vkicker$, tilt$)       = 'TILT'
    attrib_array(vkicker$, kick$)       = 'KICK'
    attrib_array(vkicker$, field_calc$) = 'FIELD_CALC'
    attrib_array(vkicker$, bl_kick$)    = 'BL_KICK'
    attrib_array(vkicker$, s_offset$)   = 'S_OFFSET'

    attrib_array(kicker$, l$)           = 'L'
    attrib_array(kicker$, tilt$)        = 'TILT'
    attrib_array(kicker$, h_displace$)  = 'H_DISPLACE'
    attrib_array(kicker$, v_displace$)  = 'V_DISPLACE'
    attrib_array(kicker$, radius$)      = 'RADIUS'
    attrib_array(kicker$, field_calc$)  = 'FIELD_CALC'
    attrib_array(kicker$, s_offset$)    = 'S_OFFSET'

    attrib_array(sbend$, l$)            = 'L'
    attrib_array(sbend$, angle$)        = 'ANGLE'
    attrib_array(sbend$, e1$)           = 'E1'
    attrib_array(sbend$, e2$)           = 'E2'
    attrib_array(sbend$, h1$)           = 'H1'
    attrib_array(sbend$, h2$)           = 'H2'
    attrib_array(sbend$, k1$)           = 'K1'
    attrib_array(sbend$, k2$)           = 'K2'
    attrib_array(sbend$, g$)            = 'G'
    attrib_array(sbend$, g_err$)        = 'G_ERR'
    attrib_array(sbend$, tilt$)         = 'TILT'
    attrib_array(sbend$, roll$)         = 'ROLL'
    attrib_array(sbend$, hgap$)         = 'HGAP'
    attrib_array(sbend$, hgapx$)        = 'HGAPX'
    attrib_array(sbend$, fint$)         = 'FINT'
    attrib_array(sbend$, fintx$)        = 'FINTX'
    attrib_array(sbend$, rho$)          = 'RHO'
    attrib_array(sbend$, l_chord$)      = 'L_CHORD'
    attrib_array(sbend$, b_field$)      = 'B_FIELD'
    attrib_array(sbend$, b_field_err$)  = 'B_FIELD_ERR'
    attrib_array(sbend$, radius$)       = 'RADIUS'
    attrib_array(sbend$, field_calc$)   = 'FIELD_CALC'

    attrib_array(rbend$, :) = attrib_array(sbend$, :)

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

    attrib_array(patch$, l$)          = 'L'
    attrib_array(patch$, x_pitch$)    = 'X_PITCH'
    attrib_array(patch$, y_pitch$)    = 'Y_PITCH'
    attrib_array(patch$, x_offset$)   = 'X_OFFSET'
    attrib_array(patch$, y_offset$)   = 'Y_OFFSET'
    attrib_array(patch$, z_offset$)   = 'Z_OFFSET'
    attrib_array(patch$, dE_offset$)  = 'DE_OFFSET'
    attrib_array(patch$, tilt$)       = 'TILT'

    attrib_array(quadrupole$, l$)           = 'L'
    attrib_array(quadrupole$, tilt$)        = 'TILT'
    attrib_array(quadrupole$, k1$)          = 'K1'
    attrib_array(quadrupole$, B1_gradient$) = 'B1_GRADIENT'
    attrib_array(quadrupole$, radius$)      = 'RADIUS'
    attrib_array(quadrupole$, field_calc$)  = 'FIELD_CALC'

    attrib_array(sextupole$, l$)           = 'L'
    attrib_array(sextupole$, tilt$)        = 'TILT'
    attrib_array(sextupole$, k2$)          = 'K2'
    attrib_array(sextupole$, B2_gradient$) = 'B2_GRADIENT'
    attrib_array(sextupole$, radius$)      = 'RADIUS'
    attrib_array(sextupole$, field_calc$)  = 'FIELD_CALC'

    attrib_array(octupole$, l$)           = 'L'
    attrib_array(octupole$, tilt$)        = 'TILT'
    attrib_array(octupole$, k3$)          = 'K3'
    attrib_array(octupole$, B3_gradient$) = 'B3_GRADIENT'
    attrib_array(octupole$, radius$)      = 'RADIUS'
    attrib_array(octupole$, field_calc$)  = 'FIELD_CALC'

    attrib_array(solenoid$, l$)          = 'L'
    attrib_array(solenoid$, ks$)         = 'KS'
    attrib_array(solenoid$, bs_field$)   = 'BS_FIELD'
    attrib_array(solenoid$, radius$)     = 'RADIUS'
    attrib_array(solenoid$, field_calc$) = 'FIELD_CALC'

    attrib_array(rfcavity$, l$)             = 'L'
    attrib_array(rfcavity$, dphi0$)         = 'DPHI0'
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

    attrib_array(wiggler$, l$)              = 'L'
    attrib_array(wiggler$, k1$)             = 'K1'
    attrib_array(wiggler$, l_pole$)         = 'L_POLE'
    attrib_array(wiggler$, b_max$)          = 'B_MAX'
    attrib_array(wiggler$, rho$)            = 'RHO'
    attrib_array(wiggler$, n_pole$)         = 'N_POLE'
    attrib_array(wiggler$, tilt$)           = 'TILT'
    attrib_array(wiggler$, radius$)         = 'RADIUS'
    attrib_array(wiggler$, term$)           = 'TERM'
    attrib_array(wiggler$, polarity$)       = 'POLARITY'
    attrib_array(wiggler$, z_patch$)        = 'Z_PATCH'
    attrib_array(wiggler$, x_patch$)        = 'X_PATCH'
    attrib_array(wiggler$, radius$)         = 'RADIUS'
    attrib_array(wiggler$, field_calc$)     = 'FIELD_CALC'
    attrib_array(wiggler$, x_ray_line_len$) = 'X_RAY_LINE_LEN'
    attrib_array(wiggler$, l_start$)        = 'L_START'
    attrib_array(wiggler$, l_end$)          = 'L_END'

    attrib_array(sol_quad$, l$)           = 'L'
    attrib_array(sol_quad$, k1$)          = 'K1'
    attrib_array(sol_quad$, ks$)          = 'KS'
    attrib_array(sol_quad$, tilt$)        = 'TILT'
    attrib_array(sol_quad$, radius$)      = 'RADIUS'
    attrib_array(sol_quad$, field_calc$)  = 'FIELD_CALC'
    attrib_array(sol_quad$, b1_gradient$) = 'B1_GRADIENT'
    attrib_array(sol_quad$, bs_field$)    = 'BS_FIELD'


    attrib_array(multipole$, l$)         = 'L'
    attrib_array(multipole$, tilt$)      = 'TILT'
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

    attrib_array(ab_multipole$, l$)        = 'L'
    attrib_array(ab_multipole$, tilt$)     = 'TILT'
    attrib_array(ab_multipole$, x_offset$) = 'X_OFFSET'
    attrib_array(ab_multipole$, y_offset$) = 'Y_OFFSET'
    attrib_array(ab_multipole$, s_offset$) = 'S_OFFSET'

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
    at_name = '!BAD ELE KEY'
  elseif (ix_att <= 0 .or. ix_att > n_attrib_special_maxx) then
    at_name = '!BAD INDEX'
  elseif (ele%control_type == overlay_lord$) then
    if (ix_att == ele%ix_value) then
      at_name = ele%attribute_name
    else
      at_name = '!INVALID INDEX'
    endif
  else
    at_name = attrib_array(key, ix_att)
  endif

end function




