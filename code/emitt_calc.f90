!+
! Subroutine EMITT_CALC (RING, WHAT, MODE)
!
! Subroutine to calculate the emittance, energy spread, and synchrotron
! integrals. This subroutine assumes that bends are in the horizontal plane.
!
! For a better, more complete calculation see the subroutine: 
!               RADIATION_INTEGRALS
! The only possible saving grace of this subroutine is that it is faster
! than RADIATION_INTEGRALS.
!
! Modules Needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!     RING -- Ring_struct: Ring to use
!     WHAT -- Integer: Which elements to use in the calculation.
!               = BENDS$      Use only the bends.
!               = WIGGLERS$   Use only wigglers.
!               = All$        Use both bends and wigglers.
!
! Output:
!     MODE -- Modes_struct: Parameters for the ("horizontal like") a-mode,
!                              ("vertical like") b-mode, and the z-mode
!       %SYNCH_INT(1:3) -- Synchrotron integrals.
!       %sigE_E         -- Sigma_E/E energy spread
!       %SIG_Z          -- Bunch Length
!       %E_LOSS         -- Energy loss in eV per turn
!       %A, %B, %Z      -- Amode_struct: Substructure
!         %EMITTANCE      -- Emittance
!         %SYNCH_INT(4:5) -- Synchrotron integrals
!         %J_DAMP         -- Damping partition factor
!         %ALPHA_DAMP     -- Exponential damping coefficient per turn
!
! Notes:
!     1) %SYNCH_INT(1) = momentum_compaction * ring_length
!     2) The calculation assumes that the Twiss parameters have been calculated.
!-


subroutine emitt_calc (ring, what, mode)

  use bmad_struct
  implicit none

  type (ring_struct)  ring
  type (modes_struct)  mode
  type (ele_struct)  ele0, ele

  real(rp) :: c_gam = 4.425e-5, c_q = 3.84e-13
  real(rp) energy_loss, k_1, k_2, k_3, k_4, k_5
  real(rp) i1, i2, i3, i4a, i4b, i4z, i5a, i5b, gamma2_factor
  real(rp) g, e1, e2, sin_kl, cos_kl, tan_e1, tan_e2
  real(rp) c11, c12, c21, c22, kl, ll, gamma_c, k1, g3_ave, energy
  real(rp) eta_x0, etap_x0, eta_a0, etap_a0, eta_b0, etap_b0
  real(rp) i4a_bend, i4z_bend, end_4a, end_4
  real(rp) gamma_a0, beta_a0, alpha_a0, gamma_b0, beta_b0, alpha_b0
  real(rp) eta_x_int, eta_x20(2), t1, t2, t3, G_max
  real(rp) eta_ax0, etap_ax0, m1, m2, n1, n2, beta_eff, alpha_eff, gamma_eff
  real(rp) arg, m65

  integer ir, what

  logical :: do_bends, do_wigs, do_type_err_message = .true. 
  logical print_wig_err_message

!----------------------------------------------------------------------
! init

  i1 = 0
  i2 = 0
  i3 = 0
  i4a = 0
  i4b = 0
  i4z = 0
  i5a = 0
  i5b = 0
  m65 = 0.0

  print_wig_err_message = .true.   ! print warning just once.

  if (what == bends$) then
    do_bends = .true.
    do_wigs  = .false.
  elseif (what == wigglers$) then
    do_bends = .false.
    do_wigs  = .true.
  elseif (what == all$) then
    do_bends = .true.
    do_wigs  = .true.
  else
    print *, 'ERROR IN EMITT_CALC: UNKNOWN "WHAT" SWITCH:', what
    call err_exit
  endif


!---------------------------------------------------------------------
! loop over all elements

  do ir = 1, ring%n_ele_use

!---------------------------------------------------------------------
! bend contribution

    if (do_bends .and. ring%ele_(ir)%key == sbend$) then

      ele  = ring%ele_(ir)

      g = ele%value(g$)
      ll = ele%value(l$)
      k1 = ele%value(k1$)
      e1 = ele%value(e1$)
      e2 = ele%value(e2$)

      k_2 = g**2 + k1

      if (k_2 > 0) then
        k_1 = sqrt(k_2)
        k_3 = k_2 * k_1
        k_4 = k_3 * k_1
        k_5 = k_4 * k_1
        kl = k_1 * ll
        sin_kl = sin(kl)
        cos_kl = cos(kl)
      else
        k_1 = sqrt(-k_2)   ! +
        k_3 = k_2 * k_1    ! -
        k_4 = k_2 * k_2    ! +
        k_5 = k_4 * k_1    ! +
        kl = k_1 * ll
        sin_kl = sinh(kl)
        cos_kl = cosh(kl)      
      endif


! if there is edge focusing then use parameters just past the edge focusing

      if (e1 == 0) then
        tan_e1 = 0
        ele0 = ring%ele_(ir-1)
      else
        tan_e1 = tan(e1)
        call mat_make_unit (ele0%mat6)
        ele0%mat6(2,1) =  tan_e1 * g
        ele0%mat6(4,3) = -tan_e1 * g
        call twiss_propagate1 (ring%ele_(ir-1), ele0)
      endif

      if (e2 == 0) then
        tan_e2 = 0
      else
        tan_e2 = tan(e2)
      endif

      gamma_c = ele0%gamma_c

      c11 = ele0%c_mat(1,1)
      c12 = ele0%c_mat(1,2)
      c21 = ele0%c_mat(2,1)
      c22 = ele0%c_mat(2,2)

      beta_a0  = ele0%x%beta
      alpha_a0 = ele0%x%alpha
      gamma_a0 = ele0%x%gamma

      beta_b0  = ele0%y%beta
      alpha_b0 = ele0%y%alpha
      gamma_b0 = ele0%y%gamma

      eta_a0  = ele0%x%eta
      etap_a0 = ele0%x%etap

      eta_b0  = ele0%y%eta
      etap_b0 = ele0%y%etap

      eta_x0  = gamma_c * eta_a0  + c11 * eta_b0 + c12 * etap_b0
      etap_x0 = gamma_c * etap_a0 + c21 * eta_b0 + c22 * etap_b0

      eta_ax0  = gamma_c * eta_a0
      etap_ax0 = gamma_c * etap_a0

! i1 calc

      eta_x_int = eta_x0 * sin_kl / k_1 + etap_x0 * (1 - cos_kl) / k_2 +  &
                                            (kl - sin_kl) * g / k_3
      i1 = i1 + eta_x_int * g

! i2 calc

      i2 = i2 + ll * g**2

! i3 calc

      i3 = i3 + ll * abs(g**3)

! i4a, i4b, i4z calcs

      t1 = eta_ax0 * sin_kl / k_1 + etap_ax0 * (1 - cos_kl) / k_2
      t2 = gamma_c**2 * (kl - sin_kl) * g / k_3
      
      end_4  = (eta_x0 * tan_e1 + ele%x%eta_lab * tan_e2) * g**2
      end_4a = (eta_ax0 * tan_e1 + gamma_c * ele%x%eta * tan_e2) * g**2

      i4a_bend = (t1 + t2) * g**3 - end_4a
      i4z_bend = eta_x_int * g**3 - end_4

      i4a = i4a + i4a_bend
      i4b = i4b + i4z_bend - i4a_bend
      i4z = i4z + i4z_bend

! i5a calc

      t1 = gamma_a0*eta_a0**2 + 2*alpha_a0*eta_a0*etap_a0 +  &
                                                    beta_a0*etap_a0**2
      t2 = (gamma_a0*eta_a0 + alpha_a0*etap_a0) * (sin_kl - kl) / k_3 +  &
                 (alpha_a0*eta_a0 + beta_a0*etap_a0) * (1 - cos_kl) / k_2
      t3 = gamma_a0 * (3*kl - 4*sin_kl + sin_kl*cos_kl) / (2*k_5) -  &
                 alpha_a0 * (1 - cos_kl)**2 / k_4 +  &
                 beta_a0 * (kl - cos_kl*sin_kl) / (2*k_3)


      i5a = i5a + (ll*t1 + 2*gamma_c * t2 * g +  &
                                   t3*(gamma_c*g)**2) * abs(g**3)

! i5b calc

      i5b = i5b + ll * (gamma_b0*eta_b0**2 + 2*alpha_b0*eta_b0*etap_b0 +  &
                                  beta_b0*etap_b0**2) * abs(g**3)

      m1 = gamma_b0 * eta_b0 + alpha_b0 * etap_b0
      m2 = alpha_b0 * eta_b0 + beta_b0  * etap_b0

      n1 =  m1 * c22 - m2 * c21
      n2 = -m1 * c12 + m2 * c11

      eta_x20(1) = (sin_kl - kl) * g / k_3
      eta_x20(2) = (1 - cos_kl) * g / k_2

      i5b = i5b + 2 * (n1 * eta_x20(1) + n2 * eta_x20(2)) * abs(g**3)

      beta_eff  = c11**2 * beta_b0 - 2*c11*c12 * alpha_b0 +  &
                                                         c12**2 * gamma_b0
      alpha_eff = -c21*c11 * beta_b0 + (c11*c22+c12*c21) * alpha_b0 -  &
                                                         c12*c22 * gamma_b0
      gamma_eff = c21**2 * beta_b0 - 2*c21*c22 * alpha_b0 +  &
                                                         c22**2 * gamma_b0

      i5b = i5b +  &
             (gamma_eff * (3*kl - 4*sin_kl + sin_kl*cos_kl) / (2*k_5) -  &
              alpha_eff * (1 - cos_kl)**3 / k_4 +  &
              beta_eff * (kl - cos_kl*sin_kl) / (2*k_3)) * abs(g**5)


!---------------------------------------------------------------------
! wiggler contribution

    elseif (do_wigs .and. ring%ele_(ir)%key == wiggler$) then

      if (ring%ele_(ir)%sub_key == map_type$) then
        if (print_wig_err_message) then
          print *, 'WARNING FROM EMITT_CALC: ',&
                        'I AM NOT UP TO HANDLING "NEW STYLE" WIGGLERS.'
          print *, '        RESULTS WILL NOT BE VALID.'
          print_wig_err_message = .false.
        endif
      endif

      ele = ring%ele_(ir)
      ele0 = ring%ele_(ir-1)
      k1 = ele%value(k1$)

      G_max = sqrt(2*abs(k1))       ! g at max B
      ll = ele%value(l$)


      beta_a0  = ele0%x%beta
      alpha_a0 = ele0%x%alpha
      gamma_a0 = ele0%x%gamma

      beta_b0  = ele0%y%beta
      alpha_b0 = ele0%y%alpha
      gamma_b0 = ele0%y%gamma

      eta_a0  = ele0%x%eta
      etap_a0 = ele0%x%etap

      eta_b0  = ele0%y%eta
      etap_b0 = ele0%y%etap

! i2 calc

      i2 = i2 + ll * G_max**2 / 2

! i3 calc

      g3_ave = 4 * G_max**3 / (3 * pi)
      i3 = i3 + ll * g3_ave

! i5a calc

      i5a = i5a + ll * g3_ave * (gamma_a0*eta_a0**2 +  &
                        2*alpha_a0*eta_a0*etap_a0 + beta_a0*etap_a0**2)

! i5b calc

      i5b = i5b + ll * g3_ave * (gamma_b0*eta_b0**2 +  &
                        2*alpha_b0*eta_b0*etap_b0 + beta_b0*etap_b0**2)


! custom contribution

    elseif (do_wigs .and. ring%ele_(ir)%key == custom$) then
       call custom_emitt_calc (ring, ir, i2, i3, i5a, i5b)

    elseif (ring%ele_(ir)%key == rfcavity$) then
      m65 = m65 + ring%ele_(ir)%mat6(6,5)   ! add up the m65s
    endif


  enddo

!---------------------------------------------------------------------
! now put everything together

  energy = ring%ele_(0)%value(beam_energy$)
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
    if(arg > 0.) mode%sigE_E = sqrt(c_q * i3 * gamma2_factor / (2*i2 + i4z))

  endif

  mode%a%alpha_damp = energy_loss * mode%a%j_damp / energy
  mode%b%alpha_damp = energy_loss * mode%b%j_damp / energy
  mode%z%alpha_damp = energy_loss * mode%z%j_damp / energy

  mode%e_loss = energy_loss

  if(abs(m65) > 0. ) then
    mode%sig_z = sqrt( mode%synch_int(1)/abs(m65) ) * mode%sigE_E
  else
    mode%sig_z = 0.
  endif

  return
  end
