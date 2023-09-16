module ibs_rates_mod

use bmad_routine_interface
use fgsl
use, intrinsic :: iso_c_binding

type ibs_struct  !these are betatron growth rates.  To get emittance growth rate use:
                 ! demit_a/dt = 2.0*emit_*inv_Ta
  real(rp) inv_Ta
  real(rp) inv_Tb
  real(rp) inv_Tz
end type

real(fgsl_double), parameter :: eps7 = 1.0d-7
integer(fgsl_size_t), parameter :: limit = 1000_fgsl_size_t

contains

!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!+
!  subroutine bjmt1(ele, coulomb_log, rates, n_part)
!
!  This is an implementation of equations 1-9 from "Intrabeam
!  scattering formulas for high energy beams" Kubo,Mtingwa,Wolski.
!
!  This formulation is one of the slowest methods for calculating IBS rates.
!
!  rates returns betatron growth rates.  Multiply by two to get transverse emittance growth rates.
!
!  Input:
!    ele               - ele_struct: contains Twiss parameters used in IBS formula
!    coulomb_log       - real(rp): Coulomb log value
!    n_part            - real(rp): number of particles in the bunch.
!  Output:
!    rates             - ibs_struct
!         %inv_Ta      - real(rp): a-mode betatron growth rate.
!         %inv_Tb      - real(rp): b-mode betatron growth rate.
!         %inv_Tz      - real(rp): energy spread growth rate.
!-

subroutine bjmt1(ele, coulomb_log, rates, n_part)

  implicit none

  type(ele_struct) :: ele
  real(rp) coulomb_log
  type(ibs_struct), intent(out) :: rates
  type(ele_struct), target :: stubele
  real(rp) n_part

  real(rp) sigma_p, emit_a, emit_b, sigma_z, energy
  real(rp) classical_radius
  real(rp) gamma, KE, rbeta, beta_a, beta_b
  real(rp) sigma_y
  real(rp) Dx, Dy, Dxp, Dyp
  real(rp) alpha_a, alpha_b
  real(rp) big_A
  real(rp) Hx, Hy
  real(rp) inv_Tz, inv_Ta, inv_Tb

  real(rp) phi_h, phi_v

  real(c_double) :: Lp(3,3), Lh(3,3), Lv(3,3), L(3,3)
  real(c_double), target :: mats(2,3,3)

  real(fgsl_double), target :: Elpha
  type(fgsl_function) :: integrand_ready
  real(fgsl_double) :: integration_result
  real(fgsl_double) :: abserr
  type(c_ptr) :: ptr
  type(fgsl_integration_workspace) :: integ_wk
  integer(fgsl_int) :: fgsl_status

  energy = ele%value(E_TOT$)
  call convert_total_energy_to(energy, ele%branch%param%particle, gamma, KE, rbeta)

  sigma_p = ele%z%sigma_p
  sigma_z = ele%z%sigma
  emit_a = ele%a%emit
  emit_b = ele%b%emit

  classical_radius = c_light*c_light*e_charge*1.0d-7*abs(charge_of(ele%branch%param%particle))/mass_of(ele%branch%param%particle)

  big_A=(classical_radius**2)*c_light*n_part/64.0/(pi**2)/(rbeta**3)/(gamma**4)/emit_a/emit_b/sigma_z/sigma_p

  alpha_a = ele%a%alpha
  alpha_b = ele%b%alpha
  beta_a = ele%a%beta
  beta_b = ele%b%beta
  Dx = ele%a%eta
  Dy = ele%b%eta
  Dxp = ele%a%etap
  Dyp = ele%b%etap
  sigma_y = sqrt(beta_b*emit_b + (Dy*sigma_p)**2)

  Hx = ( Dx**2 + (beta_a*Dxp + alpha_a*Dx)**2 ) / beta_a
  Hy = ( Dy**2 + (beta_b*Dyp + alpha_b*Dy)**2 ) / beta_b

  phi_h = Dxp + alpha_a*Dx/beta_a
  phi_v = Dyp + alpha_b*Dy/beta_b

  Lp = 0.0
  Lp(2,2) = 1.0
  Lp = (gamma**2)/(sigma_p**2)*Lp

  Lh = 0.0
  Lh(1,1) = 1.0
  Lh(1,2) = -1.0*gamma*phi_h
  Lh(2,1) = -1.0*gamma*phi_h
  Lh(2,2) = (gamma**2)*Hx/beta_a
  Lh = beta_a/emit_a*Lh

  Lv = 0.0
  Lv(2,2) = (gamma**2)*Hy/beta_b
  Lv(2,3) = -1.0*gamma*phi_v
  Lv(3,2) = -1.0*gamma*phi_v
  Lv(3,3) = 1.0
  Lv = beta_b/emit_b*Lv

  L = Lp + Lh + Lv

  mats(1,:,:) = L

  integ_wk = fgsl_integration_workspace_alloc(limit)
  ptr = c_loc(mats)
  integrand_ready = fgsl_function_init(bjmt_integrand, ptr)

  mats(2,:,:) = Lp
  fgsl_status = fgsl_integration_qag(integrand_ready, 0.0d0, 100.0d0, eps7, eps7, &
                                     limit, 3, integ_wk, integration_result, abserr)
  inv_Tz=4.0*pi*big_A*coulomb_log*integration_result

  mats(2,:,:) = Lh
  fgsl_status = fgsl_integration_qag(integrand_ready, 0.0d0, 100.0d0, eps7, eps7, &
                                     limit, 3, integ_wk, integration_result, abserr)
  inv_Ta=4.0*pi*big_A*coulomb_log*integration_result

  mats(2,:,:) = Lv
  fgsl_status = fgsl_integration_qag(integrand_ready, 0.0d0, 100.0d0, eps7, eps7, &
                                     limit, 3, integ_wk, integration_result, abserr)
  inv_Tb=4.0*pi*big_A*coulomb_log*integration_result

  call fgsl_integration_workspace_free(integ_wk)
  call fgsl_function_free(integrand_ready)

  rates%inv_Tz = inv_Tz
  rates%inv_Ta = inv_Ta
  rates%inv_Tb = inv_Tb

end subroutine bjmt1

function bjmt_integrand(xx,params) bind(c)
  real(c_double), value :: xx
  real(c_double) :: x
  type(c_ptr), value :: params
  real(c_double) :: bjmt_integrand

  real(c_double), pointer :: mats(:,:,:)

  real(c_double) :: im(3,3)
  real(c_double) :: L(3,3)
  real(c_double) :: Li(3,3)
  real(c_double) :: Li_im(3,3)

  real(c_double) :: Tr_Li
  real(c_double) :: Tr_im
  real(c_double) :: Tr_Li_im
  real(c_double) :: Det_L_l

  call c_f_pointer(params,mats,[2,3,3])

  L = mats(1,:,:)
  Li = mats(2,:,:)

  x = exp(xx) !change of variables

  Tr_Li = Li(1,1) + Li(2,2) + Li(3,3)

  im(1,1) = (L(2,2)+x)*(L(3,3)+x)-L(2,3)*L(3,2)
  im(1,2) = -L(1,2)*(L(3,3)+x)
  im(1,3) = L(1,2)*L(2,3)
  im(2,1) = -L(2,1)*(L(3,3)+x)
  im(2,2) = (L(1,1)+x)*(L(3,3)+x)
  im(2,3) = -(L(1,1)+x)*L(2,3)
  im(3,1) = L(2,3)*L(3,2)
  im(3,2) = -(L(1,1)+x)*L(3,2)
  im(3,3) = (L(1,1)+x)*(L(2,2)+x)-L(1,2)*L(2,1)

  im(:,:) = im(:,:) / ( (L(1,1)+x)*(L(2,2)+x)*(L(3,3)+x) - (L(1,1)+x)*L(3,2)*L(2,3) - L(1,2)*L(2,1)*(L(3,3)+x) )

  Tr_im = im(1,1)+im(2,2)+im(3,3)

  Li_im = matmul(Li,im)
  Tr_Li_im = Li_im(1,1)+Li_im(2,2)+Li_im(3,3)

  Det_L_l = (L(1,1)+x)*( (L(2,2)+x)*(L(3,3)+x) - L(3,2)*L(2,3) ) - L(1,2)*L(2,1)*(L(3,3)+x)

  bjmt_integrand = sqrt(x/Det_L_l) * ( Tr_Li*Tr_im - 3.0_rp*Tr_Li_im )
  bjmt_integrand = bjmt_integrand * x  !COV
end function bjmt_integrand

!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!+
!  subroutine bane1(ele, coulomb_log, rates, n_part)
!
!  This is an implementation of equations 10-15 from "Intrabeam
!  scattering formulas for high energy beams" Kubo,Mtingwa,Wolski.
!  It is a high energy approximation of the Bjorken-Mtingwa IBS
!  formulation.
!
!  Input:
!    ele               - ele_struct: contains Twiss parameters used in IBS formula
!    coulomb_log       - real(rp): Coulomb log
!    n_part            - real(rp): number of particles in the bunch.
!  Output:
!    rates             - ibs_struct
!         %inv_Ta      - real(rp): a-mode betatron growth rate.
!         %inv_Tb      - real(rp): b-mode betatron growth rate.
!         %inv_Tz      - real(rp): energy spread growth rate.
!-

subroutine bane1(ele, coulomb_log, rates, n_part)

  implicit none

  type(ele_struct) :: ele
  real(rp) coulomb_log
  type(ibs_struct), intent(out) :: rates

  real(rp) sigma_p, emit_a, emit_b, sigma_z, energy
  real(rp) classical_radius
  real(rp) gamma, KE, rbeta, beta_a, beta_b
  real(rp) sigma_b, sigma_b_beta
  real(rp) Da, Db, Dap, Dbp
  real(rp) alpha_a, alpha_b
  real(rp) a, b, g_bane
  real(rp) n_part, big_A
  real(rp) sigma_H, Ha, Hb
  real(rp) inv_Tz, inv_Ta, inv_Tb

  real(fgsl_double), target :: Elpha
  type(fgsl_function) :: integrand_ready
  real(fgsl_double) :: integration_result
  real(fgsl_double) :: abserr
  type(c_ptr) :: ptr
  type(fgsl_integration_workspace) :: integ_wk
  integer(fgsl_int) :: fgsl_status

  energy = ele%value(E_TOT$)
  call convert_total_energy_to(energy, ele%branch%param%particle, gamma, KE, rbeta)

  sigma_p = ele%z%sigma_p
  sigma_z = ele%z%sigma
  emit_a = ele%a%emit
  emit_b = ele%b%emit

  classical_radius = c_light*c_light*e_charge*1.0d-7*abs(charge_of(ele%branch%param%particle))/mass_of(ele%branch%param%particle)

  big_A=(classical_radius**2)*c_light*n_part/16.0/(gamma**3)/(emit_a**(3./4.))/(emit_b**(3./4.))/sigma_z/(sigma_p**3)

  beta_a = ele%a%beta
  beta_b = ele%b%beta
  alpha_a = ele%a%alpha
  alpha_b = ele%b%alpha
  Dap = ele%a%etap
  Dbp = ele%b%etap
  Da = ele%a%eta
  Db = ele%b%eta
  sigma_b_beta = sqrt(beta_b * emit_b)
  sigma_b = sqrt(sigma_b_beta**2 + (Db**2)*(sigma_p**2))

  Ha = ( (Da**2) + (beta_a*Dap + alpha_a*Da)**2 ) / beta_a
  Hb = ( (Db**2) + (beta_b*Dbp + alpha_b*Db)**2 ) / beta_b
  sigma_H = 1.0/sqrt(1.0/(sigma_p**2) + Ha/emit_a + Hb/emit_b)

  a = sigma_H/gamma*sqrt(beta_a/emit_a)
  b = sigma_H/gamma*sqrt(beta_b/emit_b)

  Elpha = a/b
  ptr = c_loc(Elpha)
  integ_wk = fgsl_integration_workspace_alloc(limit)
  integrand_ready = fgsl_function_init(integrand, ptr)
  fgsl_status = fgsl_integration_qagiu(integrand_ready, 0.0d0, eps7, eps7, &
                                       limit, integ_wk, integration_result, abserr)
  g_bane = 2.0d0 * sqrt(Elpha) / pi * integration_result

  call fgsl_function_free(integrand_ready)
  call fgsl_integration_workspace_free(integ_wk)

  inv_Tz = big_A*coulomb_log*sigma_H*g_bane*((beta_a*beta_b)**(-1./4.))
  inv_Ta = (sigma_p**2)*Ha/emit_a*inv_Tz
  inv_Tb = (sigma_p**2)*Hb/emit_b*inv_Tz

  rates%inv_Tz = inv_Tz
  rates%inv_Ta = inv_Ta
  rates%inv_Tb = inv_Tb

end subroutine bane1

function integrand(x, params) bind(c)
  real(c_double), value :: x
  type(c_ptr), value :: params
  real(c_double) :: integrand

  real(c_double), pointer :: alpha
  call c_f_pointer(params, alpha)

  integrand = 1.0_c_double / ( sqrt(1+x*x)*sqrt(alpha*alpha+x*x) )
end function integrand

!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!+
!  subroutine mpxx1(ele, coulomb_log, rates, n_part)
!
!  Modified Piwinski, further modified to treat Coulomb Log
!  in the same manner as Bjorken-Mtingwa, CIMP, Bane, Kubo & Oide, etc.
!  This formula is derived in Section 2.8.4 of Michael Ehrlichman's Graduate Thesis.
!
!  Input:
!    ele               - ele_struct: contains Twiss parameters used in IBS formula
!    coulomb_log       - real(rp): Coulomb log
!    n_part            - real(rp): number of particles in the bunch.
!
!  Output:
!    rates             - ibs_struct
!         %inv_Ta      - real(rp): a-mode betatron growth rate.
!         %inv_Tb      - real(rp): b-mode betatron growth rate.
!         %inv_Tz      - real(rp): energy spread growth rate.
!-

subroutine mpxx1(ele, coulomb_log, rates, n_part)

  implicit none

  type(ele_struct) :: ele
  real(rp) coulomb_log
  type(ibs_struct), intent(out) :: rates

  real(rp) sigma_p, emit_a, emit_b, sigma_z, energy
  real(rp) classical_radius
  real(rp) gamma, KE, rbeta, beta_a, beta_b
  real(rp) sigma_a, sigma_b, sigma_a_beta, sigma_b_beta
  real(rp) Da, Db, Dap, Dbp
  real(rp) alpha_a, alpha_b
  real(rp) a,b,q
  real(rp) n_part, big_A
  real(rp) sigma_H, Ha, Hb
  real(rp) inv_Tz, inv_Ta, inv_Tb
  real(rp) fab, f1b, f1a

  real(fgsl_double), target :: args(1:2)
  type(fgsl_function) :: integrand_ready
  real(fgsl_double) :: integration_result
  real(fgsl_double) :: abserr
  type(c_ptr) :: ptr
  type(fgsl_integration_workspace) :: integ_wk
  integer(fgsl_int) :: fgsl_status

  energy = ele%value(E_TOT$)
  call convert_total_energy_to(energy, ele%branch%param%particle, gamma, KE, rbeta)

  sigma_p = ele%z%sigma_p
  sigma_z = ele%z%sigma
  emit_a = ele%a%emit
  emit_b = ele%b%emit

  classical_radius = c_light*c_light*e_charge*1.0d-7*abs(charge_of(ele%branch%param%particle))/mass_of(ele%branch%param%particle)
  big_A=(classical_radius**2)*c_light*n_part/64.0/(pi**2)/(rbeta**3)/(gamma**4)/emit_a/emit_b/sigma_z/sigma_p

  alpha_a = ele%a%alpha
  alpha_b = ele%b%alpha
  beta_a = ele%a%beta
  beta_b = ele%b%beta
  sigma_a_beta = sqrt(beta_a * emit_a)
  sigma_b_beta = sqrt(beta_b * emit_b)
  Da = ele%a%eta
  Db = ele%b%eta
  Dap = ele%a%etap
  Dbp = ele%b%etap
  sigma_a = sqrt(sigma_a_beta**2 + (Da**2)*(sigma_p**2))
  sigma_b = sqrt(sigma_b_beta**2 + (Db**2)*(sigma_p**2))

  Ha = ( Da**2 + (beta_a*Dap + alpha_a*Da)**2 ) / beta_a
  Hb = ( Db**2 + (beta_b*Dbp + alpha_b*Db)**2 ) / beta_b

  sigma_H = 1.0/sqrt( 1.0/(sigma_p**2)+ Ha/emit_a + Hb/emit_b )

  a = sigma_H/gamma*sqrt(beta_a/emit_a)
  b = sigma_H/gamma*sqrt(beta_b/emit_b)

  !------------------------Begin calls to GSL integrator
  integ_wk = fgsl_integration_workspace_alloc(limit)
  ptr = c_loc(args)
  integrand_ready = fgsl_function_init(mpxx_integrand, ptr)

  args = (/a,b/)
  fgsl_status = fgsl_integration_qag(integrand_ready, 0.0d0, 1.0d0, eps7, eps7, &
                                     limit, 3, integ_wk, integration_result, abserr)
  fab = integration_result

  args = (/ 1.0_rp/a, b/a /) 
  fgsl_status = fgsl_integration_qag(integrand_ready, 0.0d0, 1.0d0, eps7, eps7, &
                                     limit, 3, integ_wk, integration_result, abserr)
  f1b = integration_result

  args = (/ 1.0_rp/b, a/b /)
  fgsl_status = fgsl_integration_qag(integrand_ready, 0.0d0, 1.0d0, eps7, eps7, &
                                     limit, 3, integ_wk, integration_result, abserr)
  f1a = integration_result

  call fgsl_integration_workspace_free(integ_wk)
  call fgsl_function_free(integrand_ready)
  !------------------------End calls to GSL integrator

  inv_Tz = coulomb_log * big_A * sigma_H**2 / sigma_p**2 * fab
  inv_Ta = coulomb_log * big_A * (f1b + Ha*sigma_H**2/emit_a*fab)
  inv_Tb = coulomb_log * big_A * (f1a + Hb*sigma_H**2/emit_b*fab)

  rates%inv_Tz = inv_Tz
  rates%inv_Ta = inv_Ta
  rates%inv_Tb = inv_Tb

end subroutine mpxx1

function mpxx_integrand(x, params) bind(c)
    real(c_double), value :: x
    type(c_ptr), value :: params
    real(c_double) :: mpxx_integrand

    real(c_double), pointer :: args(:)
    real(c_double) av, bv

    real(c_double) u

    call c_f_pointer(params,args,[2])
    av = args(1)
    bv = args(2)

    mpxx_integrand = 8.0d0*pi*(1-3.0d0*x*x) / &
                     sqrt(av*av+(1-av*av)*x*x) / &
                     sqrt(bv*bv+(1-bv*bv)*x*x)
end function mpxx_integrand

!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!+
!  subroutine mpzt1(ele, coulomb_log, rates, n_part)
!
!  Modified Piwinski with Zotter's integral.  This is Piwinski's original derivation,
!  generalized to take the derivatives of the optics functions.  Also, Piwinski's
!  original cumbersome triple integral is reaplaced by Zotter's single integral.  Zotter's
!  integral is exact, and not an approximation.
!
!  rates returns betatron growth rates.  Multiply by two to get transverse emittance growth rates.
!
!  Input:
!    ele               - ele_struct: contains Twiss parameters used in IBS formula
!    coulomb_log       - real(rp): Coulomb log
!    n_part            - real(rp): number of particles in the bunch.
!
!  Output:
!    rates             - ibs_struct
!         %inv_Ta      - real(rp): a-mode betatron growth rate.
!         %inv_Tb      - real(rp): b-mode betatron growth rate.
!         %inv_Tz      - real(rp): energy spread growth rate.
!-

subroutine mpzt1(ele, coulomb_log, rates, n_part)

  implicit none

  type(ele_struct) :: ele
  real(rp) coulomb_log
  type(ibs_struct), intent(out) :: rates

  real(rp) sigma_p, emit_a, emit_b, sigma_z, energy
  real(rp) classical_radius
  real(rp) gamma, KE, rbeta, beta_a, beta_b
  real(rp) sigma_a, sigma_b, sigma_a_beta, sigma_b_beta
  real(rp) Da, Db, Dap, Dbp
  real(rp) alpha_a, alpha_b
  real(rp) a,b,q
  real(rp) n_part, big_A
  real(rp) sigma_H, Ha, Hb
  real(rp) inv_Tz, inv_Ta, inv_Tb
  real(rp) fabq, f1bq, f1aq

  real(fgsl_double), target :: args(1:3)
  type(fgsl_function) :: integrand_ready
  real(fgsl_double) :: integration_result
  real(fgsl_double) :: abserr
  type(c_ptr) :: ptr
  type(fgsl_integration_workspace) :: integ_wk
  integer(fgsl_int) :: fgsl_status

  energy = ele%value(E_TOT$)
  call convert_total_energy_to(energy, ele%branch%param%particle, gamma, KE, rbeta)

  sigma_p = ele%z%sigma_p
  sigma_z = ele%z%sigma
  emit_a = ele%a%emit
  emit_b = ele%b%emit

  classical_radius = c_light*c_light*e_charge*1.0d-7*abs(charge_of(ele%branch%param%particle))/mass_of(ele%branch%param%particle)
  big_A=(classical_radius**2)*c_light*n_part/64.0/(pi**2)/(rbeta**3)/(gamma**4)/emit_a/emit_b/sigma_z/sigma_p

  alpha_a = ele%a%alpha
  alpha_b = ele%b%alpha
  beta_a = ele%a%beta
  beta_b = ele%b%beta
  sigma_a_beta = sqrt(beta_a * emit_a)
  sigma_b_beta = sqrt(beta_b * emit_b)
  Da = ele%a%eta
  Db = ele%b%eta
  Dap = ele%a%etap
  Dbp = ele%b%etap
  sigma_a = sqrt(sigma_a_beta**2 + (Da**2)*(sigma_p**2))
  sigma_b = sqrt(sigma_b_beta**2 + (Db**2)*(sigma_p**2))

  Ha = ( Da**2 + (beta_a*Dap + alpha_a*Da)**2 ) / beta_a
  Hb = ( Db**2 + (beta_b*Dbp + alpha_b*Db)**2 ) / beta_b

  sigma_H = 1.0/sqrt( 1.0/(sigma_p**2)+ Ha/emit_a + Hb/emit_b )

  a = sigma_H/gamma*sqrt(beta_a/emit_a)
  b = sigma_H/gamma*sqrt(beta_b/emit_b)
  q = sigma_H*rbeta*sqrt(2.0_rp*sigma_b/classical_radius)
  !---- q = (gamma**2)*sigma_b*emit_a/classical_radius/beta_a  !effective coulomb log 

  !------------------------Begin calls to GSL integrator
  integ_wk = fgsl_integration_workspace_alloc(limit)
  ptr = c_loc(args)
  integrand_ready = fgsl_function_init(zot_integrand, ptr)

  args = (/a,b,q/)
  fgsl_status = fgsl_integration_qag(integrand_ready, 0.0d0, 1.0d0, eps7, eps7, &
                                     limit, 3, integ_wk, integration_result, abserr)
  fabq = integration_result

  args = (/ 1.0_rp/a, b/a, q/a /) 
  fgsl_status = fgsl_integration_qag(integrand_ready, 0.0d0, 1.0d0, eps7, eps7, &
                                     limit, 3, integ_wk, integration_result, abserr)
  f1bq = integration_result

  args = (/ 1.0_rp/b, a/b, q/b /)
  fgsl_status = fgsl_integration_qag(integrand_ready, 0.0d0, 1.0d0, eps7, eps7, &
                                     limit, 3, integ_wk, integration_result, abserr)
  f1aq = integration_result

  call fgsl_integration_workspace_free(integ_wk)
  call fgsl_function_free(integrand_ready)
  !------------------------End calls to GSL integrator

  inv_Tz = big_A * sigma_H**2 / sigma_p**2 * fabq
  inv_Ta = big_A * (f1bq + Ha*sigma_H**2/emit_a*fabq)
  inv_Tb = big_A * (f1aq + Hb*sigma_H**2/emit_b*fabq)

  rates%inv_Tz = inv_Tz
  rates%inv_Ta = inv_Ta
  rates%inv_Tb = inv_Tb

end subroutine mpzt1

function zot_integrand(x, params) bind(c)
  real(c_double), value :: x
  type(c_ptr), value :: params
  real(c_double) :: zot_integrand

  real(c_double), pointer :: args(:)
  real(c_double) av, bv, qv, P, Q

  real(c_double) u

  call c_f_pointer(params,args,[3])
  av = args(1)
  bv = args(2)
  qv = args(3)

  !COV to remove singularity at endpoints
  u=30.*((x**5.)/5. - (x**4.)/2. + (x**3.)/3.)

  P = sqrt(av*av + (1.-av*av)*u*u)
  Q = sqrt(bv*bv + (1.-bv*bv)*u*u)
  zot_integrand = 8.0_rp*pi * (1.0_rp-3.0_rp*u*u)/P/Q * &
              (2.0_rp*log(qv/2.0_rp*(1/P+1/Q)) - 0.577215665)

  zot_integrand = zot_integrand * 30.*(x**4. - 2.*(x**3.) + x**2.)  !COV to remove singularity

  !without COV
  !P = sqrt(av*av + (1-av*av)*x*x)
  !Q = sqrt(bv*bv + (1-bv*bv)*x*x)
  !zot_integrand = 8.0_rp*pi * (1.0_rp-3.0_rp*x*x)/P/Q * &
  !            (2.0_rp*log(qv/2.0_rp*(1/P+1/Q)) - 0.577215665)

end function zot_integrand

!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!+
!  subroutine cimp1(ele, coulomb_log, rates, n_part)
!
!  This is an implementation of equations 34,38-40 from "Intrabeam
!  scattering formulas for high energy beams" Kubo,Mtingwa,Wolski.
!  It is a modified version of the Piwinski IBS formulation.
!  The integral (34) is handled with a piecewise interpolation generated
!  in mathematica.  The interpolation is accurate beyond 1% through it's
!  effective range (.0001 - 3000).
!
!  This is the quickest of the three IBS formuations in this module. 
!
!  rates returns betatron growth rates.  Multiply by two to get transverse emittance growth rates.
!
!  Input:
!    ele               - ele_struct: contains Twiss parameters used in IBS formula
!    coulomb_log       - real(rp): Coulomb log
!    n_part            - real(rp): number of particles in the bunch.
!
!  Output:
!    rates             - ibs_struct
!         %inv_Ta      - real(rp): a-mode betatron growth rate.
!         %inv_Tb      - real(rp): b-mode betatron growth rate.
!         %inv_Tz      - real(rp): energy spread growth rate.
!-

subroutine cimp1(ele, coulomb_log, rates, n_part)

  implicit none

  type(ele_struct) :: ele
  real(rp) coulomb_log
  type(ibs_struct), intent(out) :: rates
  real(rp) element_length, E_TOT, n_part

  real(rp) sigma_p, emit_a, emit_b, sigma_z
  real(rp) classical_radius
  real(rp) gamma, KE, rbeta, beta_a, beta_b
  real(rp) sigma_x, sigma_y, sigma_x_beta, sigma_y_beta
  real(rp) Dx, Dy, Dxp, Dyp
  real(rp) alpha_a, alpha_b
  real(rp) a, b
  real(rp) big_A
  real(rp) sigma_H, Hx, Hy
  real(rp) inv_Tz, inv_Ta, inv_Tb
  real(rp) g_ab,g_ba
  real(rp) bminstar, bmax
  real(rp) energy

  !- Code specific to alternitive log representation
  real(rp) q, lnqa, lnqb
  !-

  energy = ele%value(E_TOT$)

  call convert_total_energy_to(energy, ele%branch%param%particle, gamma, KE, rbeta)

  sigma_p = ele%z%sigma_p
  sigma_z = ele%z%sigma
  emit_a = ele%a%emit
  emit_b = ele%b%emit

  classical_radius = c_light*c_light*e_charge*1.0d-7*abs(charge_of(ele%branch%param%particle))/mass_of(ele%branch%param%particle)
  big_A=(classical_radius**2)*c_light*n_part/64.0/(pi**2)/(rbeta**3)/(gamma**4)/emit_a/emit_b/sigma_z/sigma_p

  alpha_a = ele%a%alpha
  alpha_b = ele%b%alpha
  beta_a = ele%a%beta
  beta_b = ele%b%beta
  sigma_x_beta = sqrt(beta_a * emit_a)
  sigma_y_beta = sqrt(beta_b * emit_b)
  Dx = ele%a%eta
  Dy = ele%b%eta
  Dxp = ele%a%etap
  Dyp = ele%b%etap
  sigma_x = sqrt(sigma_x_beta**2 + (Dx**2)*(sigma_p**2))
  sigma_y = sqrt(sigma_y_beta**2 + (Dy**2)*(sigma_p**2))

  Hx = ( Dx**2 + (beta_a*Dxp + alpha_a*Dx)**2 ) / beta_a
  Hy = ( Dy**2 + (beta_b*Dyp + alpha_b*Dy)**2 ) / beta_b

  sigma_H = 1.0/sqrt( 1.0/(sigma_p**2)+ Hx/emit_a + Hy/emit_b )

  a = sigma_H/gamma*sqrt(beta_a/emit_a)
  b = sigma_H/gamma*sqrt(beta_b/emit_b)

  g_ba = g(b/a)
  g_ab = g(a/b)

  inv_Tz = 2.*(pi**(3./2.))*big_A*(sigma_H**2)/(sigma_p**2) * &
    coulomb_log * ( g_ba/a + g_ab/b ) 
  inv_Ta = 2.*(pi**(3./2.))*big_A*coulomb_log*&
    (-a*g_ba + Hx*(sigma_H**2)/emit_a* &
    ( g_ba/a + g_ab/b ) ) 
  inv_Tb = 2.*(pi**(3./2.))*big_A*coulomb_log*&
    (-b*g_ab + Hy*(sigma_H**2)/emit_b* &
    ( g_ba/a + g_ab/b ) )

  rates%inv_Tz = inv_Tz
  rates%inv_Ta = inv_Ta
  rates%inv_Tb = inv_Tb

  contains
    !+
    !  function g(u)
    !
    !  This is an 13-degree piecewise polynomial interpolation of the
    !  integral for the CIMP ibs formulation (equation 34 in "Intrabeam 
    !  Scattering Formulas for High Energy Beams").
    !  The segments were generated in mathematica and each is accurate
    !  to beyond 1%.  The effective range of this interpolation is
    !  .0001 to 3000, and it prints an error message when that range
    !  is exceeded.
    !-

    function g(u)
                                                                                       
      real(rp), intent(in) :: u
      real(rp) :: g
      real(rp), dimension(0:13) ::  o,p,pa,pb,qb,qa,q,ra,rb,r

      o = (/ 0.27537408880308967,-0.0014411559647655005,4.760237316963749E-6,-1.056200254057877E-8, &
           1.652957676252096E-11,-1.8807607644579788E-14,1.582497859744359E-17,-9.911417850544026E-21, &
           4.606212267581107E-24,-1.5661683056801226E-27,3.7836532348260456E-31,-6.148053455724941E-35, &
           6.021954146486214E-29,-2.6856345839264216E-33 /) !last two terms have E-10 multiplied in function

      p = (/ 1.3728514796955633,-0.06329159139596323,0.0019843996282577487,-0.00004282113130209425, &
           6.588299268041797E-7,-7.411644739022654E-9,6.186473337975396E-11,-3.851874112809632E-13, &
           1.7820926145729545E-15,-6.038143148108221E-18,1.454670683108795E-20,-2.358356775328005E-23, &
           2.305700018891372E-26,-1.026694954987225E-29 /)

      pa = (/ 2.4560269972177395,-0.34483385959575796,0.04060452194930618,-0.0036548931452558006, &
            0.00024651396218017436,-0.000012205596049761074,4.224156486889531E-7,-8.746709793365587E-9, &
            2.0010958003696878E-11,5.293109492480939E-12,-1.8692549451553078E-13,3.313011636321362E-15, &
            -3.191257818162888E-17,1.333629848206322E-19 /)

      pb = (/ 0.2509305915339517, 3.655027986910339, -3.538178517122364, 2.0378515154003622, &
            -0.8131215410970869, 0.23600479506508648, -0.0508941392205539, 0.008214142113823073, &
            -0.000988882504551898, 0.00008751565040446301, -5.526659141312207E-6, &
            2.3564088895171855E-7, -6.0770418816204626E-9, 7.15762060347684E-11 /)

      qb = (/ -4.5128859976303835, 44.038623289844345, -201.50265131380843, 685.7422921643122, &
            -1733.9460775721304, 3272.7093350374225, -4637.724347309373, 4940.122242846527, &
            -3931.6350463124686, 2301.047404938234, -960.7211908043258, 270.6618330819612, &
            -46.088423620894815, 3.581345541476106 /)

      qa= (/ -7.9542527228302164, 207.51519155414007, -4406.431743148087, 75592.0699189963, &
           -992747.0161009694, 9.983494244578896E6, -7.72264540055081E7, 4.593887542691004E8, &
           -2.0856924957059484E9, 7.103526742990105E9, -1.758262228582309E10, 2.9881056304916122E10, &
           -3.11970246734184E10, 1.5092353688380434E10 /)

      q = (/ -10.427660034546474, 614.5741414073021, -37752.0301748355, 1.8467024380709291E6, &
           -6.787056202900247E7, 1.87800146779508E9, -3.9344782189063065E10, 6.245185690585099E11, &
           -7.460242970344108E12, 6.596837840619067E13, -4.186282233147678E14, 1.8023052362000775E15, &
           -4.713009373619726E15, 5.649407049035322E15 /)

      ra = (/ -13.151132987122772, 2076.0767607226853, -439312.2888781302, 7.549171067948443E7, &
            -9.918579267112797E9, 9.97612075477821E11, -7.717556597617397E13, 4.591087739048907E15, &
            -2.0844941927298595E17, 7.099636713666537E18, -1.757338267407412E20, 2.986592440773308E21, &
            -3.1181759869417024E22, 1.5085206738847978E23 /)

      rb = (/ -16.158250737444753, 7653.567961683397, -5.677248220833501E6, 3.2736621037820387E9, &
            -1.3885723919300298E12, 4.356123741560311E14, -1.018979723977809E17, 1.782083486703034E19, &
            -2.31837909173977E21, 2.209799968074478E23, -1.4978388984105023E25, 6.831861840241021E26, &
            -1.878873461985899E28, 2.3529610201488045E29 /)

      r = (/ -21.460085854969584, 80412.06002488811, -6.294777676054858E8, 3.84371292991909E12, &
            -1.7311478313121348E16, 5.7791722579515556E19, -1.4411700564165495E23, +2.6910329100895636E26, &
            -3.742607580627404E29, 3.8178627111416235E32, -2.7722017774216816E35, 1.3556698189607729E38, &
            -4.000249119089512E30, 5.378513442816346E32 /) ! E10 added to these last two in formula

      if    (u .gt. 3000.0) then
        write(*,*) "CRITICAL WARNING: interpolation range exceeded"
        g = 0.
      elseif(u .gt. 300.0) then
        g = o(0)+o(1)*u+o(2)*(u**2)+o(3)*(u**3)+o(4)*(u**4)+o(5)*(u**5)+ &
            o(6)*(u**6)+o(7)*(u**7)+o(8)*(u**8)+o(9)*(u**9)+o(10)*(u**10)+ &
            o(11)*(u**11)+o(12)*(1.0E-10)*(u**12)+o(13)*(1.0E-10)*(u**13)
      elseif(u .gt. 30.0)  then
        g = p(0)+p(1)*u+p(2)*(u**2)+p(3)*(u**3)+p(4)*(u**4)+p(5)*(u**5)+ &
            p(6)*(u**6)+p(7)*(u**7)+p(8)*(u**8)+p(9)*(u**9)+p(10)*(u**10)+ &
            p(11)*(u**11)+p(12)*(u**12)+p(13)*(u**13)
      elseif(u .gt. 10.0)  then
        g = pa(0)+pa(1)*u+pa(2)*(u**2)+pa(3)*(u**3)+pa(4)*(u**4)+pa(5)*(u**5)+ &
            pa(6)*(u**6)+pa(7)*(u**7)+pa(8)*(u**8)+pa(9)*(u**9)+pa(10)*(u**10)+ &
            pa(11)*(u**11)+pa(12)*(u**12)+pa(13)*(u**13)
      elseif(u .gt. 1.6)   then
        g = pb(0)+pb(1)*u+pb(2)*(u**2)+pb(3)*(u**3)+pb(4)*(u**4)+pb(5)*(u**5)+ &
            pb(6)*(u**6)+pb(7)*(u**7)+pb(8)*(u**8)+pb(9)*(u**9)+pb(10)*(u**10)+ &
            pb(11)*(u**11)+pb(12)*(u**12)+pb(13)*(u**13)
      elseif(u .gt. .20)   then
        g = qb(0)+qb(1)*u+qb(2)*(u**2)+qb(3)*(u**3)+qb(4)*(u**4)+qb(5)*(u**5)+ &
            qb(6)*(u**6)+qb(7)*(u**7)+qb(8)*(u**8)+qb(9)*(u**9)+qb(10)*(u**10)+ &
            qb(11)*(u**11)+qb(12)*(u**12)+qb(13)*(u**13)
      elseif(u .gt. .10)   then
        g = qa(0)+qa(1)*u+qa(2)*(u**2)+qa(3)*(u**3)+qa(4)*(u**4)+qa(5)*(u**5)+ &
            qa(6)*(u**6)+qa(7)*(u**7)+qa(8)*(u**8)+qa(9)*(u**9)+qa(10)*(u**10)+ &
            qa(11)*(u**11)+qa(12)*(u**12)+qa(13)*(u**13)
      elseif(u .gt. .02)  then
        g = q(0)+q(1)*u+q(2)*(u**2)+q(3)*(u**3)+q(4)*(u**4)+q(5)*(u**5)+ &
            q(6)*(u**6)+q(7)*(u**7)+q(8)*(u**8)+q(9)*(u**9)+q(10)*(u**10)+ &
            q(11)*(u**11)+q(12)*(u**12)+q(13)*(u**13)
      elseif(u .gt. .01)  then
        g = ra(0)+ra(1)*u+ra(2)*(u**2)+ra(3)*(u**3)+ra(4)*(u**4)+ra(5)*(u**5)+ &
            ra(6)*(u**6)+ra(7)*(u**7)+ra(8)*(u**8)+ra(9)*(u**9)+ra(10)*(u**10)+ &
            ra(11)*(u**11)+ra(12)*(u**12)+ra(13)*(u**13)
      elseif(u .gt. .001)  then
        g = rb(0)+rb(1)*u+rb(2)*(u**2)+rb(3)*(u**3)+rb(4)*(u**4)+rb(5)*(u**5)+ &
            rb(6)*(u**6)+rb(7)*(u**7)+rb(8)*(u**8)+rb(9)*(u**9)+rb(10)*(u**10)+ &
            rb(11)*(u**11)+rb(12)*(u**12)+rb(13)*(u**13)
      elseif(u .gt. .0001) then
        g = r(0)+r(1)*u+r(2)*(u**2)+r(3)*(u**3)+r(4)*(u**4)+r(5)*(u**5)+ &
            r(6)*(u**6)+r(7)*(u**7)+r(8)*(u**8)+r(9)*(u**9)+r(10)*(u**10)+ &
            r(11)*(u**11)+r(12)*(1.0E10)*(u**12)+r(13)*(1.0E10)*(u**13)
      else
        write(*,*) "CRITICAL WARNING: interpolation range exceeded"
        g = 0.
      endif
    end function g
end subroutine cimp1

end module






