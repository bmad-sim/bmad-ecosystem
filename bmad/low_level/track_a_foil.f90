!+
! Subroutine track_a_foil (orbit, ele, param, mat6, make_matrix)
!
! Bmad_standard tracking through an foil element.
!
! From Eq. (6) in 
!   Approximations to Multiple Coulomb Scattering
!   Gerald R. Lynch and Orin Dahl
!   Nuclear Insturments and methods in Physics Research B58 (7991) 6-10.
!
! Input:
!   orbit       -- Coord_struct: Starting position.
!   ele         -- ele_struct: foil element.
!   param       -- lat_param_struct: Lattice parameters.
!   make_matrix -- logical, optional: Propagate the transfer matrix? Default is False.
!
! Output:
!   orbit      -- coord_struct: End position.
!   mat6(6,6)  -- real(rp), optional: Transfer matrix through the element.
!-

subroutine track_a_foil (orbit, ele, param, mat6, make_matrix)

use bmad_interface, except_dummy => track_a_foil
use random_mod
use xraylib, dummy => r_e

implicit none

type (coord_struct) :: orbit, orb0
type (ele_struct), target :: ele
type (lat_param_struct) :: param
type (material_struct), pointer :: material

real(rp), optional :: mat6(6,6)
real(rp) xx0, sigma, rnd(2), I_excite, mass_material, dE_dx_tot, E0, E1, p2, elec_area_density
real(rp) pc_old, pc_new, chi2_c, chi2_alpha, nu, f, ln_chi_alpha_sum, zza, zza_sum, norm_sum
real(rp) atomic_mass, omega, arg, area_density

integer i, j, ns, n_step, z_material, z_particle

logical, optional :: make_matrix

character(*), parameter :: r_name = 'track_a_foil'

! On target? If not, there is nothing to be done.

call offset_particle (ele, set$, orbit, set_hvkicks = .false., mat6 = mat6, make_matrix = make_matrix)

if (orbit%vec(1) < ele%value(x1_edge$) .or. orbit%vec(1) > ele%value(x2_edge$)) return
if (orbit%vec(3) < ele%value(y1_edge$) .or. orbit%vec(3) > ele%value(y2_edge$)) return

z_particle = nint(ele%value(final_charge$))
f = ele%value(f_factor$)

!

n_step = nint(ele%value(num_steps$))
do ns = 1, n_step
  pc_old = orbit%p0c * (1.0_rp + orbit%vec(6))

  ! Angle scatter

  if (nint(ele%value(scatter_method$)) /= off$) then
    xx0 = 0; chi2_c = 0; ln_chi_alpha_sum = 0; zza_sum = 0

    do i = 1, size(ele%foil%material)
      material => ele%foil%material(i)
      area_density = this_area_density(material, ele%value(thickness$), ele%value(dthickness_dx$), orbit)
      if (orbit%state == lost$) return
      z_material = atomic_number(material%species)

      select case (nint(ele%value(scatter_method$)))
      case (highland$)
        xx0 = xx0 + area_density / material%radiation_length_used

      case (lynch_dahl$)
        atomic_mass = mass_of(material%species) / atomic_mass_unit
        zza = z_material * (z_material + 1) * area_density / atomic_mass
        zza_sum = zza_sum + zza
        chi2_c = chi2_c + 0.157e11_rp * zza * (z_particle / (pc_old * orbit%beta))**2
        ln_chi_alpha_sum = ln_chi_alpha_sum + zza * log(sqrt(2.007e7_rp * z_material**(2.0/3.0) * &
                (1.0_rp + 3.34_rp * (z_material * z_particle * fine_structure_constant / orbit%beta)**2) / pc_old**2))
      end select
    enddo

    if (is_true(ele%value(scatter_test$))) then
      rnd = 1
    else
      call ran_gauss(rnd)
    endif

    select case (nint(ele%value(scatter_method$)))
    case (highland$)
      sigma = 13.6e6_rp * z_particle * sqrt(xx0) / (pc_old * orbit%beta) * &
                                          (1.0_rp + 0.038_rp * log(xx0*z_particle**2/orbit%beta**2))
    case (lynch_dahl$)
      chi2_alpha = (exp(ln_chi_alpha_sum/zza_sum))**2
      omega = chi2_c / (1.167 * chi2_alpha) 
      nu = 0.5_rp * omega / (1 - f)
      arg = chi2_c * ((1 + nu) * log(1 + nu) / nu - 1) / (1 + f**2)
      sigma = sqrt(max(0.0_rp, arg))
    end select
    
    sigma = sigma * pc_old / orbit%p0c
    orbit%vec(2) = orbit%vec(2) + rnd(1) * sigma
    orbit%vec(4) = orbit%vec(4) + rnd(2) * sigma
  endif

  ! Energy loss

  p2 = (pc_old / mass_of(orbit%species))**2  ! beta^2 / (1 - beta^2) term

  dE_dx_tot = 0
  do j = 1, size(ele%foil%material)
    material => ele%foil%material(j)
    area_density = this_area_density(material, ele%value(thickness$), ele%value(dthickness_dx$), orbit)
    if (orbit%state == lost$) return
    z_material = atomic_number(material%species)
    I_excite = mean_excitation_energy_over_z(z_material) * z_material
    mass_material = mass_of(material%species) / atomic_mass_unit
    ! number_electrons / m^2
    elec_area_density = (1.0e3_rp * N_avogadro) * area_density * z_material / mass_material
    dE_dx_tot = dE_dx_tot - fourpi * mass_of(electron$) * elec_area_density * (z_particle * r_e / orbit%beta)**2 * &
              (log(2 * mass_of(electron$) * p2 / I_excite) - orbit%beta**2)
  enddo

  E0 = orbit%p0c * (1.0_rp + orbit%vec(6)) / orbit%beta
  E1 = E0 + dE_dx_tot / n_step
  call convert_total_energy_to(E1, orbit%species, beta = orbit%beta, pc = pc_new)
  orbit%vec(2) = orbit%vec(2) * pc_new / pc_old   ! To keep the angle x' = constant
  orbit%vec(4) = orbit%vec(4) * pc_new / pc_old   ! To keep the angle y' = constant
  orbit%vec(6) = (pc_new - orbit%p0c) / orbit%p0c
enddo

! Charge

orbit%species = set_species_charge(orbit%species, z_particle)

!

call offset_particle (ele, unset$, orbit, set_hvkicks = .false., mat6 = mat6, make_matrix = make_matrix)

!------------------------------------------------------------------------------------------------------
contains

function this_area_density(material, thickness, dthickness_dx, orbit) result (area_density)

type (material_struct), pointer :: material
type (coord_struct) orbit

real(rp) thickness, dthickness_dx, area_density

!

if (thickness == 0 .and. dthickness_dx == 0) then
  area_density = material%area_density_used
else
  area_density = material%density_used * (thickness + dthickness_dx * orbit%vec(1))
endif

if (area_density < 0) then
  call out_io (s_error$, r_name, 'FOIL THICKNESS AT PARTICLE POSITION IS NEGATIVE FOR: ' // ele%name, &
                                 'CHECK YOUR SETTING OF THICKNESS AND DTHICKNESS_DX.', &
                                 'WILL MARK PARTICLE AS LOST.')
  orbit%state = lost$
  return
endif

end function this_area_density

end subroutine track_a_foil
