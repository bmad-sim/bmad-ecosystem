!+
! Actually this module is more properly named constants since
! it evolved to define more than just physical constants.
!-

module physical_constants

use precision_def

real(rp), parameter :: pi = 3.141592653589793238462643383279d0
real(rp), parameter :: twopi = 2 * pi
real(rp), parameter :: fourpi = 4 * pi
real(rp), parameter :: sqrt_2 = 1.414213562373095048801688724209698d0
real(rp), parameter :: sqrt_3 = 1.732050807568877293527446341505872d0 

real(rp), parameter :: m_electron = 0.5109989461d6        ! Mass [eV]
real(rp), parameter :: m_proton   = 0.9382720813d9        ! Mass [eV]
real(rp), parameter :: m_muon     = 105.6583715d6        ! Mass [eV]

real(rp), parameter :: e_mass = 1d-9 * m_electron         ! [GeV] FOR MAD COMPATIBILITY USE ONLY. USE M_ELECTRON INSTEAD.
real(rp), parameter :: p_mass   = 1d-9 * m_proton         ! [GeV] FOR MAD COMPATIBILITY USE ONLY. USE M_PROTON INSTEAD.

real(rp), parameter :: m_pion_0 = 134.9766d6             ! Mass [eV]
real(rp), parameter :: m_pion_charged = 139.57018d6      ! Mass [eV]

real(rp), parameter :: m_deuteron   = 1.875612928d9      ! Mass [eV]

real(rp), parameter :: atomic_mass_unit = 931.494095d6  ! unified atomic mass unit u (or dalton) in [eV]
 
real(rp), parameter :: c_light = 2.99792458d8            ! speed of light
real(rp), parameter :: r_e = 2.8179403227d-15            ! classical electron radius
real(rp), parameter :: r_p = r_e * m_electron / m_proton ! proton radius
real(rp), parameter :: e_charge = 1.6021766208d-19          ! electron charge [Coul]
real(rp), parameter :: h_planck = 4.13566733d-15         ! Planck's constant [eV*sec]
real(rp), parameter :: h_bar_planck = 6.58211899d-16     ! h_planck/twopi [eV*sec]

real(rp), parameter :: mu_0_vac = fourpi * 1d-7                     ! Permeability of free space
real(rp), parameter :: eps_0_vac = 1 / (c_light*c_light * mu_0_vac) ! Permittivity of free space

! Radiation constants

real(rp), parameter :: classical_radius_factor = 1.439964416d-9  ! e^2 / (4 pi eps_0) [m*eV]
                                                                 !  = classical_radius * mass * c^2. 
                                                                 ! Is same for all particles of charge +/- 1.
! Chemistry

real(rp), parameter :: N_avogadro = 6.02214129d23    ! Number / mole

! Anomalous magnetic moment

real(rp), parameter :: fine_structure_constant =  7.29735257d-3
real(rp), parameter :: anomalous_mag_moment_electron = 1.159652193d-3
real(rp), parameter :: anomalous_mag_moment_proton   = 1.79285d0
real(rp), parameter :: anomalous_mag_moment_muon     = 1.1659208d-3  ! ~fine_structure_constant / twopi
real(rp), parameter :: anomalous_mag_moment_deuteron = -0.14298727047d0

! This is to suppress the ranlib "has no symbols" message

integer, private :: private_dummy

end module
