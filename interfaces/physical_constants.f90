!+
! Actually this module is more properly named constants since
! it evolved to define more than just physical constants.
!-

module physical_constants

use precision_def

type physical_const_struct
  character(40) name
  real(rp) value
end type

real(rp), parameter :: pi = 3.141592653589793238462643383279d0
real(rp), parameter :: twopi = 2 * pi
real(rp), parameter :: fourpi = 4 * pi
real(rp), parameter :: sqrt_2 = sqrt(2.0_rp)           ! 1.414213562373095048801688724209698d0
real(rp), parameter :: sqrt_3 = sqrt(3.0_rp)           ! 1.732050807568877293527446341505872d0 

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
real(rp), parameter :: anomalous_mag_moment_proton   = 1.79284735d0
real(rp), parameter :: anomalous_mag_moment_muon     = 1.1659208d-3  ! ~fine_structure_constant / twopi
real(rp), parameter :: anomalous_mag_moment_deuteron = -0.14298727047d0

! Should make physical_const_list "parameter" but there is a gcc bug (in Version 7.1 at least)
! where if you pass physical_const_list%name to a routine there will be a crash.

type (physical_const_struct) :: physical_const_list(32) = [ &
                 physical_const_struct('pi', pi), &
                 physical_const_struct('twopi', twopi), &
                 physical_const_struct('fourpi', fourpi), &
                 physical_const_struct('e_log', 2.718281828459_rp), &
                 physical_const_struct('sqrt_2', sqrt_2), &
                 physical_const_struct('degrad', 180 / pi), &
                 physical_const_struct('degrees', pi / 180), & ! From degrees to radians.
                 physical_const_struct('raddeg', pi / 180), &
                 physical_const_struct('m_electron', m_electron), &
                 physical_const_struct('m_muon', m_muon), &
                 physical_const_struct('m_pion_0', m_pion_0), &
                 physical_const_struct('m_pion_charged', m_pion_charged), &
                 physical_const_struct('m_proton', m_proton), &
                 physical_const_struct('m_deuteron', m_deuteron), &
                 physical_const_struct('c_light', c_light), &
                 physical_const_struct('r_e', r_e), &
                 physical_const_struct('r_p', r_p), &
                 physical_const_struct('e_charge', e_charge), &
                 physical_const_struct('h_planck', h_planck), &
                 physical_const_struct('h_bar_planck', h_bar_planck), &
                 physical_const_struct('pmass', p_mass), &
                 physical_const_struct('emass', e_mass), &
                 physical_const_struct('clight', c_light), &
                 physical_const_struct('fine_struct_const', fine_structure_constant), &
                 physical_const_struct('anom_mag_electron', anomalous_mag_moment_electron), &  ! Old style. Deprecated.
                 physical_const_struct('anom_mag_proton', anomalous_mag_moment_proton), &      ! Old style. Deprecated.
                 physical_const_struct('anom_mag_muon', anomalous_mag_moment_muon), &          ! Old style. Deprecated.
                 physical_const_struct('anom_mag_deuteron', anomalous_mag_moment_deuteron), &  ! Old style. Deprecated.
                 physical_const_struct('anom_moment_electron', anomalous_mag_moment_electron), &
                 physical_const_struct('anom_moment_proton', anomalous_mag_moment_proton), &
                 physical_const_struct('anom_moment_muon', anomalous_mag_moment_muon), &
                 physical_const_struct('anom_moment_deuteron', anomalous_mag_moment_deuteron)]


! This is to suppress the ranlib "has no symbols" message

integer, private :: private_dummy

end module
