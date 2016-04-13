!+
! Actually this module is more properly named constants since
! it evolved to define more than just physical constants.
!-

module physical_constants

use precision_def

real(rp), parameter :: pi = 3.14159265358979d0
real(rp), parameter :: twopi = 2 * pi
real(rp), parameter :: fourpi = 4 * pi
real(rp), parameter :: sqrt_2 = 1.41421356237310d0
real(rp), parameter :: sqrt_3 = 1.73205080757d0

real(rp), parameter :: e_mass = 0.5109989461d-3           ! [GeV] FOR MAD COMPATIBILITY USE ONLY. USE M_ELECTRON INSTEAD.
real(rp), parameter :: p_mass   = 0.9382720813d0          ! [GeV] FOR MAD COMPATIBILITY USE ONLY. USE M_PROTON INSTEAD.

real(rp), parameter :: m_electron = 0.5109989461d6        ! Mass [eV]
real(rp), parameter :: m_proton   = 0.9382720813d9        ! Mass [eV]
real(rp), parameter :: m_muon     = 105.6583745d6        ! Mass [eV]

real(rp), parameter :: m_pion_0 = 134.9766d6             ! Mass [eV]
real(rp), parameter :: m_pion_charged = 139.57018d6      ! Mass [eV]

real(rp), parameter :: m_deuteron   = 1.875612928d9      ! Mass [eV]

real(rp), parameter :: atomic_mass_unit = 931.494095d6  ! unified atomic mass unit u (or dalton) in [eV]
 
real(rp), parameter :: c_light = 2.99792458d8            ! speed of light
real(rp), parameter :: r_e = 2.8179402894d-15            ! classical electron radius
real(rp), parameter :: r_p = r_e * m_electron / m_proton ! proton radius
real(rp), parameter :: e_charge = 1.6021892d-19          ! electron charge [Coul]
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

complex(rp), parameter :: i_imaginary = (0.0d0, 1.0d0)
  
! real_garbage$ and int_garbage$ can be used, for example, to identify
! variable that have not been set.

integer, parameter :: int_garbage$ = -987654
real(rp), parameter :: real_garbage$ = -987654.3

! lf$ (the line feed or LF character) can be used to encode a multiline string.
! EG: string = 'First Line' // lf$ // 'Second Line'

character(1), parameter :: lf$ = achar(10)

! True and false

real(rp), parameter :: true$ = 1, false$ = 0
integer, parameter :: true_int$ = 1, false_int$ = 0

! Color escape sequences

character(*), parameter :: rl_prompt_start_ignore = achar(1)   ! For use with GNU readline routine.
character(*), parameter :: rl_prompt_end_ignore = achar(2)     ! For use with GNU readline routine.

character(*), parameter :: black_color = achar(27) // '[30m' 
character(*), parameter :: red_color = achar(27) // '[31m' 
character(*), parameter :: green_color = achar(27) // '[32m' 
character(*), parameter :: yellow_color = achar(27) // '[33m' 
character(*), parameter :: blue_color = achar(27) // '[34m' 
character(*), parameter :: magenta_color = achar(27) // '[35m' 
character(*), parameter :: cyan_color = achar(27) // '[36m' 
character(*), parameter :: gray_color = achar(27) // '[37m' 


character(*), parameter :: dark_gray_color = achar(27) // '[90m' 
character(*), parameter :: peach_color = achar(27) // '[91m' 
character(*), parameter :: light_green_color = achar(27) // '[92m' 
character(*), parameter :: light_yellow_color = achar(27) // '[93m' 
character(*), parameter :: light_blue_color = achar(27) // '[94m' 
character(*), parameter :: pink_color = achar(27) // '[95m' 
character(*), parameter :: aqua_color = achar(27) // '[96m' 
character(*), parameter :: white_color = achar(27) // '[97m' 

character(*), parameter :: blink_color = achar(27) // '[5m' 
character(*), parameter :: bold_color = achar(27) // '[1m' 

character(*), parameter :: reset_color = achar(27) // '[0m' 

! This is to suppress the ranlib "has no symbols" message

integer, private :: private_dummy

end module
