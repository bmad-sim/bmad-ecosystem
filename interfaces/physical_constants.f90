module physical_constants

use precision_def

real(rp), parameter :: pi = 3.14159265358979d0
real(rp), parameter :: twopi = 2 * pi
real(rp), parameter :: fourpi = 4 * pi
real(rp), parameter :: sqrt_2 = 1.41421356237310d0

real(rp), parameter :: e_mass = 0.51099906d-3   ! DO NOT USE! In GeV
real(rp), parameter :: p_mass   = 0.938271998d0   ! DO NOT USE! In GeV

real(rp), parameter :: m_electron = 0.51099906d6  ! Mass in eV
real(rp), parameter :: m_proton   = 0.938271998d9 ! Mass in eV

real(rp), parameter :: c_light = 2.99792458d8   ! speed of light
real(rp), parameter :: r_e = 2.8179380d-15      ! classical electron radius
real(rp), parameter :: r_p = r_e * m_electron / m_proton  ! proton radius
real(rp), parameter :: e_charge = 1.6021892d-19 ! electron charge

real(rp), parameter :: h_planck = 4.13566733d-15      ! eV*sec Planck's constant
real(rp), parameter :: h_bar_planck = 6.58211899d-16  ! eV*sec h_planck/twopi

! Anomalous gyro-magnetic moment

real(rp), parameter :: g_factor_electron = 0.001159652193
real(rp), parameter :: g_factor_proton   = 1.79285

complex(rp), parameter :: i_imaginary = (0.0d0, 1.0d0)
  
! real_garbage$ and int_garbage$ can be used, for example, to identify
! variable that have not been set.

integer, parameter :: int_garbage$ = -987654
real(rp), parameter :: real_garbage$ = -987654.3

! This is to suppress the ranlib "has no symbols" message

integer, private :: private_dummy

end module
