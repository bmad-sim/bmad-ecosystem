module rf_field_coef_calc_mod

use bmad

type field_cylinder_values
  complex(rp) :: rho = 0, phi = 0, z = 0
end type

!------------------------------------------------------------------
!------------------------------------------------------------------
!------------------------------------------------------------------
contains

function e_field_calc (rho, phi, z, modes, use_mode) result (E)

implicit none

type (rf_field_mode_struct), target :: modes(:)
type (rf_field_mode_struct), pointer :: mode
type (field_cylinder_values) E

complex(rp) ikz, expi
complex(rp) Im_minus, Im0, Im_plus, Im_norm, kappa_n

real(rp) rho, phi, z, k_t, kappa2_n, kap_rho, c, s, k_z, dz

integer i, imode, m, n2_z

logical :: use_mode(:)

!

E%rho = 0; E%phi = 0; E%z = 0

do imode = 1, size(modes)

  if (.not. use_mode(imode)) cycle
  mode => modes(imode)
  m = mode%m
  n2_z = size(mode%term)
  dz = mode%dz

  k_t = twopi * mode%freq / c_light

  do i = 1, n2_z
    k_z = twopi * (i - 1) / (size(mode%term) * mode%dz)
    if (2 * i > n2_z) k_z = k_z - twopi / dz
    kappa2_n = k_z**2 - k_t**2
    kappa_n = sqrt(abs(kappa2_n))
    kap_rho = kappa_n * rho
    if (kappa2_n < 0) then
      kappa_n = -i_imaginary * kappa_n
      kap_rho = -kap_rho
    endif
    expi = cmplx(cos(k_z * z), sin(k_z * z))

    if (m == 0) then
      Im0     = I_bessel(0, kap_rho)
      Im_plus = I_bessel(1, kap_rho) / kappa_n
      E%rho = E%rho + -i_imaginary * k_z * mode%term(i)%e * Im_plus * expi
      E%phi = E%phi + mode%term(i)%b * Im_plus * expi
      E%z   = E%z   + mode%term(i)%e * Im0 * expi

    else
      c = cos(m * phi - mode%phi_0)
      s = sin(m * phi - mode%phi_0)
      Im_plus  = I_bessel(m+1, kap_rho) / kappa_n**(m+1)
      Im_minus = I_bessel(m-1, kap_rho) / kappa_n**(m-1) 

      Im_norm  = (Im_minus - Im_plus * kappa_n**2) / (2 * m) ! I_m(kap_rho) / rho
      Im0      = rho * Im_norm

      E%rho = E%rho - i_imaginary * (k_z * mode%term(i)%e * Im_plus + mode%term(i)%b * Im_norm) * c * expi
      E%phi = E%phi - i_imaginary * (k_z * mode%term(i)%e * Im_plus + &
                                             mode%term(i)%b * (Im_norm - Im_minus / m)) * s * expi
      E%z   = E%z   + mode%term(i)%e * Im0 * c * expi

    endif

  enddo

enddo

end function

end module
