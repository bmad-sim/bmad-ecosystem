module crystal_param_mod

use bmad_struct

type, private :: atom_position_struct
  real(rp) vec(3)
end type

type, private :: crystal_params_struct
  type (atom_position_struct), allocatable :: r_atom(:)
  real(rp), allocatable :: a(:), b(:) 
  real(rp) c, a0
end type

type (crystal_params_struct), private, target, save :: params_si, params_b4c, params_w

private pointer_to_crystal_params

contains

!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
!+
! Subroutine multilayer_type_to_multilayer_params (ele, err_flag)
!
! Routine to set the multilayer parameters based upon the multilayer type.
!
! Multilayer types are of the form:
!   "AAA:BBB"
! Where "AAA" is the top layer crystal and "BBB" is the second layer crystal.
! the diffraction plans.
! Known elements are:
!   SI,  B4C,  W
!
! Modules needed:
!   use crystal_param_mod
!
! Input:
!   ele      -- ele_struct: Multilayer element.
!     %component_name -- Character: Multilayer type name. Assumed upper case.
!                          A blank name is not an error and results in nothing set.
!     %value(e_tot$)  -- Photon energy in eV.
!
! Output:
!   ele      -- ele_struct: Multilayer element.
!     %value(v_unitcell$)
!     %value(f0_re1$)
!     %value(f0_im1$)
!     %value(f0_re2$)
!     %value(f0_im2$)
!   err_flag -- Logical: Set True if multilayer type is unrecognized. False otherwise.
!-

subroutine multilayer_type_to_multilayer_params (ele, err_flag)

implicit none

type (ele_struct) ele

integer ix

logical err_flag

character(40) :: r_name = 'multilayer_type_to_multilayer_params'

! get types

err_flag = .true.

ix = index(ele%component_name, ':')
if (ix == 0) return

call load_layer_params (ele%component_name(1:ix-1), ele%value(f0_re1$), ele%value(f0_im1$), ele%value(v1_unitcell$))
if (err_flag) return

call load_layer_params (ele%component_name(ix+1:),  ele%value(f0_re2$), ele%value(f0_im2$), ele%value(v2_unitcell$))

!-----------------------------------------------------------------------------------------
contains

subroutine load_layer_params (material_name, f0_re, f0_im, v_unitcell)

type (crystal_params_struct), pointer :: cp

real(rp) f0_re, f0_im, v_unitcell
integer n_atom

character(*) material_name

complex f0

! Load constants for the particular crystal...
! Silicon

select case (material_name)
case ('SI') 

  f0 = atomic_f0_calc('SI', ele%value(e_tot$))
  cp => pointer_to_crystal_params('SI')

! Boron Carbide

case ('B4C')

  f0 = 4 * atomic_f0_calc('B', ele%value(e_tot$)) + atomic_f0_calc('C', ele%value(e_tot$)) 
  cp => pointer_to_crystal_params('B4C')

! Tungstan

case ('W')

  f0 = atomic_f0_calc('W', ele%value(e_tot$))
  cp => pointer_to_crystal_params('W')

case default

  call out_io (s_fatal$, r_name, 'BAD CRYSTAL TYPE: ' // ele%component_name)
  err_flag = .true.

end select

!

n_atom = size(cp%r_atom)
f0_re = n_atom * real(f0)
f0_im = n_atom * aimag(f0)

v_unitcell = cp%a0**3

err_flag = .false.

end subroutine load_layer_params 

end subroutine multilayer_type_to_multilayer_params

!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
!+
! Subroutine crystal_type_to_crystal_params (ele, err_flag)
!
! Routine to set the crystal parameters based upon the crystal type.
!
! Crystal types are of the form:
!   "ZZijk"
! Where "ZZ" is the crystal name and "ijk" is the reciprical lattice vetor specifying
! the diffraction plans.
! Known elements are:
!   SI
!
! Modules needed:
!   use crystal_param_mod
!
! Input:
!   ele      -- ele_struct: Crystal element.
!     %component_name -- Character: Crystal type name. Assumed upper case.
!                          A blank name is not an error and results in nothing set.
!     %value(e_tot$)  -- Photon energy in eV.
!
! Output:
!   ele      -- ele_struct: Crystal element.
!     %value(d_spacing$)
!     %value(v_unitcell$)
!     %value(f0_re$)
!     %value(f0_im$)
!     %value(fh_re$)
!     %value(fh_im$)
!     %value(bragg_angle$)
!     %value(alpha_angle$)
!   err_flag -- Logical: Set True if crystal type is unrecognized. False otherwise.
!-

subroutine crystal_type_to_crystal_params (ele, err_flag)

implicit none

type (ele_struct) ele
type (crystal_params_struct), pointer :: cp

real(rp) f_e, f_h, d, hkl(3), bp, r, test(3)
complex(rp) arg, atomsum, f0

integer i, ix, ios, n_coef, n_atom

logical err_flag

character(40) :: r_name = 'crystal_type_to_crystal_params'
character(10) hkl_str

! Check if component_name and e_tot are set.

err_flag = .false.

if (ele%component_name == '') return
if (ele%value(e_tot$) == 0) return

err_flag = .true.

! Load constants for the particular crystal...
! Silicon

if (ele%component_name(1:2) == 'SI') then
  hkl_str = ele%component_name(3:)
  f0 = atomic_f0_calc('SI', ele%value(e_tot$))
  cp => pointer_to_crystal_params('SI')

else

  call out_io (s_fatal$, r_name, 'BAD CRYSTAL TYPE: ' // ele%component_name)

endif

!------------------------------------------------------
! Multiply by number of atoms to get CRYSTAL structure factor, ie f0_re and f0_im

n_atom = size(cp%r_atom)
ele%value(f0_re$) = n_atom * real(f0)
ele%value(f0_im$) = n_atom * aimag(f0)

ele%value(v_unitcell$) = cp%a0**3

! Calculate HKL

read (hkl_str(1:1), *, iostat = ios) hkl(1)
if (hkl_str(1:1) == '' .or. ios /= 0) return

read (hkl_str(2:2), *, iostat = ios) hkl(2)
if (hkl_str(2:2) == '' .or. ios /= 0) return

read (hkl_str(3:3), *, iostat = ios) hkl(3)
if (hkl_str(3:3) == '' .or. ios /= 0) return

! Calculate values for the given energy

d = cp%a0 / sqrt(sum(hkl**2))
ele%value(d_spacing$) = d

! Calculate the hkl dependent part of f and
! subtract from f(hkl=0) to obtain the energy dependent part

n_coef = size(cp%a)
f_e = real(f0) - (sum(cp%a(1:n_coef)) + cp%c)

! Calculate hkl dependent part of f
! factor of 1e20 due to units of d

f_h = sum( cp%a(1:n_coef) * exp( -cp%b(1:n_coef)/(4*d*d*1e20) ) ) + cp%c

atomsum = 0
do i = 1, n_atom
  test = hkl * cp%r_atom(i)%vec
  test(1) = sum(hkl * cp%r_atom(i)%vec)
  arg = cmplx(0.0, twopi * sum(hkl * cp%r_atom(i)%vec))
  atomsum = atomsum + exp(arg)
enddo

ele%value(fh_re$) = abs(atomsum) * (f_h + f_e)
ele%value(fh_im$) = abs(atomsum) * aimag(f0)

! Set bragg and alpha angles to zero if bragg condition is not satisfied.

r = ele%value(ref_wavelength$) / (2 * d)

if (r < 1) then
  ele%value(bragg_angle$) = asin(r)
  bp = ele%value(b_param$)
  ele%value(alpha_angle$) = atan(tan(ele%value(bragg_angle$)) * (bp + 1) / (bp - 1))
else
  ele%value(bragg_angle$) = 0
  ele%value(alpha_angle$) = 0
endif


err_flag = .false.

end subroutine crystal_type_to_crystal_params

!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
!+
! Function atomic_f0_calc (atom, e_tot) result (f0)
!
! Function to calculate the F0 structure factor for a particular element at a particular energy.
! Known elements at this point are:
!     B,  C,  SI,  W
! Input:
!   atom    -- Character(*): Type of element.
!   e_tot   -- Real(rp): total energy in eV.
!
! Output:
!   f0      -- Complex(rp): Complex structure factor
!-

function atomic_f0_calc (atom, e_tot) result (f0)

use bmad_struct

implicit none

type atomic_f_struct
  real(rp) energy
  real(rp) f0_re, f0_im
end type

type (atomic_f_struct), target, save :: f0_b(94), f0_c(84), f0_si(57), f0_w(91)
type (atomic_f_struct), pointer :: f0_ptr(:)

real(rp) e_tot, r
complex f0

integer n, ix

logical :: init_needed = .true.

character(*) atom

! Init

if (init_needed) then
  init_needed = .false.
  ! si
  f0_si( 1) = atomic_f_struct(10.000000,  -0.999900,   4.0068800)
  f0_si( 2) = atomic_f_struct(12.123200,  -0.999900,   4.0001300)
  f0_si( 3) = atomic_f_struct(17.817800,  -0.999900,   2.1237800)
  f0_si( 4) = atomic_f_struct(21.600900,  -0.999900,   0.6025340)
  f0_si( 5) = atomic_f_struct(28.833700,  -0.999900,   0.3706440)
  f0_si( 6) = atomic_f_struct(29.300000,   3.838100,   0.3717420)
  f0_si( 7) = atomic_f_struct(76.728900,   2.097630,   0.4773930)
  f0_si( 8) = atomic_f_struct(85.849100,   1.136100,   0.4159400)
  f0_si( 9) = atomic_f_struct(91.539500,   0.232650,   0.4449190)
  f0_si(10) = atomic_f_struct(96.053500,  -1.001990,   0.4806420)
  f0_si(11) = atomic_f_struct(97.607100,  -1.779850,   0.4931760)
  f0_si(12) = atomic_f_struct(99.185800,  -3.569080,   0.5060360)
  f0_si(13) = atomic_f_struct(99.700000,  -6.180010,   0.5102520)
  f0_si(14) = atomic_f_struct(99.900000,  -6.058740,   4.6074200)
  f0_si(15) = atomic_f_struct(100.79000,  -2.494140,   4.3679000)
  f0_si(16) = atomic_f_struct(102.42000,  -1.410120,   3.9658400)
  f0_si(17) = atomic_f_struct(104.07700,  -1.590760,   3.6734400)
  f0_si(18) = atomic_f_struct(105.76000,  -1.863930,   4.0673800)
  f0_si(19) = atomic_f_struct(110.97500,  -1.907960,   5.4489100)
  f0_si(20) = atomic_f_struct(116.44800,  -1.786350,   7.1534100)
  f0_si(21) = atomic_f_struct(126.17500,   0.525660,   8.9108500)
  f0_si(22) = atomic_f_struct(134.53800,   2.683200,   9.5984200)
  f0_si(23) = atomic_f_struct(143.45600,   3.126990,   9.2460900)
  f0_si(24) = atomic_f_struct(152.96400,   4.181810,  10.0588000)
  f0_si(25) = atomic_f_struct(157.95300,   4.417420,   9.8635000)
  f0_si(26) = atomic_f_struct(168.42200,   6.250620,  10.9311000)
  f0_si(27) = atomic_f_struct(179.58600,   7.683370,   9.4475800)
  f0_si(28) = atomic_f_struct(232.14700,  10.217500,   8.2066400)
  f0_si(29) = atomic_f_struct(300.09200,  12.002900,   6.4372400)
  f0_si(30) = atomic_f_struct(387.92200,  13.019300,   4.7842000)
  f0_si(31) = atomic_f_struct(501.45900,  13.376000,   3.3132500)
  f0_si(32) = atomic_f_struct(648.22600,  13.372100,   2.2138500)
  f0_si(33) = atomic_f_struct(952.71500,  13.059900,   1.1738300)
  f0_si(34) = atomic_f_struct(1231.5500,  12.642000,   0.7399720)
  f0_si(35) = atomic_f_struct(1592.0100,  11.814500,   0.4776830)
  f0_si(36) = atomic_f_struct(1697.5300,  11.273600,   0.4258760)
  f0_si(37) = atomic_f_struct(1752.8900,  10.764400,   0.4016720)
  f0_si(38) = atomic_f_struct(1781.2400,  10.337500,   0.3897180)
  f0_si(39) = atomic_f_struct(1810.0500,   9.569810,   0.3780020)
  f0_si(40) = atomic_f_struct(1838.8000,   2.834190,   0.3668380)
  f0_si(41) = atomic_f_struct(1839.0000,   2.835870,   4.1573400)
  f0_si(42) = atomic_f_struct(1839.3200,   4.580550,   4.1563200)
  f0_si(43) = atomic_f_struct(1869.0700,   9.825010,   4.0645200)
  f0_si(44) = atomic_f_struct(1899.3000,  10.724100,   3.9747500)
  f0_si(45) = atomic_f_struct(1930.0200,  11.267200,   3.8868800)
  f0_si(46) = atomic_f_struct(2265.9200,  13.299100,   3.0795600)
  f0_si(47) = atomic_f_struct(2576.2600,  13.885700,   2.5186900)
  f0_si(48) = atomic_f_struct(2929.1100,  14.190600,   2.0408500)
  f0_si(49) = atomic_f_struct(4894.6000,  14.423300,   0.8266710)
  f0_si(50) = atomic_f_struct(6327.1500,  14.360900,   0.5131970)
  f0_si(51) = atomic_f_struct(10572.800,  14.209300,   0.1915190)
  f0_si(52) = atomic_f_struct(17667.400,  14.062300,   0.0687904)
  f0_si(53) = atomic_f_struct(29522.500,  14.020800,   0.0236346)
  f0_si(54) = atomic_f_struct(43000.000,  14.007800,   0.0106303)
  f0_si(55) = atomic_f_struct(64000.000,  13.996200,   0.0044384)
  f0_si(56) = atomic_f_struct(86500.000,  13.991300,   0.0022794)
  f0_si(57) = atomic_f_struct(100000.00,  13.989700,   0.00165265)

  ! Boron
  f0_b( 1) = atomic_f_struct(10.0000, -9999., 1.48933)
  f0_b( 2) = atomic_f_struct(11.0106, -9999., 1.43913)
  f0_b( 3) = atomic_f_struct(12.1232, -9999., 1.39063)
  f0_b( 4) = atomic_f_struct(13.3483, -9999., 1.34376)
  f0_b( 5) = atomic_f_struct(14.6973, -9999., 1.29847)
  f0_b( 6) = atomic_f_struct(16.1825, -9999., 1.23647)
  f0_b( 7) = atomic_f_struct(17.8178, -9999., 1.16976)
  f0_b( 8) = atomic_f_struct(19.6184, -9999., 1.10665)
  f0_b( 9) = atomic_f_struct(21.6009, -9999., 1.04694)
  f0_b(10) = atomic_f_struct(23.7838, -9999., 0.990606)
  f0_b(11) = atomic_f_struct(26.1873, -9999., 0.938720)
  f0_b(12) = atomic_f_struct(28.8337, -9999., 0.889552)
  f0_b(13) = atomic_f_struct(31.7475, 2.74123, 0.847767)
  f0_b(14) = atomic_f_struct(34.9557, 2.79220, 0.811138)
  f0_b(15) = atomic_f_struct(38.4882, 2.84076, 0.770390)
  f0_b(16) = atomic_f_struct(42.3776, 2.88508, 0.728996)
  f0_b(17) = atomic_f_struct(46.6600, 2.92380, 0.686952)
  f0_b(18) = atomic_f_struct(51.3753, 2.95768, 0.643703)
  f0_b(19) = atomic_f_struct(56.5670, 2.98242, 0.599067)
  f0_b(20) = atomic_f_struct(62.2834, 3.00149, 0.555186)
  f0_b(21) = atomic_f_struct(68.5775, 3.01030, 0.511109)
  f0_b(22) = atomic_f_struct(75.5076, 3.00995, 0.466822)
  f0_b(23) = atomic_f_struct(83.1380, 2.99744, 0.424073)
  f0_b(24) = atomic_f_struct(91.5395, 2.97030, 0.382981)
  f0_b(25) = atomic_f_struct(100.790, 2.92657, 0.343808)
  f0_b(26) = atomic_f_struct(110.975, 2.85569, 0.306228)
  f0_b(27) = atomic_f_struct(122.190, 2.75708, 0.273649)
  f0_b(28) = atomic_f_struct(134.538, 2.60723, 0.242496)
  f0_b(29) = atomic_f_struct(148.134, 2.36863, 0.214205)
  f0_b(30) = atomic_f_struct(163.103, 1.93966, 0.185553)
  f0_b(31) = atomic_f_struct(179.586, 0.782171, 0.161758)
  f0_b(32) = atomic_f_struct(191.489, -0.321484E-01, 4.16670)
  f0_b(33) = atomic_f_struct(210.840, 2.68979, 3.75945)
  f0_b(34) = atomic_f_struct(232.147, 3.67588, 3.31961)
  f0_b(35) = atomic_f_struct(255.607, 4.24126, 2.92641)
  f0_b(36) = atomic_f_struct(281.437, 4.62551, 2.57409)
  f0_b(37) = atomic_f_struct(309.878, 4.88183, 2.24035)
  f0_b(38) = atomic_f_struct(341.192, 5.05277, 1.94978)
  f0_b(39) = atomic_f_struct(375.672, 5.17136, 1.69691)
  f0_b(40) = atomic_f_struct(413.635, 5.26049, 1.46951)
  f0_b(41) = atomic_f_struct(455.435, 5.31408, 1.26256)
  f0_b(42) = atomic_f_struct(501.459, 5.34167, 1.07788)
  f0_b(43) = atomic_f_struct(552.134, 5.35135, 0.916550)
  f0_b(44) = atomic_f_struct(607.930, 5.34917, 0.779833)
  f0_b(45) = atomic_f_struct(669.365, 5.33936, 0.653250)
  f0_b(46) = atomic_f_struct(737.008, 5.31890, 0.553984)
  f0_b(47) = atomic_f_struct(811.486, 5.29797, 0.465121)
  f0_b(48) = atomic_f_struct(893.491, 5.27471, 0.388746)
  f0_b(49) = atomic_f_struct(983.783, 5.24931, 0.324999)
  f0_b(50) = atomic_f_struct(1083.20, 5.22591, 0.270598)
  f0_b(51) = atomic_f_struct(1192.66, 5.20216, 0.224555)
  f0_b(52) = atomic_f_struct(1313.19, 5.17972, 0.186228)
  f0_b(53) = atomic_f_struct(1445.89, 5.15890, 0.154338)
  f0_b(54) = atomic_f_struct(1592.01, 5.13999, 0.127556)
  f0_b(55) = atomic_f_struct(1752.89, 5.12270, 0.105266)
  f0_b(56) = atomic_f_struct(1930.02, 5.10705, 0.866612E-01)
  f0_b(57) = atomic_f_struct(2125.06, 5.09311, 0.713481E-01)
  f0_b(58) = atomic_f_struct(2339.81, 5.08070, 0.585635E-01)
  f0_b(59) = atomic_f_struct(2576.26, 5.06971, 0.479897E-01)
  f0_b(60) = atomic_f_struct(2836.61, 5.06001, 0.392606E-01)
  f0_b(61) = atomic_f_struct(3123.26, 5.05150, 0.320658E-01)
  f0_b(62) = atomic_f_struct(3438.88, 5.04406, 0.261471E-01)
  f0_b(63) = atomic_f_struct(3786.40, 5.03758, 0.212856E-01)
  f0_b(64) = atomic_f_struct(4169.03, 5.03195, 0.173072E-01)
  f0_b(65) = atomic_f_struct(4590.33, 5.02708, 0.140547E-01)
  f0_b(66) = atomic_f_struct(5054.21, 5.02287, 0.114019E-01)
  f0_b(67) = atomic_f_struct(5564.97, 5.01926, 0.923971E-02)
  f0_b(68) = atomic_f_struct(6127.33, 5.01616, 0.747830E-02)
  f0_b(69) = atomic_f_struct(6746.54, 5.01351, 0.604440E-02)
  f0_b(70) = atomic_f_struct(6855.65, 5.01311, 0.583287E-02)
  f0_b(71) = atomic_f_struct(6966.54, 5.01271, 0.562839E-02)
  f0_b(72) = atomic_f_struct(7079.22, 5.01233, 0.543100E-02)
  f0_b(73) = atomic_f_struct(7193.72, 5.01196, 0.524026E-02)
  f0_b(74) = atomic_f_struct(7310.07, 5.01160, 0.505601E-02)
  f0_b(75) = atomic_f_struct(7428.31, 5.01124, 0.487800E-02)
  f0_b(76) = atomic_f_struct(7548.45, 5.01090, 0.470595E-02)
  f0_b(77) = atomic_f_struct(7670.54, 5.01056, 0.453990E-02)
  f0_b(78) = atomic_f_struct(7794.61, 5.01024, 0.437949E-02)
  f0_b(79) = atomic_f_struct(7920.68, 5.00992, 0.422456E-02)
  f0_b(80) = atomic_f_struct(8048.79, 5.00961, 0.407489E-02)
  f0_b(81) = atomic_f_struct(8178.98, 5.00931, 0.393020E-02)
  f0_b(82) = atomic_f_struct(9005.50, 5.00766, 0.316109E-02)
  f0_b(83) = atomic_f_struct(9915.55, 5.00626, 0.253775E-02)
  f0_b(84) = atomic_f_struct(10917.6, 5.00507, 0.203335E-02)
  f0_b(85) = atomic_f_struct(12020.8, 5.00406, 0.162601E-02)
  f0_b(86) = atomic_f_struct(13235.6, 5.00321, 0.129780E-02)
  f0_b(87) = atomic_f_struct(14573.1, 5.00248, 0.103404E-02)
  f0_b(88) = atomic_f_struct(16045.8, 5.00187, 0.822718E-03)
  f0_b(89) = atomic_f_struct(17667.4, 5.00135, 0.654024E-03)
  f0_b(90) = atomic_f_struct(19452.7, 5.00092, 0.519930E-03)
  f0_b(91) = atomic_f_struct(21418.5, 5.00055, 0.413645E-03)
  f0_b(92) = atomic_f_struct(23583.0, 5.00025, 0.329702E-03)
  f0_b(93) = atomic_f_struct(25966.2, 4.99999, 0.263151E-03)
  f0_b(94) = atomic_f_struct(28590.2, 4.99977, 0.210290E-03)

  ! Carbon
  f0_c(01) = atomic_f_struct(10.0000, -9999., 0.806885)
  f0_c(02) = atomic_f_struct(11.0106, -9999., 1.12167)
  f0_c(03) = atomic_f_struct(12.1232, -9999., 1.58971)
  f0_c(04) = atomic_f_struct(13.3483, -9999., 2.20025)
  f0_c(05) = atomic_f_struct(14.6973, -9999., 2.77114)
  f0_c(06) = atomic_f_struct(16.1825, -9999., 3.24213)
  f0_c(07) = atomic_f_struct(17.8178, -9999., 3.52041)
  f0_c(08) = atomic_f_struct(19.6184, -9999., 3.65674)
  f0_c(09) = atomic_f_struct(21.6009, -9999., 3.62677)
  f0_c(10) = atomic_f_struct(23.7838, -9999., 3.42547)
  f0_c(11) = atomic_f_struct(26.1873, -9999., 3.12096)
  f0_c(12) = atomic_f_struct(28.8337, -9999., 2.79164)
  f0_c(13) = atomic_f_struct(31.7475, 3.80969, 2.48155)
  f0_c(14) = atomic_f_struct(34.9557, 3.95026, 2.19898)
  f0_c(15) = atomic_f_struct(38.4882, 4.05432, 1.96657)
  f0_c(16) = atomic_f_struct(42.3776, 4.13985, 1.72163)
  f0_c(17) = atomic_f_struct(46.6600, 4.14174, 1.51876)
  f0_c(18) = atomic_f_struct(51.3753, 4.15277, 1.37705)
  f0_c(19) = atomic_f_struct(56.5670, 4.15681, 1.27454)
  f0_c(20) = atomic_f_struct(62.2834, 4.21096, 1.18498)
  f0_c(21) = atomic_f_struct(68.5775, 4.24504, 1.06347)
  f0_c(22) = atomic_f_struct(75.5076, 4.26047, 0.954423)
  f0_c(23) = atomic_f_struct(83.1380, 4.26546, 0.856556)
  f0_c(24) = atomic_f_struct(91.5395, 4.26226, 0.768724)
  f0_c(25) = atomic_f_struct(100.790, 4.25166, 0.689898)
  f0_c(26) = atomic_f_struct(110.975, 4.23373, 0.619156)
  f0_c(27) = atomic_f_struct(122.190, 4.20806, 0.555666)
  f0_c(28) = atomic_f_struct(134.538, 4.17533, 0.497457)
  f0_c(29) = atomic_f_struct(148.134, 4.12638, 0.439868)
  f0_c(30) = atomic_f_struct(163.103, 4.05754, 0.388946)
  f0_c(31) = atomic_f_struct(179.586, 3.96398, 0.343919)
  f0_c(32) = atomic_f_struct(197.734, 3.83308, 0.295031)
  f0_c(33) = atomic_f_struct(217.716, 3.62164, 0.249231)
  f0_c(34) = atomic_f_struct(239.717, 3.26522, 0.213355)
  f0_c(35) = atomic_f_struct(263.942, 2.48138, 0.173821)
  f0_c(36) = atomic_f_struct(284.300, -4.04031, 4.18872)
  f0_c(37) = atomic_f_struct(309.878, 3.40186, 3.77095)
  f0_c(38) = atomic_f_struct(341.192, 4.53993, 3.35318)
  f0_c(39) = atomic_f_struct(375.672, 5.18910, 2.98170)
  f0_c(40) = atomic_f_struct(413.635, 5.61669, 2.59050)
  f0_c(41) = atomic_f_struct(455.435, 5.88264, 2.24703)
  f0_c(42) = atomic_f_struct(501.459, 6.06268, 1.94589)
  f0_c(43) = atomic_f_struct(552.134, 6.17864, 1.67971)
  f0_c(44) = atomic_f_struct(607.930, 6.25287, 1.44994)
  f0_c(45) = atomic_f_struct(669.365, 6.29901, 1.25160)
  f0_c(46) = atomic_f_struct(737.008, 6.33933, 1.07816)
  f0_c(47) = atomic_f_struct(811.486, 6.35140, 0.912939)
  f0_c(48) = atomic_f_struct(893.491, 6.34793, 0.772006)
  f0_c(49) = atomic_f_struct(983.783, 6.33503, 0.651378)
  f0_c(50) = atomic_f_struct(1083.20, 6.31710, 0.548198)
  f0_c(51) = atomic_f_struct(1192.66, 6.29469, 0.459933)
  f0_c(52) = atomic_f_struct(1313.19, 6.27127, 0.384872)
  f0_c(53) = atomic_f_struct(1445.89, 6.24640, 0.321159)
  f0_c(54) = atomic_f_struct(1592.01, 6.22200, 0.267545)
  f0_c(55) = atomic_f_struct(1752.89, 6.19857, 0.222665)
  f0_c(56) = atomic_f_struct(1930.02, 6.17643, 0.184850)
  f0_c(57) = atomic_f_struct(2125.06, 6.15600, 0.153384)
  f0_c(58) = atomic_f_struct(2339.81, 6.13729, 0.126998)
  f0_c(59) = atomic_f_struct(2576.26, 6.12028, 0.104990)
  f0_c(60) = atomic_f_struct(2836.61, 6.10498, 0.866617E-01)
  f0_c(61) = atomic_f_struct(3123.26, 6.09130, 0.714080E-01)
  f0_c(62) = atomic_f_struct(3438.88, 6.07912, 0.587327E-01)
  f0_c(63) = atomic_f_struct(3786.40, 6.06834, 0.482071E-01)
  f0_c(64) = atomic_f_struct(4169.03, 6.05883, 0.394919E-01)
  f0_c(65) = atomic_f_struct(4590.33, 6.05046, 0.322825E-01)
  f0_c(66) = atomic_f_struct(5054.21, 6.04313, 0.263428E-01)
  f0_c(67) = atomic_f_struct(5564.97, 6.03674, 0.214569E-01)
  f0_c(68) = atomic_f_struct(6127.33, 6.03118, 0.174453E-01)
  f0_c(69) = atomic_f_struct(6746.54, 6.02637, 0.141549E-01)
  f0_c(70) = atomic_f_struct(7428.31, 6.02220, 0.114619E-01)
  f0_c(71) = atomic_f_struct(8178.98, 6.01861, 0.926134E-02)
  f0_c(72) = atomic_f_struct(9005.50, 6.01552, 0.746686E-02)
  f0_c(73) = atomic_f_struct(9915.55, 6.01286, 0.600659E-02)
  f0_c(74) = atomic_f_struct(10917.6, 6.01059, 0.482101E-02)
  f0_c(75) = atomic_f_struct(12020.8, 6.00864, 0.386096E-02)
  f0_c(76) = atomic_f_struct(13235.6, 6.00698, 0.308580E-02)
  f0_c(77) = atomic_f_struct(14573.1, 6.00557, 0.246198E-02)
  f0_c(78) = atomic_f_struct(16045.8, 6.00436, 0.196186E-02)
  f0_c(79) = atomic_f_struct(17667.4, 6.00334, 0.156265E-02)
  f0_c(80) = atomic_f_struct(19452.7, 6.00248, 0.124561E-02)
  f0_c(81) = atomic_f_struct(21418.5, 6.00176, 0.994268E-03)
  f0_c(82) = atomic_f_struct(23583.0, 6.00115, 0.787230E-03)
  f0_c(83) = atomic_f_struct(25966.2, 6.00063, 0.624439E-03)
  f0_c(84) = atomic_f_struct(28590.2, 6.00020, 0.496161E-03)
! Tungsta
  f0_w(01) = atomic_f_struct(10.0000, -9999., 1.92551)
  f0_w(02) = atomic_f_struct(11.0105, -9999., 2.44381)
  f0_w(03) = atomic_f_struct(12.1232, -9999., 2.99166)
  f0_w(04) = atomic_f_struct(13.3483, -9999., 3.53050)
  f0_w(05) = atomic_f_struct(14.6972, -9999., 4.16638)
  f0_w(06) = atomic_f_struct(16.1825, -9999., 4.82729)
  f0_w(07) = atomic_f_struct(17.8178, -9999., 5.31009)
  f0_w(08) = atomic_f_struct(19.6184, 4.12288, 5.61231)
  f0_w(09) = atomic_f_struct(21.6009, 4.93090, 5.61526)
  f0_w(10) = atomic_f_struct(23.7838, 5.48845, 5.48096)
  f0_w(11) = atomic_f_struct(26.1872, 5.77356, 4.84732)
  f0_w(12) = atomic_f_struct(28.8336, 5.57575, 4.44943)
  f0_w(13) = atomic_f_struct(31.7474, 4.92262, 3.92770)
  f0_w(14) = atomic_f_struct(34.9556, 3.16812, 3.81788)
  f0_w(15) = atomic_f_struct(38.4880, -0.147289, 6.34746)
  f0_w(16) = atomic_f_struct(42.3774, 3.20822, 12.2742)
  f0_w(17) = atomic_f_struct(46.6599, 3.28163, 10.7655)
  f0_w(18) = atomic_f_struct(51.3751, 5.56617, 16.8755)
  f0_w(19) = atomic_f_struct(56.5668, 11.8679, 15.9366)
  f0_w(20) = atomic_f_struct(62.2832, 14.4082, 13.3481)
  f0_w(21) = atomic_f_struct(68.5772, 15.6849, 10.9127)
  f0_w(22) = atomic_f_struct(75.5072, 15.8508, 8.58572)
  f0_w(23) = atomic_f_struct(83.1376, 14.7189, 6.83406)
  f0_w(24) = atomic_f_struct(91.5391, 13.0514, 6.35435)
  f0_w(25) = atomic_f_struct(100.790, 11.8412, 7.10582)
  f0_w(26) = atomic_f_struct(110.975, 11.1605, 8.08209)
  f0_w(27) = atomic_f_struct(122.189, 10.7924, 9.33195)
  f0_w(28) = atomic_f_struct(134.537, 10.7655, 10.7258)
  f0_w(29) = atomic_f_struct(148.133, 10.9839, 12.2461)
  f0_w(30) = atomic_f_struct(163.102, 11.5855, 13.8412)
  f0_w(31) = atomic_f_struct(179.585, 12.4794, 15.2684)
  f0_w(32) = atomic_f_struct(197.733, 13.3332, 16.7162)
  f0_w(33) = atomic_f_struct(217.714, 14.4653, 18.2840)
  f0_w(34) = atomic_f_struct(239.716, 15.8455, 19.9500)
  f0_w(35) = atomic_f_struct(263.940, 17.7909, 21.7678)
  f0_w(36) = atomic_f_struct(290.612, 20.3687, 22.6118)
  f0_w(37) = atomic_f_struct(319.980, 22.5144, 23.1456)
  f0_w(38) = atomic_f_struct(352.316, 24.4799, 23.5795)
  f0_w(39) = atomic_f_struct(387.919, 26.2331, 24.1987)
  f0_w(40) = atomic_f_struct(427.120, 28.3782, 25.0905)
  f0_w(41) = atomic_f_struct(470.283, 31.4258, 25.5791)
  f0_w(42) = atomic_f_struct(517.807, 34.0772, 25.2376)
  f0_w(43) = atomic_f_struct(570.134, 36.3407, 24.4232)
  f0_w(44) = atomic_f_struct(627.749, 38.5373, 24.4150)
  f0_w(45) = atomic_f_struct(691.186, 41.1576, 23.1836)
  f0_w(46) = atomic_f_struct(761.034, 43.1546, 21.7535)
  f0_w(47) = atomic_f_struct(837.940, 44.6721, 20.1532)
  f0_w(48) = atomic_f_struct(922.618, 45.7471, 18.5711)
  f0_w(49) = atomic_f_struct(1015.85, 46.4694, 17.0243)
  f0_w(50) = atomic_f_struct(1118.51, 46.7756, 15.4670)
  f0_w(51) = atomic_f_struct(1231.54, 46.5978, 14.0353)
  f0_w(52) = atomic_f_struct(1355.99, 45.8309, 12.6248)
  f0_w(53) = atomic_f_struct(1493.02, 44.0159, 11.3266)
  f0_w(54) = atomic_f_struct(1643.90, 39.9275, 10.1275)
  f0_w(55) = atomic_f_struct(1809.10, -5.97426, 9.04256)
  f0_w(56) = atomic_f_struct(1871.70, 7.08657, 39.3036)
  f0_w(57) = atomic_f_struct(2026., 45.7756, 35.1030)
  f0_w(58) = atomic_f_struct(2190., 51.2909, 31.4014)
  f0_w(59) = atomic_f_struct(2306., 52.5275, 33.9369)
  f0_w(60) = atomic_f_struct(2493., 58.5109, 30.0011)
  f0_w(61) = atomic_f_struct(2626., 60.3896, 29.5680)
  f0_w(62) = atomic_f_struct(2819.50, 60.7116, 26.5733)
  f0_w(63) = atomic_f_struct(2990., 65.1894, 25.6028)
  f0_w(64) = atomic_f_struct(3232., 67.1358, 22.9455)
  f0_w(65) = atomic_f_struct(3494., 68.4415, 20.5202)
  f0_w(66) = atomic_f_struct(3777., 69.3333, 18.3220)
  f0_w(67) = atomic_f_struct(4083., 69.9285, 16.3379)
  f0_w(68) = atomic_f_struct(4414., 70.3006, 14.5502)
  f0_w(69) = atomic_f_struct(4771., 70.5000, 12.9478)
  f0_w(70) = atomic_f_struct(5158., 70.5631, 11.5066)
  f0_w(71) = atomic_f_struct(5576., 70.5143, 10.2171)
  f0_w(72) = atomic_f_struct(6027., 70.3696, 9.06587)
  f0_w(73) = atomic_f_struct(6516., 70.1358, 8.03454)
  f0_w(74) = atomic_f_struct(7043., 69.8112, 7.11719)
  f0_w(75) = atomic_f_struct(7614., 69.3765, 6.29741)
  f0_w(76) = atomic_f_struct(8230., 68.7860, 5.56793)
  f0_w(77) = atomic_f_struct(8897., 67.9109, 4.91726)
  f0_w(78) = atomic_f_struct(9618., 66.2737, 4.33729)
  f0_w(79) = atomic_f_struct(10223.4, 59.3366, 10.1559)
  f0_w(80) = atomic_f_struct(11256.6, 66.3010, 8.52616)
  f0_w(81) = atomic_f_struct(12002.7, 67.0878, 10.6804)
  f0_w(82) = atomic_f_struct(12798.3, 69.8814, 11.0950)
  f0_w(83) = atomic_f_struct(14091.6, 71.7428, 9.54983)
  f0_w(84) = atomic_f_struct(15515.6, 72.7063, 8.19207)
  f0_w(85) = atomic_f_struct(17083.6, 73.2506, 7.00785)
  f0_w(86) = atomic_f_struct(18809.9, 73.5487, 5.98119)
  f0_w(87) = atomic_f_struct(20710.8, 73.6901, 5.09543)
  f0_w(88) = atomic_f_struct(22803.7, 73.7287, 4.33424)
  f0_w(89) = atomic_f_struct(25108.1, 73.6993, 3.68219)
  f0_w(90) = atomic_f_struct(27645.4, 73.6248, 3.12506)
  f0_w(91) = atomic_f_struct(30000.0, 73.5377, 2.71692)

endif

! Interpolate atomic structure factors from factor struct, which are evaluated at hkl=0

select case (atom)
case ('B')
  f0_ptr => f0_w
case ('C')
  f0_ptr => f0_w
case ('SI')
  f0_ptr => f0_si
case ('W')
  f0_ptr => f0_w
end select

n = size(f0_ptr)

call bracket_index (f0_ptr%energy, 1, n, e_tot, ix)
if (ix == 0) then
  f0 = cmplx(f0_ptr(1)%f0_re, f0_ptr(1)%f0_im)
elseif (ix == n) then
  f0 = cmplx(f0_ptr(n)%f0_re, f0_ptr(n)%f0_im)
else
  r = (f0_ptr(ix+1)%energy - e_tot) / (f0_ptr(ix+1)%energy - f0_ptr(ix)%energy)
  f0 = cmplx( f0_ptr(ix)%f0_re * r + f0_ptr(ix+1)%f0_re * (1 - r), &
                f0_ptr(ix)%f0_im * r + f0_ptr(ix+1)%f0_im * (1 - r))
endif

end function atomic_f0_calc 

!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
!+
! 
! Function pointer_to_crystal_params (crystal_type) result (param_ptr)
! 
! Routine to return a pointer to the crystal structure parameters.
!
! Input:
!   crystal_type -- Character(*): Crystal type.
!
! Output:
!   param_ptr -- crystal_params_struct, pointer: Pointer to the parameters
!-

function pointer_to_crystal_params (crystal_type) result (param_ptr)

implicit none

type (crystal_params_struct), pointer :: param_ptr, cp

logical :: init_needed = .true.

character(*) crystal_type

!

if (init_needed) then
  init_needed = .false.
  ! Silicon
  cp => params_si
  allocate (cp%r_atom(8), cp%a(5), cp%b(5))
  cp%a = [4.98816795, 3.35710271, 1.50292204, 1.22172882, 2.76143663]
  cp%b = [2.53600438, 29.97580504, 0.08254945, 88.73513838, 1.16712390]
  cp%c = 0.15142442
  cp%a0 = 5.4309e-10

  cp%r_atom(1)%vec = [0.00, 0.00, 0.00]
  cp%r_atom(2)%vec = [0.00, 0.50, 0.50]
  cp%r_atom(3)%vec = [0.50, 0.00, 0.50]
  cp%r_atom(4)%vec = [0.50, 0.50, 0.00]
  cp%r_atom(5)%vec = [0.25, 0.25, 0.25]
  cp%r_atom(6)%vec = [0.25, 0.75, 0.75]
  cp%r_atom(7)%vec = [0.75, 0.25, 0.75]
  cp%r_atom(8)%vec = [0.75, 0.75, 0.25]
end if

!

select case (crystal_type)
case ('SI');   param_ptr => params_si
case ('B4C');  param_ptr => params_b4c
case ('W');    param_ptr => params_w
case default; call err_exit
end select 

end function pointer_to_crystal_params 

end module

