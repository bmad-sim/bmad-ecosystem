!+
! Subroutine crystal_type_to_crystal_params (ele, err_flag)
!
! Routine to set the crystal parameters in an element based upon the crystal type.
! Crystal types are of the form:
!   "ZZijk"
! Where "ZZ" is the element name and "ijk" is the reciprical lattice vetor specifying
! the diffraction plans.
! Known elements are:
!   SI
!
! Modules needed:
!   use bmad
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
!     %value(gamma_factor$)
!   err_flag -- Logical: Set True if crystal type is unrecognized. False otherwise.
!-

subroutine crystal_type_to_crystal_params (ele, err_flag)

use bmad_struct

implicit none

type f_factor_struct
  real(rp) energy
  real(rp) f0_re, f0_im
end type

type atom_position_struct
  real(rp) vec(3)
end type

type (ele_struct) ele
type (atom_position_struct) r_atom(8)
type (f_factor_struct) factor(100)

real(rp) f_re, f_im, f_e, f_h, fh_re, fh_im
real(rp) a(10), b(10), c, d, a0, hkl(3), bp, r, test(3)
complex(rp) arg, atomsum

integer i, ix, n, ios, n_atom, n_factor

logical err_flag

character(40) :: r_name = 'crystal_type_to_crystal_params'
character(10) hkl_str

! Check if component_name and e_tot are set.

err_flag = .false.

if (ele%component_name == '') return
if (ele%value(e_tot$) == 0) return

! Load constants for the particular crystal compound

err_flag = .true.

if (ele%component_name(1:2) == 'SI') then

  a(1:5) = [4.98816795, 3.35710271, 1.50292204, 1.22172882, 2.76143663]
  b(1:5) = [2.53600438, 29.97580504, 0.08254945, 88.73513838, 1.16712390]
  c = 0.15142442
  n = 5
  a0 = 5.4309e-10
  hkl_str = ele%component_name(3:)

  r_atom(1)%vec = [0.00, 0.00, 0.00]
  r_atom(2)%vec = [0.00, 0.50, 0.50]
  r_atom(3)%vec = [0.50, 0.00, 0.50]
  r_atom(4)%vec = [0.50, 0.50, 0.00]
  r_atom(5)%vec = [0.25, 0.25, 0.25]
  r_atom(6)%vec = [0.25, 0.75, 0.75]
  r_atom(7)%vec = [0.75, 0.25, 0.75]
  r_atom(8)%vec = [0.75, 0.75, 0.25]
  n_atom = 8

  factor( 1) = f_factor_struct(10.000000,  -0.999900,   4.0068800)
  factor( 2) = f_factor_struct(12.123200,  -0.999900,   4.0001300)
  factor( 3) = f_factor_struct(17.817800,  -0.999900,   2.1237800)
  factor( 4) = f_factor_struct(21.600900,  -0.999900,   0.6025340)
  factor( 5) = f_factor_struct(28.833700,  -0.999900,   0.3706440)
  factor( 6) = f_factor_struct(29.300000,   3.838100,   0.3717420)
  factor( 7) = f_factor_struct(76.728900,   2.097630,   0.4773930)
  factor( 8) = f_factor_struct(85.849100,   1.136100,   0.4159400)
  factor( 9) = f_factor_struct(91.539500,   0.232650,   0.4449190)
  factor(10) = f_factor_struct(96.053500,  -1.001990,   0.4806420)
  factor(11) = f_factor_struct(97.607100,  -1.779850,   0.4931760)
  factor(12) = f_factor_struct(99.185800,  -3.569080,   0.5060360)
  factor(13) = f_factor_struct(99.700000,  -6.180010,   0.5102520)
  factor(14) = f_factor_struct(99.900000,  -6.058740,   4.6074200)
  factor(15) = f_factor_struct(100.79000,  -2.494140,   4.3679000)
  factor(16) = f_factor_struct(102.42000,  -1.410120,   3.9658400)
  factor(17) = f_factor_struct(104.07700,  -1.590760,   3.6734400)
  factor(18) = f_factor_struct(105.76000,  -1.863930,   4.0673800)
  factor(19) = f_factor_struct(110.97500,  -1.907960,   5.4489100)
  factor(20) = f_factor_struct(116.44800,  -1.786350,   7.1534100)
  factor(21) = f_factor_struct(126.17500,   0.525660,   8.9108500)
  factor(22) = f_factor_struct(134.53800,   2.683200,   9.5984200)
  factor(23) = f_factor_struct(143.45600,   3.126990,   9.2460900)
  factor(24) = f_factor_struct(152.96400,   4.181810,  10.0588000)
  factor(25) = f_factor_struct(157.95300,   4.417420,   9.8635000)
  factor(26) = f_factor_struct(168.42200,   6.250620,  10.9311000)
  factor(27) = f_factor_struct(179.58600,   7.683370,   9.4475800)
  factor(28) = f_factor_struct(232.14700,  10.217500,   8.2066400)
  factor(29) = f_factor_struct(300.09200,  12.002900,   6.4372400)
  factor(30) = f_factor_struct(387.92200,  13.019300,   4.7842000)
  factor(31) = f_factor_struct(501.45900,  13.376000,   3.3132500)
  factor(32) = f_factor_struct(648.22600,  13.372100,   2.2138500)
  factor(33) = f_factor_struct(952.71500,  13.059900,   1.1738300)
  factor(34) = f_factor_struct(1231.5500,  12.642000,   0.7399720)
  factor(35) = f_factor_struct(1592.0100,  11.814500,   0.4776830)
  factor(36) = f_factor_struct(1697.5300,  11.273600,   0.4258760)
  factor(37) = f_factor_struct(1752.8900,  10.764400,   0.4016720)
  factor(38) = f_factor_struct(1781.2400,  10.337500,   0.3897180)
  factor(39) = f_factor_struct(1810.0500,   9.569810,   0.3780020)
  factor(40) = f_factor_struct(1838.8000,   2.834190,   0.3668380)
  factor(41) = f_factor_struct(1839.0000,   2.835870,   4.1573400)
  factor(42) = f_factor_struct(1839.3200,   4.580550,   4.1563200)
  factor(43) = f_factor_struct(1869.0700,   9.825010,   4.0645200)
  factor(44) = f_factor_struct(1899.3000,  10.724100,   3.9747500)
  factor(45) = f_factor_struct(1930.0200,  11.267200,   3.8868800)
  factor(46) = f_factor_struct(2265.9200,  13.299100,   3.0795600)
  factor(47) = f_factor_struct(2576.2600,  13.885700,   2.5186900)
  factor(48) = f_factor_struct(2929.1100,  14.190600,   2.0408500)
  factor(49) = f_factor_struct(4894.6000,  14.423300,   0.8266710)
  factor(50) = f_factor_struct(6327.1500,  14.360900,   0.5131970)
  factor(51) = f_factor_struct(10572.800,  14.209300,   0.1915190)
  factor(52) = f_factor_struct(17667.400,  14.062300,   0.0687904)
  factor(53) = f_factor_struct(29522.500,  14.020800,   0.0236346)
  factor(54) = f_factor_struct(43000.000,  14.007800,   0.0106303)
  factor(55) = f_factor_struct(64000.000,  13.996200,   0.0044384)
  factor(56) = f_factor_struct(86500.000,  13.991300,   0.0022794)
  factor(57) = f_factor_struct(100000.00,  13.989700,   0.00165265)
  n_factor = 57

else

  call out_io (s_fatal$, r_name, 'BAD CRYSTAL TYPE: ' // ele%component_name)

endif

! Calculate HKL

read (hkl_str(1:1), *, iostat = ios) hkl(1)
if (hkl_str(1:1) == '' .or. ios /= 0) return

read (hkl_str(2:2), *, iostat = ios) hkl(2)
if (hkl_str(2:2) == '' .or. ios /= 0) return

read (hkl_str(3:3), *, iostat = ios) hkl(3)
if (hkl_str(3:3) == '' .or. ios /= 0) return

! Calculate values for the given energy

ele%value(v_unitcell$) = a0**3

d = a0 / sqrt(sum(hkl**2))
ele%value(d_spacing$) = d

! Interpolate ATOMIC structure factors from factor struct, which are evaluated at hkl=0

call bracket_index (factor%energy, 1, n_factor, ele%value(e_tot$), ix)
if (ix == 0) then
  f_re = factor(1)%f0_re
  f_im = factor(1)%f0_im
elseif (ix == n_factor) then
  f_re = factor(n_factor)%f0_re
  f_im = factor(n_factor)%f0_im
else
  f_re = (factor(ix)%f0_re   * (factor(ix+1)%energy - ele%value(e_tot$)) + &
         factor(ix+1)%f0_re * (ele%value(e_tot$) - factor(ix)%energy)) / (factor(ix+1)%energy - factor(ix)%energy)
  f_im = (factor(ix)%f0_im   * (factor(ix+1)%energy - ele%value(e_tot$)) + &
         factor(ix+1)%f0_im * (ele%value(e_tot$) - factor(ix)%energy)) / (factor(ix+1)%energy - factor(ix)%energy)
endif

!Multiply by number of atoms to get CRYSTAL structure factor, ie f0_re and f0_im

ele%value(f0_re$) = n_atom * (f_re)
ele%value(f0_im$) = n_atom * (f_im)

! Calculate the hkl dependent part of f
f_h = sum(a(1:n)) + c

! Subtract from f(hkl=0) to obtain the energy dependent part
f_e = f_re - f_h

! Calculate hkl dependent part of f
! factor of 1e20 due to units of d
f_h = sum( a(1:n) * exp( -b(1:n)/(4*d*d*1e20) ) ) + c

atomsum = 0
do i = 1, n_atom
  test = hkl * r_atom(i)%vec
  test(1) = sum(hkl * r_atom(i)%vec)
  arg = cmplx(0.0,twopi * sum(hkl * r_atom(i)%vec))
  atomsum = atomsum + exp(arg)
enddo

fh_re = abs(atomsum)*(f_h+f_e)
fh_im = abs(atomsum)*(f_im)

ele%value(fh_re$) = fh_re
ele%value(fh_im$) = fh_im

! Set bragg and alpha angles to zero if bragg condition is not satisfied.

r = ele%value(ref_wave_length$) / (2 * d)

if (r < 1) then
  ele%value(bragg_angle$) = asin(r)
  bp = ele%value(b_param$)
  ele%value(alpha_angle$) = atan(tan(ele%value(bragg_angle$)) * (bp + 1) / (bp - 1))
else
  ele%value(bragg_angle$) = 0
  ele%value(alpha_angle$) = 0
endif


err_flag = .false.

end subroutine 
