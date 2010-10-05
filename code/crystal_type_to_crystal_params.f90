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

type atom_vec_struct
  real(rp) vec(3)
end type

type (ele_struct) ele
type (atom_vec_struct) r_atom(8), f1f2(100)

real(rp) f0_e, f0_h, f0_re, f0_im, fh_re, fh_im
real(rp) a(10), b(10), c, d, a0, hkl(3), bp, r, test(3)
complex(rp) arg, atomsum

integer i, ix, n, ios, n_atom, n_f1f2

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

  f1f2( 1)%vec = [10.000000,  -0.999900,   4.0068800]
  f1f2( 2)%vec = [12.123200,  -0.999900,   4.0001300]
  f1f2( 3)%vec = [17.817800,  -0.999900,   2.1237800]
  f1f2( 4)%vec = [21.600900,  -0.999900,   0.6025340]
  f1f2( 5)%vec = [28.833700,  -0.999900,   0.3706440]
  f1f2( 6)%vec = [29.300000,   3.838100,   0.3717420]
  f1f2( 7)%vec = [76.728900,   2.097630,   0.4773930]
  f1f2( 8)%vec = [85.849100,   1.136100,   0.4159400]
  f1f2( 9)%vec = [91.539500,   0.232650,   0.4449190]
  f1f2(10)%vec = [96.053500,  -1.001990,   0.4806420]
  f1f2(11)%vec = [97.607100,  -1.779850,   0.4931760]
  f1f2(12)%vec = [99.185800,  -3.569080,   0.5060360]
  f1f2(13)%vec = [99.700000,  -6.180010,   0.5102520]
  f1f2(14)%vec = [99.900000,  -6.058740,   4.6074200]
  f1f2(15)%vec = [100.79000,  -2.494140,   4.3679000]
  f1f2(16)%vec = [102.42000,  -1.410120,   3.9658400]
  f1f2(17)%vec = [104.07700,  -1.590760,   3.6734400]
  f1f2(18)%vec = [105.76000,  -1.863930,   4.0673800]
  f1f2(19)%vec = [110.97500,  -1.907960,   5.4489100]
  f1f2(20)%vec = [116.44800,  -1.786350,   7.1534100]
  f1f2(21)%vec = [126.17500,   0.525660,   8.9108500]
  f1f2(22)%vec = [134.53800,   2.683200,   9.5984200]
  f1f2(23)%vec = [143.45600,   3.126990,   9.2460900]
  f1f2(24)%vec = [152.96400,   4.181810,  10.0588000]
  f1f2(25)%vec = [157.95300,   4.417420,   9.8635000]
  f1f2(26)%vec = [168.42200,   6.250620,  10.9311000]
  f1f2(27)%vec = [179.58600,   7.683370,   9.4475800]
  f1f2(28)%vec = [232.14700,  10.217500,   8.2066400]
  f1f2(29)%vec = [300.09200,  12.002900,   6.4372400]
  f1f2(30)%vec = [387.92200,  13.019300,   4.7842000]
  f1f2(31)%vec = [501.45900,  13.376000,   3.3132500]
  f1f2(32)%vec = [648.22600,  13.372100,   2.2138500]
  f1f2(33)%vec = [952.71500,  13.059900,   1.1738300]
  f1f2(34)%vec = [1231.5500,  12.642000,   0.7399720]
  f1f2(35)%vec = [1592.0100,  11.814500,   0.4776830]
  f1f2(36)%vec = [1697.5300,  11.273600,   0.4258760]
  f1f2(37)%vec = [1752.8900,  10.764400,   0.4016720]
  f1f2(38)%vec = [1781.2400,  10.337500,   0.3897180]
  f1f2(39)%vec = [1810.0500,   9.569810,   0.3780020]
  f1f2(40)%vec = [1838.8000,   2.834190,   0.3668380]
  f1f2(41)%vec = [1839.0000,   2.835870,   4.1573400]
  f1f2(42)%vec = [1839.3200,   4.580550,   4.1563200]
  f1f2(43)%vec = [1869.0700,   9.825010,   4.0645200]
  f1f2(44)%vec = [1899.3000,  10.724100,   3.9747500]
  f1f2(45)%vec = [1930.0200,  11.267200,   3.8868800]
  f1f2(46)%vec = [2265.9200,  13.299100,   3.0795600]
  f1f2(47)%vec = [2576.2600,  13.885700,   2.5186900]
  f1f2(48)%vec = [2929.1100,  14.190600,   2.0408500]
  f1f2(49)%vec = [4894.6000,  14.423300,   0.8266710]
  f1f2(50)%vec = [6327.1500,  14.360900,   0.5131970]
  f1f2(51)%vec = [10572.800,  14.209300,   0.1915190]
  f1f2(52)%vec = [17667.400,  14.062300,   0.0687904]
  f1f2(53)%vec = [29522.500,  14.020800,   0.0236346]
  f1f2(54)%vec = [43000.000,  14.007800,   0.0106303]
  f1f2(55)%vec = [64000.000,  13.996200,   0.0044384]
  f1f2(56)%vec = [86500.000,  13.991300,   0.0022794]
  f1f2(57)%vec = [100000.00,  13.989700,   0.00165265]
  n_f1f2 = 57

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

! Interpolate f1 and f2 table, which should really be called the f0+f1 and f2 table
! Note throughout this code, f0_re, f0_im, fh_re, fh_im refer to ATOMIC scattering
! form factors instead of the CRYSTAL structure factors, unlike in the ele struct

call bracket_index (f1f2%vec(1), 1, n_f1f2, ele%value(e_tot$), ix)
if (ix == 0) then
  f0_re = f1f2(1)%vec(2)
  f0_im = f1f2(1)%vec(3)
elseif (ix == n_f1f2) then
  f0_re = f1f2(n_f1f2)%vec(2)
  f0_im = f1f2(n_f1f2)%vec(3)
else
  f0_re = (f1f2(ix)%vec(2)   * (f1f2(ix+1)%vec(1) - ele%value(e_tot$)) + &
         f1f2(ix+1)%vec(2) * (ele%value(e_tot$) - f1f2(ix)%vec(1))) / (f1f2(ix+1)%vec(1) - f1f2(ix)%vec(1))
  f0_im = (f1f2(ix)%vec(3)   * (f1f2(ix+1)%vec(1) - ele%value(e_tot$)) + &
         f1f2(ix+1)%vec(3) * (ele%value(e_tot$) - f1f2(ix)%vec(1))) / (f1f2(ix+1)%vec(1) - f1f2(ix)%vec(1))
endif

!Multiply by number of atoms to get CRYSTAL STRUCTURE FACTORS

ele%value(f0_re$) = n_atom * (f0_re)
ele%value(f0_im$) = n_atom * (f0_im)

! Calculate the hkl dependent part of the f evaluated at f0
f0_h = sum(a(1:n)) + c

! Subtract to obtain the energy dependent part
f0_e = f0_re - f0_h

! Calculate the new hkl dependent part
f0_h = sum( a(1:n) * exp( -b(1:n)/(4*d*d*1e20) ) ) + c

atomsum = 0
do i = 1, n_atom
  test = hkl * r_atom(i)%vec
  test(1) = sum(hkl * r_atom(i)%vec)
  arg = cmplx(0.0,twopi * sum(hkl * r_atom(i)%vec))
  atomsum = atomsum + exp(arg)
enddo

fh_re = abs(atomsum)*(f0_h+f0_e)
fh_im = abs(atomsum)*(f0_im)

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
