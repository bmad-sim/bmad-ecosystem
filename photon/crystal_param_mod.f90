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

type (atomic_f_struct), target, save :: f0_b(1), f0_c(1), f0_si(57), f0_w(1)
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
endif

! Interpolate atomic structure factors from factor struct, which are evaluated at hkl=0

select case (atom)
case ('B')
  f0_ptr => f0_c
case ('C')
  f0_ptr => f0_c
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

