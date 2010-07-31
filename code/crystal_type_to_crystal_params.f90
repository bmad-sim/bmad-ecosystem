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
type (atom_vec_struct) r_atom(8)

real(rp) f_h_factor, f_0_factor, f_0h, f_1, f_2, fh_re, fh_im
real(rp) a(10), b(10), c, d, n, a0, hkl(3), arg

integer i, ios, n_atom

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
  a0 = 5.4309d-10
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

f_0h = sum(a(1:n)*exp(-b(1:n)/(4*d*d))) + c

!  f_1 = interpolate table
!  f_2 = interpolate table

ele%value(f0_re$) = n_atom * (f_0h + f_1)
ele%value(f0_im$) = n_atom * f_2

fh_re = 0
fh_im = 0

do i = 1, n_atom
  arg = sum(hkl * r_atom(i)%vec)
  fh_re = fh_re + cos(arg) * (f_0h + f_1) - sin(arg) * f_2
  fh_im = fh_im + cos(arg) * f_2 + sin(arg) * (f_0h + f_1)
enddo

ele%value(fh_re$) = fh_re
ele%value(fh_im$) = fh_im

err_flag = .false.

end subroutine 
