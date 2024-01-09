program xraylib_test

use xraylib_interface
use xraylib, dummy => r_e

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele

real(rp) absorption, phase_shift, energy

integer i, ix, n

character (kind=c_char, len=nist_list_string_length), pointer :: nistcompounds(:)
character(60) compound, out_str

logical err_flag

!

open (1, file = 'output.now')

!

energy = 1e3

call photon_absorption_and_phase_shift ('B4C', Energy, absorption, phase_shift, err_flag)
write (1, '(a, es18.9)') '"B4C-E1"   REL 1E-8', absorption
write (1, '(a, es18.9)') '"B4C-E1"   REL 1E-8', phase_shift

call photon_absorption_and_phase_shift ('B4C', 2*Energy, absorption, phase_shift, err_flag)
write (1, '(a, es18.9)') '"B4C-E2"   REL 1E-8', absorption
write (1, '(a, es18.9)') '"B4C-E2"   REL 1E-8', phase_shift

call photon_absorption_and_phase_shift ('NaCl', Energy, absorption, phase_shift, err_flag)
write (1, '(a, es18.9)') '"NaCl-E1"   REL 1E-8', absorption
write (1, '(a, es18.9)') '"NaCl-E1"   REL 1E-8', phase_shift

call photon_absorption_and_phase_shift ('NaCl', 2*Energy, absorption, phase_shift, err_flag)
write (1, '(a, es18.9)') '"NaCl-E2"   REL 1E-8', absorption
write (1, '(a, es18.9)') '"NaCl-E2"   REL 1E-8', phase_shift

call photon_absorption_and_phase_shift ('Fe', Energy, absorption, phase_shift, err_flag)
write (1, '(a, es18.9)') '"Fe-E1"   REL 1E-8', absorption
write (1, '(a, es18.9)') '"Fe-E1"   REL 1E-8', phase_shift

call photon_absorption_and_phase_shift ('Fe', 2*Energy, absorption, phase_shift, err_flag)
write (1, '(a, es18.9)') '"Fe-E2"   REL 1E-8', absorption
write (1, '(a, es18.9)') '"Fe-E2"   REL 1E-8', phase_shift
n = atomic_number(species_id("Fe"))
write (1, '(a, es18.9)') '"Fe-Weight"  REL 1E-8', atomicweight(n)
write (1, '(a, es18.9)') '"Fe-Density" REL 1E-8', elementdensity(n)

!

out_str = '"OK"'

nistCompounds => GetCompoundDataNISTList()
do i = 1, size(nistcompounds)
  compound = nistcompounds(i)
  call upcase_string(compound)
  ix = index(compound, '(');  if (ix /= 0) compound = compound(1:ix-1) // compound(ix+1:)
  ix = index(compound, ')');  if (ix /= 0) compound = compound(1:ix-1) // compound(ix+1:)
  ix = index(compound, ',');  if (ix /= 0) compound = compound(1:ix-1) // compound(ix+1:)
  ix = index(compound, ',');  if (ix /= 0) compound = compound(1:ix-1) // compound(ix+1:)

  n = len_trim(compound)

  do
    ix = index(compound, ' ')
    if (ix > n) exit
    compound = compound(1:ix-1) // '_' // compound(ix+1:)
  enddo

  do
    ix = index(compound, '-')
    if (ix == 0) exit
    compound = compound(1:ix-1) // '_' // compound(ix+1:)
  enddo

  do
    ix = index(compound, '/')
    if (ix == 0) exit
    compound = compound(1:ix-1) // '_' // compound(ix+1:)
  enddo

  if (i /= xraylib_nist_compound(compound) + 1) then
    out_str = '"BAD: ' // trim(compound) // '"'
    exit
  endif
enddo

write (1, '(2a)') '"Compound_Names" STR ', trim(out_str)

!

call bmad_parser ('xraylib.bmad', lat)

ele => lat%ele(1)    ! crystal
write (1, '(a, f16.10)') '"Bragg_In"   REL 1E-8', ele%value(bragg_angle_in$)
write (1, '(a, f16.10)') '"Bragg_Out"  REL 1E-8', ele%value(bragg_angle_out$)
write (1, '(a, f16.10)') '"Alpha_Ang"  REL 1E-8', ele%value(alpha_angle$)
write (1, '(a, f16.10)') '"F0_Re"      REL 1E-8', real(ele%photon%material%f_0)
write (1, '(a, f16.10)') '"F0_Im"      REL 1E-8', aimag(ele%photon%material%f_0)
write (1, '(a, es20.8)') '"FH_Re"      REL 1E-8', real(ele%photon%material%f_h)
write (1, '(a, es20.8)') '"FH_Im"      REL 1E-8', aimag(ele%photon%material%f_h)
write (1, '(a, es20.8)') '"Darwin_Sig" REL 1E-8', ele%value(darwin_width_sigma$)
write (1, '(a, es20.8)') '"Darwin_Pi"  REL 1E-8', ele%value(darwin_width_pi$)
write (1, '(a, es20.8)') '"D_Spacing"  REL 1E-8', ele%value(d_spacing$)

ele => lat%ele(2)    ! multilayer_mirror
write (1, '(a, es20.8)') '"Graze_Ang"  REL 1E-8', ele%value(graze_angle$)
write (1, '(a, es20.8)') '"V1_Cell"    REL 1E-8', ele%value(v1_unitcell$)
write (1, '(a, es20.8)') '"F0_Re1"     REL 1E-8', real(ele%photon%material%f0_m1)
write (1, '(a, es20.8)') '"F0_Im2"     REL 1E-8', aimag(ele%photon%material%f0_m2)

!

end program
