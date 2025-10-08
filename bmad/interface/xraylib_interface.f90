!+
! Interface to the XRAYLIB library
!
! This replaces the routines in crystal_param_mod.f90
!-

module xraylib_interface

use bmad_routine_interface

integer, parameter :: xraylib_z_max$ = 98  ! Maximum atomic Z value of tables.

contains

!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
!+
! Subroutine photon_absorption_and_phase_shift (material, Energy, absorption, phase_shift, err_flag)
!
! Routine to calcualte the absorption and phase shift values for a photon with a given 
! energy going through a particular material.
!
! Input:
!   material -- character(*): Material name.
!   Energy   -- real(rp): Photon energy (eV).
!
! Output:
!   absorption    -- real(rp): E_field ~ Exp(-absorption * length)
!   phase_shift   -- real(rp): E_field Phase shift (radians) per unit length relative to vacuum.
!   err_flag      -- logical: Set true if material not recognized.
!-

subroutine photon_absorption_and_phase_shift (material, Energy, absorption, phase_shift, err_flag)

use xraylib, dummy => r_e

implicit none

type (crystal_struct), pointer :: cryst
type (compoundDataNIST), pointer :: compound

character(*) material

real(rp) Energy, absorption, phase_shift
real(rp) wavelength, factor, volume, number_fraction
real(c_double) debye_temp_factor, f0, fp, fpp, rel_angle, q, E_kev

complex(rp) f0_tot

integer n, ix, int_err

logical err_flag

character(*), parameter :: r_name = 'photon_absorption_and_phase_shift'

!

err_flag = .true.
E_kev = Energy * 1d-3
debye_temp_factor = 1.0          
q = 0   
wavelength = c_light * h_planck / energy

! Is this a crystal?

cryst => Crystal_GetCrystal (material)
if (associated(cryst)) then
  f0_tot = Crystal_F_H_StructureFactor (cryst, E_kev, 0, 0, 0, debye_temp_factor, rel_angle)
  factor = r_e * wavelength  / (1d-30 * cryst%volume)
  phase_shift = factor * real(f0_tot)
  absorption = factor * aimag(f0_tot)
  return
endif

! Is this an element?

do n = 1, xraylib_z_max$
  if (material /= AtomicNumberToSymbol(n)) cycle
  volume = 1d-6 * AtomicWeight(n) / (N_avogadro * ElementDensity(n))
  factor = r_e * wavelength  / volume
  int_err =  Atomic_Factors (n, E_kev, q, debye_temp_factor, f0, fp, fpp)
  if (int_err == 0) then
    call out_io (s_error$, r_name, 'PROBLEM RETRIEVING FORM FACTORS FOR ELEMENT: ' // material)
    return
  endif
  phase_shift = factor * (f0 + fp)
  absorption = factor * fpp
  err_flag = .false.
  return
enddo

! Is this a NIST material?

ix = xraylib_nist_compound(material)
if (ix > -1) then
  compound => GetCompoundDataNISTByIndex(ix)

  f0_tot = 0
  do n = 1, compound%nElements
    number_fraction = compound%massFractions(n) / AtomicWeight(compound%elements(n))
    int_err = Atomic_Factors (compound%elements(n), E_kev, q, debye_temp_factor, f0, fp, fpp)
    if (int_err == 0) then
      call out_io (s_error$, r_name, 'PROBLEM RETRIEVING FORM FACTORS FOR ELEMENT IN COMPOUND: ' // material)
      return
    endif
    f0_tot = f0_tot + cmplx(f0 + fp, fpp, rp) * number_fraction
  enddo

  volume = 1d-6 / (N_avogadro * compound%density)
  factor = r_e * wavelength  / volume
  phase_shift = factor * real(f0_tot)
  absorption = factor * aimag(f0_tot)

  !! call FreeCompoundDataNIST(compound)  Used in Version 3.1 but removed for 3.3
  err_flag = .false.
  return
endif

! Fail

call out_io (s_error$, r_name, 'MATERIAL NOT RECOGNIZED: ' // trim(material))

end subroutine photon_absorption_and_phase_shift

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
! Where "AAA" is the atomic formula for the top layer crystal and "BBB" is the second layer atomic formula.
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
!     %photon%material%f0_m1
!     %photon%material%f0_m2
!   err_flag -- Logical: Set True if multilayer type is unrecognized. False otherwise.
!-

subroutine multilayer_type_to_multilayer_params (ele, err_flag)

use xraylib, dummy => r_e

implicit none

type (ele_struct), target :: ele

real(rp) D, sin_theta0, del
real(rp), pointer :: v(:)

complex(rp) xi1, xi2, f0

integer ix

logical err_flag

character(*), parameter :: r_name = 'multilayer_type_to_multilayer_params'

! get types

err_flag = .false.
if (ele%value(E_tot$) == 0) return ! Can happen when lattice parsing

ix = index(ele%component_name, ':')
if (ix == 0) return
v => ele%value

call load_layer_params (ele%component_name(1:ix-1), ele%photon%material%f0_m1, v(v1_unitcell$), err_flag)
if (err_flag) return

call load_layer_params (ele%component_name(ix+1:), ele%photon%material%f0_m2, v(v2_unitcell$), err_flag)
if (err_flag) return

! Calc graze angle. See Kohn Eq 36.

D = v(d1_thickness$) + v(d2_thickness$)
if (D == 0) then
  call out_io (s_error$, r_name, 'MULTILAYER_MIRROR HAS ZERO THICKNESS: ' // ele%name)
  err_flag = .true.
  return
endif

xi1 = -conjg(ele%photon%material%f0_m1) * v(ref_wavelength$)**2 * r_e / (pi * ele%value(v1_unitcell$))
xi2 = -conjg(ele%photon%material%f0_m2) * v(ref_wavelength$)**2 * r_e / (pi * ele%value(v2_unitcell$)) 

sin_theta0 = v(ref_wavelength$) / (2 * D)
del = (v(d1_thickness$) * real(sqrt(sin_theta0**2 + xi1) - sin_theta0) + &
       v(d2_thickness$) * real(sqrt(sin_theta0**2 + xi2) - sin_theta0)) / D
v(graze_angle$) = asin(sin_theta0 - del)

!-----------------------------------------------------------------------------------------
contains

subroutine load_layer_params (material_name, f0, v_unitcell, err_flag)

type (crystal_struct), pointer :: cryst
type (compoundDataNIST), pointer :: compound

real(rp) v_unitcell
real(rp)  density, number_fraction_min, number_fraction
real(c_double) energy, debye_temp_factor, f0_atom, fp, fpp, rel_angle, q

complex(rp) f0

integer n, n_atom, int_err
logical err_flag

character(*) material_name

!

energy = ele%value(E_tot$) * 1d-3
debye_temp_factor = 1.0
q = 0
err_flag = .true.

! Is this a crystal?

cryst => Crystal_GetCrystal (material_name)
if (associated(cryst)) then
  f0 = Crystal_F_H_StructureFactor (cryst, energy, 0, 0, 0, debye_temp_factor, rel_angle)
  v_unitcell = 1d-30 * cryst%volume
  return
endif

! Is this an element?

do n = 1, xraylib_z_max$
  if (material_name /= AtomicNumberToSymbol(n)) cycle
  v_unitcell = 1d-6 * AtomicWeight(n) / (N_avogadro * ElementDensity(n))
  int_err = Atomic_Factors (n, energy, q, debye_temp_factor, f0_atom, fp, fpp)
  if (int_err == 0) then
    call out_io (s_error$, r_name, 'PROBLEM RETRIEVING FORM FACTORS FOR ELEMENT: ' // material_name)
    return
  endif
  f0 = cmplx(f0_atom + fp, fpp, rp)
  err_flag = .false.
  return
enddo

! Is this a NIST material?
! The size of the "unit cell" is arbitrarily taken to be 1/number_fraction_min.
! The reason why this scaling is used is that f0 and v_unitcell are stored in the ele structure
! and this scaling makes the numbers look reasonable. This does not affect the physics.

ix = xraylib_nist_compound(material_name)
if (ix > -1) then
  compound => GetCompoundDataNISTByIndex(ix)

  f0 = 0
  number_fraction_min = 1d10  ! something larget

  do n = 1, compound%nElements
    number_fraction = compound%massFractions(n) / AtomicWeight(compound%elements(n))
    int_err = Atomic_Factors (compound%elements(n), energy, q, debye_temp_factor, f0_atom, fp, fpp)
    if (int_err == 0) then
      call out_io (s_error$, r_name, 'PROBLEM RETRIEVING FORM FACTORS FOR ELEMENT IN COMPOUND: ' // material_name)
      return
    endif
    f0 = f0 + cmplx(f0_atom + fp, fpp, rp) * number_fraction
    number_fraction_min = min(number_fraction_min, number_fraction)
  enddo
  f0 = f0 / number_fraction_min
  v_unitcell = 1d-6 / (N_avogadro * compound%density)  / number_fraction_min
  !! call FreeCompoundDataNIST(compound)  Used in Version 3.1 but removed for 3.3
  err_flag = .false.
  return
endif

! Fail

call out_io (s_error$, r_name, 'MATERIAL NOT RECOGNIZED: ' // trim(material_name))

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
!   "ZZZ(ijk)"
! Where "ZZZ" is the atomic formula of the crystal material and "ijk" is the reciprical lattice 
! vetor specifying the diffraction plans.
!
! Input:
!   ele      -- ele_struct: Crystal element.
!     %component_name -- Character: Crystal type name. Assumed upper case.
!                          A blank name is not an error and results in nothing set.
!     %value(e_tot$)  -- Photon energy in eV.
!
! Output:
!   ele      -- ele_struct: Crystal element with computed parameter..
!   err_flag -- Logical: Set True if crystal type is unrecognized. False otherwise.
!-

subroutine crystal_type_to_crystal_params (ele, err_flag)

use xraylib

implicit none

type (ele_struct) ele
type (crystal_struct), pointer :: cryst

real(rp) bp, r
real(c_double) energy, debye_temp_factor, fp, fpp, rel_angle

complex(rp) f0, fh

integer i, ix, ix1, ix2, ios, hkl(3), ios1, ios2, ios3

logical err_flag

character(*), parameter :: r_name = 'crystal_type_to_crystal_params'
character(16) hkl_str, cyrstal_type

! Check if component_name and e_tot are set.

err_flag = .false.

if (ele%component_name == '') return
if (ele%value(e_tot$) == 0) return

err_flag = .true.

! Separate Miller indices

ix1 = index(ele%component_name, '(')
ix2 = len_trim(ele%component_name)
if (ix1 == 0 .or. ele%component_name(ix2:ix2) /= ')') then
  call out_io (s_fatal$, r_name, 'MALFORMED CRYSTAL_TYPE: ' // ele%component_name, 'FOR ELEMENT: ' // ele%name)
  if (global_com%exit_on_error) call err_exit
endif

cyrstal_type = ele%component_name(1:ix1-1)
hkl_str = ele%component_name(ix1+1:ix2-1)

! Load constants for the particular crystal...

! Calculate HKL

ix = index(hkl_str, ',')
if (ix == 0) then
  if (.not. is_integer(hkl_str) .or. len_trim(hkl_str) /= 3) then
    call out_io (s_fatal$, r_name, 'MALFORMED CRYSTAL_TYPE: ' // ele%component_name, 'FOR ELEMENT: ' // ele%name)
    if (global_com%exit_on_error) call err_exit
    return
  endif
  read (hkl_str(1:1), *, iostat = ios1) hkl(1)
  read (hkl_str(2:2), *, iostat = ios2) hkl(2)
  read (hkl_str(3:3), *, iostat = ios3) hkl(3)
  if (ios1 /= 0 .or. ios2 /= 0 .or. ios3 /= 0) then
    call out_io (s_fatal$, r_name, 'MALFORMED CRYSTAL_TYPE: ' // ele%component_name, 'FOR ELEMENT: ' // ele%name)
    if (global_com%exit_on_error) call err_exit
    return
  endif

else
  do i = 1, 3
    if (ix == 0) then
      call out_io (s_fatal$, r_name, 'MALFORMED CRYSTAL_TYPE: ' // ele%component_name, 'FOR ELEMENT: ' // ele%name)
      if (global_com%exit_on_error) call err_exit
      return
    endif

    read(hkl_str(1:ix-1), *, iostat = ios) hkl(i)
    if (ios /= 0) then
      call out_io (s_fatal$, r_name, 'MALFORMED CRYSTAL_TYPE: ' // ele%component_name, 'FOR ELEMENT: ' // ele%name)
      if (global_com%exit_on_error) call err_exit
      return
    endif

    if (i == 3) exit
    hkl_str = hkl_str(ix+1:)
    ix = index(hkl_str, ',')
    if (i == 2) ix = len(hkl_str)
  enddo
endif

! Calculate values for the given energy

energy = ele%value(E_tot$) * 1d-3
debye_temp_factor = 1
rel_angle = 1

cryst => Crystal_GetCrystal (cyrstal_type)
if (.not. associated(cryst)) then
  call out_io(s_error$, r_name, 'Crystal_type not recognized by xraylib: ' // cyrstal_type)
  return
endif

ele%value(v_unitcell$) = 1d-30 * cryst%volume
ele%value(d_spacing$) = 1d-10 * Crystal_dSpacing(cryst, hkl(1), hkl(2), hkl(3))

ele%photon%material%f_0    = Crystal_F_H_StructureFactor (cryst, energy, 0, 0, 0, debye_temp_factor, rel_angle)
ele%photon%material%f_h    = Crystal_F_H_StructureFactor (cryst, energy, hkl(1), hkl(2), hkl(3), debye_temp_factor, rel_angle)
ele%photon%material%f_hbar = Crystal_F_H_StructureFactor (cryst, energy, -hkl(1), -hkl(2), -hkl(3), debye_temp_factor, rel_angle)
ele%photon%material%f_hkl = sqrt(ele%photon%material%f_h * ele%photon%material%f_hbar)

! Set bragg and alpha angles to zero if bragg condition is not satisfied.

r = ele%value(ref_wavelength$) / (2 * ele%value(d_spacing$))

if (r < 1) then
  ele%value(bragg_angle$) = asin(r)
  bp = ele%value(b_param$)
  ! These two formulas for alpha_angle are the same except where the branch cut is.
  if (bp < 0) then  ! Bragg
    ele%value(alpha_angle$) = atan2(-sin(ele%value(bragg_angle$)) * (1 + bp), cos(ele%value(bragg_angle$)) * (1 - bp))
  else
    ele%value(alpha_angle$) = pi/2 + atan2(cos(ele%value(bragg_angle$)) * (1 - bp), sin(ele%value(bragg_angle$)) * (1 + bp))
  endif
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
! Function xraylib_nist_compound (name) result (indx)
!
! Routine to return the xraylib index for a given NIST compound.
! Taken from file xraylib/include/xraylib-nist_compounds.h
!
! Input:
!   name  -- character(*): Name of compound
!
! Output:
!   indx  -- integer: Compound index. -1 if not found.
!-

function xraylib_nist_compound (name) result (indx)

implicit none

character(*) name
integer indx

! Full name is NIST_COMPOUND_XXX

select case (name)
case default;                                        indx = -1
case ('A_150_TISSUE_EQUIVALENT_PLASTIC');            indx = 0
case ('ACETONE');                                    indx = 1
case ('ACETYLENE');                                  indx = 2
case ('ADENINE');                                    indx = 3
case ('ADIPOSE_TISSUE_ICRP');                        indx = 4
case ('AIR_DRY_NEAR_SEA_LEVEL');                     indx = 5
case ('ALANINE');                                    indx = 6
case ('ALUMINUM_OXIDE', 'Al2O3');                    indx = 7
case ('AMBER');                                      indx = 8
case ('AMMONIA', 'NH3');                             indx = 9
case ('ANILINE');                                    indx = 10
case ('ANTHRACENE');                                 indx = 11
case ('B_100_BONE_EQUIVALENT_PLASTIC');              indx = 12
case ('BAKELITE');                                   indx = 13
case ('BARIUM_FLUORIDE');                            indx = 14
case ('BARIUM_SULFATE');                             indx = 15
case ('BENZENE', 'C6H6');                            indx = 16
case ('BERYLLIUM_OXIDE');                            indx = 17
case ('BISMUTH_GERMANIUM_OXIDE');                    indx = 18
case ('BLOOD_ICRP');                                 indx = 19
case ('BONE_COMPACT_ICRU');                          indx = 20
case ('BONE_CORTICAL_ICRP');                         indx = 21
case ('BORON_CARBIDE', 'B4C');                       indx = 22
case ('BORON_OXIDE', 'B2O3');                        indx = 23
case ('BRAIN_ICRP');                                 indx = 24
case ('BUTANE');                                     indx = 25
case ('N_BUTYL_ALCOHOL');                            indx = 26
case ('C_552_AIR_EQUIVALENT_PLASTIC');               indx = 27
case ('CADMIUM_TELLURIDE');                          indx = 28
case ('CADMIUM_TUNGSTATE');                          indx = 29
case ('CALCIUM_CARBONATE');                          indx = 30
case ('CALCIUM_FLUORIDE');                           indx = 31
case ('CALCIUM_OXIDE');                              indx = 32
case ('CALCIUM_SULFATE');                            indx = 33
case ('CALCIUM_TUNGSTATE');                          indx = 34
case ('CARBON_DIOXIDE');                             indx = 35
case ('CARBON_TETRACHLORIDE');                       indx = 36
case ('CELLULOSE_ACETATE_CELLOPHANE');               indx = 37
case ('CELLULOSE_ACETATE_BUTYRATE');                 indx = 38
case ('CELLULOSE_NITRATE');                          indx = 39
case ('CERIC_SULFATE_DOSIMETER_SOLUTION');           indx = 40
case ('CESIUM_FLUORIDE');                            indx = 41
case ('CESIUM_IODIDE');                              indx = 42
case ('CHLOROBENZENE');                              indx = 43
case ('CHLOROFORM');                                 indx = 44
case ('CONCRETE_PORTLAND');                          indx = 45
case ('CYCLOHEXANE');                                indx = 46
case ('12_DDIHLOROBENZENE');                         indx = 47
case ('DICHLORODIETHYL_ETHER');                      indx = 48
case ('12_DICHLOROETHANE');                          indx = 49
case ('DIETHYL_ETHER');                              indx = 50
case ('NN_DIMETHYL_FORMAMIDE');                      indx = 51
case ('DIMETHYL_SULFOXIDE');                         indx = 52
case ('ETHANE');                                     indx = 53
case ('ETHYL_ALCOHOL');                              indx = 54
case ('ETHYL_CELLULOSE');                            indx = 55
case ('ETHYLENE');                                   indx = 56
case ('EYE_LENS_ICRP');                              indx = 57
case ('FERRIC_OXIDE');                               indx = 58
case ('FERROBORIDE');                                indx = 59
case ('FERROUS_OXIDE');                              indx = 60
case ('FERROUS_SULFATE_DOSIMETER_SOLUTION');         indx = 61
case ('FREON_12');                                   indx = 62
case ('FREON_12B2');                                 indx = 63
case ('FREON_13');                                   indx = 64
case ('FREON_13B1');                                 indx = 65
case ('FREON_13I1');                                 indx = 66
case ('GADOLINIUM_OXYSULFIDE');                      indx = 67
case ('GALLIUM_ARSENIDE');                           indx = 68
case ('GEL_IN_PHOTOGRAPHIC_EMULSION');               indx = 69
case ('GLASS_PYREX');                                indx = 70
case ('GLASS_LEAD');                                 indx = 71
case ('GLASS_PLATE');                                indx = 72
case ('GLUCOSE');                                    indx = 73
case ('GLUTAMINE');                                  indx = 74
case ('GLYCEROL');                                   indx = 75
case ('GUANINE');                                    indx = 76
case ('GYPSUM_PLASTER_OF_PARIS');                    indx = 77
case ('N_HEPTANE');                                  indx = 78
case ('N_HEXANE');                                   indx = 79
case ('KAPTON_POLYIMIDE_FILM');                      indx = 80
case ('LANTHANUM_OXYBROMIDE');                       indx = 81
case ('LANTHANUM_OXYSULFIDE');                       indx = 82
case ('LEAD_OXIDE');                                 indx = 83
case ('LITHIUM_AMIDE');                              indx = 84
case ('LITHIUM_CARBONATE');                          indx = 85
case ('LITHIUM_FLUORIDE');                           indx = 86
case ('LITHIUM_HYDRIDE');                            indx = 87
case ('LITHIUM_IODIDE');                             indx = 88
case ('LITHIUM_OXIDE');                              indx = 89
case ('LITHIUM_TETRABORATE');                        indx = 90
case ('LUNG_ICRP');                                  indx = 91
case ('M3_WAX');                                     indx = 92
case ('MAGNESIUM_CARBONATE');                        indx = 93
case ('MAGNESIUM_FLUORIDE');                         indx = 94
case ('MAGNESIUM_OXIDE');                            indx = 95
case ('MAGNESIUM_TETRABORATE');                      indx = 96
case ('MERCURIC_IODIDE');                            indx = 97
case ('METHANE');                                    indx = 98
case ('METHANOL');                                   indx = 99
case ('MIX_D_WAX');                                  indx = 100
case ('MS20_TISSUE_SUBSTITUTE');                     indx = 101
case ('MUSCLE_SKELETAL');                            indx = 102
case ('MUSCLE_STRIATED');                            indx = 103
case ('MUSCLE_EQUIVALENT_LIQUID_WITH_SUCROSE');      indx = 104
case ('MUSCLE_EQUIVALENT_LIQUID_WITHOUT_SUCROSE');   indx = 105
case ('NAPHTHALENE');                                indx = 106
case ('NITROBENZENE');                               indx = 107
case ('NITROUS_OXIDE');                              indx = 108
case ('NYLON_DU_PONT_ELVAMIDE_8062');                indx = 109
case ('NYLON_TYPE_6_AND_TYPE_6_6');                  indx = 110
case ('NYLON_TYPE_6_10');                            indx = 111
case ('NYLON_TYPE_11_RILSAN');                       indx = 112
case ('OCTANE_LIQUID');                              indx = 113
case ('PARAFFIN_WAX');                               indx = 114
case ('N_PENTANE');                                  indx = 115
case ('PHOTOGRAPHIC_EMULSION');                      indx = 116
case ('PLASTIC_SCINTILLATOR_VINYLTOLUENE_BASED');    indx = 117
case ('PLUTONIUM_DIOXIDE');                          indx = 118
case ('POLYACRYLONITRILE');                          indx = 119
case ('POLYCARBONATE_MAKROLON_LEXAN');               indx = 120
case ('POLYCHLOROSTYRENE');                          indx = 121
case ('POLYETHYLENE');                               indx = 122
case ('POLYETHYLENE_TEREPHTHALATE_MYLAR');           indx = 123
case ('POLYMETHYL_METHACRALATE_LUCITE_PERSPEX');     indx = 124
case ('POLYOXYMETHYLENE');                           indx = 125
case ('POLYPROPYLENE');                              indx = 126
case ('POLYSTYRENE');                                indx = 127
case ('POLYTETRAFLUOROETHYLENE_TEFLON');             indx = 128
case ('POLYTRIFLUOROCHLOROETHYLENE');                indx = 129
case ('POLYVINYL_ACETATE');                          indx = 130
case ('POLYVINYL_ALCOHOL');                          indx = 131
case ('POLYVINYL_BUTYRAL');                          indx = 132
case ('POLYVINYL_CHLORIDE');                         indx = 133
case ('POLYVINYLIDENE_CHLORIDE_SARAN');              indx = 134
case ('POLYVINYLIDENE_FLUORIDE');                    indx = 135
case ('POLYVINYL_PYRROLIDONE');                      indx = 136
case ('POTASSIUM_IODIDE');                           indx = 137
case ('POTASSIUM_OXIDE');                            indx = 138
case ('PROPANE');                                    indx = 139
case ('PROPANE_LIQUID');                             indx = 140
case ('N_PROPYL_ALCOHOL');                           indx = 141
case ('PYRIDINE');                                   indx = 142
case ('RUBBER_BUTYL');                               indx = 143
case ('RUBBER_NATURAL');                             indx = 144
case ('RUBBER_NEOPRENE');                            indx = 145
case ('SILICON_DIOXIDE');                            indx = 146
case ('SILVER_BROMIDE');                             indx = 147
case ('SILVER_CHLORIDE');                            indx = 148
case ('SILVER_HALIDES_IN_PHOTOGRAPHIC_EMULSION');    indx = 149
case ('SILVER_IODIDE');                              indx = 150
case ('SKIN_ICRP');                                  indx = 151
case ('SODIUM_CARBONATE');                           indx = 152
case ('SODIUM_IODIDE');                              indx = 153
case ('SODIUM_MONOXIDE');                            indx = 154
case ('SODIUM_NITRATE');                             indx = 155
case ('STILBENE');                                   indx = 156
case ('SUCROSE');                                    indx = 157
case ('TERPHENYL');                                  indx = 158
case ('TESTES_ICRP');                                indx = 159
case ('TETRACHLOROETHYLENE');                        indx = 160
case ('THALLIUM_CHLORIDE');                          indx = 161
case ('TISSUE_SOFT_ICRP');                           indx = 162
case ('TISSUE_SOFT_ICRU_FOUR_COMPONENT');            indx = 163
case ('TISSUE_EQUIVALENT_GAS_METHANE_BASED');        indx = 164
case ('TISSUE_EQUIVALENT_GAS_PROPANE_BASED');        indx = 165
case ('TITANIUM_DIOXIDE');                           indx = 166
case ('TOLUENE');                                    indx = 167
case ('TRICHLOROETHYLENE');                          indx = 168
case ('TRIETHYL_PHOSPHATE');                         indx = 169
case ('TUNGSTEN_HEXAFLUORIDE');                      indx = 170
case ('URANIUM_DICARBIDE');                          indx = 171
case ('URANIUM_MONOCARBIDE');                        indx = 172
case ('URANIUM_OXIDE');                              indx = 173
case ('UREA');                                       indx = 174
case ('VALINE');                                     indx = 175
case ('VITON_FLUOROELASTOMER');                      indx = 176
case ('WATER_LIQUID', 'H2O');                        indx = 177
case ('WATER_VAPOR');                                indx = 178
case ('XYLENE');                                     indx = 179
end select

end function xraylib_nist_compound

end module


