!+
! Module particle_species_mod
!
! This module defines the differnet types of particles that Bmad knows about along
! with masses, etc.
!
! IMPORTANT: DO NOT USE HARD CODED ID NUMBERS IN YOUR CODE!! 
! For example, the association of positrons with ID = 1 is not assured.
! In general, use species_id(name) to get the species ID number or you can use positron$, 
! proton$, etc, for particles that have named paramters (see below for a list).
!
! Particles are divided into several categories:
!   1) Subatomic particles.
!   2) Atoms
!   3) "Known" molecules which are listed in the molecular_name() array.
!   4) "Unknown" molecules where the mass and charge are specified.
! Decoding of species integer:
! If |species| < 1000: Use elementary particle mapping (electron$ = -1, etc.)
! Else if |species| > 1000 mapping is:
! Write in hex: species = CCPPMMMM (Hex)
! Where:
!   CC   = Charge (2 Hex digits with range [-127, 127]). Set to 0 for subatomic particles.
!   PP   = Particle ID (2 Hex digits with range [0, 255]).
!           if PP = 0                --> Used for subatomic particles.
!           if 0 < PP < 199 (C8 Hex) --> Atom with PP = # Protons 
!           if PP = 199 (C8 Hex)     --> Anti atom.
!           if PP = 200 (C8 Hex)     --> Molecule of unknown type.
!           if PP > 200 (C8 Hex)     --> "Named" molecule. See molecular_name array below for a list. 
!                                        In this case PP = Species ID. EG: nh2$ = 201, etc.
!   MMMM (4 Hex digits):
!          For subatomic particles (where CC = PP = 0): Particle integer ID. 
!          For atoms: Number of nucleons. If zero then number of nucleons is unknown (EG: "C+")
!          For anti atoms: Split MMMM = Xb13b4YY where Xb13 is number of protons (7 bits) and b4YY = number of nucleons (9bits).
!          For Molecules: 100*Mass (That is, resolution is hundredths of an AMU). 0 = Use default (only valid for "Named" molecules).
!
! Example external input names:
!   NH3+            Molecule                           01 201 00000
!   CH3++ or CH3+2  Molecule                           02 204 00000
!   CH2@M37.5-      Molecule With specified mass      -01 201 03750
!   @M37.5+         Unknown Molecule with given mass   01 200 03750
!   C+              Atom:                              01 006 00000
!   #12C+           Atom: Carbon-12                    01 006 00012 
!-

module particle_species_mod

use output_mod
use sim_utils_struct
use sign_of_mod

implicit none

!----------------------
! Subatomic particles.
! Note: It is convenient for debugging to define an "antiparticle" with reversed sign even though it does not exit in practice. 
! Example: anti_deuteron$.

integer, parameter :: pion_0$            = +9
integer, parameter :: helion$            = +8
integer, parameter :: ref_particle$      = +7
integer, parameter :: neutron$           = +6
integer, parameter :: deuteron$          = +5
integer, parameter :: pion_plus$         = +4
integer, parameter :: antimuon$          = +3
integer, parameter :: proton$            = +2
integer, parameter :: positron$          = +1
integer, parameter :: photon$            =  0
integer, parameter :: electron$          = -1
integer, parameter :: antiproton$        = -2
integer, parameter :: muon$              = -3
integer, parameter :: pion_minus$        = -4
integer, parameter :: anti_deuteron$     = -5
integer, parameter :: anti_neutron$      = -6
integer, parameter :: anti_ref_particle$ = -7
integer, parameter :: anti_helion$       = -8

integer, parameter :: lb_subatomic = -8, ub_subatomic = 9
integer, parameter :: anti_atom$ = 199

character(20), parameter:: subatomic_species_name(lb_subatomic:ub_subatomic) = [character(20):: 'Anti_Helion', &
              'Anti_Ref_Particle', 'Anti_Neutron', 'Anti_Deuteron', 'Pion-', 'Muon', 'Antiproton', 'Electron', &
              'Photon', 'Positron', 'Proton', 'Antimuon', 'Pion+', 'Deuteron', 'Neutron', 'Ref_Particle', 'Helion', 'Pion0']

character(20), parameter:: openPMD_subatomic_species_name(lb_subatomic:ub_subatomic) = [character(20):: 'Garbage!', 'Garbage!', &
                      'anti-neutron', 'anti-deuteron', 'pion-', 'muon', 'anti-proton', 'electron', &
                      'photon', 'positron', 'proton', 'anti-muon', 'pion+', 'deuteron', 'neutron', 'Garbage!', 'Garbage!', 'pion0']

integer, parameter :: charge_of_subatomic(lb_subatomic:ub_subatomic) = [-2, 0, 0, -1, -1, -1, -1, -1, 0, 1, 1, 1, 1, 1, 0, 0, 2, 0]

real(rp), parameter :: mass_of_subatomic(lb_subatomic:ub_subatomic) = [m_helion, 0.0_rp, m_neutron, m_deuteron, &
                                m_pion_charged, m_muon, m_proton, m_electron, 0.0_rp, m_electron, m_proton, m_muon, &
                                m_pion_charged, m_deuteron, m_neutron, 0.0_rp, m_helion, m_pion_0]

real(rp), parameter :: anomalous_moment_of_subatomic(lb_subatomic:ub_subatomic) = [anomalous_mag_moment_He3, 0.0_rp, &
                        anomalous_mag_moment_neutron, anomalous_mag_moment_deuteron,0.0_rp,  &
                        anomalous_mag_moment_muon, anomalous_mag_moment_proton, anomalous_mag_moment_electron, 0.0_rp, &
                        anomalous_mag_moment_electron, anomalous_mag_moment_proton, anomalous_mag_moment_muon, &
                        0.0_rp, anomalous_mag_moment_deuteron, anomalous_mag_moment_neutron, 0.0_rp, anomalous_mag_moment_He3, 0.0_rp]

integer, parameter :: antiparticle_of_subatomic(lb_subatomic:ub_subatomic) = [helion$, &
                                                                   7, 6, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5, -6, -7, anti_helion$, 9]

real(rp), parameter :: spin_of_subatomic(lb_subatomic:ub_subatomic) = [0.0_rp, real_garbage$, 0.5_rp, 1.0_rp, 0.0_rp, &
                          0.5_rp, 0.5_rp, 0.5_rp, 0.0_rp, 0.5_rp, 0.5_rp, 0.5_rp, 0.0_rp, 1.0_rp, 0.5_rp, real_garbage$, 0.0_rp, 0.0_rp]

!----------------------
! Atoms

character(2), parameter :: atomic_name(118) = [character(2) :: &
                    'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', &
                    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', &
                    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', &
                    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', &
                    'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', &
                    'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', &
                    'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', &
                    'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', &
                    'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', &
                    'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', &
                    'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', & 
                    'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']

real(rp), parameter, private :: no_iso = -1 ! mass if there is no known isotope

type atom_struct
  integer :: z                    ! Number of protons
  character(2) :: name            ! Element name
  integer :: i_offset             ! isotope number offset 
  real(rp) :: mass(0:46) = no_iso ! isotope masses in units of the unified atomic mass unit.
                                  ! mass(0) is the standard atomic weight
                                  ! mass(n) is the mass of isotope n + i_offset
end type

! Radiation length X_0 from Y.S. Tsai, Rev. Mod. Phys. 46, 815 (1974). Units are gm/cm^2

real(rp) :: x0_rad_length(92) = [63.0470, 94.3221, 82.7559, 65.1899, 52.6868, 42.6983, 37.9879, 34.2381, 32.9303, 28.9367, &
                                       27.7362, 25.0387, 24.0111, 21.8234, 21.2053, 19.4953, 19.2783, 19.5489, 17.3167, 16.1442, &
                                       16.5455, 16.1745, 15.8425, 14.9444, 14.6398, 13.8389, 13.6174, 12.6820, 12.8616, 12.4269, &
                                       12.4734, 12.2459, 11.9401, 11.9082, 11.4230, 11.3722, 11.0272, 10.7623, 10.4101, 10.1949, &
                                       9.9225, 9.8029, 9.6881, 9.4825, 9.2654, 9.2025, 8.9701, 8.9945, 8.8491, 8.8170, &
                                       8.7244, 8.8267, 8.4803, 8.4819, 8.3052, 8.3073, 8.1381, 7.9557, 7.7579, 7.7051, &
                                       7.5193, 7.5727, 7.4377, 7.4830, 7.3563, 7.3199, 7.2332, 7.1448, 7.0318, 7.0214, &
                                       6.9237, 6.8907, 6.8177, 6.7630, 6.6897, 6.6763, 6.5936, 6.5433, 6.4608, 6.4368, &
                                       6.4176, 6.3688, 6.2899, 6.1907, 6.0651, 6.2833, 6.1868, 6.1477, 6.0560, 6.0726, &
                                       5.9319, 5.9990]

! Mean excitation energy (in eV) normalized by atomic Z. 
! This is used in the Bethe-Bloch formula for calculating dE/dx of charged particles through matter.
! Data from: https://physics.nist.gov/PhysRefData/XrayMassCoef/tab1.html
!
! This is suplemented with values from:
!     Hans Bichsel,
!    "Stopping Power of Fast Charged Particles in Heavy Elements",
!     NIST NISTIR 4550, April 1991.
! Values in this paper are from the b_e column in table 1 (pg 34) and are for 19 element in the range Z = 57 to 92.

real(rp), parameter :: mean_excitation_energy_over_z(92) = [ &
                19.20, 20.90, 13.33, 15.93, 15.20, 13.00, 11.71, 11.88, 12.78, 13.70, &
                13.55, 13.00, 12.77, 12.36, 11.53, 11.25, 10.24, 10.44, 10.00,  9.55, &
                10.29, 10.59, 10.65, 10.71, 10.88, 11.00, 11.00, 11.11, 11.10, 11.00, &
                10.77, 10.94, 10.52, 10.24,  9.80,  9.78,  9.81,  9.63,  9.72,  9.82, &
                10.17, 10.10,  9.95, 10.02,  9.98, 10.22, 10.00,  9.77,  9.96,  9.76, &
                 9.55,  9.33,  9.26,  8.93,  8.87,  8.77,  8.32,  8.76,  8.64,  9.10, &
                 9.18,  9.05,  9.21,  8.83,  9.45,  9.17,  9.55,  9.56,  9.77,  9.66, &
                 9.77,  9.32, 10.05, 10.53,  9.81,  9.82, 10.23, 10.08, 10.00, 10.00, &
                10.00,  9.50,  8.98,  9.88,  9.71,  9.23,  9.51,  9.39,  9.45,  8.51, &
                 9.65,  9.09]

! Table from NIST before Bichsel substitution
! 19.20, 20.90, 13.33, 15.93, 15.20, 13.00, 11.71, 11.88, 12.78, 13.70, &
! 13.55, 13.00, 12.77, 12.36, 11.53, 11.25, 10.24, 10.44, 10.00,  9.55, &
! 10.29, 10.59, 10.65, 10.71, 10.88, 11.00, 11.00, 11.11, 11.10, 11.00, &
! 10.77, 10.94, 10.52, 10.24,  9.80,  9.78,  9.81,  9.63,  9.72,  9.82, &
! 10.17, 10.10,  9.95, 10.02,  9.98, 10.22, 10.00,  9.77,  9.96,  9.76, &
!  9.55,  9.33,  9.26,  8.93,  8.87,  8.77,  8.79,  9.02,  9.07,  9.10, &
!  9.18,  9.26,  9.21,  9.23,  9.45,  9.52,  9.70,  9.68,  9.77,  9.77, &
!  9.77,  9.79,  9.84,  9.82,  9.81,  9.82,  9.83, 10.13, 10.00, 10.00, &
! 10.00, 10.04,  9.92,  9.88,  9.71,  9.23,  9.51,  9.39,  9.45,  9.41, &
!  9.65,  9.67]

! Isotopes from NIST data 2020/Jan.
! The first number in a row is the average mass for natually occuring isotope mixtures.

type (atom_struct), parameter, private :: atm1 = atom_struct(1, "H", 0, [1.007975_rp, 1.00782503223_rp, 2.01410177812_rp, 3.0160492779_rp, 4.02643_rp, &
      5.035311_rp, 6.04496_rp, 7.0527_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm2 = atom_struct(2, "He", 2, [4.002602_rp, 3.0160293201_rp, 4.00260325413_rp, 5.012057_rp, 6.018885891_rp, &
      7.0279907_rp, 8.033934390_rp, 9.043946_rp, 10.05279_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm3 = atom_struct(3, "Li", 2, [6.9675_rp, 3.0308_rp, 4.02719_rp, 5.012538_rp, 6.0151228874_rp, 7.0160034366_rp, &
      8.022486246_rp, 9.02679019_rp, 10.035483_rp, 11.04372358_rp, 12.052517_rp, 13.06263_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm4 = atom_struct(4, "Be", 4, [9.0121831_rp, 5.0399_rp, 6.0197264_rp, 7.016928717_rp, 8.005305102_rp, 9.012183065_rp, &
      10.013534695_rp, 11.02166108_rp, 12.0269221_rp, 13.036135_rp, 14.04289_rp, 15.05342_rp, 16.06167_rp, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm5 = atom_struct(5, "B", 5, [10.8135_rp, 6.0508_rp, 7.029712_rp, 8.0246073_rp, 9.01332965_rp, 10.01293695_rp, &
      11.00930536_rp, 12.0143527_rp, 13.0177802_rp, 14.025404_rp, 15.031088_rp, 16.039842_rp, 17.04699_rp, 18.05566_rp, 19.06310_rp, 20.07207_rp, &
      21.08129_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm6 = atom_struct(6, "C", 7, [12.0106_rp, 8.037643_rp, 9.0310372_rp, 10.01685331_rp, 11.0114336_rp, 12.0000000_rp, &
      13.00335483507_rp, 14.0032419884_rp, 15.01059926_rp, 16.0147013_rp, 17.022577_rp, 18.026751_rp, 19.03480_rp, 20.04032_rp, 21.04900_rp, 22.05753_rp, &
      23.0689_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm7 = atom_struct(7, "N", 9, [14.006855_rp, 10.04165_rp, 11.026091_rp, 12.0186132_rp, 13.00573861_rp, 14.00307400443_rp, &
      15.00010889888_rp, 16.0061019_rp, 17.008449_rp, 18.014078_rp, 19.017022_rp, 20.023366_rp, 21.02711_rp, 22.03439_rp, 23.04114_rp, 24.05039_rp, &
      25.06010_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm8 = atom_struct(8, "O", 11, [15.9994_rp, 12.034262_rp, 13.024815_rp, 14.00859636_rp, 15.00306562_rp, &
      15.99491461957_rp, 16.99913175650_rp, 17.99915961286_rp, 19.0035780_rp, 20.00407535_rp, 21.008655_rp, 22.009966_rp, 23.015696_rp, 24.01986_rp, &
      25.02936_rp, 26.03729_rp, 27.04772_rp, 28.05591_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm9 = atom_struct(9, "F", 13, [18.998403163_rp, 14.034315_rp, 15.018043_rp, 16.0114657_rp, 17.00209524_rp, &
      18.00093733_rp, 18.99840316273_rp, 19.999981252_rp, 20.9999489_rp, 22.002999_rp, 23.003557_rp, 24.008115_rp, 25.012199_rp, 26.020038_rp, 27.02644_rp, &
      28.03534_rp, 29.04254_rp, 30.05165_rp, 31.05971_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm10 = atom_struct(10, "Ne", 15, [20.1797_rp, 16.025750_rp, 17.01771396_rp, 18.00570870_rp, 19.00188091_rp, &
      19.9924401762_rp, 20.993846685_rp, 21.991385114_rp, 22.99446691_rp, 23.99361065_rp, 24.997789_rp, 26.000515_rp, 27.007553_rp, 28.01212_rp, 29.01975_rp, &
      30.02473_rp, 31.0331_rp, 32.03972_rp, 33.04938_rp, 34.05673_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm11 = atom_struct(11, "Na", 17, [22.98976928_rp, 18.02688_rp, 19.013880_rp, 20.0073544_rp, 20.99765469_rp, &
      21.99443741_rp, 22.9897692820_rp, 23.990962950_rp, 24.9899540_rp, 25.9926346_rp, 26.9940765_rp, 27.998939_rp, 29.0028771_rp, 30.0090979_rp, &
      31.013163_rp, 32.02019_rp, 33.02573_rp, 34.03359_rp, 35.04062_rp, 36.04929_rp, 37.05705_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso])
type (atom_struct), parameter, private :: atm12 = atom_struct(12, "Mg", 18, [24.3055_rp, 19.034169_rp, 20.018850_rp, 21.011716_rp, 21.99957065_rp, &
      22.99412421_rp, 23.985041697_rp, 24.985836976_rp, 25.982592968_rp, 26.984340624_rp, 27.9838767_rp, 28.988617_rp, 29.9904629_rp, 30.9966480_rp, &
      31.9991102_rp, 33.0053271_rp, 34.008935_rp, 35.01679_rp, 36.02188_rp, 37.03037_rp, 38.03658_rp, 39.04538_rp, 40.05218_rp, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm13 = atom_struct(13, "Al", 20, [26.9815385_rp, 21.02897_rp, 22.01954_rp, 23.00724435_rp, 23.9999489_rp, &
      24.99042810_rp, 25.986891904_rp, 26.98153853_rp, 27.98191021_rp, 28.9804565_rp, 29.982960_rp, 30.983945_rp, 31.988085_rp, 32.990909_rp, 33.996705_rp, &
      34.999764_rp, 36.00639_rp, 37.01053_rp, 38.01740_rp, 39.02254_rp, 40.03003_rp, 41.03638_rp, 42.04384_rp, 43.05147_rp, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso])
type (atom_struct), parameter, private :: atm14 = atom_struct(14, "Si", 21, [28.085_rp, 22.03579_rp, 23.02544_rp, 24.011535_rp, 25.004109_rp, 25.99233384_rp, &
      26.98670481_rp, 27.97692653465_rp, 28.97649466490_rp, 29.973770136_rp, 30.975363194_rp, 31.97415154_rp, 32.97797696_rp, 33.978576_rp, 34.984583_rp, &
      35.986695_rp, 36.992921_rp, 37.995523_rp, 39.002491_rp, 40.00583_rp, 41.01301_rp, 42.01778_rp, 43.02480_rp, 44.03061_rp, 45.03995_rp, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm15 = atom_struct(15, "P", 23, [30.973761998_rp, 24.03577_rp, 25.02119_rp, 26.01178_rp, 26.999224_rp, &
      27.9923266_rp, 28.98180079_rp, 29.97831375_rp, 30.97376199842_rp, 31.973907643_rp, 32.9717257_rp, 33.97364589_rp, 34.9733141_rp, 35.978260_rp, &
      36.979607_rp, 37.984252_rp, 38.986227_rp, 39.99133_rp, 40.994654_rp, 42.00108_rp, 43.00502_rp, 44.01121_rp, 45.01645_rp, 46.02446_rp, 47.03139_rp, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm16 = atom_struct(16, "S", 25, [32.0675_rp, 26.02907_rp, 27.01828_rp, 28.00437_rp, 28.996611_rp, 29.98490703_rp, &
      30.97955701_rp, 31.9720711744_rp, 32.9714589098_rp, 33.967867004_rp, 34.969032310_rp, 35.96708071_rp, 36.97112551_rp, 37.9711633_rp, 38.975134_rp, &
      39.9754826_rp, 40.9795935_rp, 41.9810651_rp, 42.9869076_rp, 43.9901188_rp, 44.99572_rp, 46.00004_rp, 47.00795_rp, 48.01370_rp, 49.02276_rp, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm17 = atom_struct(17, "Cl", 27, [35.4515_rp, 28.02954_rp, 29.01478_rp, 30.00477_rp, 30.992414_rp, 31.98568464_rp, &
      32.97745199_rp, 33.973762485_rp, 34.968852682_rp, 35.968306809_rp, 36.965902602_rp, 37.96801044_rp, 38.9680082_rp, 39.970415_rp, 40.970685_rp, &
      41.97325_rp, 42.97389_rp, 43.97787_rp, 44.98029_rp, 45.98517_rp, 46.98916_rp, 47.99564_rp, 49.00123_rp, 50.00905_rp, 51.01554_rp, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso])
type (atom_struct), parameter, private :: atm18 = atom_struct(18, "Ar", 29, [39.948_rp, 30.02307_rp, 31.01212_rp, 31.9976378_rp, 32.98992555_rp, 33.980270090_rp, &
      34.97525759_rp, 35.967545105_rp, 36.96677633_rp, 37.96273211_rp, 38.9643130_rp, 39.9623831237_rp, 40.96450057_rp, 41.9630457_rp, 42.9656361_rp, &
      43.9649238_rp, 44.96803973_rp, 45.968083_rp, 46.972935_rp, 47.97591_rp, 48.98190_rp, 49.98613_rp, 50.99370_rp, 51.99896_rp, 53.00729_rp, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm19 = atom_struct(19, "K", 31, [39.0983_rp, 32.02265_rp, 33.00756_rp, 33.99869_rp, 34.98800541_rp, 35.98130201_rp, &
      36.97337589_rp, 37.96908112_rp, 38.9637064864_rp, 39.963998166_rp, 40.9618252579_rp, 41.96240231_rp, 42.96073470_rp, 43.96158699_rp, 44.96069149_rp, &
      45.96198159_rp, 46.9616616_rp, 47.96534119_rp, 48.96821075_rp, 49.9723800_rp, 50.975828_rp, 51.98224_rp, 52.98746_rp, 53.99463_rp, 55.00076_rp, &
      56.00851_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm20 = atom_struct(20, "Ca", 33, [40.078_rp, 34.01487_rp, 35.00514_rp, 35.993074_rp, 36.98589785_rp, 37.97631922_rp, &
      38.97071081_rp, 39.962590863_rp, 40.96227792_rp, 41.95861783_rp, 42.95876644_rp, 43.95548156_rp, 44.95618635_rp, 45.9536890_rp, 46.9545424_rp, &
      47.95252276_rp, 48.95566274_rp, 49.9574992_rp, 50.960989_rp, 51.963217_rp, 52.96945_rp, 53.97340_rp, 54.98030_rp, 55.98508_rp, 56.99262_rp, &
      57.99794_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm21 = atom_struct(21, "Sc", 35, [44.955908_rp, 36.01648_rp, 37.00374_rp, 37.99512_rp, 38.984785_rp, 39.9779673_rp, &
      40.969251105_rp, 41.96551653_rp, 42.9611505_rp, 43.9594029_rp, 44.95590828_rp, 45.95516826_rp, 46.9524037_rp, 47.9522236_rp, 48.9500146_rp, &
      49.952176_rp, 50.953592_rp, 51.95688_rp, 52.95909_rp, 53.96393_rp, 54.96782_rp, 55.97345_rp, 56.97777_rp, 57.98403_rp, 58.98894_rp, 59.99565_rp, &
      61.00100_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm22 = atom_struct(22, "Ti", 37, [47.867_rp, 38.01145_rp, 39.00236_rp, 39.99050_rp, 40.983148_rp, 41.97304903_rp, &
      42.9685225_rp, 43.95968995_rp, 44.95812198_rp, 45.95262772_rp, 46.95175879_rp, 47.94794198_rp, 48.94786568_rp, 49.94478689_rp, 50.94661065_rp, &
      51.9468930_rp, 52.94973_rp, 53.95105_rp, 54.95527_rp, 55.95791_rp, 56.96364_rp, 57.96660_rp, 58.97247_rp, 59.97603_rp, 60.98245_rp, 61.98651_rp, &
      62.99375_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm23 = atom_struct(23, "V", 39, [50.9415_rp, 40.01276_rp, 41.00021_rp, 41.99182_rp, 42.980766_rp, 43.97411_rp, &
      44.9657748_rp, 45.96019878_rp, 46.95490491_rp, 47.9522522_rp, 48.94851180_rp, 49.94715601_rp, 50.94395704_rp, 51.94477301_rp, 52.9443367_rp, &
      53.946439_rp, 54.94724_rp, 55.95048_rp, 56.95252_rp, 57.95672_rp, 58.95939_rp, 59.96431_rp, 60.96725_rp, 61.97265_rp, 62.97639_rp, 63.98264_rp, &
      64.98750_rp, 65.99398_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm24 = atom_struct(24, "Cr", 41, [51.9961_rp, 42.00670_rp, 42.99753_rp, 43.98536_rp, 44.979050_rp, 45.968359_rp, &
      46.9628974_rp, 47.9540291_rp, 48.9513333_rp, 49.94604183_rp, 50.94476502_rp, 51.94050623_rp, 52.94064815_rp, 53.93887916_rp, 54.94083843_rp, &
      55.9406531_rp, 56.9436130_rp, 57.94435_rp, 58.94859_rp, 59.95008_rp, 60.95442_rp, 61.95610_rp, 62.96165_rp, 63.96408_rp, 64.96996_rp, 65.97366_rp, &
      66.98016_rp, 67.98403_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm25 = atom_struct(25, "Mn", 43, [54.938044_rp, 44.00715_rp, 44.99449_rp, 45.98609_rp, 46.975775_rp, 47.96852_rp, &
      48.959595_rp, 49.95423778_rp, 50.94820847_rp, 51.9455639_rp, 52.94128889_rp, 53.9403576_rp, 54.93804391_rp, 55.93890369_rp, 56.9382861_rp, 57.9400666_rp, &
      58.9403911_rp, 59.9431366_rp, 60.9444525_rp, 61.94795_rp, 62.9496647_rp, 63.9538494_rp, 64.9560198_rp, 65.960547_rp, 66.96424_rp, 67.96962_rp, &
      68.97366_rp, 69.97937_rp, 70.98368_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm26 = atom_struct(26, "Fe", 44, [55.845_rp, 45.01442_rp, 46.00063_rp, 46.99185_rp, 47.98023_rp, 48.973429_rp, &
      49.962975_rp, 50.9568410_rp, 51.9481131_rp, 52.9453064_rp, 53.93960899_rp, 54.93829199_rp, 55.93493633_rp, 56.93539284_rp, 57.93327443_rp, 58.93487434_rp, &
      59.9340711_rp, 60.9367462_rp, 61.9367918_rp, 62.9402727_rp, 63.9409878_rp, 64.9450115_rp, 65.9462500_rp, 66.95054_rp, 67.95295_rp, 68.95807_rp, &
      69.96102_rp, 70.96672_rp, 71.96983_rp, 72.97572_rp, 73.97935_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm27 = atom_struct(27, "Co", 46, [58.933194_rp, 47.01057_rp, 48.00093_rp, 48.98891_rp, 49.98091_rp, 50.970647_rp, &
      51.96351_rp, 52.9542041_rp, 53.94845987_rp, 54.94199720_rp, 55.93983880_rp, 56.93629057_rp, 57.9357521_rp, 58.93319429_rp, 59.93381630_rp, 60.93247662_rp, &
      61.934059_rp, 62.933600_rp, 63.935811_rp, 64.9364621_rp, 65.939443_rp, 66.9406096_rp, 67.94426_rp, 68.94614_rp, 69.94963_rp, 70.95237_rp, 71.95729_rp, &
      72.96039_rp, 73.96515_rp, 74.96876_rp, 75.97413_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm28 = atom_struct(28, "Ni", 47, [58.6934_rp, 48.01769_rp, 49.00770_rp, 49.99474_rp, 50.98611_rp, 51.97480_rp, &
      52.968190_rp, 53.957892_rp, 54.95133063_rp, 55.94212855_rp, 56.93979218_rp, 57.93534241_rp, 58.93434620_rp, 59.93078588_rp, 60.93105557_rp, &
      61.92834537_rp, 62.92966963_rp, 63.92796682_rp, 64.93008517_rp, 65.9291393_rp, 66.9315694_rp, 67.9318688_rp, 68.9356103_rp, 69.9364313_rp, 70.9405190_rp, &
      71.9417859_rp, 72.9462067_rp, 73.94798_rp, 74.95250_rp, 75.95533_rp, 76.96055_rp, 77.96336_rp, 78.97025_rp, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm29 = atom_struct(29, "Cu", 51, [63.546_rp, 51.99671_rp, 52.98459_rp, 53.97666_rp, 54.96604_rp, 55.95895_rp, &
      56.94921250_rp, 57.94453305_rp, 58.93949748_rp, 59.9373645_rp, 60.9334576_rp, 61.93259541_rp, 62.92959772_rp, 63.92976434_rp, 64.92778970_rp, &
      65.92886903_rp, 66.9277303_rp, 67.9296109_rp, 68.9294293_rp, 69.9323921_rp, 70.9326768_rp, 71.9358203_rp, 72.9366744_rp, 73.9398749_rp, 74.9415226_rp, &
      75.9452750_rp, 76.94792_rp, 77.95223_rp, 78.95502_rp, 79.96089_rp, 80.96587_rp, 81.97244_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm30 = atom_struct(30, "Zn", 53, [65.38_rp, 53.99204_rp, 54.98398_rp, 55.97254_rp, 56.96506_rp, 57.954591_rp, &
      58.94931266_rp, 59.94184210_rp, 60.939507_rp, 61.93433397_rp, 62.9332115_rp, 63.92914201_rp, 64.92924077_rp, 65.92603381_rp, 66.92712775_rp, &
      67.92484455_rp, 68.9265507_rp, 69.9253192_rp, 70.9277196_rp, 71.9268428_rp, 72.9295826_rp, 73.9294073_rp, 74.9328402_rp, 75.9331150_rp, 76.9368872_rp, &
      77.9382892_rp, 78.9426381_rp, 79.9445529_rp, 80.9504026_rp, 81.95426_rp, 82.96056_rp, 83.96521_rp, 84.97226_rp, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm31 = atom_struct(31, "Ga", 55, [69.723_rp, 55.99536_rp, 56.98320_rp, 57.97478_rp, 58.96353_rp, 59.95729_rp, &
      60.949399_rp, 61.94419025_rp, 62.9392942_rp, 63.9368404_rp, 64.93273459_rp, 65.9315894_rp, 66.9282025_rp, 67.9279805_rp, 68.9255735_rp, 69.9260219_rp, &
      70.92470258_rp, 71.92636747_rp, 72.9251747_rp, 73.9269457_rp, 74.9265002_rp, 75.9288276_rp, 76.9291543_rp, 77.9316088_rp, 78.9328523_rp, 79.9364208_rp, &
      80.9381338_rp, 81.9431765_rp, 82.9471203_rp, 83.95246_rp, 84.95699_rp, 85.96301_rp, 86.96824_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm32 = atom_struct(32, "Ge", 57, [72.630_rp, 57.99172_rp, 58.98249_rp, 59.97036_rp, 60.96379_rp, 61.95502_rp, &
      62.949628_rp, 63.9416899_rp, 64.9393681_rp, 65.9338621_rp, 66.9327339_rp, 67.9280953_rp, 68.9279645_rp, 69.92424875_rp, 70.92495233_rp, 71.922075826_rp, &
      72.923458956_rp, 73.921177761_rp, 74.922858370_rp, 75.921402726_rp, 76.923549843_rp, 77.9228529_rp, 78.925360_rp, 79.9253508_rp, 80.9288329_rp, &
      81.9297740_rp, 82.9345391_rp, 83.9375751_rp, 84.9429697_rp, 85.94658_rp, 86.95268_rp, 87.95691_rp, 88.96379_rp, 89.96863_rp, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm33 = atom_struct(33, "As", 59, [74.921595_rp, 59.99388_rp, 60.98112_rp, 61.97361_rp, 62.96390_rp, 63.95743_rp, &
      64.949611_rp, 65.9441488_rp, 66.93925111_rp, 67.9367741_rp, 68.932246_rp, 69.930926_rp, 70.9271138_rp, 71.9267523_rp, 72.9238291_rp, 73.9239286_rp, &
      74.92159457_rp, 75.92239202_rp, 76.9206476_rp, 77.921828_rp, 78.9209484_rp, 79.9224746_rp, 80.9221323_rp, 81.9247412_rp, 82.9252069_rp, 83.9293033_rp, &
      84.9321637_rp, 85.9367015_rp, 86.9402917_rp, 87.94555_rp, 88.94976_rp, 89.95563_rp, 90.96039_rp, 91.96674_rp, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm34 = atom_struct(34, "Se", 63, [78.971_rp, 63.97109_rp, 64.96440_rp, 65.95559_rp, 66.949994_rp, 67.94182524_rp, &
      68.9394148_rp, 69.9335155_rp, 70.9322094_rp, 71.9271405_rp, 72.9267549_rp, 73.922475934_rp, 74.922522870_rp, 75.919213704_rp, 76.919914154_rp, &
      77.91730928_rp, 78.91849929_rp, 79.9165218_rp, 80.9179930_rp, 81.9166995_rp, 82.9191186_rp, 83.9184668_rp, 84.9222608_rp, 85.9243117_rp, 86.9286886_rp, &
      87.9314175_rp, 88.9366691_rp, 89.94010_rp, 90.94596_rp, 91.94984_rp, 92.95629_rp, 93.96049_rp, 94.96730_rp, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm35 = atom_struct(35, "Br", 66, [79.904_rp, 66.96465_rp, 67.95873_rp, 68.950497_rp, 69.944792_rp, 70.9393422_rp, &
      71.9365886_rp, 72.9316715_rp, 73.9299102_rp, 74.9258105_rp, 75.924542_rp, 76.9213792_rp, 77.9211459_rp, 78.9183376_rp, 79.9185298_rp, 80.9162897_rp, &
      81.9168032_rp, 82.9151756_rp, 83.916496_rp, 84.9156458_rp, 85.9188054_rp, 86.9206740_rp, 87.9240833_rp, 88.9267046_rp, 89.9312928_rp, 90.9343986_rp, &
      91.9396316_rp, 92.94313_rp, 93.94890_rp, 94.95301_rp, 95.95903_rp, 96.96344_rp, 97.96946_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm36 = atom_struct(36, "Kr", 68, [83.798_rp, 68.96518_rp, 69.95604_rp, 70.95027_rp, 71.9420924_rp, 72.9392892_rp, &
      73.9330840_rp, 74.9309457_rp, 75.9259103_rp, 76.9246700_rp, 77.92036494_rp, 78.9200829_rp, 79.91637808_rp, 80.9165912_rp, 81.91348273_rp, 82.91412716_rp, &
      83.9114977282_rp, 84.9125273_rp, 85.9106106269_rp, 86.91335476_rp, 87.9144479_rp, 88.9178355_rp, 89.9195279_rp, 90.9238063_rp, 91.9261731_rp, &
      92.9311472_rp, 93.934140_rp, 94.939711_rp, 95.943017_rp, 96.94909_rp, 97.95243_rp, 98.95839_rp, 99.96237_rp, 100.96873_rp, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm37 = atom_struct(37, "Rb", 70, [85.4678_rp, 70.96532_rp, 71.95908_rp, 72.95053_rp, 73.9442659_rp, 74.9385732_rp, &
      75.9350730_rp, 76.9304016_rp, 77.9281419_rp, 78.9239899_rp, 79.9225164_rp, 80.9189939_rp, 81.9182090_rp, 82.9151142_rp, 83.9143752_rp, 84.9117897379_rp, &
      85.91116743_rp, 86.9091805310_rp, 87.91131559_rp, 88.9122783_rp, 89.9147985_rp, 90.9165372_rp, 91.9197284_rp, 92.9220393_rp, 93.9263948_rp, &
      94.929260_rp, 95.9341334_rp, 96.9371771_rp, 97.9416869_rp, 98.94503_rp, 99.95003_rp, 100.95404_rp, 101.95952_rp, 102.96392_rp, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm38 = atom_struct(38, "Sr", 72, [87.62_rp, 72.96570_rp, 73.95617_rp, 74.94995_rp, 75.941763_rp, 76.9379455_rp, &
      77.9321800_rp, 78.9297077_rp, 79.9245175_rp, 80.9232114_rp, 81.9183999_rp, 82.9175544_rp, 83.9134191_rp, 84.9129320_rp, 85.9092606_rp, 86.9088775_rp, &
      87.9056125_rp, 88.9074511_rp, 89.9077300_rp, 90.9101954_rp, 91.9110382_rp, 92.9140242_rp, 93.9153556_rp, 94.9193529_rp, 95.9217066_rp, 96.9263740_rp, &
      97.9286888_rp, 98.9328907_rp, 99.935770_rp, 100.940352_rp, 101.943791_rp, 102.94909_rp, 103.95265_rp, 104.95855_rp, 105.96265_rp, 106.96897_rp, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm39 = atom_struct(39, "Y", 75, [88.90584_rp, 75.95856_rp, 76.949781_rp, 77.94361_rp, 78.93735_rp, 79.9343561_rp, &
      80.9294556_rp, 81.9269314_rp, 82.922485_rp, 83.9206721_rp, 84.916433_rp, 85.914886_rp, 86.9108761_rp, 87.9095016_rp, 88.9058403_rp, 89.9071439_rp, &
      90.9072974_rp, 91.9089451_rp, 92.909578_rp, 93.9115906_rp, 94.9128161_rp, 95.9158968_rp, 96.9182741_rp, 97.9223821_rp, 98.9241480_rp, 99.927715_rp, &
      100.9301477_rp, 101.9343277_rp, 102.937243_rp, 103.94196_rp, 104.94544_rp, 105.95056_rp, 106.95452_rp, 107.95996_rp, 108.96436_rp, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm40 = atom_struct(40, "Zr", 77, [91.224_rp, 77.95566_rp, 78.94948_rp, 79.9404_rp, 80.93731_rp, 81.93135_rp, &
      82.9292421_rp, 83.9233269_rp, 84.9214444_rp, 85.9162972_rp, 86.9148180_rp, 87.9102213_rp, 88.9088814_rp, 89.9046977_rp, 90.9056396_rp, 91.9050347_rp, &
      92.9064699_rp, 93.9063108_rp, 94.9080385_rp, 95.9082714_rp, 96.9109512_rp, 97.9127289_rp, 98.916667_rp, 99.9180006_rp, 100.9214480_rp, 101.9231409_rp, &
      102.927191_rp, 103.929436_rp, 104.934008_rp, 105.93676_rp, 106.94174_rp, 107.94487_rp, 108.95041_rp, 109.95396_rp, 110.95968_rp, 111.96370_rp, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm41 = atom_struct(41, "Nb", 80, [92.90637_rp, 80.94960_rp, 81.94396_rp, 82.93729_rp, 83.93449_rp, 84.9288458_rp, &
      85.9257828_rp, 86.9206937_rp, 87.918222_rp, 88.913445_rp, 89.9112584_rp, 90.9069897_rp, 91.9071881_rp, 92.9063730_rp, 93.9072788_rp, 94.90683240_rp, &
      95.9080973_rp, 96.9080959_rp, 97.9103265_rp, 98.911613_rp, 99.9143276_rp, 100.9153103_rp, 101.9180772_rp, 102.9194572_rp, 103.9228925_rp, 104.9249465_rp, &
      105.9289317_rp, 106.9315937_rp, 107.9360748_rp, 108.93922_rp, 109.94403_rp, 110.94753_rp, 111.95247_rp, 112.95651_rp, 113.96201_rp, 114.96634_rp, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm42 = atom_struct(42, "Mo", 82, [95.95_rp, 82.94988_rp, 83.94149_rp, 84.938261_rp, 85.9311748_rp, 86.9281962_rp, &
      87.9219678_rp, 88.9194682_rp, 89.9139309_rp, 90.9117453_rp, 91.90680796_rp, 92.90680958_rp, 93.90508490_rp, 94.90583877_rp, 95.90467612_rp, &
      96.90601812_rp, 97.90540482_rp, 98.90770851_rp, 99.9074718_rp, 100.9103414_rp, 101.9102834_rp, 102.913079_rp, 103.9137344_rp, 104.916969_rp, &
      105.918259_rp, 106.922106_rp, 107.924033_rp, 108.928424_rp, 109.930704_rp, 110.935654_rp, 111.93831_rp, 112.94335_rp, 113.94653_rp, 114.95196_rp, &
      115.95545_rp, 116.96117_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm43 = atom_struct(43, "Tc", 84, [98.0_rp, 84.95058_rp, 85.94493_rp, 86.9380672_rp, 87.93378_rp, 88.9276487_rp, &
      89.9240739_rp, 90.9184254_rp, 91.9152698_rp, 92.9102460_rp, 93.9096536_rp, 94.9076536_rp, 95.9078680_rp, 96.9063667_rp, 97.9072124_rp, 98.9062508_rp, &
      99.9076539_rp, 100.907309_rp, 101.9092097_rp, 102.909176_rp, 103.911425_rp, 104.911655_rp, 105.914358_rp, 106.9154606_rp, 107.9184957_rp, 108.920256_rp, &
      109.923744_rp, 110.925901_rp, 111.9299458_rp, 112.9325690_rp, 113.93691_rp, 114.93998_rp, 115.94476_rp, 116.94806_rp, 117.95299_rp, 118.95666_rp, &
      119.96187_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm44 = atom_struct(44, "Ru", 86, [101.07_rp, 86.95069_rp, 87.94160_rp, 88.93762_rp, 89.9303444_rp, 90.9267419_rp, &
      91.9202344_rp, 92.9171044_rp, 93.9113429_rp, 94.910406_rp, 95.90759025_rp, 96.9075471_rp, 97.9052868_rp, 98.9059341_rp, 99.9042143_rp, 100.9055769_rp, &
      101.9043441_rp, 102.9063186_rp, 103.9054275_rp, 104.9077476_rp, 105.9073291_rp, 106.9099720_rp, 107.9101880_rp, 108.9133260_rp, 109.9140407_rp, &
      110.917570_rp, 111.918809_rp, 112.922844_rp, 113.9246136_rp, 114.928820_rp, 115.9312192_rp, 116.93610_rp, 117.93853_rp, 118.94357_rp, 119.94631_rp, &
      120.95164_rp, 121.95447_rp, 122.95989_rp, 123.96305_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm45 = atom_struct(45, "Rh", 88, [102.90550_rp, 88.95058_rp, 89.94422_rp, 90.93688_rp, 91.9323677_rp, 92.9259128_rp, &
      93.9217305_rp, 94.9158979_rp, 95.914453_rp, 96.911329_rp, 97.910708_rp, 98.9081282_rp, 99.908117_rp, 100.9061606_rp, 101.9068374_rp, 102.9054980_rp, &
      103.9066492_rp, 104.9056885_rp, 105.9072868_rp, 106.906748_rp, 107.908714_rp, 108.9087488_rp, 109.911079_rp, 110.9116423_rp, 111.914403_rp, &
      112.9154393_rp, 113.918718_rp, 114.9203116_rp, 115.924059_rp, 116.9260354_rp, 117.930340_rp, 118.932557_rp, 119.93686_rp, 120.93942_rp, 121.94399_rp, &
      122.94685_rp, 123.95151_rp, 124.95469_rp, 125.95946_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm46 = atom_struct(46, "Pd", 90, [106.42_rp, 90.95032_rp, 91.94088_rp, 92.93651_rp, 93.9290376_rp, 94.9248898_rp, &
      95.9182151_rp, 96.9164720_rp, 97.9126983_rp, 98.9117748_rp, 99.908505_rp, 100.9082864_rp, 101.9056022_rp, 102.9060809_rp, 103.9040305_rp, 104.9050796_rp, &
      105.9034804_rp, 106.9051282_rp, 107.9038916_rp, 108.9059504_rp, 109.90517220_rp, 110.90768968_rp, 111.9073297_rp, 112.9102610_rp, 113.9103686_rp, &
      114.913659_rp, 115.9142970_rp, 116.9179547_rp, 117.9190667_rp, 118.9233402_rp, 119.9245511_rp, 120.9289503_rp, 121.930632_rp, 122.93514_rp, &
      123.93714_rp, 124.94179_rp, 125.94416_rp, 126.94907_rp, 127.95183_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm47 = atom_struct(47, "Ag", 92, [107.8682_rp, 92.95033_rp, 93.94373_rp, 94.93602_rp, 95.930744_rp, 96.92397_rp, &
      97.921560_rp, 98.9176458_rp, 99.9161154_rp, 100.9126840_rp, 101.9117047_rp, 102.9089631_rp, 103.9086239_rp, 104.9065256_rp, 105.9066636_rp, &
      106.9050916_rp, 107.9059503_rp, 108.9047553_rp, 109.9061102_rp, 110.9052959_rp, 111.9070486_rp, 112.906573_rp, 113.9088230_rp, 114.908767_rp, &
      115.9113868_rp, 116.911774_rp, 117.9145955_rp, 118.915570_rp, 119.9187848_rp, 120.920125_rp, 121.923664_rp, 122.925337_rp, 123.92893_rp, 124.93105_rp, &
      125.93475_rp, 126.93711_rp, 127.94106_rp, 128.94395_rp, 129.95070_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm48 = atom_struct(48, "Cd", 94, [112.414_rp, 94.94994_rp, 95.94034_rp, 96.93510_rp, 97.927389_rp, 98.9249258_rp, &
      99.9203488_rp, 100.9185862_rp, 101.9144820_rp, 102.9134165_rp, 103.9098564_rp, 104.9094639_rp, 105.9064599_rp, 106.9066121_rp, 107.9041834_rp, &
      108.9049867_rp, 109.90300661_rp, 110.90418287_rp, 111.90276287_rp, 112.90440813_rp, 113.90336509_rp, 114.90543751_rp, 115.90476315_rp, 116.9072260_rp, &
      117.906922_rp, 118.909847_rp, 119.9098681_rp, 120.9129637_rp, 121.9134591_rp, 122.9168925_rp, 123.9176574_rp, 124.9212576_rp, 125.9224291_rp, &
      126.926472_rp, 127.9278129_rp, 128.93182_rp, 129.93394_rp, 130.94060_rp, 131.94604_rp, 132.95285_rp, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso])
type (atom_struct), parameter, private :: atm49 = atom_struct(49, "In", 96, [114.818_rp, 96.94934_rp, 97.94214_rp, 98.93411_rp, 99.93096_rp, 100.92634_rp, &
      101.9241071_rp, 102.9198819_rp, 103.9182145_rp, 104.914502_rp, 105.913464_rp, 106.910290_rp, 107.9096935_rp, 108.9071514_rp, 109.907170_rp, &
      110.9051085_rp, 111.9055377_rp, 112.90406184_rp, 113.90491791_rp, 114.903878776_rp, 115.90525999_rp, 116.9045157_rp, 117.9063566_rp, 118.9058507_rp, &
      119.907967_rp, 120.907851_rp, 121.910281_rp, 122.910434_rp, 123.913182_rp, 124.913605_rp, 125.916507_rp, 126.917446_rp, 127.92040_rp, 128.9218053_rp, &
      129.924977_rp, 130.9269715_rp, 131.933001_rp, 132.93831_rp, 133.94454_rp, 134.95005_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm50 = atom_struct(50, "Sn", 98, [118.710_rp, 98.94853_rp, 99.93850_rp, 100.93526_rp, 101.93029_rp, 102.928105_rp, &
      103.9231052_rp, 104.9212684_rp, 105.9169574_rp, 106.9157137_rp, 107.9118943_rp, 108.9112921_rp, 109.907845_rp, 110.9077401_rp, 111.90482387_rp, &
      112.9051757_rp, 113.9027827_rp, 114.903344699_rp, 115.90174280_rp, 116.90295398_rp, 117.90160657_rp, 118.90331117_rp, 119.90220163_rp, 120.9042426_rp, &
      121.9034438_rp, 122.9057252_rp, 123.9052766_rp, 124.9077864_rp, 125.907659_rp, 126.910390_rp, 127.910507_rp, 128.913465_rp, 129.9139738_rp, &
      130.9170450_rp, 131.9178267_rp, 132.9239134_rp, 133.9286821_rp, 134.9349086_rp, 135.93999_rp, 136.94655_rp, 137.95184_rp, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm51 = atom_struct(51, "Sb", 102, [121.760_rp, 102.93969_rp, 103.93648_rp, 104.931276_rp, 105.9286380_rp, &
      106.9241506_rp, 107.9222267_rp, 108.9181411_rp, 109.9168543_rp, 110.9132182_rp, 111.912400_rp, 112.909375_rp, 113.909290_rp, 114.906598_rp, &
      115.9067931_rp, 116.9048415_rp, 117.9055321_rp, 118.9039455_rp, 119.9050794_rp, 120.9038120_rp, 121.9051699_rp, 122.9042132_rp, 123.9059350_rp, &
      124.9052530_rp, 125.907253_rp, 126.9069243_rp, 127.909146_rp, 128.909147_rp, 129.911662_rp, 130.9119888_rp, 131.9145077_rp, 132.9152732_rp, &
      133.9205357_rp, 134.9251851_rp, 135.9307459_rp, 136.93555_rp, 137.94145_rp, 138.94655_rp, 139.95283_rp, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm52 = atom_struct(52, "Te", 104, [127.60_rp, 104.94330_rp, 105.93750_rp, 106.935012_rp, 107.9293805_rp, &
      108.9273045_rp, 109.9224581_rp, 110.9210006_rp, 111.9167279_rp, 112.915891_rp, 113.912089_rp, 114.911902_rp, 115.908460_rp, 116.908646_rp, 117.905854_rp, &
      118.9064071_rp, 119.9040593_rp, 120.904944_rp, 121.9030435_rp, 122.9042698_rp, 123.9028171_rp, 124.9044299_rp, 125.9033109_rp, 126.9052257_rp, &
      127.90446128_rp, 128.90659646_rp, 129.906222748_rp, 130.908522213_rp, 131.9085467_rp, 132.9109688_rp, 133.9113940_rp, 134.9165557_rp, 135.9201006_rp, &
      136.9255989_rp, 137.9294722_rp, 138.9353672_rp, 139.939499_rp, 140.94580_rp, 141.95022_rp, 142.95676_rp, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso])
type (atom_struct), parameter, private :: atm53 = atom_struct(53, "I", 106, [126.90447_rp, 106.94678_rp, 107.94348_rp, 108.9380853_rp, 109.935089_rp, &
      110.9302692_rp, 111.928005_rp, 112.9236501_rp, 113.92185_rp, 114.918048_rp, 115.91681_rp, 116.913648_rp, 117.913074_rp, 118.910074_rp, 119.910087_rp, &
      120.9074051_rp, 121.9075888_rp, 122.9055885_rp, 123.9062090_rp, 124.9046294_rp, 125.9056233_rp, 126.9044719_rp, 127.9058086_rp, 128.9049837_rp, &
      129.9066702_rp, 130.90612630_rp, 131.9079935_rp, 132.9077970_rp, 133.9097588_rp, 134.9100488_rp, 135.914604_rp, 136.9180282_rp, 137.9227264_rp, &
      138.926506_rp, 139.93173_rp, 140.93569_rp, 141.94120_rp, 142.94565_rp, 143.95139_rp, 144.95605_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso])
type (atom_struct), parameter, private :: atm54 = atom_struct(54, "Xe", 108, [131.293_rp, 108.95043_rp, 109.94426_rp, 110.941607_rp, 111.9355590_rp, &
      112.9332217_rp, 113.927980_rp, 114.926294_rp, 115.921581_rp, 116.920359_rp, 117.916179_rp, 118.915411_rp, 119.911784_rp, 120.911453_rp, 121.908368_rp, &
      122.908482_rp, 123.9058920_rp, 124.9063944_rp, 125.9042983_rp, 126.9051829_rp, 127.9035310_rp, 128.9047808611_rp, 129.903509349_rp, 130.90508406_rp, &
      131.9041550856_rp, 132.9059108_rp, 133.90539466_rp, 134.9072278_rp, 135.907214484_rp, 136.91155778_rp, 137.9141463_rp, 138.9187922_rp, 139.9216458_rp, &
      140.9267872_rp, 141.9299731_rp, 142.9353696_rp, 143.9389451_rp, 144.944720_rp, 145.948518_rp, 146.95426_rp, 147.95813_rp, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm55 = atom_struct(55, "Cs", 111, [132.90545196_rp, 111.950309_rp, 112.9444291_rp, 113.941296_rp, 114.93591_rp, &
      115.93337_rp, 116.928617_rp, 117.926560_rp, 118.922377_rp, 119.920677_rp, 120.917227_rp, 121.916108_rp, 122.912996_rp, 123.9122578_rp, 124.9097280_rp, &
      125.909446_rp, 126.9074174_rp, 127.9077487_rp, 128.9060657_rp, 129.9067093_rp, 130.9054649_rp, 131.9064339_rp, 132.9054519610_rp, 133.906718503_rp, &
      134.9059770_rp, 135.9073114_rp, 136.90708923_rp, 137.9110171_rp, 138.9133638_rp, 139.9172831_rp, 140.9200455_rp, 141.9242960_rp, 142.927349_rp, &
      143.932076_rp, 144.935527_rp, 145.940344_rp, 146.944156_rp, 147.94923_rp, 148.95302_rp, 149.95833_rp, 150.96258_rp, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm56 = atom_struct(56, "Ba", 113, [137.327_rp, 113.95066_rp, 114.94737_rp, 115.94128_rp, 116.93814_rp, &
      117.93306_rp, 118.93066_rp, 119.92605_rp, 120.92405_rp, 121.919904_rp, 122.918781_rp, 123.915094_rp, 124.914472_rp, 125.911250_rp, 126.911091_rp, &
      127.9083420_rp, 128.908681_rp, 129.9063207_rp, 130.9069410_rp, 131.9050611_rp, 132.9060074_rp, 133.90450818_rp, 134.90568838_rp, 135.90457573_rp, &
      136.90582714_rp, 137.90524700_rp, 138.90884110_rp, 139.9106057_rp, 140.9144033_rp, 141.9164324_rp, 142.9206253_rp, 143.9229549_rp, 144.9275184_rp, &
      145.930284_rp, 146.935304_rp, 147.938171_rp, 148.94308_rp, 149.94605_rp, 150.95127_rp, 151.95481_rp, 152.96036_rp, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso])
type (atom_struct), parameter, private :: atm57 = atom_struct(57, "La", 115, [138.90547_rp, 115.95630_rp, 116.94999_rp, 117.94673_rp, 118.94099_rp, &
      119.93807_rp, 120.93315_rp, 121.93071_rp, 122.92630_rp, 123.924574_rp, 124.920816_rp, 125.919513_rp, 126.916375_rp, 127.915592_rp, 128.912694_rp, &
      129.912369_rp, 130.910070_rp, 131.910119_rp, 132.908218_rp, 133.908514_rp, 134.906984_rp, 135.907635_rp, 136.9064504_rp, 137.9071149_rp, 138.9063563_rp, &
      139.9094806_rp, 140.9109660_rp, 141.9140909_rp, 142.9160795_rp, 143.919646_rp, 144.921808_rp, 145.925875_rp, 146.928418_rp, 147.932679_rp, 148.93535_rp, &
      149.93947_rp, 150.94232_rp, 151.94682_rp, 152.95036_rp, 153.95517_rp, 154.95901_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm58 = atom_struct(58, "Ce", 118, [140.116_rp, 118.95271_rp, 119.94654_rp, 120.94335_rp, 121.93787_rp, &
      122.93528_rp, 123.93031_rp, 124.92844_rp, 125.923971_rp, 126.922727_rp, 127.918911_rp, 128.918102_rp, 129.914736_rp, 130.914429_rp, 131.911464_rp, &
      132.911520_rp, 133.908928_rp, 134.909161_rp, 135.90712921_rp, 136.90776236_rp, 137.905991_rp, 138.9066551_rp, 139.9054431_rp, 140.9082807_rp, &
      141.9092504_rp, 142.9123921_rp, 143.9136529_rp, 144.917265_rp, 145.918802_rp, 146.9226899_rp, 147.924424_rp, 148.928427_rp, 149.930384_rp, 150.934272_rp, &
      151.93660_rp, 152.94093_rp, 153.94380_rp, 154.94855_rp, 155.95183_rp, 156.95705_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm59 = atom_struct(59, "Pr", 120, [140.90766_rp, 120.95532_rp, 121.95175_rp, 122.94596_rp, 123.94294_rp, &
      124.93770_rp, 125.93524_rp, 126.93071_rp, 127.928791_rp, 128.925095_rp, 129.923590_rp, 130.920235_rp, 131.919255_rp, 132.916331_rp, 133.915697_rp, &
      134.913112_rp, 135.912677_rp, 136.9106792_rp, 137.910754_rp, 138.9089408_rp, 139.9090803_rp, 140.9076576_rp, 141.9100496_rp, 142.9108228_rp, &
      143.9133109_rp, 144.9145182_rp, 145.917680_rp, 146.919008_rp, 147.922130_rp, 148.923736_rp, 149.9266765_rp, 150.928309_rp, 151.931553_rp, 152.933904_rp, &
      153.93753_rp, 154.940509_rp, 155.94464_rp, 156.94789_rp, 157.95241_rp, 158.95589_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm60 = atom_struct(60, "Nd", 123, [144.242_rp, 123.95220_rp, 124.94890_rp, 125.94311_rp, 126.94038_rp, &
      127.93525_rp, 128.93310_rp, 129.928506_rp, 130.927248_rp, 131.923321_rp, 132.922348_rp, 133.918790_rp, 134.918181_rp, 135.914976_rp, 136.914562_rp, &
      137.911950_rp, 138.911954_rp, 139.909550_rp, 140.9096147_rp, 141.9077290_rp, 142.9098200_rp, 143.9100930_rp, 144.9125793_rp, 145.9131226_rp, &
      146.9161061_rp, 147.9168993_rp, 148.9201548_rp, 149.9209022_rp, 150.9238403_rp, 151.924692_rp, 152.9277180_rp, 153.92948_rp, 154.9331357_rp, &
      155.93508_rp, 156.939386_rp, 157.94197_rp, 158.94653_rp, 159.94940_rp, 160.95428_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso])
type (atom_struct), parameter, private :: atm61 = atom_struct(61, "Pm", 125, [145.0_rp, 125.95792_rp, 126.95192_rp, 127.94870_rp, 128.94323_rp, 129.94053_rp, &
      130.93567_rp, 131.93384_rp, 132.929782_rp, 133.928353_rp, 134.924823_rp, 135.923585_rp, 136.920480_rp, 137.919548_rp, 138.916800_rp, 139.916040_rp, &
      140.913555_rp, 141.912890_rp, 142.9109383_rp, 143.9125964_rp, 144.9127559_rp, 145.9147024_rp, 146.9151450_rp, 147.9174819_rp, 148.9183423_rp, &
      149.920991_rp, 150.9212175_rp, 151.923506_rp, 152.9241567_rp, 153.926472_rp, 154.9281370_rp, 155.9311175_rp, 156.9331214_rp, 157.936565_rp, &
      158.939287_rp, 159.94310_rp, 160.94607_rp, 161.95022_rp, 162.95357_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm62 = atom_struct(62, "Sm", 127, [150.36_rp, 127.95842_rp, 128.95476_rp, 129.94900_rp, 130.94618_rp, 131.94087_rp, &
      132.93856_rp, 133.93411_rp, 134.93252_rp, 135.928276_rp, 136.926971_rp, 137.923244_rp, 138.922297_rp, 139.918995_rp, 140.9184816_rp, 141.9152044_rp, &
      142.9146353_rp, 143.9120065_rp, 144.9134173_rp, 145.9130470_rp, 146.9149044_rp, 147.9148292_rp, 148.9171921_rp, 149.9172829_rp, 150.9199398_rp, &
      151.9197397_rp, 152.9221047_rp, 153.9222169_rp, 154.9246477_rp, 155.925536_rp, 156.9284187_rp, 157.9299510_rp, 158.9332172_rp, 159.9353353_rp, &
      160.9391602_rp, 161.94146_rp, 162.94555_rp, 163.94836_rp, 164.95297_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm63 = atom_struct(63, "Eu", 129, [151.964_rp, 129.96369_rp, 130.95784_rp, 131.95467_rp, 132.94929_rp, &
      133.94640_rp, 134.94187_rp, 135.93962_rp, 136.93546_rp, 137.933709_rp, 138.929792_rp, 139.928088_rp, 140.924932_rp, 141.923442_rp, 142.920299_rp, &
      143.918820_rp, 144.9162726_rp, 145.9172110_rp, 146.9167527_rp, 147.918089_rp, 148.9179378_rp, 149.9197077_rp, 150.9198578_rp, 151.9217522_rp, &
      152.9212380_rp, 153.9229870_rp, 154.9229011_rp, 155.9247605_rp, 156.9254334_rp, 157.927799_rp, 158.9291001_rp, 159.931851_rp, 160.933664_rp, &
      161.936989_rp, 162.939196_rp, 163.94274_rp, 164.94559_rp, 165.94962_rp, 166.95289_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso])
type (atom_struct), parameter, private :: atm64 = atom_struct(64, "Gd", 132, [157.25_rp, 132.96133_rp, 133.95566_rp, 134.95245_rp, 135.94730_rp, 136.94502_rp, &
      137.94025_rp, 138.93813_rp, 139.933674_rp, 140.932126_rp, 141.928116_rp, 142.92675_rp, 143.922963_rp, 144.921713_rp, 145.9183188_rp, 146.9191014_rp, &
      147.9181215_rp, 148.9193481_rp, 149.9186644_rp, 150.9203560_rp, 151.9197995_rp, 152.9217580_rp, 153.9208741_rp, 154.9226305_rp, 155.9221312_rp, &
      156.9239686_rp, 157.9241123_rp, 158.9263970_rp, 159.9270624_rp, 160.9296775_rp, 161.9309930_rp, 162.9341769_rp, 163.93583_rp, 164.93936_rp, &
      165.94146_rp, 166.94545_rp, 167.94808_rp, 168.95260_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm65 = atom_struct(65, "Tb", 134, [158.92535_rp, 134.96476_rp, 135.96129_rp, 136.95602_rp, 137.95312_rp, &
      138.94833_rp, 139.94581_rp, 140.94145_rp, 141.93928_rp, 142.935137_rp, 143.933045_rp, 144.92882_rp, 145.927253_rp, 146.9240548_rp, 147.924282_rp, &
      148.9232535_rp, 149.9236649_rp, 150.9231096_rp, 151.924083_rp, 152.9234424_rp, 153.924685_rp, 154.923511_rp, 155.9247552_rp, 156.9240330_rp, &
      157.9254209_rp, 158.9253547_rp, 159.9271756_rp, 160.9275778_rp, 161.929495_rp, 162.9306547_rp, 163.93336_rp, 164.93498_rp, 165.937860_rp, 166.93996_rp, &
      167.94340_rp, 168.94597_rp, 169.94984_rp, 170.95273_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm66 = atom_struct(66, "Dy", 137, [162.500_rp, 137.96250_rp, 138.95959_rp, 139.95402_rp, 140.95128_rp, &
      141.94619_rp, 142.943994_rp, 143.9392695_rp, 144.9374740_rp, 145.9328445_rp, 146.9310827_rp, 147.927157_rp, 148.927322_rp, 149.9255933_rp, 150.9261916_rp, &
      151.9247253_rp, 152.9257724_rp, 153.9244293_rp, 154.925759_rp, 155.9242847_rp, 156.9254707_rp, 157.9244159_rp, 158.9257470_rp, 159.9252046_rp, &
      160.9269405_rp, 161.9268056_rp, 162.9287383_rp, 163.9291819_rp, 164.9317105_rp, 165.9328139_rp, 166.935661_rp, 167.93713_rp, 168.94031_rp, 169.94239_rp, &
      170.94612_rp, 171.94846_rp, 172.95283_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm67 = atom_struct(67, "Ho", 139, [164.93033_rp, 139.96859_rp, 140.96311_rp, 141.96001_rp, 142.95486_rp, &
      143.9521097_rp, 144.9472674_rp, 145.9449935_rp, 146.9401423_rp, 147.937744_rp, 148.933803_rp, 149.933498_rp, 150.9316983_rp, 151.931724_rp, &
      152.9302064_rp, 153.9306068_rp, 154.929104_rp, 155.929706_rp, 156.928254_rp, 157.928946_rp, 158.9277197_rp, 159.928737_rp, 160.9278615_rp, 161.9291023_rp, &
      162.9287410_rp, 163.9302403_rp, 164.9303288_rp, 165.9322909_rp, 166.9331385_rp, 167.935522_rp, 168.936878_rp, 169.939625_rp, 170.94147_rp, 171.94473_rp, &
      172.94702_rp, 173.95095_rp, 174.95362_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm68 = atom_struct(68, "Er", 141, [167.259_rp, 141.97010_rp, 142.96662_rp, 143.96070_rp, 144.95805_rp, &
      145.9524184_rp, 146.949964_rp, 147.944735_rp, 148.942306_rp, 149.937916_rp, 150.937449_rp, 151.935057_rp, 152.935080_rp, 153.9327908_rp, 154.9332159_rp, &
      155.931067_rp, 156.931949_rp, 157.929893_rp, 158.9306918_rp, 159.929077_rp, 160.9300046_rp, 161.9287884_rp, 162.9300408_rp, 163.9292088_rp, &
      164.9307345_rp, 165.9302995_rp, 166.9320546_rp, 167.9323767_rp, 168.9345968_rp, 169.9354702_rp, 170.9380357_rp, 171.9393619_rp, 172.94240_rp, &
      173.94423_rp, 174.94777_rp, 175.94994_rp, 176.95399_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm69 = atom_struct(69, "Tm", 143, [168.93422_rp, 143.97628_rp, 144.97039_rp, 145.96684_rp, 146.9613799_rp, &
      147.958384_rp, 148.95289_rp, 149.95009_rp, 150.945488_rp, 151.944422_rp, 152.942040_rp, 153.941570_rp, 154.939210_rp, 155.938992_rp, 156.936944_rp, &
      157.936980_rp, 158.934975_rp, 159.935263_rp, 160.933549_rp, 161.934002_rp, 162.9326592_rp, 163.933544_rp, 164.9324431_rp, 165.933561_rp, 166.9328562_rp, &
      167.9341774_rp, 168.9342179_rp, 169.9358060_rp, 170.9364339_rp, 171.9384055_rp, 172.9396084_rp, 173.942173_rp, 174.943841_rp, 175.94700_rp, &
      176.94904_rp, 177.95264_rp, 178.95534_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm70 = atom_struct(70, "Yb", 147, [173.054_rp, 147.96758_rp, 148.96436_rp, 149.95852_rp, 150.95540_rp, &
      151.95027_rp, 152.94932_rp, 153.946396_rp, 154.945783_rp, 155.942825_rp, 156.942645_rp, 157.9398705_rp, 158.940055_rp, 159.937557_rp, 160.937907_rp, &
      161.935774_rp, 162.936340_rp, 163.934495_rp, 164.935270_rp, 165.9338747_rp, 166.9349530_rp, 167.9338896_rp, 168.9351825_rp, 169.9347664_rp, &
      170.9363302_rp, 171.9363859_rp, 172.9382151_rp, 173.9388664_rp, 174.9412808_rp, 175.9425764_rp, 176.9452656_rp, 177.946651_rp, 178.95004_rp, &
      179.95212_rp, 180.95589_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm71 = atom_struct(71, "Lu", 149, [174.9668_rp, 149.97355_rp, 150.96768_rp, 151.96412_rp, 152.95875_rp, &
      153.95736_rp, 154.954321_rp, 155.953033_rp, 156.950127_rp, 157.949316_rp, 158.946636_rp, 159.946033_rp, 160.943572_rp, 161.943283_rp, 162.941179_rp, &
      163.941339_rp, 164.939407_rp, 165.939859_rp, 166.938270_rp, 167.938736_rp, 168.9376441_rp, 169.938478_rp, 170.9379170_rp, 171.9390891_rp, 172.9389340_rp, &
      173.9403409_rp, 174.9407752_rp, 175.9426897_rp, 176.9437615_rp, 177.9459580_rp, 178.9473309_rp, 179.949888_rp, 180.95191_rp, 181.95504_rp, 182.957363_rp, &
      183.96091_rp, 184.96362_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm72 = atom_struct(72, "Hf", 152, [178.49_rp, 152.97069_rp, 153.96486_rp, 154.96311_rp, 155.95935_rp, 156.95824_rp, &
      157.954801_rp, 158.953996_rp, 159.950691_rp, 160.950278_rp, 161.9472148_rp, 162.947113_rp, 163.944371_rp, 164.944567_rp, 165.942180_rp, 166.942600_rp, &
      167.940568_rp, 168.941259_rp, 169.939609_rp, 170.940492_rp, 171.939450_rp, 172.940513_rp, 173.9400461_rp, 174.9415092_rp, 175.9414076_rp, 176.9432277_rp, &
      177.9437058_rp, 178.9458232_rp, 179.9465570_rp, 180.9491083_rp, 181.9505612_rp, 182.953530_rp, 183.955446_rp, 184.958862_rp, 185.960897_rp, &
      186.96477_rp, 187.96685_rp, 188.97084_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm73 = atom_struct(73, "Ta", 154, [180.94788_rp, 154.97424_rp, 155.97203_rp, 156.96818_rp, 157.96654_rp, &
      158.963023_rp, 159.961488_rp, 160.958452_rp, 161.957294_rp, 162.954337_rp, 163.953534_rp, 164.950781_rp, 165.950512_rp, 166.948093_rp, 167.948047_rp, &
      168.946011_rp, 169.946175_rp, 170.944476_rp, 171.944895_rp, 172.943750_rp, 173.944454_rp, 174.943737_rp, 175.944857_rp, 176.9444795_rp, 177.945678_rp, &
      178.9459366_rp, 179.9474648_rp, 180.9479958_rp, 181.9501519_rp, 182.9513726_rp, 183.954008_rp, 184.955559_rp, 185.958551_rp, 186.960386_rp, &
      187.963916_rp, 188.96583_rp, 189.96939_rp, 190.97156_rp, 191.97514_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm74 = atom_struct(74, "W", 156, [183.84_rp, 156.97884_rp, 157.97456_rp, 158.97264_rp, 159.96846_rp, 160.96720_rp, &
      161.963499_rp, 162.962524_rp, 163.958961_rp, 164.958281_rp, 165.955031_rp, 166.954805_rp, 167.951806_rp, 168.951779_rp, 169.949232_rp, 170.949451_rp, &
      171.947292_rp, 172.947689_rp, 173.946079_rp, 174.946717_rp, 175.945634_rp, 176.946643_rp, 177.945883_rp, 178.947077_rp, 179.9467108_rp, 180.9481978_rp, &
      181.94820394_rp, 182.95022275_rp, 183.95093092_rp, 184.95341897_rp, 185.9543628_rp, 186.9571588_rp, 187.9584862_rp, 188.961763_rp, 189.963091_rp, &
      190.966531_rp, 191.96817_rp, 192.97178_rp, 193.97367_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm75 = atom_struct(75, "Re", 158, [186.207_rp, 158.98418_rp, 159.98182_rp, 160.97757_rp, 161.97584_rp, &
      162.972080_rp, 163.970453_rp, 164.967103_rp, 165.965761_rp, 166.962595_rp, 167.961573_rp, 168.958766_rp, 169.958220_rp, 170.955716_rp, 171.955420_rp, &
      172.953243_rp, 173.953115_rp, 174.951381_rp, 175.951623_rp, 176.950328_rp, 177.950989_rp, 178.949989_rp, 179.950792_rp, 180.950058_rp, 181.95121_rp, &
      182.9508196_rp, 183.9525228_rp, 184.9529545_rp, 185.9549856_rp, 186.9557501_rp, 187.9581115_rp, 188.9592260_rp, 189.961744_rp, 190.963122_rp, &
      191.966088_rp, 192.967541_rp, 193.97076_rp, 194.97254_rp, 195.97580_rp, 196.97799_rp, 197.98160_rp, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso])
type (atom_struct), parameter, private :: atm76 = atom_struct(76, "Os", 160, [190.23_rp, 160.98903_rp, 161.98443_rp, 162.98241_rp, 163.97802_rp, 164.97660_rp, &
      165.972692_rp, 166.971549_rp, 167.967808_rp, 168.967018_rp, 169.963578_rp, 170.963174_rp, 171.960017_rp, 172.959808_rp, 173.957064_rp, 174.956945_rp, &
      175.954806_rp, 176.954966_rp, 177.953254_rp, 178.953817_rp, 179.952375_rp, 180.953247_rp, 181.952110_rp, 182.953125_rp, 183.9524885_rp, 184.9540417_rp, &
      185.9538350_rp, 186.9557474_rp, 187.9558352_rp, 188.9581442_rp, 189.9584437_rp, 190.9609264_rp, 191.9614770_rp, 192.9641479_rp, 193.9651772_rp, &
      194.968318_rp, 195.969641_rp, 196.97283_rp, 197.97441_rp, 198.97801_rp, 199.97984_rp, 200.98364_rp, 201.98595_rp, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm77 = atom_struct(77, "Ir", 163, [192.217_rp, 163.99191_rp, 164.98750_rp, 165.98566_rp, 166.981666_rp, &
      167.979907_rp, 168.976298_rp, 169.974922_rp, 170.971640_rp, 171.970607_rp, 172.967506_rp, 173.966861_rp, 174.964150_rp, 175.963650_rp, 176.961301_rp, &
      177.961082_rp, 178.959120_rp, 179.959229_rp, 180.957625_rp, 181.958076_rp, 182.956840_rp, 183.957476_rp, 184.956698_rp, 185.957944_rp, 186.957542_rp, &
      187.958828_rp, 188.958715_rp, 189.9605412_rp, 190.9605893_rp, 191.9626002_rp, 192.9629216_rp, 193.9650735_rp, 194.9659747_rp, 195.968397_rp, &
      196.969655_rp, 197.97228_rp, 198.973805_rp, 199.97680_rp, 200.97864_rp, 201.98199_rp, 202.98423_rp, 203.98960_rp, no_iso, no_iso, no_iso, no_iso, &
      no_iso])
type (atom_struct), parameter, private :: atm78 = atom_struct(78, "Pt", 165, [195.084_rp, 165.99486_rp, 166.99269_rp, 167.98813_rp, 168.98657_rp, &
      169.982496_rp, 170.981245_rp, 171.977351_rp, 172.976443_rp, 173.972820_rp, 174.972410_rp, 175.968938_rp, 176.968470_rp, 177.965650_rp, 178.9653590_rp, &
      179.963032_rp, 180.963098_rp, 181.961172_rp, 182.961597_rp, 183.959915_rp, 184.960614_rp, 185.959351_rp, 186.960617_rp, 187.9593889_rp, 188.960831_rp, &
      189.9599297_rp, 190.9616729_rp, 191.9610387_rp, 192.9629824_rp, 193.9626809_rp, 194.9647917_rp, 195.96495209_rp, 196.96734069_rp, 197.9678949_rp, &
      198.9705952_rp, 199.971443_rp, 200.974513_rp, 201.975639_rp, 202.97893_rp, 203.98076_rp, 204.98608_rp, 205.98966_rp, no_iso, no_iso, no_iso, &
      no_iso, no_iso])
type (atom_struct), parameter, private :: atm79 = atom_struct(79, "Au", 168, [196.966569_rp, 168.99808_rp, 169.99597_rp, 170.991876_rp, 171.989942_rp, &
      172.986241_rp, 173.984717_rp, 174.981304_rp, 175.980250_rp, 176.976870_rp, 177.976032_rp, 178.973174_rp, 179.972523_rp, 180.970079_rp, 181.969618_rp, &
      182.967591_rp, 183.967452_rp, 184.965790_rp, 185.965953_rp, 186.964543_rp, 187.965349_rp, 188.963948_rp, 189.964698_rp, 190.963702_rp, 191.964814_rp, &
      192.9641373_rp, 193.9654178_rp, 194.9650352_rp, 195.9665699_rp, 196.96656879_rp, 197.96824242_rp, 198.96876528_rp, 199.970756_rp, 200.9716575_rp, &
      201.973856_rp, 202.9751544_rp, 203.97783_rp, 204.97985_rp, 205.98474_rp, 206.98840_rp, 207.99345_rp, 208.99735_rp, 210.00250_rp, no_iso, no_iso, &
      no_iso, no_iso])
type (atom_struct), parameter, private :: atm80 = atom_struct(80, "Hg", 170, [200.592_rp, 171.00353_rp, 171.99881_rp, 172.99709_rp, 173.992865_rp, &
      174.991441_rp, 175.987361_rp, 176.986277_rp, 177.982484_rp, 178.981831_rp, 179.978260_rp, 180.977819_rp, 181.974689_rp, 182.9744448_rp, 183.971714_rp, &
      184.971899_rp, 185.969362_rp, 186.969814_rp, 187.967567_rp, 188.968195_rp, 189.966323_rp, 190.967157_rp, 191.965635_rp, 192.966653_rp, 193.9654491_rp, &
      194.966721_rp, 195.9658326_rp, 196.9672128_rp, 197.96676860_rp, 198.96828064_rp, 199.96832659_rp, 200.97030284_rp, 201.97064340_rp, 202.9728728_rp, &
      203.97349398_rp, 204.9760734_rp, 205.977514_rp, 206.982300_rp, 207.985759_rp, 208.99072_rp, 209.99424_rp, 210.99933_rp, 212.00296_rp, 213.00823_rp, &
      214.01200_rp, 215.01740_rp, 216.02132_rp])
type (atom_struct), parameter, private :: atm81 = atom_struct(81, "Tl", 175, [204.3835_rp, 176.000624_rp, 176.996431_rp, 177.99485_rp, 178.991111_rp, &
      179.990057_rp, 180.9862600_rp, 181.985713_rp, 182.982193_rp, 183.981886_rp, 184.978789_rp, 185.978651_rp, 186.9759063_rp, 187.976021_rp, 188.973588_rp, &
      189.973828_rp, 190.9717842_rp, 191.972225_rp, 192.9705020_rp, 193.971081_rp, 194.969774_rp, 195.970481_rp, 196.969576_rp, 197.970483_rp, 198.969877_rp, &
      199.9709633_rp, 200.970822_rp, 201.972102_rp, 202.9723446_rp, 203.9738639_rp, 204.9744278_rp, 205.9761106_rp, 206.9774197_rp, 207.9820190_rp, &
      208.9853594_rp, 209.990074_rp, 210.993475_rp, 211.99834_rp, 213.001915_rp, 214.00694_rp, 215.01064_rp, 216.01580_rp, 217.01966_rp, 218.02479_rp, &
      no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm82 = atom_struct(82, "Pb", 177, [207.2_rp, 178.003831_rp, 179.002201_rp, 179.997928_rp, 180.996653_rp, &
      181.992672_rp, 182.991872_rp, 183.988136_rp, 184.987610_rp, 185.984238_rp, 186.9839109_rp, 187.980875_rp, 188.980807_rp, 189.978082_rp, 190.978276_rp, &
      191.975775_rp, 192.976173_rp, 193.974012_rp, 194.974543_rp, 195.972774_rp, 196.9734312_rp, 197.972034_rp, 198.972913_rp, 199.971819_rp, 200.972883_rp, &
      201.9721520_rp, 202.9733911_rp, 203.9730440_rp, 204.9744822_rp, 205.9744657_rp, 206.9758973_rp, 207.9766525_rp, 208.9810905_rp, 209.9841889_rp, &
      210.9887371_rp, 211.9918977_rp, 212.9965629_rp, 213.9998059_rp, 215.00474_rp, 216.00803_rp, 217.01314_rp, 218.01659_rp, 219.02177_rp, 220.02541_rp, &
      no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm83 = atom_struct(83, "Bi", 183, [208.98040_rp, 184.001275_rp, 184.997600_rp, 185.996644_rp, 186.993147_rp, &
      187.992287_rp, 188.989195_rp, 189.988622_rp, 190.9857866_rp, 191.985469_rp, 192.982960_rp, 193.982785_rp, 194.9806488_rp, 195.980667_rp, 196.9788651_rp, &
      197.979206_rp, 198.977673_rp, 199.978131_rp, 200.977010_rp, 201.977734_rp, 202.976893_rp, 203.9778361_rp, 204.9773867_rp, 205.9784993_rp, 206.9784710_rp, &
      207.9797425_rp, 208.9803991_rp, 209.9841207_rp, 210.9872697_rp, 211.9912860_rp, 212.9943851_rp, 213.998712_rp, 215.001770_rp, 216.006306_rp, &
      217.009372_rp, 218.014188_rp, 219.01748_rp, 220.02235_rp, 221.02587_rp, 222.03078_rp, 223.03450_rp, 224.03947_rp, no_iso, no_iso, no_iso, no_iso, &
      no_iso])
type (atom_struct), parameter, private :: atm84 = atom_struct(84, "Po", 185, [209.0_rp, 186.004393_rp, 187.003041_rp, 187.999416_rp, 188.998473_rp, &
      189.995101_rp, 190.9945585_rp, 191.991336_rp, 192.991026_rp, 193.988186_rp, 194.988126_rp, 195.985526_rp, 196.985660_rp, 197.983389_rp, 198.983667_rp, &
      199.981799_rp, 200.9822598_rp, 201.980758_rp, 202.9814161_rp, 203.980310_rp, 204.981203_rp, 205.9804740_rp, 206.9815938_rp, 207.9812461_rp, &
      208.9824308_rp, 209.9828741_rp, 210.9866536_rp, 211.9888684_rp, 212.9928576_rp, 213.9952017_rp, 214.9994201_rp, 216.0019152_rp, 217.0063182_rp, &
      218.0089735_rp, 219.013614_rp, 220.016386_rp, 221.021228_rp, 222.024140_rp, 223.02907_rp, 224.03211_rp, 225.03707_rp, 226.04031_rp, 227.04539_rp, &
      no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm85 = atom_struct(85, "At", 190, [210.0_rp, 191.004148_rp, 192.003152_rp, 192.999927_rp, 193.999236_rp, &
      194.9962685_rp, 195.995800_rp, 196.993189_rp, 197.992784_rp, 198.9905277_rp, 199.990351_rp, 200.9884171_rp, 201.988630_rp, 202.986943_rp, 203.987251_rp, &
      204.986076_rp, 205.986657_rp, 206.985800_rp, 207.9866133_rp, 208.9861702_rp, 209.9871479_rp, 210.9874966_rp, 211.9907377_rp, 212.9929370_rp, &
      213.9963721_rp, 214.9986528_rp, 216.0024236_rp, 217.0047192_rp, 218.008695_rp, 219.0111618_rp, 220.015433_rp, 221.018017_rp, 222.022494_rp, &
      223.025151_rp, 224.029749_rp, 225.03263_rp, 226.03716_rp, 227.04024_rp, 228.04475_rp, 229.04812_rp, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso])
type (atom_struct), parameter, private :: atm86 = atom_struct(86, "Rn", 192, [222.0_rp, 193.009708_rp, 194.006144_rp, 195.005422_rp, 196.002116_rp, &
      197.001585_rp, 197.998679_rp, 198.998390_rp, 199.995690_rp, 200.995628_rp, 201.993264_rp, 202.993388_rp, 203.991430_rp, 204.991719_rp, 205.990214_rp, &
      206.9907303_rp, 207.989635_rp, 208.990415_rp, 209.9896891_rp, 210.9906011_rp, 211.9907039_rp, 212.9938831_rp, 213.9953630_rp, 214.9987459_rp, &
      216.0002719_rp, 217.0039280_rp, 218.0056016_rp, 219.0094804_rp, 220.0113941_rp, 221.0155371_rp, 222.0175782_rp, 223.0218893_rp, 224.024096_rp, &
      225.028486_rp, 226.030861_rp, 227.035304_rp, 228.037835_rp, 229.042257_rp, 230.04514_rp, 231.04987_rp, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso])
type (atom_struct), parameter, private :: atm87 = atom_struct(87, "Fr", 198, [223.0_rp, 199.007259_rp, 200.006586_rp, 201.003867_rp, 202.003320_rp, &
      203.0009407_rp, 204.000652_rp, 204.9985939_rp, 205.998666_rp, 206.996946_rp, 207.997138_rp, 208.995955_rp, 209.996422_rp, 210.995556_rp, 211.9962257_rp, &
      212.9961860_rp, 213.9989713_rp, 215.0003418_rp, 216.0031899_rp, 217.0046323_rp, 218.0075787_rp, 219.0092524_rp, 220.0123277_rp, 221.0142552_rp, &
      222.017552_rp, 223.0197360_rp, 224.023398_rp, 225.025573_rp, 226.029566_rp, 227.031869_rp, 228.035823_rp, 229.038298_rp, 230.042416_rp, 231.045158_rp, &
      232.04937_rp, 233.05264_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm88 = atom_struct(88, "Ra", 200, [226.0_rp, 201.01271_rp, 202.009760_rp, 203.009304_rp, 204.006492_rp, &
      205.006268_rp, 206.003828_rp, 207.003799_rp, 208.001841_rp, 209.001990_rp, 210.000494_rp, 211.0008932_rp, 211.999787_rp, 213.000384_rp, 214.0000997_rp, &
      215.0027204_rp, 216.0035334_rp, 217.0063207_rp, 218.007141_rp, 219.0100855_rp, 220.0110259_rp, 221.0139177_rp, 222.0153748_rp, 223.0185023_rp, &
      224.0202120_rp, 225.0236119_rp, 226.0254103_rp, 227.0291783_rp, 228.0310707_rp, 229.034942_rp, 230.037055_rp, 231.041027_rp, 232.0434753_rp, &
      233.047582_rp, 234.050342_rp, 235.05497_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm89 = atom_struct(89, "Ac", 205, [227.0_rp, 206.014452_rp, 207.011966_rp, 208.011550_rp, 209.009495_rp, &
      210.009436_rp, 211.007732_rp, 212.007813_rp, 213.006609_rp, 214.006918_rp, 215.006475_rp, 216.008743_rp, 217.009344_rp, 218.011642_rp, 219.012421_rp, &
      220.0147549_rp, 221.015592_rp, 222.0178442_rp, 223.0191377_rp, 224.0217232_rp, 225.0232300_rp, 226.0260984_rp, 227.0277523_rp, 228.0310215_rp, &
      229.032956_rp, 230.036327_rp, 231.038393_rp, 232.042034_rp, 233.044346_rp, 234.048139_rp, 235.050840_rp, 236.054988_rp, 237.05827_rp, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm90 = atom_struct(90, "Th", 207, [232.0377_rp, 208.017900_rp, 209.017753_rp, 210.015094_rp, 211.014929_rp, &
      212.012988_rp, 213.013009_rp, 214.011500_rp, 215.0117248_rp, 216.011056_rp, 217.013117_rp, 218.013276_rp, 219.015537_rp, 220.015748_rp, 221.018184_rp, &
      222.018469_rp, 223.0208119_rp, 224.021464_rp, 225.0239514_rp, 226.0249034_rp, 227.0277042_rp, 228.0287413_rp, 229.0317627_rp, 230.0331341_rp, &
      231.0363046_rp, 232.0380558_rp, 233.0415823_rp, 234.0436014_rp, 235.047255_rp, 236.049657_rp, 237.053629_rp, 238.05650_rp, 239.06077_rp, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm91 = atom_struct(91, "Pa", 211, [231.03588_rp, 212.023203_rp, 213.021109_rp, 214.020918_rp, 215.019183_rp, &
      216.019109_rp, 217.018325_rp, 218.020059_rp, 219.019904_rp, 220.021705_rp, 221.021875_rp, 222.023784_rp, 223.023963_rp, 224.0256176_rp, 225.026131_rp, &
      226.027948_rp, 227.0288054_rp, 228.0310517_rp, 229.0320972_rp, 230.0345410_rp, 231.0358842_rp, 232.0385917_rp, 233.0402472_rp, 234.0433072_rp, &
      235.045399_rp, 236.048668_rp, 237.051023_rp, 238.054637_rp, 239.05726_rp, 240.06098_rp, 241.06408_rp, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm92 = atom_struct(92, "U", 216, [238.02891_rp, 217.02466_rp, 218.023523_rp, 219.024999_rp, 220.02462_rp, &
      221.02628_rp, 222.02600_rp, 223.027739_rp, 224.027605_rp, 225.029391_rp, 226.029339_rp, 227.031157_rp, 228.031371_rp, 229.0335063_rp, 230.0339401_rp, &
      231.0362939_rp, 232.0371563_rp, 233.0396355_rp, 234.0409523_rp, 235.0439301_rp, 236.0455682_rp, 237.0487304_rp, 238.0507884_rp, 239.0542935_rp, &
      240.0565934_rp, 241.06033_rp, 242.06293_rp, 243.06699_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm93 = atom_struct(93, "Np", 218, [237.0_rp, 219.03143_rp, 220.03254_rp, 221.03204_rp, 222.03330_rp, 223.03285_rp, &
      224.03422_rp, 225.033911_rp, 226.035188_rp, 227.034957_rp, 228.036067_rp, 229.036264_rp, 230.037828_rp, 231.038245_rp, 232.04011_rp, 233.040741_rp, &
      234.0428953_rp, 235.0440635_rp, 236.046570_rp, 237.0481736_rp, 238.0509466_rp, 239.0529392_rp, 240.056165_rp, 241.058253_rp, 242.06164_rp, 243.064280_rp, &
      244.06785_rp, 245.07080_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm94 = atom_struct(94, "Pu", 227, [244.0_rp, 228.038732_rp, 229.040144_rp, 230.039650_rp, 231.041102_rp, &
      232.041185_rp, 233.042998_rp, 234.0433174_rp, 235.045286_rp, 236.0460581_rp, 237.0484098_rp, 238.0495601_rp, 239.0521636_rp, 240.0538138_rp, &
      241.0568517_rp, 242.0587428_rp, 243.0620036_rp, 244.0642053_rp, 245.067826_rp, 246.070205_rp, 247.07419_rp, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm95 = atom_struct(95, "Am", 229, [243.0_rp, 230.04609_rp, 231.04556_rp, 232.04645_rp, 233.04644_rp, 234.04773_rp, &
      235.047908_rp, 236.04943_rp, 237.049996_rp, 238.051985_rp, 239.0530247_rp, 240.055300_rp, 241.0568293_rp, 242.0595494_rp, 243.0613813_rp, 244.0642851_rp, &
      245.0664548_rp, 246.069775_rp, 247.07209_rp, 248.07575_rp, 249.07848_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm96 = atom_struct(96, "Cm", 231, [247.0_rp, 232.04982_rp, 233.050770_rp, 234.050160_rp, 235.05154_rp, &
      236.051374_rp, 237.052869_rp, 238.053081_rp, 239.054910_rp, 240.0555297_rp, 241.0576532_rp, 242.0588360_rp, 243.0613893_rp, 244.0627528_rp, &
      245.0654915_rp, 246.0672238_rp, 247.0703541_rp, 248.0723499_rp, 249.0759548_rp, 250.078358_rp, 251.082286_rp, 252.08487_rp, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm97 = atom_struct(97, "Bk", 233, [247.0_rp, 234.05727_rp, 235.05658_rp, 236.05748_rp, 237.05710_rp, 238.05820_rp, &
      239.05824_rp, 240.05976_rp, 241.06016_rp, 242.06198_rp, 243.0630078_rp, 244.065181_rp, 245.0663618_rp, 246.068673_rp, 247.0703073_rp, 248.073088_rp, &
      249.0749877_rp, 250.0783167_rp, 251.080762_rp, 252.08431_rp, 253.08688_rp, 254.09060_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm98 = atom_struct(98, "Cf", 236, [251.0_rp, 237.062198_rp, 238.06149_rp, 239.06253_rp, 240.062256_rp, &
      241.06369_rp, 242.063754_rp, 243.06548_rp, 244.0660008_rp, 245.0680487_rp, 246.0688055_rp, 247.070965_rp, 248.0721851_rp, 249.0748539_rp, 250.0764062_rp, &
      251.0795886_rp, 252.0816272_rp, 253.0851345_rp, 254.087324_rp, 255.09105_rp, 256.09344_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso])
type (atom_struct), parameter, private :: atm99 = atom_struct(99, "Es", 238, [252.0_rp, 239.06823_rp, 240.06892_rp, 241.06856_rp, 242.06957_rp, 243.06951_rp, &
      244.07088_rp, 245.07125_rp, 246.07290_rp, 247.073622_rp, 248.075471_rp, 249.076411_rp, 250.07861_rp, 251.0799936_rp, 252.082980_rp, 253.0848257_rp, &
      254.0880222_rp, 255.090275_rp, 256.09360_rp, 257.09598_rp, 258.09952_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm100 = atom_struct(100, "Fm", 240, [257.0_rp, 241.07421_rp, 242.07343_rp, 243.07446_rp, 244.07404_rp, &
      245.07535_rp, 246.075350_rp, 247.07694_rp, 248.0771865_rp, 249.0789275_rp, 250.0795210_rp, 251.081540_rp, 252.0824671_rp, 253.0851846_rp, 254.0868544_rp, &
      255.0899640_rp, 256.0917745_rp, 257.0951061_rp, 258.09708_rp, 259.10060_rp, 260.10281_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso])
type (atom_struct), parameter, private :: atm101 = atom_struct(101, "Md", 244, [258.0_rp, 245.08081_rp, 246.08171_rp, 247.08152_rp, 248.08282_rp, &
      249.08291_rp, 250.08441_rp, 251.084774_rp, 252.08643_rp, 253.087144_rp, 254.08959_rp, 255.0910841_rp, 256.09389_rp, 257.0955424_rp, 258.0984315_rp, &
      259.10051_rp, 260.10365_rp, 261.10583_rp, 262.10910_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm102 = atom_struct(102, "No", 247, [259.0_rp, 248.08655_rp, 249.08780_rp, 250.08756_rp, 251.08894_rp, &
      252.088967_rp, 253.0905641_rp, 254.090956_rp, 255.093191_rp, 256.0942829_rp, 257.0968878_rp, 258.09821_rp, 259.10103_rp, 260.10264_rp, 261.10570_rp, &
      262.10746_rp, 263.11071_rp, 264.11273_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm103 = atom_struct(103, "Lr", 250, [262.0_rp, 251.09418_rp, 252.09526_rp, 253.09509_rp, 254.09648_rp, &
      255.096562_rp, 256.098494_rp, 257.099418_rp, 258.10176_rp, 259.102902_rp, 260.10550_rp, 261.10688_rp, 262.10961_rp, 263.11136_rp, 264.11420_rp, &
      265.11619_rp, 266.11983_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm104 = atom_struct(104, "Rf", 252, [261.0_rp, 253.10044_rp, 254.10005_rp, 255.10127_rp, 256.101152_rp, &
      257.102918_rp, 258.103428_rp, 259.105596_rp, 260.10644_rp, 261.108773_rp, 262.10992_rp, 263.11249_rp, 264.11388_rp, 265.11668_rp, 266.11817_rp, &
      267.12179_rp, 268.12397_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm105 = atom_struct(105, "Db", 254, [262.0_rp, 255.10707_rp, 256.10789_rp, 257.10758_rp, 258.10928_rp, &
      259.109492_rp, 260.11130_rp, 261.11192_rp, 262.11407_rp, 263.11499_rp, 264.11741_rp, 265.11861_rp, 266.12103_rp, 267.12247_rp, 268.12567_rp, &
      269.12791_rp, 270.13136_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm106 = atom_struct(106, "Sg", 257, [266.0_rp, 258.11298_rp, 259.11440_rp, 260.114384_rp, 261.115949_rp, &
      262.116337_rp, 263.11829_rp, 264.11893_rp, 265.12109_rp, 266.12198_rp, 267.12436_rp, 268.12539_rp, 269.12863_rp, 270.13043_rp, 271.13393_rp, &
      272.13589_rp, 273.13958_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm107 = atom_struct(107, "Bh", 259, [264.0_rp, 260.12166_rp, 261.12145_rp, 262.12297_rp, 263.12292_rp, &
      264.12459_rp, 265.12491_rp, 266.12679_rp, 267.12750_rp, 268.12969_rp, 269.13042_rp, 270.13336_rp, 271.13526_rp, 272.13826_rp, 273.14024_rp, &
      274.14355_rp, 275.14567_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm108 = atom_struct(108, "Hs", 262, [277.0_rp, 263.12852_rp, 264.128357_rp, 265.129793_rp, 266.130046_rp, &
      267.13167_rp, 268.13186_rp, 269.13375_rp, 270.13429_rp, 271.13717_rp, 272.13850_rp, 273.14168_rp, 274.14330_rp, 275.14667_rp, 276.14846_rp, &
      277.15190_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm109 = atom_struct(109, "Mt", 264, [268.0_rp, 265.13600_rp, 266.13737_rp, 267.13719_rp, 268.13865_rp, &
      269.13882_rp, 270.14033_rp, 271.14074_rp, 272.14341_rp, 273.14440_rp, 274.14724_rp, 275.14882_rp, 276.15159_rp, 277.15327_rp, 278.15631_rp, &
      279.15808_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm110 = atom_struct(110, "Ds", 266, [274.0_rp, 267.14377_rp, 268.14348_rp, 269.144752_rp, 270.144584_rp, &
      271.14595_rp, 272.14602_rp, 273.14856_rp, 274.14941_rp, 275.15203_rp, 276.15303_rp, 277.15591_rp, 278.15704_rp, 279.16010_rp, 280.16131_rp, &
      281.16451_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm111 = atom_struct(111, "Rg", 271, [278.0_rp, 272.15327_rp, 273.15313_rp, 274.15525_rp, 275.15594_rp, &
      276.15833_rp, 277.15907_rp, 278.16149_rp, 279.16272_rp, 280.16514_rp, 281.16636_rp, 282.16912_rp, 283.17054_rp, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm112 = atom_struct(112, "Cn", 275, [281.0_rp, 276.16141_rp, 277.16364_rp, 278.16416_rp, 279.16654_rp, &
      280.16715_rp, 281.16975_rp, 282.17050_rp, 283.17327_rp, 284.17416_rp, 285.17712_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm113 = atom_struct(113, "Nh", 277, [283.0_rp, 278.17058_rp, 279.17095_rp, 280.17293_rp, 281.17348_rp, &
      282.17567_rp, 283.17657_rp, 284.17873_rp, 285.17973_rp, 286.18221_rp, 287.18339_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm114 = atom_struct(114, "Fl", 284, [287.0_rp, 285.18364_rp, 286.18423_rp, 287.18678_rp, 288.18757_rp, &
      289.19042_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm115 = atom_struct(115, "Mc", 286, [289.0_rp, 287.19070_rp, 288.19274_rp, 289.19363_rp, 290.19598_rp, &
      291.19707_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm116 = atom_struct(116, "Lv", 288, [291.0_rp, 289.19816_rp, 290.19864_rp, 291.20108_rp, 292.20174_rp, &
      293.20449_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm117 = atom_struct(117, "Ts", 290, [293.0_rp, 291.20553_rp, 292.20746_rp, 293.20824_rp, 294.21046_rp, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso])
type (atom_struct), parameter, private :: atm118 = atom_struct(118, "Og", 292, [294.0_rp, 293.21356_rp, 294.21392_rp, 295.21624_rp, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, &
      no_iso, no_iso, no_iso, no_iso, no_iso])
type(atom_struct), parameter, private :: atom(1:118) = [atm1, atm2, atm3, atm4, atm5, atm6, atm7, atm8, atm9, atm10, atm11, atm12, atm13, atm14, atm15, atm16, &
      atm17, atm18, atm19, atm20, atm21, atm22, atm23, atm24, atm25, atm26, atm27, atm28, atm29, atm30, atm31, atm32, atm33, atm34, atm35, atm36, atm37, &
      atm38, atm39, atm40, atm41, atm42, atm43, atm44, atm45, atm46, atm47, atm48, atm49, atm50, atm51, atm52, atm53, atm54, atm55, atm56, atm57, atm58, &
      atm59, atm60, atm61, atm62, atm63, atm64, atm65, atm66, atm67, atm68, atm69, atm70, atm71, atm72, atm73, atm74, atm75, atm76, atm77, atm78, atm79, &
      atm80, atm81, atm82, atm83, atm84, atm85, atm86, atm87, atm88, atm89, atm90, atm91, atm92, atm93, atm94, atm95, atm96, atm97, atm98, atm99, atm100, &
      atm101, atm102, atm103, atm104, atm105, atm106, atm107, atm108, atm109, atm110, atm111, atm112, atm113, atm114, atm115, atm116, atm117, atm118]

!----------------------
! Known Molecules

character(8), parameter :: molecular_name(18) = [character(8) :: &
                        'CH2', 'CH3', 'CH4', 'CO', 'CO2', 'D2', 'D2O', 'H2', 'H2O', 'N2', 'HF',  &
                        'OH', 'O2', 'NH2', 'NH3', 'C2H3', 'C2H4', 'C2H5']

real(rp), parameter, private ::  molecular_mass(18) = [ &
  14.026579004642441_rp, & ! for CH2
  15.034498933118178_rp, & ! for CH3
  16.04249886158824_rp, &  ! for CH4
  28.01009801234052_rp, &  ! for CO
  44.009496876987235_rp, & ! for CO2
  4.028203269749646_rp, &  ! for D2
  20.02759857879661_rp, &  ! for D2O
  2.015879856948637_rp, &  ! for H2
  18.01527872159535_rp, &  ! for H2O
  28.013398012106343_rp, & ! for N2
  20.006341580305058_rp, & ! for HF
  17.00733879312103_rp, &  ! for OH
  31.998797729293425_rp, & ! for O2
  16.0228_rp, 17.030518791476126_rp, &       ! for NH2, NH3
  27.0457_rp, 28.0536_rp, 29.0615_rp]        ! For C2H3, C2H4, C2H5

contains

!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!+
! Function antiparticle (species) result (anti_species)
!
! Routine to return the antiparticle ID given the particle ID.
! For atoms, the "antiparticle" is just the atom with the charge negated.
!
! Input:
!   species       -- integer: Particle ID.
!
! Output:
!   anti_species  -- integer: Antiparticle ID.
!-

function antiparticle (species) result (anti_species)

integer species, anti_species

!

if (lb_subatomic <= species .and. species <= ub_subatomic) then
  anti_species = antiparticle_of_subatomic(species)
  return

else
  anti_species = -species
endif

end function antiparticle

!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!+
! Function species_of (mass, charge) result (species)
!
! Routine to return the integer ID index of a particle species given the mass and charge.
! Note: Currently this routine only works for subatomic particles and is used for decoding PTC flat files.
!
! Input:
!   mass    -- real(rp): Mass of the particle
!   charge  -- integer: Charge of the particle.
!
! Output:
!   species -- integer: Species ID. Will return invalid$ if name is not valid.
!-

function species_of (mass, charge) result(species)

real(rp) mass
integer charge, species

!

do species = lb_subatomic, ub_subatomic
  if (charge == charge_of_subatomic(species) .and. abs(mass - mass_of_subatomic(species)) <= 1d-6 * mass) return
enddo

species = invalid$

end function species_of

!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!+
! Function species_id (name, default, print_err) result(species)
!
! Routine to return the integer ID index of a particle species given the name.
!
! For subatomic particles, the case does not matter. 
! For all other types of particles, the case does matter.
!
! Input:
!   name        -- character(20): Name of the species.
!   default     -- integer, optional: Default species to use if name is blank or 'ref_species'.
!                   If not present, a blank name is an error.
!   print_err   -- logical, optional: Print error message? Default is True. If False, return species = invalid$,
!
! Output:
!   species     -- integer: Species ID. Will return invalid$ if name is not valid.
!                       Will return not_set$ if name is blank
!-

function species_id (name, default, print_err) result(species)

real(rp) :: mol_mass
integer, optional :: default
integer ::  species, charge, i, ix, ix1, ix2, iso, ios, n_nuc
character(*) :: name
character(20) :: nam
character(*), parameter :: r_name = 'species_id'
logical, optional :: print_err
logical anti, do_print

! Init

if (name == '' .or. upcase(name) == 'REF_SPECIES') then
  if (present(default)) then
    species = default
  else
    species = not_set$
  endif
  return
endif

do_print = logic_option(.true., print_err)
species = invalid$
iso = 0
nam = name
mol_mass = 0
n_nuc = 0

if (upcase(nam) == 'PION_PLUS')  nam = 'pion+'  ! Old style
if (upcase(nam) == 'PION_0')     nam = 'pion0'  ! Old style
if (upcase(nam) == 'PION_MINUS') nam = 'pion-'  ! Old style

! Subatomic particle?

call match_word(nam, subatomic_species_name, ix, .false., .false.)
if (ix > 0) then
  species = ix + lb_subatomic - 1
  return
endif

! Get mass and charge

ix1 = max(index(nam, '+'), index(nam, '-'))
ix2 = index(nam, '@M')

if (ix1 > ix2) then
  if (.not. get_this_charge()) return
  if (.not. get_this_mass()) return
else
  if (.not. get_this_mass()) return
  if (.not. get_this_charge()) return
endif

! If unknown molecule

if (nam == '') then
  if (mol_mass == 0) return
  species = abs(charge * int(z'1000000')) + 200 * int(z'10000') + nint(100 * mol_mass)
  if (charge < 0) species = -species
  return
endif

! Isotopic number (e.g. #12C)

if (nam(1:1) == '#') then
 do i = 2, 10
   if (index('0123456789', nam(i:i)) == 0) exit
 enddo
 if (i == 2) then 
   n_nuc = 0 ! use the standard if no number was given, i.e. #C 
 else
   read (nam(2:i-1), *) n_nuc
 endif
 
 nam = nam(i:)
endif

! Is molecule?

call match_word(nam, molecular_name, ix, .true., .false.)
if (ix > 0) then
  if (n_nuc /= 0) return  ! Isotopic number not allowed with molecules
  ! Only 3+2 digits are allowed
  if (mol_mass*100 > 99999) then
    if (do_print) then
      call out_io (s_abort$, r_name, 'SPECIFIED MOLECULE MASS TOO LARGE FOR ' // name)
      if (global_com%exit_on_error) call err_exit
    endif
    return
  endif

  if (n_nuc /= 0) then
    if (do_print) then
      call out_io (s_abort$, r_name, 'SPECIFYING THE NUMBER OF NUCLEONS FOR A MOLECULE IS INVALID: ' // name)
      if (global_com%exit_on_error) call err_exit
    endif
    return
  endif  

  species = abs(charge * int(z'1000000')) + (ix+200) * int(z'10000') + nint(mol_mass * 100)
  if (charge < 0) species = -species
  return  
endif

! Is atom or anti-atom

if (mol_mass /= 0) then
  if (do_print) then
    call out_io (s_abort$, r_name, 'SPECIFYING A MASS FOR AN ATOM IS INVALID: ' // name)
    if (global_com%exit_on_error) call err_exit
  endif
  return
endif  

anti = .false.
if (nam(1:4) == 'anti') then
  anti = .true.
  nam = nam(5:)
endif

select case (nam)
case ('Uut');  ix = 113
case ('Uup');  ix = 115
case ('Uus');  ix = 117
case ('Uuo');  ix = 118
case default
  call match_word(nam, atomic_name, ix, .true., .false.)
  if (ix <= 0) then
    if (do_print) then
      call out_io (s_error$, r_name, 'CANNOT DECODE ATOM NAME: ' // name)
      if (global_com%exit_on_error) call err_exit
    endif
    return
  endif
end select

if (anti) then
  species = (abs(charge * int(z'1000000')) + anti_atom$ * int(z'10000') + ix * 2 * int(z'100') + n_nuc)
else
  species = (abs(charge * int(z'1000000')) + ix * int(z'10000') + n_nuc)
endif

if (charge < 0) species = -species

!----------------------------------------------------------------------------------------
contains

! Strip off charge.

function get_this_charge () result (is_ok)

logical is_ok

!

is_ok = .false.
charge = 0

ix = index(nam, '+')
if (ix /= 0) then
  select case (nam(ix+1:ix+1))
  case ('+')
    do i = ix+1, len(nam)
      if (nam(i:i) == ' ') exit
      if (nam(i:i) /= '+') return
    enddo
    charge = i - ix
  case (' ')
    charge = +1
  case default
    read (nam(ix+1:), '(i4)', iostat = ios) charge
    if (ios /= 0) return
  end select
  nam = nam(1:ix-1)
endif
  
ix = index(nam, '-')
if (ix /= 0) then
  if (charge /= 0) return
  select case (nam(ix+1:ix+1))
  case ('-')
    do i = ix+1, len(nam)
      if (nam(i:i) == ' ') exit
      if (nam(i:i) /= '-') return
    enddo
    charge = -(i - ix)
  case (' ')
    charge = -1
  case default
    read (nam(ix+1:), '(i4)', iostat = ios) charge
    if (ios /= 0) return
    charge = -charge
  end select
  nam = nam(1:ix-1)
endif

if (abs(charge) > 127) then
  if (do_print) then
    call out_io (s_abort$, r_name, 'CHARGE > 127 not allowed: ' // name)
    if (global_com%exit_on_error) call err_exit
  endif
  return
endif

is_ok = .true.

end function get_this_charge

!----------------------------------------------------------------------------------------
! contains

! Strip off mass @M

function get_this_mass () result (is_ok)

logical is_ok

!

is_ok = .false.

ix = index(nam, '@M')
if (ix /= 0) then
  read (nam(ix+2:), *, iostat = ios) mol_mass
  if (ios /= 0) return 
  nam = nam(1:ix-1)
endif

is_ok = .true.

end function get_this_mass

end function species_id

!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!+
! Function species_name (species) result(name)
!
! Routine to return the name of a particle species given the integer index.
!
! Input:
!   species -- integer: Species ID.
!
! Output:
!   name    -- Character(20): Name of the species.
!               Will return 'INVALID!' (= invalid_name) if index is not valid.
!-

function species_name(species) result(name)

integer :: charge, species, ia, im, pp, mmmm
character(20) :: name
character(6)  :: extra

!

select case (species)
case (invalid$)
  name = invalid_name
  return
case (not_set$)
  name = "Not_Set!"
  return
end select

!

if (lb_subatomic <= species .and. species <= ub_subatomic) then
  name = subatomic_species_name(species)
  return
endif

pp = mod(abs(species), int(z'1000000')) / int(z'10000') 
if (pp == 0) then
  name = invalid_name
  return
endif

mmmm = mod(abs(species), int(z'10000'))
name = ''

! Atom?

if (pp < 200) then
  if (pp == anti_atom$) then
    pp = mmmm / 512
    mmmm = mmmm - pp * 512
    name = 'anti' // atomic_name(pp)
  else
    name = atomic_name(pp) 
  endif

  ! Isotope?
  if (mmmm > 0) then
    write(extra, '(i0)') mmmm
    name = '#' // trim(extra) // trim(name)
  endif

! molecule, possibly with specified mass
else
  if (mmmm > 0) then
    if (mod(mmmm, 100) == 0) then
      write (name, '(i0)') mmmm / 100
    elseif (mod(mmmm, 10) == 0) then
      write (name, '(f6.1)') mmmm / 100.0_rp
    else
      write (name, '(f6.2)') mmmm / 100.0_rp
    endif

    name = '@M' // trim(adjustl(name))
  endif

  if (pp > 200) then
    name = trim(molecular_name(pp-200)) // trim(name)
  endif
endif


! Add on charge
charge = charge_of(species)
select case (charge)
case (0)
  ! Nothing to do
case (1)
  name = trim(name) // '+'
case (-1)
  name = trim(name) // '-'
case (2)
  name = trim(name) // '++'
case (-2)
  name = trim(name) // '--'
case default
  if (charge > 0) then
    write (name, '(2a,i0)') trim(name), '+', charge
  else
    write (name, '(2a,i0)') trim(name), '-', abs(charge)
  endif
end select

end function species_name

!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!+
! Function species_id_from_openpmd (pmd_name, charge) result(species)
!
! Routine to return the Bmad species ID given the openPMD species name and given particle charge.
! Note: If pmd_name corresponds to a subatomic particle, the charge argument is ignored.
!
! Input:
!   pmd_name      -- character(*): OpenPMD species name.
!   charge        -- integer: Species charge. Ignored for subatomic particles.
!
! Output:
!   species       -- integer: Bmad spicies ID number.
!-

function species_id_from_openpmd (pmd_name, charge) result(species)

integer charge, species
integer i
character(*) pmd_name

! Subatomic particle

do i = lb_subatomic, ub_subatomic
  if (pmd_name /= openPMD_subatomic_species_name(i)) cycle
  species = i
  return
enddo

! All else. In this case pmd_name corresponds to the Bmad name of a neutral particle.

species = set_species_charge(species_id(pmd_name), charge)

end function species_id_from_openpmd

!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!+
! Function openpmd_species_name (species) result(pmd_name)
!
! Routine to return the openPMD name of a particle species given the Bmad species ID.
! Note: the pmd_name does not include the particle charge. For example, if species
! corresponds to He+ then the pmd_name will be "He".
!
! Input:
!   species   -- integer: Bmad species ID number.
!
! Output:
!   pmd_name  -- Character(20): Name of the species.
!                Will return 'INVALID!' (= invalid_name) if index is not valid.
!-

function openpmd_species_name(species) result(pmd_name)

integer :: species, ix
character(20) :: pmd_name

! Subatomic particles

if (lb_subatomic <= species .and. species <= ub_subatomic) then
  pmd_name = openpmd_subatomic_species_name(species)

! All else just remove any charge suffix. EG: "H-" -> "H".
else
  pmd_name = species_name(species)
  ix = index(pmd_name, '+')
  if (ix /= 0) pmd_name = pmd_name(1:ix-1)
  ix = index(pmd_name, '-')
  if (ix /= 0) pmd_name = pmd_name(1:ix-1)
endif

end function openpmd_species_name

!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!+
! Function anomalous_moment_of (species) result (moment)
!
! Routine to return the anomolous moment for subatomic species type. Otherwise returns 0.
!
! Input:
!   species -- integer: Species ID.
!
! Output:
!   moment  -- real(rp): Anomalous moment.
!-

function anomalous_moment_of (species) result (moment)

integer :: species, sp
integer, parameter :: He3$ = 2 * int(z'10000') + 3
real(rp) :: moment

!

if (lb_subatomic <= species .and. species <= ub_subatomic) then
  moment = anomalous_moment_of_subatomic(species)
  return
endif

sp = abs(species)
sp = sp - int(z'1000000') * (sp / int(z'1000000'))  ! Subtract off charge

if (sp == He3$) then
  moment = anomalous_mag_moment_He3
else
  moment = 0
endif

end function anomalous_moment_of

!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!+
! Function spin_of (species, non_subatomic_default) result (spin)
!
! Routine to return the spin, in units of hbar, of a particle.
! This routine is only valid for subatomic particles. 
! For all other particles, the returned spin value will be the value of non_subatomic_default.
!
! Input:
!   species               -- integer: Species ID.
!   non_subatomic_default -- real(rp), optional: Default value to be used for non-subatomic species.
!                             Default value of this argument is zero.
!
! Output:
!   spin    -- real(rp): Particle spin.
!-

function spin_of (species, non_subatomic_default) result (spin)

integer :: species
real(rp) spin
real(rp), optional :: non_subatomic_default

!

if (lb_subatomic <= species .and. species <= ub_subatomic) then
  spin = spin_of_subatomic(species)
else
  spin = real_option(0.0_rp, non_subatomic_default)
endif

end function spin_of

!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!+
! Function charge_of (species, default) result (charge)
!
! Routine to return the charge, in units of e+, of a particle.
!
! Input:
!   species -- integer: Species ID.
!   default -- integer, optional: If present then use default value if species = not_set$.
!
! Output:
!   charge -- integer: particle charge.
!-

function charge_of (species, default) result (charge)

integer :: charge, species
integer, optional :: default
character(*), parameter :: r_name = 'charge_of'


!

if (lb_subatomic <= species .and. species <= ub_subatomic) then
  charge = charge_of_subatomic(species)
  return
endif

! Invalid
if (abs(species) < 1000) then
  if (present(default) .and. species == not_set$) then
    charge = default
    return
  endif
  call out_io (s_abort$, r_name, 'CHARGE NOT KNOWN FOR SPECIES: \i0\ ', species)
  if (global_com%exit_on_error) call err_exit
  charge = int_garbage$
  return
endif

! |species| > 1000, decode CC PP MMMM
! Charge encoded in first two hex digits of species but if negative must avoid two's complement representation.
charge = sign_of(species) * abs(species) / int(z'1000000')  

end function charge_of

!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!+
! Function mass_of (species) result (mass)
!
! Routine to return the mass, in units of eV/c^2, of a particle.
! To convert to AMU divide mass_of value by the constant atomic_mass_unit.
!
! Note: For atoms where the isotopic number is given, the mass is calculated using the neutral atomic mass
! adjusted by the weight of any added or missing electrons. The calculated mass is off very slightly due to 
! binding energy effects. Exception: For #1H+ (proton) and #2H+ (deuteron) the exact mass is used since it is known.
!
! Input:
!   species -- integer: Species ID.
!
! Output:
!   mass -- real(rp): particle mass. Set to real_garbage$ if species value is invalid.
!-

function mass_of (species) result (mass)

real(rp) mass
integer n, ix, species, n_nuc, pp, charge, sp
character(*), parameter :: r_name = 'mass_of'
logical anti

! Subatomic particle

mass = real_garbage$

if (lb_subatomic <= species .and. species <= ub_subatomic) then
  mass = mass_of_subatomic(species)
  return
endif

! Invalid
if (abs(species) < 1000) then
  call out_io (s_abort$, r_name, 'MASS NOT KNOWN FOR SPECIES: \i0\ ', species)
  if (global_com%exit_on_error) call err_exit
  return
endif


! |species| > 1000, decode CC PP MMMM
sp = abs(species)
pp = mod(sp, int(z'1000000')) / int(z'10000')
charge = sign_of(species) * sp / int(z'1000000')  ! Charge encoded in first two hex digits of species.

! Atom?

if (pp < 200) then
  anti = (pp == anti_atom$) 
  if (anti) then
    pp = mod(sp, int(z'10000')) / 512
    n_nuc = mod(sp, int(z'10000')) - pp * 512
    charge = -charge
  else
    n_nuc = mod(sp, int(z'10000'))
  endif

  if (n_nuc == 0) then
    ! Average naturally occuring mass
    mass = atom(pp)%mass(0) * atomic_mass_unit - charge * m_electron
  else
    ! Isotope. In general the mass is calculated using the neutral atom plus/minus the weight of any added/missing electrons.
    ! This will be off very slightly due to binding energy effects. 

    ! For #1H+ (proton) and #2H+ (deuteron) use the exact mass.
    if (pp == 1 .and. n_nuc == 1 .and. charge == 1) then 
      mass = m_proton
      return
    elseif (pp == 1 .and. n_nuc == 2 .and. charge == 1) then 
      mass = m_deuteron
      return
    endif

    ix = n_nuc - atom(pp)%i_offset
    if (ix < 1) then ! Extrapolate
      mass = atom(pp)%mass(1) + (ix - 1) * (atom(pp)%mass(2) - atom(pp)%mass(1))
    elseif (ix > ubound(atom(pp)%mass, 1)) then ! Extrapolate
      do n = ubound(atom(pp)%mass, 1), 1, -1
        if (atom(pp)%mass(n) == no_iso) cycle
        mass = atom(pp)%mass(n) + (ix - n) * (atom(pp)%mass(n) - atom(pp)%mass(n-1))
        exit
      enddo      
    elseif (atom(pp)%mass(ix) == no_iso) then  ! Extrapolate
      do n = ix-1, 1, -1
        if (atom(pp)%mass(n) == no_iso) cycle
        mass = atom(pp)%mass(n) + (ix - n) * (atom(pp)%mass(n) - atom(pp)%mass(n-1))
        exit
      enddo      
    else
      mass = atom(pp)%mass(ix)
    endif
    mass = mass * atomic_mass_unit - charge * m_electron
  endif
  return
endif

! Molecule
if (pp == 200) then
  ! unknown, mass is specified directly in units of u/100
  mass = mod(sp, int(z'10000')) * atomic_mass_unit / 100
  return
else
  ! known molecule
  mass = molecular_mass(pp-200) * atomic_mass_unit - charge * m_electron
  return
endif

! ERROR
call out_io (s_abort$, r_name, 'ERROR: CANNOT DECODE MASS FOR SPECIES: \i0\ ', species)
if (global_com%exit_on_error) call err_exit

end function mass_of

!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!+
! Function charge_to_mass_of (species) result (charge_mass_ratio)
!
! Routine to return the charge (in units of e+) to mass (in units of eV) ratio of a particle.
!
! Input:
!   species -- integer: Species ID.
!
! Output:
!   charge_mass_ratio -- real(rp): particle charge to mass ratio. (1/eV)
!-

function charge_to_mass_of (species) result (charge_mass_ratio)

integer :: species
real(rp) charge_mass_ratio

character(*), parameter :: r_name = 'charge_to_mass_of'

!

charge_mass_ratio = charge_of(species) / mass_of(species)

end function charge_to_mass_of

!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!+
! Function set_species_charge(species_in, charge) result(species_charged)
!
! Routine to return the ID for a particle of the same type as species_in but with a different charge.
! Exception: If species_in corresponds to a subatomic particle, the charge argument is ignored and
! species_charged will be set equal to species_in.
!
! Input:
!   species_in      -- integer: Input species.
!   charge          -- integer: Charge to set species_charged to.
!
! Output:
!   species_charged -- integer: Species of the same type as species_in but with different charge.
!-

function set_species_charge(species_in, charge) result(species_charged)

integer species_in, charge, species_charged, sp
character(*), parameter :: r_name = 'set_species_charge'

!

if (lb_subatomic <= species_in .and. species_in <= ub_subatomic) then
  species_charged = species_in
  return
endif

if (charge < -127 .or. charge > 127) then
  call out_io (s_error$, r_name, 'CHARGE TO SET TO DOES NOT MAKE SENSE: ' // int_str(charge))
  return
endif

sp = abs(species_in)
species_charged = sign_of(charge) * (sp + int(z'1000000') * (abs(charge) - sp / int(z'1000000')))

end function set_species_charge

!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!+
! Function x0_radiation_length(species) result (x0)
!
! Routine to return the X0 raidation length for atomes.
!
! Input:
!   Species       -- integer: Species ID.
!
! Output:
!   x0            -- real(rp): Radiation length in kg/m^2. Set to real_garbage$ if species is not atomic
!                     or has atomic index greater than 92.
!-

function x0_radiation_length(species) result (x0)

integer species, ixa
real(rp) x0

!

ixa = atomic_number(species)
if (ixa < 1 .or. ixa > size(x0_rad_length)) then
  x0 =real_garbage$
  return
endif

x0 = x0_rad_length(ixa) * 10.0_rp

end function x0_radiation_length

!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!+
! Function atomic_number(species) result (atomic_num)
!
! Routine to return the atomic number Z if species argument corresponds to an atomic particle or is a proton.
! Set to zero otherwise.
!
! Input:
!   species       -- integer: Spicies ID.
!
! Output:
!   atomic_num    -- Integer: Atomic index. Set to zero if species is not atomic
!-

elemental function atomic_number(species) result (atomic_num)

integer, intent(in) :: species
integer atomic_num

!

if (species == proton$) then
  atomic_num = 1
  return
endif

atomic_num = mod(abs(species), int(z'1000000')) / int(z'10000')
if (atomic_num == anti_atom$) atomic_num = mod(abs(species), int(z'10000')) / 512
if (atomic_num > 199) atomic_num = 0

end function atomic_number

!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!+
! Function is_subatomic_species(species) result (is_subatomic)
!
! Routine to return True if species argument corresponds to a subatomic particle.
!
! Input:
!   species     -- integer: Spicies ID.
!
! Output:
!   is_subatomic  -- logical: Set True if species corresponds to a subatomic particle.
!-

elemental function is_subatomic_species(species) result (is_subatomic)

integer, intent(in) :: species
logical is_subatomic

!

is_subatomic = (lb_subatomic <= species .and. species <= ub_subatomic)

end function is_subatomic_species

end module



