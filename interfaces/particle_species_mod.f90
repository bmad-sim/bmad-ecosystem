!+
! Module particle_species_mod
!
! This module defines the differnet types of particles that Bmad knows about along
! with masses, etc.
!
! IMPORTANT: DO USE HARD CODED ID NUMBERS IN YOUR CODE!! 
! For example, the association of positrons with ID = 1 is not assured.
! In general, use species_id(name) to get the species ID number or you can use positron$, 
! proton$, etc, for particles that have named paramters (see below for a list).
!
! Decoding of species integer:
! If |species| < 1000: Use elementary particle mapping (electron$ = -1, etc.)
! Else if |species| > 1000 mapping is:
! Write in hex: species = CCPPMMMM (Hex)
! Where:
!   CC   = Charge (2 Hex digits with range [-127, 127]). Set to 0 for fundamental particles.
!   PP   = Particle ID (2 Hex digits with range [0, 255]).
!           if PP = 0                --> Used for fundamental particles.
!           if 0 < PP < 200 (C8 Hex) --> Atom with PP = # Protons 
!           if PP = 200 (C8 Hex)     --> Molecule of unknown type.
!           if PP > 200 (C8 Hex)     --> "Named" molecule. See molecular_name array below for a list. 
!                                        In this case PP = Species ID. EG: nh2$ = 201, etc.
!   MMMM (4 Hex digits):
!          For fundamental particles (where CC = PP = 0): Particle integer ID. 
!          For atoms: Number of nucleons .
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

implicit none


integer, parameter :: invalid$ = -666
integer, parameter :: not_set$ = -999

character(*), parameter :: invalid_name = 'INVALID!'

!----------------------
! Fundamental particles.
! Note: It is convenient for debugging to define an "antiparticle" with reversed sign even though it does not exit in practice. 
! Example: anti_deuteron$.

integer, parameter :: pion_0$            = +8
integer, parameter :: ref_particle$      = +7
integer, parameter :: neutron            = +6
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
integer, parameter :: anti_neutron       = -6
integer, parameter :: anti_ref_particle$ = -7


character(20), parameter:: fundamental_species_name(-7:8) = [character(20):: 'Anti_Ref_Particle', 'Anti_Neutron     ', &
              'Anti_Deuteron    ', 'Pion-            ', 'Muon             ', 'Antiproton       ', 'Electron         ', &
              'Photon           ', 'Positron         ', 'Proton           ', 'Antimuon         ', 'Pion+            ', &
              'Deuteron         ', 'Neutron          ', 'Ref_Particle     ', 'Pion0            ']

character(20), parameter:: openPMD_fundamental_species_name(-7:8) = [character(20):: 'Garbage!         ', 'anti-neutron     ', &
                      'anti-deuteron    ', 'pion-            ', 'muon             ', 'anti-proton      ', 'electron         ', &
                      'photon           ', 'positron         ', 'proton           ', 'anti-muon        ', 'pion+            ', &
                      'deuteron         ', 'neutron          ', 'Garbage!         ', 'pion0            ']

integer, parameter :: charge_of_fundamental(-7:8) = [0, 0, -1, -1, -1, -1, -1, 0, 1, 1, 1, 1, 1, 0, 0, 0]

real(rp), parameter :: mass_of_fundamental(-7:8) = [0.0_rp, m_neutron, m_deuteron, m_pion_charged, m_muon, m_proton, m_electron, 0.0_rp, &
                                m_electron, m_proton, m_muon, m_pion_charged, m_deuteron, m_neutron, 0.0_rp, m_pion_0]

real(rp), parameter :: anomalous_moment_of_fundamental(-7:8) = [0.0_rp, anomalous_mag_moment_neutron, anomalous_mag_moment_deuteron, &
                        0.0_rp, anomalous_mag_moment_muon, anomalous_mag_moment_proton, anomalous_mag_moment_electron, 0.0_rp, &
                        anomalous_mag_moment_electron, anomalous_mag_moment_proton, anomalous_mag_moment_muon, &
                        0.0_rp, anomalous_mag_moment_deuteron, anomalous_mag_moment_neutron, 0.0_rp, 0.0_rp]

integer, parameter :: antiparticle(-7:8) = [7, 6, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5, -6, -7, 8]

!----------------------
! Atoms

character(3), parameter :: atomic_name(118) = [character(3) :: &
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
                    'Rg', 'Cn', 'Uut', 'Fl', 'Uup', 'Lv', 'Uus', 'Uuo']

real(rp), parameter, private :: no_iso = -1 ! mass if there is no known isotope

type atom_struct
  integer :: i_offset             ! isotope number offset 
  real(rp) :: mass(0:13) = no_iso ! isotope masses in units of the unified atomic mass unit.
                                  ! mass(0) is the standard atomic weight
                                  ! mass(n) is the mass of isotope n + i_offset
end type


! Atoms from NIST data.
! The first number in a row is the average mass for natually occuring isotope mixtures.

type(atom_struct), parameter, private :: atom(1:118) = [ &
! 1: <H>, H, D, T
atom_struct(0, [1.00784_rp, 1.00782503223_rp, 2.01410177812_rp, 3.0160492779_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 2: He
atom_struct(2, [4.002602_rp, 3.0160293201_rp, 4.00260325413_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 3: Li
atom_struct(5, [6.938_rp, 6.0151228874_rp, 7.0160034366_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 4: Be
atom_struct(8, [9.0121831_rp, 9.012183065_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 5: B
atom_struct(9, [10.806_rp, 10.01293695_rp, 11.00930536_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 6: C
atom_struct(11, [12.0096_rp, 12.0000000_rp, 13.00335483507_rp, 14.0032419884_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 7: N
atom_struct(13, [14.00643_rp, 14.00307400443_rp, 15.00010889888_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 8: O
atom_struct(15, [15.99903_rp, 15.99491461957_rp, 16.99913175650_rp, 17.99915961286_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 9: F
atom_struct(18, [18.998403163_rp, 18.99840316273_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 10: Ne
atom_struct(19, [20.1797_rp, 19.9924401762_rp, 20.993846685_rp, 21.991385114_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 11: Na
atom_struct(22, [22.98976928_rp, 22.9897692820_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 12: Mg
atom_struct(23, [24.304_rp, 23.985041697_rp, 24.985836976_rp, 25.982592968_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 13: Al
atom_struct(26, [26.9815385_rp, 26.98153853_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 14: Si
atom_struct(27, [28.084_rp, 27.97692653465_rp, 28.97649466490_rp, 29.973770136_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 15: P
atom_struct(30, [30.973761998_rp, 30.97376199842_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 16: S
atom_struct(31, [32.059_rp, 31.9720711744_rp, 32.9714589098_rp, 33.967867004_rp, no_iso, 35.96708071_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 17: Cl
atom_struct(34, [35.446_rp, 34.968852682_rp, no_iso, 36.965902602_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 18: Ar
atom_struct(35, [39.948_rp, 35.967545105_rp, no_iso, 37.96273211_rp, no_iso, 39.9623831237_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 19: K
atom_struct(38, [39.0983_rp, 38.9637064864_rp, 39.963998166_rp, 40.9618252579_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 20: Ca
atom_struct(39, [40.078_rp, 39.962590863_rp, no_iso, 41.95861783_rp, 42.95876644_rp, 43.95548156_rp, no_iso, 45.9536890_rp, no_iso, 47.95252276_rp, no_iso, no_iso, no_iso, no_iso]), &
! 21: Sc
atom_struct(44, [44.955908_rp, 44.95590828_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 22: Ti
atom_struct(45, [47.867_rp, 45.95262772_rp, 46.95175879_rp, 47.94794198_rp, 48.94786568_rp, 49.94478689_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 23: V
atom_struct(49, [50.9415_rp, 49.94715601_rp, 50.94395704_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 24: Cr
atom_struct(49, [51.9961_rp, 49.94604183_rp, no_iso, 51.94050623_rp, 52.94064815_rp, 53.93887916_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 25: Mn
atom_struct(54, [54.938044_rp, 54.93804391_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 26: Fe
atom_struct(53, [55.845_rp, 53.93960899_rp, no_iso, 55.93493633_rp, 56.93539284_rp, 57.93327443_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 27: Co
atom_struct(58, [58.933194_rp, 58.93319429_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 28: Ni
atom_struct(57, [58.6934_rp, 57.93534241_rp, no_iso, 59.93078588_rp, 60.93105557_rp, 61.92834537_rp, no_iso, 63.92796682_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 29: Cu
atom_struct(62, [63.546_rp, 62.92959772_rp, no_iso, 64.92778970_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 30: Zn
atom_struct(63, [65.38_rp, 63.92914201_rp, no_iso, 65.92603381_rp, 66.92712775_rp, 67.92484455_rp, no_iso, 69.9253192_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 31: Ga
atom_struct(68, [69.723_rp, 68.9255735_rp, no_iso, 70.92470258_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 32: Ge
atom_struct(69, [72.630_rp, 69.92424875_rp, no_iso, 71.922075826_rp, 72.923458956_rp, 73.921177761_rp, no_iso, 75.921402726_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 33: As
atom_struct(74, [74.921595_rp, 74.92159457_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 34: Se
atom_struct(73, [78.971_rp, 73.922475934_rp, no_iso, 75.919213704_rp, 76.919914154_rp, 77.91730928_rp, no_iso, 79.9165218_rp, no_iso, 81.9166995_rp, no_iso, no_iso, no_iso, no_iso]), &
! 35: Br
atom_struct(78, [79.901_rp, 78.9183376_rp, no_iso, 80.9162897_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 36: Kr
atom_struct(77, [83.798_rp, 77.92036494_rp, no_iso, 79.91637808_rp, no_iso, 81.91348273_rp, 82.91412716_rp, 83.9114977282_rp, no_iso, 85.9106106269_rp, no_iso, no_iso, no_iso, no_iso]), &
! 37: Rb
atom_struct(84, [85.4678_rp, 84.9117897379_rp, no_iso, 86.9091805310_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 38: Sr
atom_struct(83, [87.62_rp, 83.9134191_rp, no_iso, 85.9092606_rp, 86.9088775_rp, 87.9056125_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 39: Y
atom_struct(88, [88.90584_rp, 88.9058403_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 40: Zr
atom_struct(89, [91.224_rp, 89.9046977_rp, 90.9056396_rp, 91.9050347_rp, no_iso, 93.9063108_rp, no_iso, 95.9082714_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 41: Nb
atom_struct(92, [92.90637_rp, 92.9063730_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 42: Mo
atom_struct(91, [95.95_rp, 91.90680796_rp, no_iso, 93.90508490_rp, 94.90583877_rp, 95.90467612_rp, 96.90601812_rp, 97.90540482_rp, no_iso, 99.9074718_rp, no_iso, no_iso, no_iso, no_iso]), &
! 43: Tc
atom_struct(96, [98.0_rp, 96.9063667_rp, 97.9072124_rp, 98.9062508_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 44: Ru
atom_struct(95, [101.07_rp, 95.90759025_rp, no_iso, 97.9052868_rp, 98.9059341_rp, 99.9042143_rp, 100.9055769_rp, 101.9043441_rp, no_iso, 103.9054275_rp, no_iso, no_iso, no_iso, no_iso]), &
! 45: Rh
atom_struct(102, [102.90550_rp, 102.9054980_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 46: Pd
atom_struct(101, [106.42_rp, 101.9056022_rp, no_iso, 103.9040305_rp, 104.9050796_rp, 105.9034804_rp, no_iso, 107.9038916_rp, no_iso, 109.90517220_rp, no_iso, no_iso, no_iso, no_iso]), &
! 47: Ag
atom_struct(106, [107.8682_rp, 106.9050916_rp, no_iso, 108.9047553_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 48: Cd
atom_struct(105, [112.414_rp, 105.9064599_rp, no_iso, 107.9041834_rp, no_iso, 109.90300661_rp, 110.90418287_rp, 111.90276287_rp, 112.90440813_rp, 113.90336509_rp, no_iso, 115.90476315_rp, no_iso, no_iso]), &
! 49: In
atom_struct(112, [114.818_rp, 112.90406184_rp, no_iso, 114.903878776_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 50: Sn
atom_struct(111, [118.710_rp, 111.90482387_rp, no_iso, 113.9027827_rp, 114.903344699_rp, 115.90174280_rp, 116.90295398_rp, 117.90160657_rp, 118.90331117_rp, 119.90220163_rp, no_iso, 121.9034438_rp, no_iso, 123.9052766_rp]), &
! 51: Sb
atom_struct(120, [121.760_rp, 120.9038120_rp, no_iso, 122.9042132_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 52: Te
atom_struct(119, [127.60_rp, 119.9040593_rp, no_iso, 121.9030435_rp, 122.9042698_rp, 123.9028171_rp, 124.9044299_rp, 125.9033109_rp, no_iso, 127.90446128_rp, no_iso, 129.906222748_rp, no_iso, no_iso]), &
! 53: I
atom_struct(126, [126.90447_rp, 126.9044719_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 54: Xe
atom_struct(123, [131.293_rp, 123.9058920_rp, no_iso, 125.9042983_rp, no_iso, 127.9035310_rp, 128.9047808611_rp, 129.903509349_rp, 130.90508406_rp, 131.9041550856_rp, no_iso, 133.90539466_rp, no_iso, 135.907214484_rp]), &
! 55: Cs
atom_struct(132, [132.90545196_rp, 132.9054519610_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 56: Ba
atom_struct(129, [137.327_rp, 129.9063207_rp, no_iso, 131.9050611_rp, no_iso, 133.90450818_rp, 134.90568838_rp, 135.90457573_rp, 136.90582714_rp, 137.90524700_rp, no_iso, no_iso, no_iso, no_iso]), &
! 57: La
atom_struct(137, [138.90547_rp, 137.9071149_rp, 138.9063563_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 58: Ce
atom_struct(135, [140.116_rp, 135.90712921_rp, no_iso, 137.905991_rp, no_iso, 139.9054431_rp, no_iso, 141.9092504_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 59: Pr
atom_struct(140, [140.90766_rp, 140.9076576_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 60: Nd
atom_struct(141, [144.242_rp, 141.9077290_rp, 142.9098200_rp, 143.9100930_rp, 144.9125793_rp, 145.9131226_rp, no_iso, 147.9168993_rp, no_iso, 149.9209022_rp, no_iso, no_iso, no_iso, no_iso]), &
! 61: Pm
atom_struct(144, [145.0_rp, 144.9127559_rp, no_iso, 146.9151450_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 62: Sm
atom_struct(143, [150.36_rp, 143.9120065_rp, no_iso, no_iso, 146.9149044_rp, 147.9148292_rp, 148.9171921_rp, 149.9172829_rp, no_iso, 151.9197397_rp, no_iso, 153.9222169_rp, no_iso, no_iso]), &
! 63: Eu
atom_struct(150, [151.964_rp, 150.9198578_rp, no_iso, 152.9212380_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 64: Gd
atom_struct(151, [157.25_rp, 151.9197995_rp, no_iso, 153.9208741_rp, 154.9226305_rp, 155.9221312_rp, 156.9239686_rp, 157.9241123_rp, no_iso, 159.9270624_rp, no_iso, no_iso, no_iso, no_iso]), &
! 65: Tb
atom_struct(158, [158.92535_rp, 158.9253547_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 66: Dy
atom_struct(155, [162.500_rp, 155.9242847_rp, no_iso, 157.9244159_rp, no_iso, 159.9252046_rp, 160.9269405_rp, 161.9268056_rp, 162.9287383_rp, 163.9291819_rp, no_iso, no_iso, no_iso, no_iso]), &
! 67: Ho
atom_struct(164, [164.93033_rp, 164.9303288_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 68: Er
atom_struct(161, [167.259_rp, 161.9287884_rp, no_iso, 163.9292088_rp, no_iso, 165.9302995_rp, 166.9320546_rp, 167.9323767_rp, no_iso, 169.9354702_rp, no_iso, no_iso, no_iso, no_iso]), &
! 69: Tm
atom_struct(168, [168.93422_rp, 168.9342179_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 70: Yb
atom_struct(167, [173.054_rp, 167.9338896_rp, no_iso, 169.9347664_rp, 170.9363302_rp, 171.9363859_rp, 172.9382151_rp, 173.9388664_rp, no_iso, 175.9425764_rp, no_iso, no_iso, no_iso, no_iso]), &
! 71: Lu
atom_struct(174, [174.9668_rp, 174.9407752_rp, 175.9426897_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 72: Hf
atom_struct(173, [178.49_rp, 173.9400461_rp, no_iso, 175.9414076_rp, 176.9432277_rp, 177.9437058_rp, 178.9458232_rp, 179.9465570_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 73: Ta
atom_struct(179, [180.94788_rp, 179.9474648_rp, 180.9479958_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 74: W
atom_struct(179, [183.84_rp, 179.9467108_rp, no_iso, 181.94820394_rp, 182.95022275_rp, 183.95093092_rp, no_iso, 185.9543628_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 75: Re
atom_struct(184, [186.207_rp, 184.9529545_rp, no_iso, 186.9557501_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 76: Os
atom_struct(183, [190.23_rp, 183.9524885_rp, no_iso, 185.9538350_rp, 186.9557474_rp, 187.9558352_rp, 188.9581442_rp, 189.9584437_rp, no_iso, 191.9614770_rp, no_iso, no_iso, no_iso, no_iso]), &
! 77: Ir
atom_struct(190, [192.217_rp, 190.9605893_rp, no_iso, 192.9629216_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 78: Pt
atom_struct(189, [195.084_rp, 189.9599297_rp, no_iso, 191.9610387_rp, no_iso, 193.9626809_rp, 194.9647917_rp, 195.96495209_rp, no_iso, 197.9678949_rp, no_iso, no_iso, no_iso, no_iso]), &
! 79: Au
atom_struct(196, [196.966569_rp, 196.96656879_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 80: Hg
atom_struct(195, [200.592_rp, 195.9658326_rp, no_iso, 197.96676860_rp, 198.96828064_rp, 199.96832659_rp, 200.97030284_rp, 201.97064340_rp, no_iso, 203.97349398_rp, no_iso, no_iso, no_iso, no_iso]), &
! 81: Tl
atom_struct(202, [204.382_rp, 202.9723446_rp, no_iso, 204.9744278_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 82: Pb
atom_struct(203, [207.2_rp, 203.9730440_rp, no_iso, 205.9744657_rp, 206.9758973_rp, 207.9766525_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 83: Bi
atom_struct(208, [208.98040_rp, 208.9803991_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 84: Po
atom_struct(208, [209.0_rp, 208.9824308_rp, 209.9828741_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 85: At
atom_struct(209, [210.0_rp, 209.9871479_rp, 210.9874966_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 86: Rn
atom_struct(210, [222.0_rp, 210.9906011_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, 220.0113941_rp, no_iso, 222.0175782_rp, no_iso]), &
! 87: Fr
atom_struct(222, [223.0_rp, 223.0197360_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 88: Ra
atom_struct(222, [226.0_rp, 223.0185023_rp, 224.0202120_rp, no_iso, 226.0254103_rp, no_iso, 228.0310707_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 89: Ac
atom_struct(226, [227.0_rp, 227.0277523_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 90: Th
atom_struct(229, [232.0377_rp, 230.0331341_rp, no_iso, 232.0380558_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 91: Pa
atom_struct(230, [231.03588_rp, 231.0358842_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 92: U
atom_struct(232, [238.02891_rp, 233.0396355_rp, 234.0409523_rp, 235.0439301_rp, 236.0455682_rp, no_iso, 238.0507884_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 93: Np
atom_struct(235, [237.0_rp, 236.046570_rp, 237.0481736_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 94: Pu
atom_struct(237, [244.0_rp, 238.0495601_rp, 239.0521636_rp, 240.0538138_rp, 241.0568517_rp, 242.0587428_rp, no_iso, 244.0642053_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 95: Am
atom_struct(240, [241.0568293_rp, 241.0568293_rp, no_iso, 243.0613813_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 96: Cm
atom_struct(242, [243.0613893_rp, 243.0613893_rp, 244.0627528_rp, 245.0654915_rp, 246.0672238_rp, 247.0703541_rp, 248.0723499_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 97: Bk
atom_struct(246, [247.0703073_rp, 247.0703073_rp, no_iso, 249.0749877_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 98: Cf
atom_struct(248, [249.0748539_rp, 249.0748539_rp, 250.0764062_rp, 251.0795886_rp, 252.0816272_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 99: Es
atom_struct(251, [252.082980_rp, 252.082980_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 100: Fm
atom_struct(256, [257.0951061_rp, 257.0951061_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 101: Md
atom_struct(257, [258.0984315_rp, 258.0984315_rp, no_iso, 260.10365_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 102: No
atom_struct(258, [259.10103_rp, 259.10103_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 103: Lr
atom_struct(261, [262.10961_rp, 262.10961_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 104: Rf
atom_struct(266, [267.12179_rp, 267.12179_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 105: Db
atom_struct(267, [268.12567_rp, 268.12567_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 106: Sg
atom_struct(270, [271.13393_rp, 271.13393_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 107: Bh
atom_struct(271, [272.13826_rp, 272.13826_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 108: Hs
atom_struct(269, [270.13429_rp, 270.13429_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 109: Mt
atom_struct(275, [276.15159_rp, 276.15159_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 110: Ds
atom_struct(280, [281.16451_rp, 281.16451_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 111: Rg
atom_struct(279, [280.16514_rp, 280.16514_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 112: Cn
atom_struct(284, [285.17712_rp, 285.17712_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 113: Uut
atom_struct(283, [284.17873_rp, 284.17873_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 114: Fl
atom_struct(288, [289.19042_rp, 289.19042_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 115: Uup
atom_struct(287, [288.19274_rp, 288.19274_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 116: Lv
atom_struct(292, [293.20449_rp, 293.20449_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 117: Uus
atom_struct(291, [292.20746_rp, 292.20746_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]), &
! 118: Uuo
atom_struct(293, [294.21392_rp, 294.21392_rp, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso, no_iso]) ]


!----------------------
! Molecules

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
! Function species_of (mass, charge) result(species)
!
! Routine to return the integer ID index of a particle species given the mass and charge.
! Note: Currently this routine only works for fundamental particles.
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

do species = lbound(mass_of_fundamental, 1), ubound(mass_of_fundamental, 1)
  if (charge == charge_of_fundamental(species) .and. abs(mass - mass_of_fundamental(species)) <= 1d-6 * mass) return
enddo

species = invalid$

end function species_of

!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!+
! Function species_id (name) result(species)
!
! Routine to return the integer ID index of a particle species given the name.
!
! For fundamental particles, the case does not matter. 
! For all other types of particles, the case does matter.
!
! Input:
!   name    -- Character(20): Name of the species.
!
! Output:
!   species -- integer: Species ID. Will return invalid$ if name is not valid.
!                       Will return not_set$ if name is blank
!-

function species_id (name) result(species)

integer ::  species, charge, i, ix, ix1, ix2, iso, ios, n_nuc
real(rp) :: mol_mass
character(*) :: name
character(20) :: nam
character(40), parameter :: r_name = 'species_id'

! Init

if (name == '') then
  species = not_set$
  return
endif

species = invalid$
iso = 0
nam = name
mol_mass = 0
n_nuc = 0

if (upcase(nam) == 'PION_PLUS')  nam = 'pion+'  ! Old style
if (upcase(nam) == 'PION_0')     nam = 'pion0'  ! Old style
if (upcase(nam) == 'PION_MINUS') nam = 'pion-'  ! Old style

! Fundamental particle?

call match_word(nam, fundamental_species_name, ix, .false., .false.)
if (ix > 0) then
  species = ix + lbound(fundamental_species_name, 1) - 1
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
  species = abs(charge * z'1000000') + 200 * z'10000' + nint(100 * mol_mass)
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
    call out_io (s_abort$, r_name, 'SPECIFIED MOLECULE MASS TOO LARGE FOR ' // name)
    if (global_com%exit_on_error) call err_exit
    return
  endif

  if (n_nuc /= 0) then
    call out_io (s_abort$, r_name, 'SPECIFYING THE NUMBER OF NUCLEONS FOR A MOLECULE IS INVALID: ' // name)
    if (global_com%exit_on_error) call err_exit
    return
  endif  

  species = abs(charge * z'1000000') + (ix+200) * z'10000' + nint(mol_mass * 100)
  if (charge < 0) species = -species
  return  
endif


! Is Atom?

call match_word(nam, atomic_name, ix, .true., .false.)
if (ix > 0) then
  if (mol_mass /= 0) then
    call out_io (s_abort$, r_name, 'SPECIFYING A MASS FOR AN ATOM IS INVALID: ' // name)
    if (global_com%exit_on_error) call err_exit
    return
  endif  

  species =  abs(charge * z'1000000') + ix * z'10000' + n_nuc
  if (charge < 0) species = -species
  return  
endif

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
  call out_io (s_abort$, r_name, 'CHARGE > 127 not allowed: ' // name)
  if (global_com%exit_on_error) call err_exit
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
  name = "Invalid!"
  return
case (not_set$)
  name = "Not_Set!"
  return
end select

!

if (is_fundamental_species(species)) then
  name = fundamental_species_name(species)
  return
endif

pp = mod(abs(species), z'1000000') / z'10000' 
mmmm = mod(abs(species), z'10000')
name = ''

! Atom
if (pp < 200) then
  name = atomic_name(pp) 
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
! Note: If pmd_name corresponds to a fundamental particle, the charge argument is ignored.
!
! Input:
!   pmd_name      -- character(*): OpenPMD species name.
!   charge        -- integer: Species charge. Ignored for fundamental particles.
!
! Output:
!   species       -- integer: Bmad spicies ID number.
!-

function species_id_from_openpmd (pmd_name, charge) result(species)

integer charge, species
integer i
character(*) pmd_name

! Fundamental particle

do i = lbound(fundamental_species_name, 1), ubound(fundamental_species_name, 1) 
  if (pmd_name /= openPMD_fundamental_species_name(i)) cycle
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

! Fundamental particles

if (is_fundamental_species(species)) then
  pmd_name = openpmd_fundamental_species_name(species)

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
! Routine to return the anomolous moment for fundamental species type. Otherwise returns 0.
!
! Input:
!   species -- integer: Species ID.
!
! Output:
!   moment  -- real(rp) 
!-

function anomalous_moment_of (species) result (moment)
integer :: species, sp
integer, parameter :: He3$ = 2 * z'10000' + 3
real(rp) :: moment

!

if (is_fundamental_species(species)) then
	moment = anomalous_moment_of_fundamental(species)
	return
endif

sp = abs(species - z'1000000' * (species / z'1000000'))  ! Subtract off charge

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
! Function charge_of (species) result (charge)
!
! Routine to return the charge, in units of e+, of a particle.
!
! Input:
!   species -- integer: Species ID.
!
! Output:
!   charge -- integer: particle charge.
!-

function charge_of (species) result (charge)
integer :: charge, species
character(*), parameter :: r_name = 'charge_of'
!

if (is_fundamental_species(species)) then
	charge = charge_of_fundamental(species)
	return
endif

! Invalid
if (abs(species) < 1000) then
  call out_io (s_abort$, r_name, 'CHARGE NOT KNOWN FOR SPECIES: \i0\ ', species)
  if (global_com%exit_on_error) call err_exit
  charge = int_garbage$
  return
endif

! |species| > 1000, decode CC PP MMMM
charge = species / z'1000000'  ! Charge encoded in first two hex digits of species.

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
integer species, n_nuc, pp, charge
character(*), parameter :: r_name = 'mass_of'

! Fundamental particle

mass = real_garbage$

if (is_fundamental_species(species)) then
	mass = mass_of_fundamental(species)
	return
endif

! Invalid
if (abs(species) < 1000) then
  call out_io (s_abort$, r_name, 'MASS NOT KNOWN FOR SPECIES: \i0\ ', species)
  if (global_com%exit_on_error) call err_exit
  return
endif


! |species| > 1000, decode CC PP MMMM
pp = mod(abs(species), z'1000000') / z'10000'

! Atom
if (pp<200) then
  n_nuc = mod(abs(species), z'10000')
  if (n_nuc == 0) then
    ! Average naturally occuring mass
    mass = atom(pp)%mass(0) * atomic_mass_unit
  else
    ! Isotope. In general the mass is calculated using the neutral atom plus/minus the weight of any added/missing electrons.
    ! This will be off very slightly due to binding energy effects. 
    charge = species / z'1000000'  ! Charge encoded in first two hex digits of species.

    ! For #1H+ (proton) and #2H+ (deuteron) use the exact mass.
    if (pp == 1 .and. n_nuc == 1 .and. charge == 1) then 
      mass = m_proton
      return
    elseif (pp == 1 .and. n_nuc == 2 .and. charge == 1) then 
      mass = m_deuteron
      return
    endif

    mass = atom(pp)%mass(n_nuc - atom(pp)%i_offset)  
    if (mass == no_iso) then
      call out_io (s_abort$, r_name, 'ISOTOPE MASS NOT KNOWN FOR SPECIES: \i0\ ', species)
      if (global_com%exit_on_error) call err_exit
      return
    endif
    mass = mass * atomic_mass_unit - charge * m_electron
  endif
  return
endif

! Molecule
if (pp == 200) then
  ! unknown, mass is specified directly in units of u/100
  mass = mod(abs(species), z'10000') * atomic_mass_unit / 100
  return
else
  ! known molecule
  mass = molecular_mass(pp-200) * atomic_mass_unit
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
! Routine to return the charge (in units of e+) to mass (units of eV/c^2) ration of a particle.
!
! Input:
!   species -- integer: Species ID.
!
! Output:
!   charge_mass_ratio -- real(rp): particle charge to mass ratio.
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
! Exception: If species_in corresponds to a fundamental particle, the charge argument is ignored and
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

integer species_in, charge, species_charged
character(*), parameter :: r_name = 'set_species_charge'

!

if (is_fundamental_species(species_in)) then
  species_charged = species_in
  return
endif

species_charged = species_in +  z'1000000' * (charge - species_in / z'1000000')

end function set_species_charge

!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!+
! Function is_fundamental_species(species) result (is_fundamental)
!
! Routine to return True if species argument corresponds to a fundamental particle.
!
! Input:
!   species     -- integer: Spicies ID.
!
! Output:
!   is_fundamental  -- logical: Set True if species corresponds to a fundamental particle.
!-

elemental function is_fundamental_species(species) result (is_fundamental)

integer, intent(in) :: species
logical is_fundamental

!

is_fundamental = (lbound(fundamental_species_name, 1) <= species .and. &
                                species <= ubound(fundamental_species_name, 1)) 

end function is_fundamental_species

end module



