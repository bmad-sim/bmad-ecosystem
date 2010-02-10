#ifndef BMAD_PARAMETERS

namespace Bmad {
  const int BMAD_INC_VERSION = 91;
  const int N_ATTRIB_MAXX = 60;
  const int HYPER_Y = 1, HYPER_XY = 2, HYPER_X = 3;
  const int DRIFT = 1, SBEND = 2, QUADRUPOLE = 3, GROUP = 4;
  const int SEXTUPOLE = 5, OVERLAY = 6, CUSTOM = 7, TAYLOR = 8;
  const int RFCAVITY = 9;
  const int ELSEPARATOR = 10, BEAMBEAM = 11, WIGGLER = 12;
  const int SOL_QUAD = 13, MARKER = 14, KICKER = 15;
  const int HYBRID = 16, OCTUPOLE = 17, RBEND = 18;
  const int MULTIPOLE = 19, ACCEL_SOL = 20;
  const int DEF_BEAM = 21, AB_MULTIPOLE = 22, SOLENOID = 23;
  const int PATCH = 24, LCAVITY = 25, DEF_PARAMETER = 26;
  const int NULL_ELE = 27, INIT_ELE = 28, HOM = 29;
  const int MATCH = 30, MONITOR = 31, INSTRUMENT = 32;
  const int HKICKER = 33, VKICKER = 34, RCOLLIMATOR = 35;
  const int ECOLLIMATOR = 36, GIRDER = 37, BEND_SOL_QUAD = 38;
  const int DEF_BEAM_START = 39, PHOTON_BRANCH = 40;
  const int BRANCH = 41, MIRROR = 42, CRYSTAL = 43;
  const int PIPE = 44;
  const int N_KEY = 44;
  const int PARTICLE = 1, N_PART   = 2, TAYLOR_ORDER = 3;
  const int LATTICE_TYPE = 5, SYMMETRY = 6;
  const int X_BEG_LIMIT=2, Y_BEG_LIMIT=3, B_X2=4,
          B_Y2=5, L_ST2=9, B_Z=10, L_ST1=11, S_ST2=12, S_ST1=13,
          B_X1=14, B_Y1=15;
  const int VAL1=3, VAL2=4, VAL3=5, VAL4=6, VAL5=7,
          VAL6=8, VAL7=9, VAL8=10, VAL9=11, VAL10=12, VAL11=13,
          VAL12=15;
  const int BETA_A0 = 2, ALPHA_A0 = 3, BETA_B0 = 4,
          ALPHA_B0 = 5, BETA_A1 = 6, ALPHA_A1 = 7, BETA_B1 = 8,
          ALPHA_B1 = 9, DPHI_A = 10, DPHI_B = 11,
          ETA_X0 = 12, ETAP_X0 = 13, ETA_Y0 = 14, ETAP_Y0 = 15,
          ETA_X1 = 16, ETAP_X1 = 17, ETA_Y1 = 18, ETAP_Y1 = 19,
          MATCH_END = 20,
          X0 = 21, PX0 = 22, Y0 = 23, PY0 = 24, Z0 = 25, PZ0 = 26,
          X1 = 34, PX1 = 35, Y1 = 36, PY1 = 37, Z1 = 38, PZ1 = 39,
          MATCH_END_ORBIT = 40;
  const int X = 1, PX = 2, Y = 3, PY = 4, Z = 5, PZ = 6;
  const int L=1;
  const int TILT=2, COMMAND=2, IX_BRANCH_TO=2;
  const int OLD_COMMAND=3, ANGLE=3, KICK=3, GRADIENT_ERR=3, X_GAIN_ERR=3;
  const int DIRECTION=3, GRAZE_ANGLE=3;
  const int K1=4, SIG_X=4, HARMON=4, H_DISPLACE=4, E_LOSS=4, Y_GAIN_ERR=4;
  const int       GRAZE_ANGLE_ERR = 4;
  const int K2=5, SIG_Y=5, B_MAX=5, V_DISPLACE=5, PHI0_ERR=5, CRUNCH=5;
  const int       CRITICAL_ANGLE = 5;
  const int K3=6, SIG_Z=6, RF_WAVELENGTH=6, G_ERR=6, NOISE=6;
  const int       DKS_DS=6, LRAD=6;
  const int G=7, KS=7, VOLTAGE=7, N_POLE=7, BBI_CONST=7, OSC_AMPLITUDE=7;
  const int       G_GRAZE = 7;
  const int E1=8, CHARGE=8, GAP=8, DPHI0=8, X_GAIN_CALIB=8, G_TRANS=8;
  const int N_SLICE=9, E2=9, RF_FREQUENCY=9, Y_GAIN_CALIB=9;
  const int FINT=10, POLARITY=10, GRADIENT=10, CRUNCH_CALIB=10;
  const int FINTX=11, Z_PATCH=11, PHI0=11, X_OFFSET_CALIB=11;
  const int RHO=12, S_CENTER=12, P0C_START=12, Y_OFFSET_CALIB=12;
  const int HGAP=13, E_TOT_START=13, TILT_CALIB=13;
  const int COEF=14, CURRENT=14, HGAPX=14, DELTA_E=14, L_POLE=14;
  const int       DE_ETA_MEAS=14;
  const int ROLL=15, QUAD_TILT=15, LR_FREQ_SPREAD=15, X_RAY_LINE_LEN=15;
  const int N_SAMPLE=15, DELTA_REF_TIME=15;
  const int L_ORIGINAL=16, L_CHORD=16, BEND_TILT=16;
  const int L_START=17, H1=17, X_QUAD=17;
  const int L_END=18, H2=18, Y_QUAD=18;
  const int X_PITCH=19;
  const int Y_PITCH=20;
  const int HKICK=21;
  const int VKICK=22;
  const int BL_HKICK=23;
  const int BL_VKICK=24;
  const int X_OFFSET=25;
  const int Y_OFFSET=26;
  const int S_OFFSET=27, Z_OFFSET=27;
  const int B_FIELD_ERR=28, BL_KICK = 28;
  const int RADIUS=29;
  const int N_REF_PASS=30;
  const int TILT_ERR=31;
  const int P0C=32;
  const int E_TOT=33;
  const int BS_FIELD=34;
  const int B_FIELD=35, E_FIELD=35;
  const int B_GRADIENT=36, E_GRADIENT=36;
  const int B1_GRADIENT=37, E1_GRADIENT=37;
  const int B2_GRADIENT=38, E2_GRADIENT=38, PATCH_END = 38;
  const int B3_GRADIENT=39, E3_GRADIENT=39, TRANSLATE_AFTER=39;
  const int TILT_TOT=40;
  const int X_PITCH_TOT=41;
  const int Y_PITCH_TOT=42;
  const int X_OFFSET_TOT=43;
  const int Y_OFFSET_TOT=44;
  const int S_OFFSET_TOT=45;
  const int COUPLER_STRENGTH = 46, PZ_OFFSET = 46, C_11 = 46;
  const int COUPLER_PHASE = 47, C_12 = 47;
  const int COUPLER_ANGLE = 48, C_21 = 48;
  const int C_22 = 49, POLE_RADIUS = 49, COUPLER_AT = 49;
  const int DS_STEP = 50, GAMMA_C = 50;
  const int GENERAL1 = 51;
  const int GENERAL2 = 52;
  const int GENERAL3 = 53;
  const int GENERAL4 = 54;
  const int GENERAL5 = 55;
  const int X1_LIMIT = 56;
  const int X2_LIMIT = 57;
  const int Y1_LIMIT = 58;
  const int Y2_LIMIT = 59;
  const int CHECK_SUM = 60;
  const int REF_ORBIT = 61, TERM = 61;
  const int X_POSITION = 62;
  const int SYMPLECTIFY = 63, Y_POSITION = 63;
  const int DESCRIP = 64, Z_POSITION = 64;
  const int IS_ON = 65, THETA_POSITION = 65;
  const int FIELD_CALC = 66, PHI_POSITION = 66;
  const int TYPE = 67, PSI_POSITION = 67;
  const int APERTURE_AT = 68, BETA_A = 68;
  const int RAN_SEED = 69, BETA_B = 69;
  const int SR_WAKE_FILE = 70, ALPHA_A = 70, REF_PATCH = 70;
  const int LR_WAKE_FILE = 71, ALPHA_B = 71;
  const int ALIAS =72, ETA_X = 72;
  const int START_EDGE =73, ETA_Y = 73;
  const int END_EDGE =74, ETAP_X = 74;
  const int ACCORDION_EDGE =75, ETAP_Y = 75;
  const int LATTICE = 76, PHI_A = 76;
  const int APERTURE_TYPE = 77, PHI_B = 77;
  const int MAP_WITH_OFFSETS = 78, CMAT_11 = 78;
  const int CSR_CALC_ON = 79, CMAT_12 = 79;
  const int SYMMETRIC_EDGE = 80, CMAT_21 = 80;
  const int MAT6_CALC_METHOD = 81, CMAT_22 = 81;
  const int TRACKING_METHOD  = 82, S_LONG = 82;
  const int NUM_STEPS = 83, REF_TIME = 83;
  const int INTEGRATOR_ORDER = 84;
  const int APERTURE = 85;
  const int X_LIMIT = 86;
  const int Y_LIMIT = 87;
  const int OFFSET_MOVES_APERTURE = 88;
  const int APERTURE_LIMIT_ON = 89;
  const int SUPERIMPOSE    = 90;
  const int OFFSET         = 91;
  const int REFERENCE      = 92;
  const int ELE_BEGINNING  = 93;
  const int ELE_CENTER     = 94;
  const int ELE_END        = 95;
  const int REF_BEGINNING  = 96;
  const int REF_CENTER     = 97;
  const int REF_END        = 98;
  const int COMMON_LORD    = 99;
  const int TO = 100;
  const int FIELD_MASTER = 101;
  const int A0  = 110, K0L  = 110;
  const int A20 = 130, K20L = 130;
  const int B0  = 140, T0  = 140;
  const int B20 = 160, T20 = 160;
  const int N_ATTRIB_SPECIAL_MAXX = T20;
  const int PROTON     = +2;
  const int POSITRON   = +1;
  const int PHOTON     = 0;
  const int ELECTRON   = -1;
  const int ANTIPROTON = -2;
  const int LINEAR_LATTICE = 10;
  const int CIRCULAR_LATTICE = 12;
  const int FREE = 1, SUPER_SLAVE = 2, OVERLAY_SLAVE = 3;
  const int GROUP_LORD = 4, SUPER_LORD = 5, OVERLAY_LORD = 6;
  const int GIRDER_LORD = 7, MULTIPASS_LORD = 8, MULTIPASS_SLAVE = 9;
  const int NOT_A_LORD = 10, GROUP_SLAVE = 11, PATCH_IN_SLAVE = 12;
  const int X_PLANE = 1, Y_PLANE = 2;
  const int Z_PLANE = 3, N_PLANE = 4, S_PLANE = 5;
  const int BMAD_STANDARD = 1, SYMP_LIE_PTC = 2;
  const int RUNGE_KUTTA = 3;
  const int LINEAR = 4, TRACKING = 5, SYMP_MAP = 6;
  const int SYMP_LIE_BMAD = 10, NO_METHOD = 11;
  const int BORIS = 12, ADAPTIVE_BORIS = 13, MAD = 14;
  const int MAP_TYPE = 1, PERIODIC_TYPE = 3;
  const int ENTRANCE_END = 1, EXIT_END = 2, BOTH_ENDS = 3;
  const int NO_END = 4;
  const int SINGLE_REF = 1, MATCH_AT_ENTRANCE = 2, MATCH_AT_EXIT = 3;
  const int MATCH_GLOBAL_COORDS = 4, PATCH_IN = 5, PATCH_OUT = 6;
  const int BENDS = 201;
  const int WIGGLERS = 202;
  const int ALL = 203;
  const int RADIANS = 1, DEGREES = 2, CYCLES = 3, KHZ = 4;
  const int KICK_FIELD = 1, EM_FIELD = 2;
  const int NOT_LOST = -1;
  const int IS_LOGICAL = 1, IS_INTEGER = 2, IS_REAL = 3, IS_NAME = 4;
  const int RECTANGULAR = 1, ELLIPTICAL = 2;
    const int N_POLE_MAXX = 20;
}

#define BMAD_PARAMETERS
#endif
