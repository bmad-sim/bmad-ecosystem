!+
! Subroutine bend_exact_multipole_field (ele, param, orbit, local_ref_frame, field, calc_dfield, calc_potential)
!
! Routine to calculate the electric and magnetic field at a given point in a bend element due to any multipoles.
! The field due to a multipole in a bend is different from a straight element since Maxwell's equations
! are modified due to the curvature of the reference orbit.
!
! Note: The returned field does not include the bend field (g) itself or the bend field error (dg) but does
! include contributions from h/vkick, k1, and k2.
!
! Input:
!   ele               -- ele_stuct: Bend element.
!   param             -- lat_param_struct: Lattice branch parameters.
!   orbit             -- coord_struct: particle position.
!   local_ref_frame   -- logical: Is the particle position in the local element ref 
!                         frame (as opposed to the lab frame)?
!   calc_dfield       -- logical, optional: If present and True then calculate the field derivatives.
!   calc_potential    -- logical, optional: Calc electric and magnetic potentials? Default is false. 
!
! Output:
!   field             -- em_field_struct: Field
!-

subroutine bend_exact_multipole_field (ele, param, orbit, local_ref_frame, field, calc_dfield, calc_potential)

use bmad_interface, except_dummy => bend_exact_multipole_field

implicit none

type bend_exact_coef_struct
	integer :: order = 0
	real(rp) :: cutoff_minus = 0    ! Crossover between exact formula and pade approximant
	real(rp) :: cutoff_plus = 0     ! Crossover between exact formula and pade approximant
  ! Exact formula coefs
	integer :: n_exact_non = 0
	real(rp) :: exact_non_coef(0:11) = 0  ! Non-zero coefs are in the range [0:n_exact_non]
	integer :: n_exact_log = 0
	real(rp) :: exact_log_coef(0:10) = 0  ! Non-zero coefs are in the range [0:n_exact_log]
  ! pade approximant
	integer :: n_pade_numer = 0
	real(rp) :: pade_numer_coef(0:7) = 0  ! Non-zero coefs are in the range [0:n_pade_numer]
	integer :: n_pade_denom = 0
	real(rp) :: pade_denom_coef(0:8) = 0  ! Non-zero coefs are in the range [0:n_pade_denom]
  ! 2nd derivative pade approximant.
	integer :: n_d2_pade_numer = 0
	real(rp) :: d2_pade_numer_coef(0:7) = 0  ! Non-zero coefs are in the range [0:n_d2_pade_numer]
	integer :: n_d2_pade_denom = 0
	real(rp) :: d2_pade_denom_coef(0:8) = 0  ! Non-zero coefs are in the range [0:n_d2_pade_denom]
end type

type (bend_exact_coef_struct), parameter :: F_coef(0:22) = [&
  bend_exact_coef_struct(), &
  bend_exact_coef_struct(1, -0.01, 0.01, &
    -1, [0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    0, [1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    3, [1.0_rp, 0.83333333333333_rp, 0.066666666666667_rp, -0.0055555555555556_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    2, [1.0_rp, 1.3333333333333_rp, 0.4_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    0, [-1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    2, [1.0_rp, 2.0_rp, 1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp]), &
  bend_exact_coef_struct(2, -0.0224, 0.0316, &
    1, [-0.5_rp, 0.5_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    0, [-1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    4, [1.0_rp, 1.1666666666667_rp, 0.28571428571429_rp, -0.0035714285714286_rp, 0.0005952380952381_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    2, [1.0_rp, 1.5_rp, 0.53571428571429_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    2, [2.0_rp, 2.0_rp, 1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    2, [1.0_rp, 2.0_rp, 1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp]), &
  bend_exact_coef_struct(3, -0.0794, 0.1, &
    1, [1.5_rp, -1.5_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    1, [1.5_rp, 1.5_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    4, [1.0_rp, 1.3834196891192_rp, 0.48255613126079_rp, 0.022366148531952_rp, -0.0010270170244264_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    3, [1.0_rp, 1.8834196891192_rp, 1.0742659758204_rp, 0.17530224525043_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    5, [6.0_rp, 8.8009313154831_rp, 3.7459254947614_rp, 0.21162481290537_rp, -0.020676035256943_rp, 0.0020784827318588_rp, 0.0_rp, 0.0_rp], &
    3, [1.0_rp, 2.4668218859139_rp, 1.9244761350407_rp, 0.45678807029214_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp]), &
  bend_exact_coef_struct(4, -0.15, 0.178, &
    2, [-1.875_rp, 1.5_rp, 0.375_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    1, [-1.5_rp, -3.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    4, [1.0_rp, 1.9565979090138_rp, 1.267864573968_rp, 0.30992013543582_rp, 0.023174683703788_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    4, [1.0_rp, 2.3565979090138_rp, 1.9105037375735_rp, 0.60999940061822_rp, 0.060982814868093_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    4, [12.0_rp, 25.856179775281_rp, 19.732584269663_rp, 6.0797752808989_rp, 0.62005350454789_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    4, [1.0_rp, 2.8213483146067_rp, 2.7752808988764_rp, 1.0908239700375_rp, 0.13723916532905_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp]), &
  bend_exact_coef_struct(5, -0.219, 0.251, &
    2, [2.8125_rp, 0.0_rp, -2.8125_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    2, [1.875_rp, 7.5_rp, 1.875_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    4, [1.0_rp, 2.0_rp, 1.3133394383394_rp, 0.31333943833944_rp, 0.019614644614645_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    5, [1.0_rp, 2.5_rp, 2.2061965811966_rp, 0.80929487179487_rp, 0.10954901579902_rp, 0.0032253626003626_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    4, [20.0_rp, 42.871328274455_rp, 30.882576910024_rp, 8.1715064226414_rp, 0.56499719380436_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    5, [1.0_rp, 2.8935664137227_rp, 2.9643036557932_rp, 1.2616282526849_rp, 0.19823929548004_rp, 0.0068953548215405_rp, 0.0_rp, 0.0_rp, 0.0_rp]), &
  bend_exact_coef_struct(6, -0.282, 0.316, &
    3, [-3.125_rp, -2.8125_rp, 5.625_rp, 0.3125_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    2, [-1.875_rp, -11.25_rp, -5.625_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    5, [1.0_rp, 2.4187964980921_rp, 2.1040382083956_rp, 0.79642117021259_rp, 0.1255086337457_rp, 0.0062963159736057_rp, 0.0_rp, 0.0_rp], &
    5, [1.0_rp, 2.8473679266635_rp, 3.0029101769656_rp, 1.4300620315322_rp, 0.29569616640861_rp, 0.019880608705292_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    5, [30.0_rp, 78.027333646035_rp, 74.54248410026_rp, 31.68672408437_rp, 5.7744847541202_rp, 0.34895663426502_rp, 0.0_rp, 0.0_rp], &
    5, [1.0_rp, 3.2009111215345_rp, 3.8052961429294_rp, 2.047426577554_rp, 0.4801195051624_rp, 0.036942712949507_rp, 0.0_rp, 0.0_rp, 0.0_rp]), &
  bend_exact_coef_struct(7, -0.338, 0.373, &
    3, [4.0104166666667_rp, 9.84375_rp, -9.84375_rp, -4.0104166666667_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    3, [2.1875_rp, 19.6875_rp, 19.6875_rp, 2.1875_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    5, [1.0_rp, 2.5_rp, 2.2514276625317_rp, 0.87714149379762_rp, 0.13753081472819_rp, 0.0059084917311593_rp, 0.0_rp, 0.0_rp], &
    6, [1.0_rp, 3.0_rp, 3.3903165514206_rp, 1.7806331028413_rp, 0.43046032864087_rp, 0.040143777220235_rp, 0.00082665505743365_rp, 0.0_rp, 0.0_rp], &
    5, [42.0_rp, 111.0687306549_rp, 106.71844005378_rp, 44.594530783095_rp, 7.5149065501092_rp, 0.34550740757978_rp, 0.0_rp, 0.0_rp], &
    6, [1.0_rp, 3.311160253688_rp, 4.1293077894533_rp, 2.3898805307595_rp, 0.63561544948173_rp, 0.065272073243703_rp, 0.0014988941223356_rp, 0.0_rp, 0.0_rp]), &
  bend_exact_coef_struct(8, -0.387, 0.422, &
    4, [-4.2838541666667_rp, -17.5_rp, 9.84375_rp, 11.666666666667_rp, 0.2734375_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    3, [-2.1875_rp, -26.25_rp, -39.375_rp, -8.75_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    5, [1.0_rp, 2.5610637208739_rp, 2.3919099978286_rp, 0.98937170028964_rp, 0.17357822943969_rp, 0.0096995892929036_rp, 0.0_rp, 0.0_rp], &
    6, [1.0_rp, 3.0055081653184_rp, 3.3943580713034_rp, 1.7688665051567_rp, 0.41438661668697_rp, 0.033650020420609_rp, 0.000050495737648433_rp, 0.0_rp, 0.0_rp], &
    5, [56.0_rp, 154.27768697566_rp, 156.50816672089_rp, 71.038574426579_rp, 13.876393477537_rp, 0.88729180719484_rp, 0.0_rp, 0.0_rp], &
    6, [1.0_rp, 3.3263872674226_rp, 4.1598671299717_rp, 2.3993340100535_rp, 0.62173843909719_rp, 0.055825305974345_rp, 0.000067884096308701_rp, 0.0_rp, 0.0_rp]), &
  bend_exact_coef_struct(9, -0.43, 0.464, &
    4, [5.126953125_rp, 32.8125_rp, 0.0_rp, -32.8125_rp, -5.126953125_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    4, [2.4609375_rp, 39.375_rp, 88.59375_rp, 39.375_rp, 2.4609375_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    6, [1.0_rp, 2.733043364401_rp, 2.7718291902232_rp, 1.2772433471491_rp, 0.25997188487676_rp, 0.018326394926351_rp, 0.00013278743887644_rp, 0.0_rp], &
    6, [1.0_rp, 3.233043364401_rp, 4.0247145087873_rp, 2.4094030144879_rp, 0.70373582351232_rp, 0.088893885031519_rp, 0.0033132610003509_rp, 0.0_rp, 0.0_rp], &
    6, [72.0_rp, 208.47164020037_rp, 224.93003875114_rp, 110.54311991155_rp, 24.017935388763_rp, 1.8035342017279_rp, 0.013500562516323_rp, 0.0_rp], &
    6, [1.0_rp, 3.5204394472274_rp, 4.7687474149497_rp, 3.1016552179887_rp, 0.98236904351452_rp, 0.13439444126767_rp, 0.005443115713299_rp, 0.0_rp, 0.0_rp]), &
  bend_exact_coef_struct(10, -0.468, 0.501, &
    5, [-5.373046875_rp, -47.16796875_rp, -24.609375_rp, 57.421875_rp, 19.482421875_rp, 0.24609375_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    4, [-2.4609375_rp, -49.21875_rp, -147.65625_rp, -98.4375_rp, -12.3046875_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    6, [1.0_rp, 2.8752528324977_rp, 3.1372522022764_rp, 1.6195427416905_rp, 0.40170620422812_rp, 0.04303954725897_rp, 0.0014669363563495_rp, 0.0_rp], &
    6, [1.0_rp, 3.3297982870431_rp, 4.309887787296_rp, 2.7231444177806_rp, 0.86124668679305_rp, 0.12371395040216_rp, 0.0058577600405354_rp, 0.0_rp, 0.0_rp], &
    6, [90.0_rp, 276.26611885817_rp, 323.47593022248_rp, 180.17720474476_rp, 48.57986931215_rp, 5.7293537920329_rp, 0.21975002984642_rp, 0.0_rp], &
    6, [1.0_rp, 3.6251790984241_rp, 5.1081653904854_rp, 3.5120975386257_rp, 1.2084104117485_rp, 0.18902561595337_rp, 0.0097883838077107_rp, 0.0_rp, 0.0_rp]), &
  bend_exact_coef_struct(11, -0.501, 0.534, &
    5, [6.1810546875_rp, 73.3154296875_rp, 90.234375_rp, -90.234375_rp, -73.3154296875_rp, -6.1810546875_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    5, [2.70703125_rp, 67.67578125_rp, 270.703125_rp, 270.703125_rp, 67.67578125_rp, 2.70703125_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    6, [1.0_rp, 2.7333931468221_rp, 2.7730902668328_rp, 1.2786959487213_rp, 0.26066077963333_rp, 0.018454839938474_rp, 0.00013857880170088_rp, 0.0_rp], &
    6, [1.0_rp, 3.2333931468221_rp, 4.0244022248592_rp, 2.407541872889_rp, 0.70200770652817_rp, 0.088320424888959_rp, 0.0032543708103201_rp, 0.0_rp, 0.0_rp], &
    6, [110.0_rp, 323.77181590469_rp, 354.24422343341_rp, 176.22933901445_rp, 38.732075353657_rp, 2.9491593740496_rp, 0.02328200154313_rp, 0.0_rp], &
    6, [1.0_rp, 3.5433801445881_rp, 4.8282482997839_rp, 3.1561005323515_rp, 1.0031035616685_rp, 0.13726458794807_rp, 0.0055066226115528_rp, 0.0_rp, 0.0_rp]), &
  bend_exact_coef_struct(12, -0.531, 0.562, &
    6, [-6.406640625_rp, -96.099609375_rp, -186.1083984375_rp, 90.234375_rp, 169.189453125_rp, 28.965234375_rp, 0.2255859375_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    5, [-2.70703125_rp, -81.2109375_rp, -406.0546875_rp, -541.40625_rp, -203.02734375_rp, -16.2421875_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    6, [1.0_rp, 3.0448349228732_rp, 3.5652325648319_rp, 2.0107624226303_rp, 0.5581463905672_rp, 0.068912599775289_rp, 0.0027503351125044_rp, 0.0_rp], &
    7, [1.0_rp, 3.5063733844117_rp, 4.8374048960989_rp, 3.3142815877642_rp, 1.1661040072182_rp, 0.19412177384366_rp, 0.011475769024855_rp, 0.000014175308484462_rp, 0.0_rp], &
    6, [132.0_rp, 427.16040102496_rp, 533.03535391816_rp, 321.29547888758_rp, 95.686069229407_rp, 12.763674635952_rp, 0.5589152377982_rp, 0.0_rp], &
    7, [1.0_rp, 3.781518189583_rp, 5.6235201785465_rp, 4.1494130043063_rp, 1.5707144983039_rp, 0.28109618522875_rp, 0.017854714867322_rp, 0.000020024640468448_rp, 0.0_rp]), &
  bend_exact_coef_struct(13, -0.557, 0.588, &
    6, [7.184912109375_rp, 135.4869140625_rp, 384.90600585937_rp, 0.0_rp, -384.90600585937_rp, -135.4869140625_rp, -7.184912109375_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    6, [2.9326171875_rp, 105.57421875_rp, 659.8388671875_rp, 1173.046875_rp, 659.8388671875_rp, 105.57421875_rp, 2.9326171875_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    6, [1.0_rp, 3.0_rp, 3.4398700386372_rp, 1.8797400772743_rp, 0.49529660107961_rp, 0.055426562442451_rp, 0.0017537166619788_rp, 0.0_rp], &
    7, [1.0_rp, 3.5_rp, 4.8232033719705_rp, 3.3080084299262_rp, 1.1727046776928_rp, 0.20104858661291_rp, 0.013546152700783_rp, 0.00019859291244789_rp, 0.0_rp], &
    6, [156.0_rp, 499.61786495051_rp, 612.01027032877_rp, 357.32351842723_rp, 100.55321903531_rp, 12.004833475776_rp, 0.40405958757359_rp, 0.0_rp], &
    7, [1.0_rp, 3.7860119548109_rp, 5.6380599885677_rp, 4.1722054927337_rp, 1.5928600032099_rp, 0.2935580177983_rp, 0.021248055653563_rp, 0.0003365447018836_rp, 0.0_rp]), &
  bend_exact_coef_struct(14, -0.518, 0.611, &
    7, [-7.394384765625_rp, -168.3322265625_rp, -631.24584960937_rp, -256.60400390625_rp, 641.51000976562_rp, 381.8267578125_rp, 40.030224609375_rp, 0.20947265625_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    6, [-2.9326171875_rp, -123.169921875_rp, -923.7744140625_rp, -2052.83203125_rp, -1539.6240234375_rp, -369.509765625_rp, -20.5283203125_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    6, [1.0_rp, 3.034629705342_rp, 3.5369070660493_rp, 1.9813638048949_rp, 0.54414752125926_rp, 0.065935628026009_rp, 0.0025331166641352_rp, 0.0_rp], &
    7, [1.0_rp, 3.5012963720087_rp, 4.8208453729867_rp, 3.2938732095366_rp, 1.154415151606_rp, 0.19109330453762_rp, 0.011210062927602_rp, 0.000016581712694799_rp, 0.0_rp], &
    6, [182.0_rp, 593.4520238077_rp, 743.7913255258_rp, 448.4350991714_rp, 132.7488657447_rp, 17.404307385932_rp, 0.73106807871057_rp, 0.0_rp], &
    7, [1.0_rp, 3.7991869439983_rp, 5.6709431100969_rp, 4.1948122014563_rp, 1.5890796519162_rp, 0.28387953303457_rp, 0.017944734057527_rp, 0.000025465677304369_rp, 0.0_rp]), &
  bend_exact_coef_struct(15, -0.541, 0.631, &
    7, [8.1469900948661_rp, 223.24548339844_rp, 1085.4349365234_rp, 962.26501464844_rp, -962.26501464844_rp, -1085.4349365234_rp, -223.24548339844_rp, -8.1469900948661_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    7, [3.14208984375_rp, 153.96240234375_rp, 1385.6616210938_rp, 3849.0600585937_rp, 3849.0600585937_rp, 1385.6616210938_rp, 153.96240234375_rp, 3.14208984375_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    6, [1.0_rp, 3.0_rp, 3.4400585392107_rp, 1.8801170784215_rp, 0.49556813038953_rp, 0.055509591178792_rp, 0.0017626354878635_rp, 0.0_rp], &
    7, [1.0_rp, 3.5_rp, 4.8224114803872_rp, 3.306028700968_rp, 1.1709157214433_rp, 0.20034488119687_rp, 0.013433848489043_rp, 0.00019373407731112_rp, 0.0_rp], &
    6, [210.0_rp, 681.6635152695_rp, 844.90272570991_rp, 498.45355598581_rp, 141.58804772541_rp, 17.057362707399_rp, 0.58074667960623_rp, 0.0_rp], &
    7, [1.0_rp, 3.8174453108071_rp, 5.7285531571751_rp, 4.2685019226957_rp, 1.6392252031112_rp, 0.30337174133444_rp, 0.021966095923949_rp, 0.0003432421486395_rp, 0.0_rp]), &
  bend_exact_coef_struct(16, -0.562, 0.649, &
    8, [-8.3433707101004_rp, -267.7060546875_rp, -1601.208984375_rp, -2155.4736328125_rp, 962.26501464844_rp, 2278.6435546875_rp, 739.01953125_rp, 52.607561383929_rp, 0.19638061523438_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    7, [-3.14208984375_rp, -175.95703125_rp, -1847.548828125_rp, -6158.49609375_rp, -7698.1201171875_rp, -3695.09765625_rp, -615.849609375_rp, -25.13671875_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    6, [1.0_rp, 3.0267474784445_rp, 3.5150700003594_rp, 1.9587569664169_rp, 0.53342251019288_rp, 0.063667813336172_rp, 0.0023691677567435_rp, 0.0_rp], &
    7, [1.0_rp, 3.4973357137386_rp, 4.8079338656481_rp, 3.2779819701702_rp, 1.1453369246616_rp, 0.18875391331486_rp, 0.011008376720253_rp, 0.000018907401960289_rp, 0.0_rp], &
    6, [240.0_rp, 792.01643279814_rp, 1001.4014310909_rp, 606.72552157455_rp, 179.50984240347_rp, 23.297394195417_rp, 0.94806295689106_rp, 0.0_rp], &
    7, [1.0_rp, 3.8334018033256_rp, 5.7669869246522_rp, 4.2934244372316_rp, 1.6340898906471_rp, 0.29261955200841_rp, 0.018495995257021_rp, 0.00003187163381539_rp, 0.0_rp]), &
  bend_exact_coef_struct(17, -0.582, 0.666, &
    8, [9.0734857831682_rp, 340.33321707589_rp, 2486.4927978516_rp, 4711.2495117188_rp, 0.0_rp, -4711.2495117188_rp, -2486.4927978516_rp, -340.33321707589_rp, -9.0734857831682_rp, 0.0_rp, 0.0_rp, 0.0_rp], &
    8, [3.3384704589844_rp, 213.662109375_rp, 2617.3608398437_rp, 10469.443359375_rp, 16358.505249023_rp, 10469.443359375_rp, 2617.3608398437_rp, 213.662109375_rp, 3.3384704589844_rp, 0.0_rp, 0.0_rp], &
    6, [1.0_rp, 3.2157192672142_rp, 4.0246058913212_rp, 2.4656382544766_rp, 0.7600290519396_rp, 0.10732582666566_rp, 0.0050545105774104_rp, 0.0_rp], &
    8, [1.0_rp, 3.7157192672142_rp, 5.5140444722967_rp, 4.1563428658619_rp, 1.6699262094031_rp, 0.34103513472824_rp, 0.02989625911893_rp, 0.00066188011699055_rp, -3.8321173171098e-6_rp], &
    6, [272.0_rp, 937.11399143158_rp, 1255.5205301028_rp, 822.55656066273_rp, 270.82539830172_rp, 40.79710536667_rp, 2.0467627828159_rp, 0.0_rp], &
    8, [1.0_rp, 4.007772027322_rp, 6.4070207730997_rp, 5.1943068511192_rp, 2.2405392284895_rp, 0.49034590748123_rp, 0.046007827720805_rp, 0.0010935402932571_rp, -6.6496478399894e-6_rp]), &
  bend_exact_coef_struct(18, -0.599, 0.681, &
    9, [-9.2589563642229_rp, -397.89798627581_rp, -3437.2891845703_rp, -8375.5546875_rp, -2944.5309448242_rp, 7655.780456543_rp, 6150.7979736328_rp, 1291.1296037946_rp, 66.638254983085_rp, 0.18547058105469_rp, 0.0_rp, 0.0_rp], &
    8, [-3.3384704589844_rp, -240.36987304687_rp, -3365.1782226562_rp, -15704.165039063_rp, -29445.309448242_rp, -23556.247558594_rp, -7852.0825195312_rp, -961.4794921875_rp, -30.046234130859_rp, 0.0_rp, 0.0_rp], &
    7, [1.0_rp, 3.5306836929127_rp, 4.9787705805848_rp, 3.5748288749893_rp, 1.3787221678681_rp, 0.27591307847442_rp, 0.025151905086668_rp, 0.00074062043812306_rp], &
    8, [1.0_rp, 4.004367903439_rp, 6.5203132716875_rp, 5.5340270154316_rp, 2.6041300484626_rp, 0.66504710948496_rp, 0.082751975852322_rp, 0.0037091269405145_rp, 4.1162181059434e-6_rp], &
    7, [306.0_rp, 1140.6790109412_rp, 1698.7240636962_rp, 1288.4739697103_rp, 525.2265240594_rp, 111.23014189022_rp, 10.762624082882_rp, 0.33914776274957_rp], &
    8, [1.0_rp, 4.2571209507881_rp, 7.3639792736075_rp, 6.6336116941671_rp, 3.3097896353312_rp, 0.89535506565443_rp, 0.11791802055689_rp, 0.005590233444132_rp, 6.0115785877696e-6_rp]), &
  bend_exact_coef_struct(19, -0.616, 0.695, &
    9, [9.9691173311264_rp, 490.34381021772_rp, 4991.1087210519_rp, 15333.372253418_rp, 11189.217590332_rp, -11189.217590332_rp, -15333.372253418_rp, -4991.1087210519_rp, -490.34381021772_rp, -9.9691173311264_rp, 0.0_rp, 0.0_rp], &
    9, [3.5239410400391_rp, 285.43922424316_rp, 4567.0275878906_rp, 24864.927978516_rp, 55946.08795166_rp, 55946.08795166_rp, 24864.927978516_rp, 4567.0275878906_rp, 285.43922424316_rp, 3.5239410400391_rp, 0.0_rp], &
    7, [1.0_rp, 3.5_rp, 4.8775459718386_rp, 3.4438649295966_rp, 1.294825988365_rp, 0.24837405295099_rp, 0.020877112471247_rp, 0.00050504506367791_rp], &
    8, [1.0_rp, 4.0_rp, 6.508498352791_rp, 5.525495058373_rp, 2.6074910509159_rp, 0.67249033787685_rp, 0.086419935430204_rp, 0.0044239428872841_rp, 0.000049085206240032_rp], &
    7, [342.0_rp, 1268.294295087_rp, 1871.7680196343_rp, 1398.6279551451_rp, 556.06487759462_rp, 112.68765968972_rp, 9.9953901402249_rp, 0.25464964602582_rp], &
    8, [1.0_rp, 4.2640184066872_rp, 7.3886886459206_rp, 6.6719472071155_rp, 3.3442426142721_rp, 0.91480097920062_rp, 0.1245219370701_rp, 0.006747813323692_rp, 0.000079518432011526_rp]), &
  bend_exact_coef_struct(20, -0.631, 0.708, &
    10, [-10.145314383128_rp, -562.44616099766_rp, -6595.6849316188_rp, -24442.055053711_rp, -26418.985977173_rp, 11189.217590332_rp, 30563.140640259_rp, 14099.791521345_rp, 2094.9200207847_rp, 82.071468111068_rp, 0.17619705200195_rp, 0.0_rp], &
    9, [-3.5239410400391_rp, -317.15469360352_rp, -5708.7844848633_rp, -35521.325683594_rp, -93243.479919434_rp, -111892.17590332_rp, -62162.319946289_rp, -15223.425292969_rp, -1427.1961212158_rp, -35.239410400391_rp, 0.0_rp], &
    7, [1.0_rp, 3.5251396690978_rp, 4.9605987099841_rp, 3.5514756331536_rp, 1.3638765739909_rp, 0.27108740107008_rp, 0.024413032466054_rp, 0.00070074454442973_rp], &
    8, [1.0_rp, 4.0013301452883_rp, 6.5088511601214_rp, 5.5169131419853_rp, 2.5913408457078_rp, 0.66008942597404_rp, 0.081834003934574_rp, 0.0036501401808332_rp, 4.5632350175799e-6_rp], &
    7, [380.0_rp, 1424.9935396178_rp, 2131.8507602905_rp, 1621.6929490686_rp, 661.47149290292_rp, 139.6665019865_rp, 13.381952622433_rp, 0.41096212893076_rp], &
    8, [1.0_rp, 4.276298788468_rp, 7.426606626274_rp, 6.7123965487154_rp, 3.3575778944078_rp, 0.90958736078044_rp, 0.11977948728242_rp, 0.0056681330420259_rp, 7.0484394263371e-6_rp]), &
  bend_exact_coef_struct(21, -0.645, 0.72, &
    10, [10.837587006887_rp, 676.74351056417_rp, 9125.1352000237_rp, 40468.938903809_rp, 59831.232948303_rp, 0.0_rp, -59831.232948303_rp, -40468.938903809_rp, -9125.1352000237_rp, -676.74351056417_rp, -10.837587006887_rp, 0.0_rp], &
    10, [3.700138092041_rp, 370.0138092041_rp, 7492.7796363831_rp, 53281.988525391_rp, 163176.08985901_rp, 234973.56939697_rp, 163176.08985901_rp, 53281.988525391_rp, 7492.7796363831_rp, 370.0138092041_rp, 3.700138092041_rp], &
    7, [1.0_rp, 3.5_rp, 4.877620852413_rp, 3.4440521310325_rp, 1.2950082042224_rp, 0.24846017530121_rp, 0.020896613402104_rp, 0.00050668185193194_rp], &
    8, [1.0_rp, 4.0_rp, 6.5080556350217_rp, 5.524166905065_rp, 2.6059624003339_rp, 0.67164662555942_rp, 0.0861932030478_rp, 0.004397707778928_rp, 0.000048223327017921_rp], &
    7, [420.0_rp, 1571.6018012631_rp, 2338.2364313731_rp, 1760.0102394784_rp, 704.41667268967_rp, 143.63481807921_rp, 12.817978876798_rp, 0.32906711217672_rp], &
    8, [1.0_rp, 4.2919090506264_rp, 7.4825414811137_rp, 6.7949769741632_rp, 3.4234287574474_rp, 0.94064349656167_rp, 0.12846989936403_rp, 0.0069674870175204_rp, 0.000081373008986209_rp]), &
  bend_exact_coef_struct(22, -0.658, 0.731, &
    11, [-11.00577510198_rp, -764.76862112681_rp, -11661.712009907_rp, -60223.711881638_rp, -112336.19247437_rp, -35898.739768982_rp, 95729.972717285_rp, 92973.898429871_rp, 28879.908177853_rp, 3213.3203204473_rp, 98.862697569529_rp, 0.16818809509277_rp], &
    10, [-3.700138092041_rp, -407.01519012451_rp, -9157.8417778015_rp, -73262.734222412_rp, -256419.56977844_rp, -430784.87722778_rp, -358987.39768982_rp, -146525.46844482_rp, -27473.525333405_rp, -2035.0759506226_rp, -40.701519012451_rp], &
    7, [1.0_rp, 3.5205293526005_rp, 4.9455011102732_rp, 3.5320978036904_rp, 1.3515804842329_rp, 0.26710119359406_rp, 0.023805221453277_rp, 0.00066816326407054_rp], &
    8, [1.0_rp, 3.9987902221657_rp, 6.4992703469612_rp, 5.5026175638515_rp, 2.5806701431507_rp, 0.65596147985212_rp, 0.081072889433702_rp, 0.0036018948868731_rp, 4.9980796385363e-6_rp], &
    7, [462.0_rp, 1747.6654079823_rp, 2633.3075427931_rp, 2013.9436889295_rp, 824.05725763474_rp, 173.96415989744_rp, 16.562835537541_rp, 0.49811550114185_rp], &
    8, [1.0_rp, 4.3066350822128_rp, 7.5270849192393_rp, 6.8413366592761_rp, 3.4381213820185_rp, 0.93471971417383_rp, 0.12334485091487_rp, 0.005840413754959_rp, 8.2586867428032e-6_rp]) &
  ]

type pole_coef_struct
  real(rp) :: value = 0
  real(rp) :: derivative = 0
  real(rp) :: derivative2 = 0
end type pole_coef_struct

type (pole_coef_struct) pole_coef(0:n_pole_maxx+1), pole_elec_coef(0:n_pole_maxx+1)

type (ele_struct), target :: ele
type (lat_param_struct) param
type (coord_struct) orbit
type (em_field_struct) field

real(rp) b0, f, x, y, g, rho, rho_n, yg, xg, fact_n
real(rp) a_pole(0:n_pole_maxx), b_pole(0:n_pole_maxx), a_pole_elec(0:n_pole_maxx), b_pole_elec(0:n_pole_maxx)
real(rp) yg_n(0:n_pole_maxx+1)

integer n, i, j, sgn, ix_mag_max, ix_elec_max

logical local_ref_frame
logical, optional :: calc_dfield, calc_potential
logical do_dfield_calc

character(*), parameter :: r_name = 'bend_exact_multipole_field'

!

g = ele%value(g$)
rho = ele%value(rho$)

if (ele%value(g$) == 0) then
  call out_io (s_fatal$, r_name, 'this routine called with g = 0!')
  if (global_com%exit_on_error) call err_exit
  return
endif

x = orbit%vec(1)
y = orbit%vec(3)
yg = y * g
xg = x * g

field = em_field_struct()
do_dfield_calc = logic_option(.false., calc_dfield)

!

call multipole_ele_to_ab (ele, .not. local_ref_frame, ix_mag_max, a_pole, b_pole, magnetic$, include_kicks$)
call multipole_ele_to_ab (ele, .not. local_ref_frame, ix_elec_max, a_pole_elec, b_pole_elec, electric$)

if (nint(ele%value(exact_multipoles$)) == horizontally_pure$) then
  if (ix_mag_max /= -1) then
    call convert_bend_exact_multipole(g, vertically_pure$, a_pole, b_pole)
    ix_mag_max = n_pole_maxx
  endif
  if (ix_elec_max /= -1) then
    ! Notice that a_pole and b_pole are reversed for electric fields.
    call convert_bend_exact_multipole(g, vertically_pure$, b_pole_elec, a_pole_elec)
    ix_elec_max = n_pole_maxx
  endif
endif

! Calculate y-dependent coefs

yg_n(0) = 1
do n = 0, max(ix_mag_max, ix_elec_max)
  yg_n(n+1) = yg_n(n) * yg
enddo

pole_coef = pole_coef_struct()
pole_elec_coef = pole_coef_struct()

do n = 0, max(ix_mag_max, ix_elec_max)
  fact_n = factorial(n)

  if (n == 0) then
    rho_n = 1
  else
    rho_n = rho * rho_n
  endif

  if (a_pole(n) /= 0) then
    sgn = 1
    do j = n+1, 0, -2
      call add_to_this_coef (a_pole(n) * fact_n, pole_coef(j))
    enddo
  endif

  if (b_pole(n) /= 0) then
    sgn = 1
    do j = n, 0, -2
      call add_to_this_coef (b_pole(n) * fact_n, pole_coef(j))
    enddo
  endif

  if (a_pole_elec(n) /= 0) then
    sgn = 1
    do j = n, 0, -2
      call add_to_this_coef (a_pole_elec(n) * fact_n, pole_elec_coef(j))
    enddo
  endif

  if (b_pole_elec(n) /= 0) then
    sgn = 1
    do j = n+1, 0, -2
      call add_to_this_coef (b_pole_elec(n) * fact_n, pole_elec_coef(j))
    enddo
  endif
enddo

! Combine y-dependent coefs with F factors to get the field.

do j = 0, max(ix_mag_max, ix_elec_max)+1
  call add_this_field (pole_coef(j), field%B, field%dB, magnetic$)
  call add_this_field (pole_elec_coef(j), field%E, field%dE, electric$)
enddo

if (ix_mag_max >= 0) then
  f = ele%value(p0c$) / (c_light * charge_of(param%particle) * ele%value(l$))
  field%B = field%B * f
  if (do_dfield_calc) field%dB = field%dB * f
  if (logic_option(.false., calc_potential)) then
    field%phi_B = field%phi_B * f
    field%a = field%a * f
  endif
endif

!------------------------------------------------------------------------
contains

subroutine add_to_this_coef (coef, pole_coef)

type (pole_coef_struct) pole_coef
real(rp) coef, Z(3), dZ(3,3), crs

!

crs = coef * rho_n * sgn 

              pole_coef%value = pole_coef%value             + crs * yg_n(n+1-j) / (factorial(j) * factorial(n+1-j))
if (n-j > -1) pole_coef%derivative = pole_coef%derivative   + crs * yg_n(n-j) / (factorial(j) * factorial(n-j))
if (n-j>0)    pole_coef%derivative2 = pole_coef%derivative2 + crs * yg_n(n-1-j) / (factorial(j) * factorial(n-1-j))

sgn = -sgn

end subroutine add_to_this_coef

!------------------------------------------------------------------------
! contains

subroutine add_this_field (pole_coef, Z, dZ, field_type)

type (pole_coef_struct) pole_coef
real(rp) Z(3), dZ(3,3), crs, value, derivative
integer field_type

!

if (pole_coef%value /= 0) then
      Z(1) = Z(1) + F_derivative(F_coef(j), xg) * pole_coef%value
endif

if (pole_coef%derivative /= 0 .and. n-j > -1) then
     Z(2) = Z(2) + F_value(F_coef(j), xg) * pole_coef%derivative
endif

if (do_dfield_calc) then
  if (pole_coef%value /= 0) then
     dZ(1,1) = dZ(1,1) + g * F_derivative2(F_coef(j), xg) * pole_coef%value
  endif
  if (pole_coef%derivative /= 0 .and. n-j > -1) then
     dZ(1,2) = dZ(1,2) + g * F_derivative(F_coef(j), xg) * pole_coef%derivative
     dZ(2,1) = dZ(1,2)
  endif
  if (pole_coef%derivative2 /= 0 .and. n-j > 0) then
     dZ(2,2) = dZ(2,2) + g * F_value(F_coef(j), xg) * pole_coef%derivative2
  endif
endif

if (logic_option(.false., calc_potential) .and. pole_coef%value /= 0 .and. n-j > -1) then
  if (field_type == magnetic$) then
    field%phi_B = field%phi_B - rho * F_value(F_coef(j), xg) * pole_coef%value
  else
    field%phi = field%phi - rho * F_value(F_coef(j), xg) * pole_coef%value
  endif
endif

end subroutine add_this_field

!------------------------------------------------------------------------
! contains

function F_value(coef, xg) result (value)

type (bend_exact_coef_struct) coef
real(rp) xg, value
real(rp) r, numer, denom
integer i

! Order 0 case

if (coef%order == 0) then
  value = 1
  return
endif

! Use exact log formula?

if (xg < coef%cutoff_minus .or. xg > coef%cutoff_plus) then
  r = 1 + xg

  value = 0
  do i = 0, coef%n_exact_log
    value = value + coef%exact_log_coef(i) * r**(2*i)
  enddo
  value = value * log(r)

  do i = 0, coef%n_exact_non
    value = value + coef%exact_non_coef(i) * r**(2*i)
  enddo  

! Else use Taylor expansion.

else
  numer = 0
  do i = 0, coef%n_pade_numer
    numer = numer + coef%pade_numer_coef(i) * xg**(i+coef%order)
  enddo

  denom = 0
  do i = 0, coef%n_pade_denom
    denom = denom + coef%pade_denom_coef(i) * xg**i
  enddo

  value = numer / denom

endif

end function F_value

!------------------------------------------------------------------------
! contains


function F_derivative(coef, xg) result (value)

type (bend_exact_coef_struct) coef
real(rp) xg, value, v0, f
real(rp) r, numer, denom, d_numer, d_denom
integer i

! Order 0 case

if (coef%order == 0) then
  value = 0
  return
endif

! Use exact log formula?

if (xg < coef%cutoff_minus .or. xg > coef%cutoff_plus) then
  r = 1 + xg

  value = 0
  v0 = 0
  do i = 0, coef%n_exact_log
    f = coef%exact_log_coef(i) * r**(2*i-1)
    v0    = v0 + f
    value = value + 2 * i * f
  enddo
  value = value * log(r) + v0

  do i = 1, coef%n_exact_non
    value = value + 2 * i * coef%exact_non_coef(i) * r**(2*i-1)
  enddo

! Else use Taylor expansion.

else
  numer = 0
  d_numer = 0
  do i = 0, coef%n_pade_numer
    numer = numer + coef%pade_numer_coef(i) * xg**(i+coef%order)
    if (i+coef%order < 1) cycle
    d_numer = d_numer + (i+coef%order) * coef%pade_numer_coef(i) * xg**(i-1+coef%order)
  enddo

  denom = 0
  d_denom = 0
  do i = 0, coef%n_pade_denom
    denom = denom + coef%pade_denom_coef(i) * xg**i
    if (i == 0) cycle
    d_denom = d_denom + i * coef%pade_denom_coef(i) * xg**(i-1)
  enddo

  value = (denom * d_numer - numer * d_denom) / denom**2

endif

end function F_derivative

!------------------------------------------------------------------------
! contains

function F_derivative2(coef, xg) result (value)

type (bend_exact_coef_struct) coef
real(rp) xg, value
real(rp) r, numer, denom
integer i, min_o

! Order 0 case

if (coef%order == 0) then
  value = 0
  return
endif

! Use exact log formula?

if (xg < coef%cutoff_minus .or. xg > coef%cutoff_plus) then
  r = 1 + xg

  value = 0
  do i = 1, coef%n_exact_log
    value = value + 2 * i * (2*i-1) * coef%exact_log_coef(i) * r**(2*i-2)
  enddo
  value = value * log(r)

  do i = 0, coef%n_exact_log
    value = value + (4*i-1) * coef%exact_log_coef(i) * r**(2*i-2)
  enddo

  do i = 1, coef%n_exact_non
    value = value + 2 * i * (2*i-1) * coef%exact_non_coef(i) * r**(2*i-2)
  enddo  

! Else use Taylor expansion.

else
  min_o = max(0, coef%order-2)
  numer = 0
  do i = 0, coef%n_d2_pade_numer
    numer = numer + coef%d2_pade_numer_coef(i) * xg**(i+min_o)
  enddo

  denom = 0
  do i = 0, coef%n_d2_pade_denom
    denom = denom + coef%d2_pade_denom_coef(i) * xg**i
  enddo

  value = numer / denom

endif

end function F_derivative2

end subroutine bend_exact_multipole_field

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!+
! Subroutine convert_bend_exact_multipole (g, out_type, an, bn)
!
! Routine to convert exact bend multipole coefficients between horizontally pure and vertically pure bases.
! Note: Scaling by r0 multipole radius not handled here.
! 
! Input:
!   g                 -- real(rp): 1/rho bending strength.
!   out_type          -- integer: Output type: horizontally_pure$ or vertically_pure$.
!   an(0:n_pole_maxx) -- real(rp): Skew multipoles. 
!   bn(0:n_pole_maxx) -- real(rp): Non-skew multipoles.
!
! Output:
!   an(0:n_pole_maxx) -- real(rp): Converted skew multipoles.
!   bn(0:n_pole_maxx) -- real(rp): Converted Non-skew multipoles.
!-

subroutine convert_bend_exact_multipole (g, out_type, an, bn)

use bmad_interface, except_dummy => convert_bend_exact_multipole

implicit none

! Horizontal to Vertical Imaginary Potential

type pure_bend_multipole_struct
  real(rp) convert(0:n_pole_maxx)
end type

type (pure_bend_multipole_struct), parameter :: v_to_h_imag(0:n_pole_maxx) = [ &
     pure_bend_multipole_struct([1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
        0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp]), &
     pure_bend_multipole_struct([0.0_rp, 1.0_rp, -0.5_rp, 0.33333333333333_rp, -0.25_rp, 0.2_rp, -0.16666666666667_rp, 0.14285714285714_rp, -0.125_rp, 0.11111111111111_rp, -0.1_rp, 0.090909090909091_rp, &
        -0.083333333333333_rp, 0.076923076923077_rp, -0.071428571428571_rp, 0.066666666666667_rp, -0.0625_rp, 0.058823529411765_rp, -0.055555555555556_rp, 0.052631578947368_rp, -0.05_rp, 0.047619047619048_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 1.0_rp, -0.33333333333333_rp, 0.25_rp, -0.2_rp, 0.16666666666667_rp, -0.14285714285714_rp, 0.125_rp, -0.11111111111111_rp, 0.1_rp, -0.090909090909091_rp, &
        0.083333333333333_rp, -0.076923076923077_rp, 0.071428571428571_rp, -0.066666666666667_rp, 0.0625_rp, -0.058823529411765_rp, 0.055555555555556_rp, -0.052631578947368_rp, 0.05_rp, -0.047619047619048_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, -0.5_rp, 0.35_rp, -0.275_rp, 0.22857142857143_rp, -0.19642857142857_rp, 0.17261904761905_rp, -0.15416666666667_rp, 0.13939393939394_rp, &
        -0.12727272727273_rp, 0.11713286713287_rp, -0.10851648351648_rp, 0.1010989010989_rp, -0.094642857142857_rp, 0.088970588235294_rp, -0.083946078431373_rp, 0.079463364293086_rp, -0.075438596491228_rp, 0.071804511278195_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, -0.4_rp, 0.3_rp, -0.24285714285714_rp, 0.20535714285714_rp, -0.17857142857143_rp, 0.15833333333333_rp, -0.14242424242424_rp, &
        0.12954545454545_rp, -0.11888111888112_rp, 0.10989010989011_rp, -0.1021978021978_rp, 0.095535714285714_rp, -0.089705882352941_rp, 0.084558823529412_rp, -0.079979360165119_rp, 0.075877192982456_rp, -0.07218045112782_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, -0.5_rp, 0.35714285714286_rp, -0.28571428571429_rp, 0.24107142857143_rp, -0.20982142857143_rp, 0.18641774891775_rp, &
        -0.16808712121212_rp, 0.1532634032634_rp, -0.14098401598402_rp, 0.13061938061938_rp, -0.12173763736264_rp, 0.11403118939884_rp, -0.10727415966387_rp, 0.10129643962848_rp, -0.095967169762642_rp, 0.091183841957836_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, -0.42857142857143_rp, 0.32142857142857_rp, -0.26190476190476_rp, 0.22321428571429_rp, -0.19561688311688_rp, &
        0.17471590909091_rp, -0.15821678321678_rp, 0.14479270729271_rp, -0.13361638361638_rp, 0.12414148351648_rp, -0.11599062702004_rp, 0.10889355742297_rp, -0.10265092879257_rp, 0.097112035603715_rp, -0.092160548429898_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, -0.5_rp, 0.36111111111111_rp, -0.29166666666667_rp, 0.24810606060606_rp, &
        -0.21748737373737_rp, 0.19445658508159_rp, -0.17633668414918_rp, 0.16161616161616_rp, -0.14936625874126_rp, 0.13897958658988_rp, -0.13003943564605_rp, 0.12224884893229_rp, -0.11538957688338_rp, 0.1092970631235_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, -0.44444444444444_rp, 0.33333333333333_rp, -0.27272727272727_rp, &
        0.23358585858586_rp, -0.20571095571096_rp, 0.18458624708625_rp, -0.16788073038073_rp, 0.15425590034965_rp, -0.14288101604278_rp, 0.13320921317245_rp, -0.12486355878384_rp, 0.11757449690402_rp, -0.11114336085311_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, -0.5_rp, 0.36363636363636_rp, &
        -0.29545454545455_rp, 0.25262237762238_rp, -0.22246503496503_rp, 0.19973776223776_rp, -0.18181818181818_rp, 0.16722638574661_rp, -0.15505370089469_rp, 0.14470631197363_rp, -0.13577707303822_rp, 0.12797596898071_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, -0.45454545454545_rp, &
        0.34090909090909_rp, -0.27972027972028_rp, 0.24038461538462_rp, -0.21241258741259_rp, 0.19121503496503_rp, -0.17443953105718_rp, 0.16074114304813_rp, -0.14928769228063_rp, 0.1395326175131_rp, -0.13109985264349_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, &
        -0.5_rp, 0.36538461538462_rp, -0.29807692307692_rp, 0.25576923076923_rp, -0.22596153846154_rp, 0.20347850678733_rp, -0.18573246606335_rp, 0.17126380537033_rp, -0.15917694317099_rp, 0.14888705830257_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
        1.0_rp, -0.46153846153846_rp, 0.34615384615385_rp, -0.28461538461538_rp, 0.24519230769231_rp, -0.21719457013575_rp, 0.19598416289593_rp, -0.17919445105978_rp, 0.16547634109312_rp, -0.15399685937128_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
        0.0_rp, 1.0_rp, -0.5_rp, 0.36666666666667_rp, -0.3_rp, 0.25808823529412_rp, -0.22855392156863_rp, 0.2062693498452_rp, -0.18867066563467_rp, 0.17431227425181_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
        0.0_rp, 0.0_rp, 1.0_rp, -0.46666666666667_rp, 0.35_rp, -0.28823529411765_rp, 0.24877450980392_rp, -0.22078173374613_rp, 0.19958397832817_rp, -0.1828044375645_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
        0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, -0.5_rp, 0.36764705882353_rp, -0.30147058823529_rp, 0.25986842105263_rp, -0.23055340557276_rp, 0.20843238611234_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
        0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, -0.47058823529412_rp, 0.35294117647059_rp, -0.29102167182663_rp, 0.2515479876161_rp, -0.22357363998231_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
        0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, -0.5_rp, 0.36842105263158_rp, -0.30263157894737_rp, 0.26127819548872_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
        0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, -0.47368421052632_rp, 0.35526315789474_rp, -0.29323308270677_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
        0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, -0.5_rp, 0.36904761904762_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
        0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, -0.47619047619048_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
        0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp]) &
]

type (pure_bend_multipole_struct), parameter :: h_to_v_imag(0:n_pole_maxx) = [ &
     pure_bend_multipole_struct([1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
        0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp]), &
     pure_bend_multipole_struct([0.0_rp, 1.0_rp, 0.5_rp, -0.16666666666667_rp, 0.041666666666667_rp, -0.025_rp, 0.0125_rp, -0.0089285714285714_rp, 0.0055803571428571_rp, -0.0043402777777778_rp, 0.0030381944444444_rp, -0.0024857954545455_rp, &
        0.0018643465909091_rp, -0.0015775240384615_rp, 0.0012394831730769_rp, -0.00107421875_rp, 0.000872802734375_rp, -0.00077012005974265_rp, 0.00064176671645221_rp, -0.00057421232524671_rp, 0.0004880804764597_rp, -0.00044159662155878_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 1.0_rp, 0.33333333333333_rp, -0.083333333333333_rp, 0.05_rp, -0.025_rp, 0.017857142857143_rp, -0.011160714285714_rp, 0.0086805555555556_rp, -0.0060763888888889_rp, 0.0049715909090909_rp, &
        -0.0037286931818182_rp, 0.0031550480769231_rp, -0.0024789663461538_rp, 0.0021484375_rp, -0.00174560546875_rp, 0.0015402401194853_rp, -0.0012835334329044_rp, 0.0011484246504934_rp, -0.00097616095291941_rp, 0.00088319324311756_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 0.5_rp, -0.15_rp, 0.05_rp, -0.032142857142857_rp, 0.01875_rp, -0.014136904761905_rp, 0.009672619047619_rp, -0.0078125_rp, &
        0.0058001893939394_rp, -0.0048759833916084_rp, 0.0038106424825175_rp, -0.0032902644230769_rp, 0.0026648888221154_rp, -0.0023459041819853_rp, 0.0019509708180147_rp, -0.0017429032931018_rp, 0.0014794411674003_rp, -0.0013370944145031_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 0.4_rp, -0.1_rp, 0.057142857142857_rp, -0.030357142857143_rp, 0.021825396825397_rp, -0.014384920634921_rp, 0.011363636363636_rp, &
        -0.0082859848484848_rp, 0.0068837412587413_rp, -0.0053267045454545_rp, 0.0045673076923077_rp, -0.0036771334134615_rp, 0.00322265625_rp, -0.0026697495404412_rp, 0.0023779145704334_rp, -0.0020131208579238_rp, 0.001815604685542_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 0.5_rp, -0.14285714285714_rp, 0.053571428571429_rp, -0.034722222222222_rp, 0.021329365079365_rp, -0.016233766233766_rp, &
        0.011498917748918_rp, -0.0093786421911422_rp, 0.0071478001165501_rp, -0.0060642482517483_rp, 0.0048384232954545_rp, -0.0042122543127828_rp, 0.0034693377050339_rp, -0.0030763230456656_rp, 0.002594088622291_rp, -0.0023322349344754_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 0.42857142857143_rp, -0.10714285714286_rp, 0.05952380952381_rp, -0.032738095238095_rp, 0.023538961038961_rp, &
        -0.015963203463203_rp, 0.012674825174825_rp, -0.009451486013986_rp, 0.0078999125874126_rp, -0.0062240493881119_rp, 0.0053688843325792_rp, -0.0043867982890271_rp, 0.0038661897736068_rp, -0.0032426107778638_rp, 0.0029029086963733_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 0.5_rp, -0.13888888888889_rp, 0.055555555555556_rp, -0.035984848484848_rp, &
        0.022727272727273_rp, -0.01733682983683_rp, 0.012529137529138_rp, -0.010257624320124_rp, 0.0079439223970474_rp, -0.0067684498406006_rp, 0.00547212438297_rp, -0.0047842608545189_rp, 0.0039845098907478_rp, -0.003547471620227_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 0.44444444444444_rp, -0.11111111111111_rp, 0.060606060606061_rp, &
        -0.034090909090909_rp, 0.024475524475524_rp, -0.016899766899767_rp, 0.013442113442113_rp, -0.010169604700855_rp, 0.0085240127519539_rp, -0.0067966500411353_rp, 0.0058812366039533_rp, -0.0048543226437842_rp, 0.0042918071833997_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 0.5_rp, -0.13636363636364_rp, &
        0.056818181818182_rp, -0.036713286713287_rp, 0.023601398601399_rp, -0.018006993006993_rp, 0.013188374125874_rp, -0.010814415364048_rp, 0.0084676123508844_rp, -0.0072304400804304_rp, 0.0059002545079727_rp, -0.0051710334943694_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 0.45454545454545_rp, &
        -0.11363636363636_rp, 0.061188811188811_rp, -0.034965034965035_rp, 0.025058275058275_rp, -0.017518939393939_rp, 0.013941022213081_rp, -0.010655015939942_rp, 0.0089429605804412_rp, -0.0071924042723917_rp, 0.0062345994192211_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, &
        0.5_rp, -0.13461538461538_rp, 0.057692307692308_rp, -0.037179487179487_rp, 0.024198717948718_rp, -0.018453054298643_rp, 0.01364536199095_rp, -0.011195932067159_rp, 0.0088372045427483_rp, -0.0075548232324306_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
        1.0_rp, 0.46153846153846_rp, -0.11538461538462_rp, 0.061538461538462_rp, -0.035576923076923_rp, 0.025452488687783_rp, -0.01795814479638_rp, 0.014289116456299_rp, -0.011004294177185_rp, 0.0092419432347838_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
        0.0_rp, 1.0_rp, 0.5_rp, -0.13333333333333_rp, 0.058333333333333_rp, -0.0375_rp, 0.024632352941176_rp, -0.018769349845201_rp, 0.013980263157895_rp, -0.011472269828984_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
        0.0_rp, 0.0_rp, 1.0_rp, 0.46666666666667_rp, -0.11666666666667_rp, 0.061764705882353_rp, -0.036029411764706_rp, 0.025735294117647_rp, -0.01828560371517_rp, 0.014544633642931_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
        0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 0.5_rp, -0.13235294117647_rp, 0.058823529411765_rp, -0.037732198142415_rp, 0.024961300309598_rp, -0.019004312251216_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
        0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 0.47058823529412_rp, -0.11764705882353_rp, 0.061919504643963_rp, -0.036377708978328_rp, 0.025947220993661_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
        0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 0.5_rp, -0.13157894736842_rp, 0.059210526315789_rp, -0.037907268170426_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
        0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 0.47368421052632_rp, -0.11842105263158_rp, 0.06203007518797_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
        0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 0.5_rp, -0.13095238095238_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
        0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 0.47619047619048_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
        0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp]) &
]

type (pure_bend_multipole_struct), parameter :: v_to_h_real(0:n_pole_maxx) = [ &
     pure_bend_multipole_struct([1.0_rp, -1.0_rp, 1.0_rp, -1.0_rp, 1.0_rp, -1.0_rp, 1.0_rp, -1.0_rp, 1.0_rp, -1.0_rp, 1.0_rp, -1.0_rp, &
        1.0_rp, -1.0_rp, 1.0_rp, -1.0_rp, 1.0_rp, -1.0_rp, 1.0_rp, -1.0_rp, 1.0_rp, -1.0_rp]), &
     pure_bend_multipole_struct([0.0_rp, 1.0_rp, -0.5_rp, 0.5_rp, -0.5_rp, 0.5_rp, -0.5_rp, 0.5_rp, -0.5_rp, 0.5_rp, -0.5_rp, 0.5_rp, &
        -0.5_rp, 0.5_rp, -0.5_rp, 0.5_rp, -0.5_rp, 0.5_rp, -0.5_rp, 0.5_rp, -0.5_rp, 0.5_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 1.0_rp, -0.66666666666667_rp, 0.58333333333333_rp, -0.55_rp, 0.53333333333333_rp, -0.52380952380952_rp, 0.51785714285714_rp, -0.51388888888889_rp, 0.51111111111111_rp, -0.50909090909091_rp, &
        0.50757575757576_rp, -0.50641025641026_rp, 0.50549450549451_rp, -0.5047619047619_rp, 0.50416666666667_rp, -0.50367647058824_rp, 0.50326797385621_rp, -0.50292397660819_rp, 0.50263157894737_rp, -0.50238095238095_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, -0.5_rp, 0.45_rp, -0.425_rp, 0.41071428571429_rp, -0.40178571428571_rp, 0.39583333333333_rp, -0.39166666666667_rp, 0.38863636363636_rp, &
        -0.38636363636364_rp, 0.38461538461538_rp, -0.38324175824176_rp, 0.38214285714286_rp, -0.38125_rp, 0.38051470588235_rp, -0.37990196078431_rp, 0.37938596491228_rp, -0.37894736842105_rp, 0.37857142857143_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, -0.6_rp, 0.5_rp, -0.45714285714286_rp, 0.43392857142857_rp, -0.41964285714286_rp, 0.41011904761905_rp, -0.40340909090909_rp, &
        0.39848484848485_rp, -0.39475524475524_rp, 0.39185814185814_rp, -0.38956043956044_rp, 0.38770604395604_rp, -0.38618697478992_rp, 0.38492647058824_rp, -0.38386867905057_rp, 0.38297213622291_rp, -0.38220551378446_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, -0.5_rp, 0.42857142857143_rp, -0.39285714285714_rp, 0.37202380952381_rp, -0.35863095238095_rp, 0.34943181818182_rp, &
        -0.34280303030303_rp, 0.33784965034965_rp, -0.33404095904096_rp, 0.33104395604396_rp, -0.32864010989011_rp, 0.32668067226891_rp, -0.3250612745098_rp, 0.32370678534572_rp, -0.32256191950464_rp, 0.32158521303258_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, -0.57142857142857_rp, 0.46428571428571_rp, -0.41666666666667_rp, 0.38988095238095_rp, -0.3728354978355_rp, &
        0.36113365800866_rp, -0.35267336829837_rp, 0.34632034632035_rp, -0.34140859140859_rp, 0.33752185314685_rp, -0.33438712023271_rp, 0.33181830424477_rp, -0.3296845053811_rp, 0.32789118937049_rp, -0.32636854083739_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, -0.5_rp, 0.41666666666667_rp, -0.375_rp, 0.35037878787879_rp, &
        -0.3342803030303_rp, 0.32302593240093_rp, -0.31477636946387_rp, 0.3085118006993_rp, -0.30362215909091_rp, 0.29972072963801_rp, -0.29655095211161_rp, 0.29393624226006_rp, -0.29175132223942_rp, 0.2899050245098_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, -0.55555555555556_rp, 0.44444444444444_rp, -0.39393939393939_rp, &
        0.36489898989899_rp, -0.34605672105672_rp, 0.33289627039627_rp, -0.32323232323232_rp, 0.31587206196581_rp, -0.31010740178939_rp, 0.30549110305544_rp, -0.30172682897383_rp, 0.29861059428832_rp, -0.29599753826969_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, -0.5_rp, 0.40909090909091_rp, &
        -0.36363636363636_rp, 0.33653846153846_rp, -0.31861888111888_rp, 0.30594405594406_rp, -0.2965472027972_rp, 0.28933405748663_rp, -0.2836466153332_rp, 0.2790652350262_rp, -0.27530969055132_rp, 0.27218580688854_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, -0.54545454545455_rp, &
        0.43181818181818_rp, -0.37937062937063_rp, 0.34877622377622_rp, -0.32867132867133_rp, 0.31446678321678_rp, -0.30392585355821_rp, 0.29581930018511_rp, -0.28941262394726_rp, 0.28423892948673_rp, -0.27998691094606_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, &
        -0.5_rp, 0.40384615384615_rp, -0.35576923076923_rp, 0.32692307692308_rp, -0.30769230769231_rp, 0.29397624434389_rp, -0.28372454751131_rp, 0.27579390182186_rp, -0.26949450389974_rp, 0.26438470283103_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
        1.0_rp, -0.53846153846154_rp, 0.42307692307692_rp, -0.36923076923077_rp, 0.3375_rp, -0.3164592760181_rp, 0.30147058823529_rp, -0.29026256251488_rp, 0.28158136609907_rp, -0.27467458769945_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
        0.0_rp, 1.0_rp, -0.5_rp, 0.4_rp, -0.35_rp, 0.31985294117647_rp, -0.29963235294118_rp, 0.28511996904025_rp, -0.27420665634675_rp, 0.26571449303406_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
        0.0_rp, 0.0_rp, 1.0_rp, -0.53333333333333_rp, 0.41666666666667_rp, -0.36176470588235_rp, 0.32916666666667_rp, -0.30740454076367_rp, 0.29180534055728_rp, -0.28007288441692_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
        0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, -0.5_rp, 0.39705882352941_rp, -0.34558823529412_rp, 0.31443498452012_rp, -0.29344040247678_rp, 0.27829914860681_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
        0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, -0.52941176470588_rp, 0.41176470588235_rp, -0.35603715170279_rp, 0.32275541795666_rp, -0.30042016806723_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
        0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, -0.5_rp, 0.39473684210526_rp, -0.34210526315789_rp, 0.31015037593985_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
        0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, -0.52631578947368_rp, 0.40789473684211_rp, -0.3515037593985_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
        0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, -0.5_rp, 0.39285714285714_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
        0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, -0.52380952380952_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
        0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp]) &
]

type (pure_bend_multipole_struct), parameter :: h_to_v_real(0:n_pole_maxx) = [ &
     pure_bend_multipole_struct([1.0_rp, 1.0_rp, -0.5_rp, 0.16666666666667_rp, -0.125_rp, 0.075_rp, -0.0625_rp, 0.044642857142857_rp, -0.0390625_rp, 0.030381944444444_rp, -0.02734375_rp, 0.022372159090909_rp, &
        -0.0205078125_rp, 0.017352764423077_rp, -0.01611328125_rp, 0.01396484375_rp, -0.013092041015625_rp, 0.01155180089614_rp, -0.010910034179688_rp, 0.0097616095291941_rp, -0.0092735290527344_rp, 0.0083903358096168_rp]), &
     pure_bend_multipole_struct([0.0_rp, 1.0_rp, 0.5_rp, -0.16666666666667_rp, 0.125_rp, -0.075_rp, 0.0625_rp, -0.044642857142857_rp, 0.0390625_rp, -0.030381944444444_rp, 0.02734375_rp, -0.022372159090909_rp, &
        0.0205078125_rp, -0.017352764423077_rp, 0.01611328125_rp, -0.01396484375_rp, 0.013092041015625_rp, -0.01155180089614_rp, 0.010910034179688_rp, -0.0097616095291941_rp, 0.0092735290527344_rp, -0.0083903358096168_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 1.0_rp, 0.66666666666667_rp, -0.25_rp, 0.1_rp, -0.075_rp, 0.05_rp, -0.042410714285714_rp, 0.032242063492063_rp, -0.028645833333333_rp, 0.023200757575758_rp, &
        -0.021129261363636_rp, 0.017782998251748_rp, -0.016451322115385_rp, 0.014212740384615_rp, -0.01329345703125_rp, 0.011705824908088_rp, -0.011038387522978_rp, 0.0098629411160023_rp, -0.0093596609015214_rp, 0.0084600615919682_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 0.5_rp, -0.15_rp, 0.1_rp, -0.060714285714286_rp, 0.049107142857143_rp, -0.035962301587302_rp, 0.03125_rp, -0.024857954545455_rp, &
        0.022372159090909_rp, -0.018643465909091_rp, 0.017127403846154_rp, -0.014708533653846_rp, 0.0136962890625_rp, -0.012013872931985_rp, 0.011295094209559_rp, -0.010065604289619_rp, 0.0095319245990954_rp, -0.008599513156671_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 0.6_rp, -0.2_rp, 0.085714285714286_rp, -0.0625_rp, 0.04265873015873_rp, -0.035714285714286_rp, 0.027597402597403_rp, &
        -0.02438446969697_rp, 0.02001384032634_rp, -0.018192744755245_rp, 0.015482954545455_rp, -0.014321664663462_rp, 0.012489615738122_rp, -0.011690027573529_rp, 0.010376354489164_rp, -0.0097953867247968_rp, 0.0088123094889683_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 0.5_rp, -0.14285714285714_rp, 0.089285714285714_rp, -0.054563492063492_rp, 0.043154761904762_rp, -0.031926406926407_rp, &
        0.027462121212121_rp, -0.022053467365967_rp, 0.019749781468531_rp, -0.016597465034965_rp, 0.015211838942308_rp, -0.013160394867081_rp, 0.012242934283088_rp, -0.010808702592879_rp, 0.010160180437307_rp, -0.0091056885593464_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 0.57142857142857_rp, -0.17857142857143_rp, 0.079365079365079_rp, -0.056547619047619_rp, 0.038961038961039_rp, &
        -0.03219696969697_rp, 0.025058275058275_rp, -0.021980623543124_rp, 0.018157536907537_rp, -0.016437663898601_rp, 0.01407117698478_rp, -0.012985850890837_rp, 0.011384313973565_rp, -0.010642414860681_rp, 0.0094910702767949_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 0.5_rp, -0.13888888888889_rp, 0.083333333333333_rp, -0.051136363636364_rp, &
        0.039772727272727_rp, -0.029574592074592_rp, 0.025203962703963_rp, -0.020339209401709_rp, 0.018113527097902_rp, -0.015292462592555_rp, 0.013967936934389_rp, -0.012135806609461_rp, 0.011265993856424_rp, -0.0099851823953266_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 0.55555555555556_rp, -0.16666666666667_rp, 0.075757575757576_rp, &
        -0.053030303030303_rp, 0.036713286713287_rp, -0.030011655011655_rp, 0.023445998445998_rp, -0.020427229020979_rp, 0.016935224701769_rp, -0.01526426239202_rp, 0.013111676684384_rp, -0.012065744820195_rp, 0.010611959232073_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 0.5_rp, -0.13636363636364_rp, &
        0.079545454545455_rp, -0.048951048951049_rp, 0.037587412587413_rp, -0.028030303030303_rp, 0.023699737762238_rp, -0.019179028691896_rp, 0.016991625102838_rp, -0.014384808544783_rp, 0.013092658780364_rp, -0.011405632913591_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 0.54545454545455_rp, &
        -0.15909090909091_rp, 0.073426573426573_rp, -0.050699300699301_rp, 0.035198135198135_rp, -0.028518356643357_rp, 0.022328774167009_rp, -0.019338428116002_rp, 0.016067644623179_rp, -0.014422844352822_rp, 0.012415818427573_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, &
        0.5_rp, -0.13461538461538_rp, 0.076923076923077_rp, -0.047435897435897_rp, 0.036057692307692_rp, -0.02693721719457_rp, 0.02262443438914_rp, -0.018340490295308_rp, 0.016173400660872_rp, -0.013716118722286_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
        1.0_rp, 0.53846153846154_rp, -0.15384615384615_rp, 0.071794871794872_rp, -0.049038461538462_rp, 0.034106334841629_rp, -0.027432126696833_rp, 0.021508097165992_rp, -0.018532128185282_rp, 0.015419319296545_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
        0.0_rp, 1.0_rp, 0.5_rp, -0.13333333333333_rp, 0.075_rp, -0.046323529411765_rp, 0.034926470588235_rp, -0.026122291021672_rp, 0.021816950464396_rp, -0.017705684247383_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
        0.0_rp, 0.0_rp, 1.0_rp, 0.53333333333333_rp, -0.15_rp, 0.070588235294118_rp, -0.047794117647059_rp, 0.03328173374613_rp, -0.026606037151703_rp, 0.020879404393336_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
        0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 0.5_rp, -0.13235294117647_rp, 0.073529411764706_rp, -0.04547213622291_rp, 0.03405572755418_rp, -0.025491117499631_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
        0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 0.52941176470588_rp, -0.14705882352941_rp, 0.069659442724458_rp, -0.046826625386997_rp, 0.032636738906089_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
        0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 0.5_rp, -0.13157894736842_rp, 0.072368421052632_rp, -0.044799498746867_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
        0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 0.52631578947368_rp, -0.14473684210526_rp, 0.068922305764411_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
        0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 0.5_rp, -0.13095238095238_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
        0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 0.52380952380952_rp]), &
     pure_bend_multipole_struct([0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
        0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp]) &
]

real(rp) g_n(0:n_pole_maxx), g
real(rp) an(0:n_pole_maxx), bn(0:n_pole_maxx), a_temp(0:n_pole_maxx), b_temp(0:n_pole_maxx)
integer out_type, i, j, np, np0

!

np = n_pole_maxx
do np0 = n_pole_maxx, 0, -1
  if (an(np0) /= 0 .or. bn(np0) /= 0) exit
enddo

g_n(0) = 1
do i = 1, np
  g_n(i) = g * g_n(i-1)
enddo

!

a_temp(0:np0) = an(0:np0)
b_temp(0:np0) = bn(0:np0)

an = 0
bn = 0

select case (out_type)
case (vertically_pure$)
  do i = 0, np0
    if (a_temp(i) /= 0) an(i:np) = an(i:np) + a_temp(i) * (h_to_v_real(i)%convert(i:np) * g_n(0:np-i))
    if (b_temp(i) /= 0) bn(i:np) = bn(i:np) + b_temp(i) * (h_to_v_imag(i)%convert(i:np) * g_n(0:np-i))
  enddo

case (horizontally_pure$)
  do i = 0, np0
    if (a_temp(i) /= 0) an(i:np) = an(i:np) + a_temp(i) * (v_to_h_real(i)%convert(i:np) * g_n(0:np-i))
    if (b_temp(i) /= 0) bn(i:np) = bn(i:np) + b_temp(i) * (v_to_h_imag(i)%convert(i:np) * g_n(0:np-i))
  enddo

case default
  if (global_com%exit_on_error) call err_exit
end select

end subroutine convert_bend_exact_multipole
