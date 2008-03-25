module mia_types

  !
  !Module contains global variables, constants, and type definitions.
  !It contains one method to initialize several global variables.
  !

  use nr
  use precision_def
  use physical_constants

  type cbpm_data
     real(rp), allocatable :: x(:), y(:)  !X and Y vectors from raw data
  end type cbpm_data

  type twiss_parameters
     real(rp) :: gam2_beta_ratio   !from preceding CBPM to the present CBPM
     real(rp) :: d_phase_adv 	   !from preceding CBPM to the present CBPM
     real(rp) :: d_delta           !Delta_delta = out-of-plane phase shifts
     !from preceding CBPM to the present CBPM
     real(rp) :: beta              !Beta
     real(rp) :: alpha             !Alpha
     real(rp) :: phi               !Phase advance
     real(rp) :: ratio(2),&        !numer/denom
          magnitude2(2), &         !numer^2+denom^2
          numer(2), denom(2)       !lambda*pi of positive and negative cols
     real(rp) :: inv_gamma_cbar_sqrt_betas(2,2) 
                          !(C bar of i,j * sqrt(beta-B/beta-A))/gamma (pg 16)
     real(rp) :: gam2_beta         !gamma^2*beta
  end type twiss_parameters

  type processor_analysis
     type (twiss_parameters) a, b         !A and B modes
     real(rp) :: inv_gamma_cbar(2,2)      !(1/gamma)*Cbar
     real(rp) :: inv_gamma_cbar_check(2)  !Consistancy check value
     real(rp) :: gamma                    !Gamma (gamma^2 + det_Cbar = 1)
     real(rp) :: cbar(2,2)                !Cbar
     real(rp) :: sqrt_beta_cbar(2,2)      !sqrt(beta)*cbar
  end type processor_analysis

  type active_processors
     character (13) :: label      !Full BPM label (ex. BPM09E)
     real(rp) :: number               !BPM number (negative for east BPMs)
     logical :: is_west           !If the BPM is in the west   !Not needed now
  end type active_processors

  type cbpm_analysis
     integer :: col_a_p, col_a_n, &  !Horizontal positive and negative
          col_b_p, col_b_n           !Vertical positive and negative
     integer :: set_num_a, set_num_b !Which of the two input files (1 or 2) is 
     !horizontal (a) or vertical(b)
     type (processor_analysis), allocatable :: loc(:)
     !Information for each BPM location
     type (active_processors), allocatable :: proc(:)
     !Array of processors 
     real(rp) :: tune(2)                     !Tune of the machine
     real(rp) :: phi_t(2)                    !phi(t)
     real(rp) :: j_amp_ave(2,2)           !Average amplitude (A and B modes)
                                          !Second index is the file
  end type cbpm_analysis

  type data_set
     character (60) :: filename    !Full filename
     character (20) :: shortname   !Filename without path (used for printing)
     real(rp), allocatable :: tau_mat(:,:), &   !Position matrix
          pi_mat(:,:), &           !Pi matrix
          poshis(:,:), &           !Position history
          phi_spec(:,:), &         !Takes data from polar FFT
          spectrum(:,:), &         !Spectrum (FFT)
          lambda(:)                !Eigenvalues from SVD
     integer, allocatable :: fr_peak(:)  !Frequency peak (was nmax)
     integer :: bpmproc            !Number of active BPM processors
     integer :: numturns           !Number of turns
     type (cbpm_data), allocatable :: cdata(:)     !Temporarily holds data
  end type data_set                                !from file being read in.

  type data_file
     type (active_processors), allocatable :: proc(:) !Data for BPM processors
  end type data_file

  type bpm_point
     integer :: bpm_pntr(2) !Points to a BPM
  end type bpm_point


  type known_spacings
     character (13) :: bpm_name(2) !Name of the BPM
     logical :: in_use,&           !If the BPM is in use?
          has_one                  !True if a BPM has a known spacing
     real(rp) :: length            !Distance between BPMs
     type (bpm_point) :: file(2)   !File pointer
     integer :: number             !Number of BPM pairs
  end type known_spacings

  !Global variables
!  type (data_file), allocatable :: file(:) !File being used
  type(cbpm_analysis) data_struc
  !BPM pairs with known spacing:
  type(known_spacings), allocatable:: bpm_pairs(:) 
  integer:: NUM_BPMS, &                  !Number of BPMs
       NUM_TURNS                         !Number of turns
  !Change FREQ for use with machines other than CESR
  real(rp), parameter :: FREQ=390.12         !Frequency of the machine (in MHz)

contains

  subroutine initialize_structures(data, nset)
    !
    !Allocate global arrays, assign global constants
    !
    type(data_set) data
    integer :: nset, i
    NUM_BPMS = data%bpmproc
    NUM_TURNS = data%numturns
    allocate(data_struc%loc(NUM_BPMS))
    allocate(data_struc%proc(NUM_BPMS))
!    allocate(file(nset))
!    do i=1, nset
!       allocate(file(i)%proc(1:data%bpmproc))
!    enddo
  end subroutine initialize_structures

  subroutine allocate_bpm_pairs(num)
    integer:: num
    allocate(bpm_pairs(num))
  end subroutine allocate_bpm_pairs
end module mia_types
