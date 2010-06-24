module mia_types

  !
  !Module contains global variables, constants, and type definitions.
  !It contains one method to initialize several global variables.
  !

  use nr
  use precision_def
  use physical_constants

  type cbpm_data
     !
     !Contains X and Y position data from input file for a single BPM
     !
     real(rp), allocatable :: x(:), y(:)  !X and Y vectors from raw data
  end type cbpm_data

  type twiss_parameters
     !
     !Contains the twiss parameters for a particular detector.
     !
     real(rp) :: gam2_beta_ratio  !from preceding CBPM to the present CBPM
     real(rp) :: d_phase_adv 	  !phase advance from preceding CBPM
                                  !to the present CBPM
     real(rp) :: d_delta          !change in delta = out-of-plane phase shifts
                                  !from preceding CBPM to the present CBPM
     real(rp) :: beta             !Beta
     real(rp) :: alpha            !Alpha
     real(rp) :: phi              !Phase advance (cumulative)
   !The following variables are arrays because they are calculated for both A and B modes.
     real(rp) :: ratio(2),&       !numer/denom
          magnitude2(2), &        !numer^2+denom^2
          numer(2), denom(2)      !lambda*pi of positive and negative cols
     real(rp) :: inv_gamma_cbar_sqrt_betas(2,2) 
                          !(C bar of i,j * sqrt(beta-B/beta-A))/gamma (pg 16)
                          !This is an intermediate variable
     real(rp) :: gam2_beta        !gamma^2*beta
     real(rp) :: j_amp(2)         !Amplitude (A and B modes)
  end type twiss_parameters

  type processor_analysis
     !
     !Contains data and A and B mode twiss parameters for a single detector.
     !
     type (twiss_parameters) :: a, b      !A and B mode twiss parameters
     real(rp) :: inv_gamma_cbar(2,2)      !(1/gamma)*Cbar
     real(rp) :: inv_gamma_cbar_check(2)  !Consistancy check value
     real(rp) :: gamma                    !Gamma (gamma^2 + det_Cbar = 1)
     real(rp) :: cbar(2,2)                !Cbar
     real(rp) :: sqrt_beta_cbar(2,2)      !sqrt(beta)*cbar
  end type processor_analysis

  type active_processors
     !
     !Contains identifying information for a single detector
     !
     character (13) :: label      !Full BPM label (ex. BPM09E)
     real(rp) :: number, &        !BPM number (negative for east BPMs)
          sPos                    !S position (for sorting and CESRV)
     logical :: is_west           !If the BPM is in the west   !Not needed now
     integer :: eleNum            !Element number (for CESRV)
  end type active_processors

  type cbpm_analysis
     !
     !Contains information for the full set of BPMs including their identification and
     !data. Also contains several derived values for the entire machine, such as average
     !amplitude and tune.
     !
     integer :: col_a_p, col_a_n, &  !A mode positive and negative
          col_b_p, col_b_n           !B mode positive and negative
     integer :: set_num_a, set_num_b !File number (1 or 2) of the data in that mode
     type (processor_analysis), allocatable :: loc(:) !Twiss, etc calculated at BPM
     type (active_processors), allocatable :: proc(:) !BPM data (name, element #, etc)
     real(rp) :: tune(2)                     !Tune of the machine (A&B modes)
     real(rp) :: phi_t(2)                    !Total phase
     real(rp) :: j_amp_ave(2)              !Average amplitude (A/B mode, File #) (no more file #)
     integer :: intTune(2)                   !Integer tunes
  end type cbpm_analysis

  type data_set
     !
     !Contains raw data from one input file as well as
     !results from SVD and FFT.
     !
     character (120) :: filename   !Full filename
     character (40) :: shortname   !Filename without path (used for display)
     integer :: bunch              !Bunch # of interest (default is 1)
     real(rp), allocatable :: tau_mat(:,:), &   !Position matrix
          pi_mat(:,:), &           !Pi matrix
          poshis(:,:), &           !Position history
          phi_spec(:,:), &         !Takes data from polar FFT
          spectrum(:,:), &         !Spectrum (FFT)
          lambda(:)                !Eigenvalues from SVD
     integer, allocatable :: fr_peak(:)  !Frequency peak (was nmax)
     integer :: bpmproc            !Number of active BPM processors
     integer :: numturns           !Number of turns
     type (active_processors), allocatable :: proc(:) !Data for BPM processors
     type (cbpm_data), allocatable :: cdata(:)     !Raw data from the file being read in
     real(rp) :: noise           !Relative noise level from higher eigenvalues
  end type data_set                                


  type known_spacings
     !
     !Data for BPMs separated by a drift space; used
     !to find them and set length scale for beta ratios
     !
     character (13) :: bpm_name(2) !Name of the BPM
     logical :: in_use,&           !If the BPM is in use?
          has_one                  !True if a BPM has a known spacing
     real(rp) :: length            !Distance between BPMs
     integer :: number             !Number of BPM pairs
     integer :: bpm_pntr(2)        !Points to the pair's position in BPM array
  end type known_spacings


  !****************
  !Global variables
  !****************
  type(cbpm_analysis), target :: data_struc !Contains data from calculations
  !BPM pairs with known spacing:
  type(known_spacings), allocatable:: bpm_pairs(:) 

  integer:: NUM_BPMS, &                  !Number of BPMs
       NUM_TURNS, &                      !Number of turns
       nset                              !Number of files
  real(rp), parameter :: FREQ=390.12     !Frequency of the machine (in MHz)
  real(rp), parameter :: ip_L3=384.213   !Location (m) of IP_L3
  logical :: outfile                     !True if user specified an output file
  character(100) :: outname              !Name of output file
  logical :: fileread(2)                 !True if filename has been read
  character (120) :: xfile,yfile         !Filenames for X and Y modes
  logical :: debugMode, &                !Enables special plots for diagnosing
                                         !problems with input files
       silentMode,&                      !MIA runs without plotting
       phase,&                           !Plot phase advance instead of beta
       postScript                        !Print postscript plot only
  logical :: vetoBPM                     !Remove detectors; specified
                                         !as an option in the plotting routine
  real(rp) :: endLoc                     !Location of the end of the machine
                                         !Used in plotting east before west
  character(120) :: inputPasser          !Not very elegant, but passes
                                         !input between modules


contains

  subroutine initialize_struct(data)
    !
    !Allocate global arrays, assign global constants
    !
    type(data_set) data
    integer :: i

    vetoBPM = 0
    NUM_BPMS = data%bpmproc
    NUM_TURNS = data%numturns
    allocate(data_struc%loc(NUM_BPMS))
    allocate(data_struc%proc(NUM_BPMS))
    data_struc%col_a_p = 0
    data_struc%col_a_n = 0
    data_struc%col_b_p = 0
    data_struc%col_b_n = 0
    data_struc%set_num_a = 0
    data_struc%set_num_b = 0
    data_struc%tune(1:2) = 0.
    data_struc%phi_t(1:2) = 0.
    data_struc%j_amp_ave(1:2) = 0.

    do i=1, NUM_BPMS
       data_struc%loc(i)%inv_gamma_cbar(1:2,1:2) = 0.
       data_struc%loc(i)%inv_gamma_cbar_check(1:2) = 0.
       data_struc%loc(i)%gamma = 0.
       data_struc%loc(i)%cbar(1:2,1:2) = 0.
       data_struc%loc(i)%sqrt_beta_cbar(1:2,1:2) = 0.
       data_struc%loc(i)%gamma = 0.

       data_struc%loc(i)%a%gam2_beta_ratio = 0.
       data_struc%loc(i)%a%d_phase_adv = 0.
       data_struc%loc(i)%a%d_delta = 0.
       data_struc%loc(i)%a%beta = 0.
       data_struc%loc(i)%a%alpha = 0.
       data_struc%loc(i)%a%phi = 0.
       data_struc%loc(i)%a%ratio(1:2) = 0.
       data_struc%loc(i)%a%magnitude2(1:2) = 0.
       data_struc%loc(i)%a%numer(1:2) = 0.
       data_struc%loc(i)%a%denom(1:2) = 0.
       data_struc%loc(i)%a%inv_gamma_cbar_sqrt_betas(1:2,1:2) = 0.
       data_struc%loc(i)%a%gam2_beta = 0.

       data_struc%loc(i)%b%gam2_beta_ratio = 0.
       data_struc%loc(i)%b%d_phase_adv = 0.
       data_struc%loc(i)%b%d_delta = 0.
       data_struc%loc(i)%b%beta = 0.
       data_struc%loc(i)%b%alpha = 0.
       data_struc%loc(i)%b%phi = 0.
       data_struc%loc(i)%b%ratio(1:2) = 0.
       data_struc%loc(i)%b%magnitude2(1:2) = 0.
       data_struc%loc(i)%b%numer(1:2) = 0.
       data_struc%loc(i)%b%denom(1:2) = 0.
       data_struc%loc(i)%b%inv_gamma_cbar_sqrt_betas(1:2,1:2) = 0.
       data_struc%loc(i)%b%gam2_beta = 0.

       data_struc%proc(i)%label = ""
       data_struc%proc(i)%is_west = .false.
       data_struc%proc(i)%number = 0
       data_struc%proc(i)%sPos = -9999.0
       data_struc%proc(i)%eleNum = -999
    enddo

  end subroutine initialize_struct

  subroutine allocate_bpm_pairs(length)
    integer :: length
    allocate (bpm_pairs(length))
  end subroutine allocate_bpm_pairs

  subroutine clean(keepPairs)
    !
    !Arg is true if new files will be read in;
    !don't reread knownl.inp to get paired detectors.
    !
    !Deallocates data_struc and bpm_pairs either when the program
    !exits or when a new analysis is being done.
    !

    logical :: keepPairs
    if (.not. keepPairs) then
       deallocate (bpm_pairs)
    endif
    deallocate (data_struc%loc)
    deallocate (data_struc%proc)
  end subroutine clean

  subroutine sort_l(list, length)
    real (rp) :: list(:), temp
    integer:: length, i, j

    do i=1, length
       do j=1, length
          if (list(i) > list(j)) then
             temp = list(i)
             list(i) = list(j)
             list(j) = temp
          endif
       enddo
    enddo

  end subroutine sort_l

end module mia_types
