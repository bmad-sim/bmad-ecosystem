module orbit_mia

  !this module is linked to mia.f90 and contains subroutines needed to read,
  !structure, process, and display the data taken from CESR concerning the 
  !motion of a beam with transverse coupling

  use mia_types
  use mia_input
  use mia_matrixops

  implicit none

  integer :: iset, icolumn, &      !Current set number, column to be used
       arr_length                  !Length of arrays
  !*Make ints local variables

contains 

  subroutine convert_data_from_pi_matrix(data)

    !
    !can find beta ratio and magnitude^2 from paper page 10
    !

    type(data_set) data(*)        !Data from file
    integer ::  n_bpm, &          !BPM #
         i, count, &              !Counters
         count_x, count_y, &      !Counters
         mag, &                   !Counter for magnitude^2; either 1 or 2
         numWest
    real(rp) ::  sum_a, sum_b     !Sum of phase advances (phi) for A, B     
    integer :: east, west, bpm, n_bpm_e
    real(rp) :: dS

    do n_bpm = 1, NUM_BPMS

       DO count = 1, 2

          call pi_calc(data(data_struc%set_num_a),&
               data_struc%loc(n_bpm)%a, &
               data_struc%col_a_p, &
               data_struc%col_a_n, count, n_bpm)
          call pi_calc(data(data_struc%set_num_b),&
               data_struc%loc(n_bpm)%b, &
               data_struc%col_b_p, &
               data_struc%col_b_n,count, n_bpm)
       ENDDO
!!$       If (n_bpm == 10 .or. n_bpm == 11) then
!!$          Print *, data_struc%proc(n_bpm)%label
!!$          Print *, "numer:", data_struc%loc(n_bpm)%a%numer(1), &
!!$               data_struc%loc(n_bpm)%b%numer(2)
!!$          Print *, "denom:", data_struc%loc(n_bpm)%a%denom(1),&
!!$               data_struc%loc(n_bpm)%b%denom(2)
!!$
!!$          Print *, data(data_struc%set_num_a)%lambda(data_struc%col_a_p), &
!!$               data(data_struc%set_num_a)%pi_mat(2*n_bpm-1, data_struc%col_a_p)
!!$          Print *, data(data_struc%set_num_a)%lambda(data_struc%col_a_n), &
!!$               data(data_struc%set_num_a)%pi_mat(2*n_bpm-1, data_struc%col_a_n)
!!$
!!$          Print *, data(data_struc%set_num_b)%lambda(data_struc%col_b_p), &
!!$               data(data_struc%set_num_b)%pi_mat(2*n_bpm, data_struc%col_b_p)
!!$          Print *, data(data_struc%set_num_b)%lambda(data_struc%col_b_n), &
!!$               data(data_struc%set_num_b)%pi_mat(2*n_bpm, data_struc%col_b_n)
!!$
!!$          Print *, ""
!!$       endif

    end do

    !
    !Calculates phase adv around ring in clockwise direction(pg15 Delta Phi)
    !Phase is defined as 0 at BPM 0W.
    !
    !Find BPMs 0E and 0W
!    do bpm = 1, NUM_BPMS
!       if (data_struc%proc(bpm)%number == 0) then
!          if (data_struc%proc(bpm)%is_west) then
!             west = bpm
!          else
!             east = bpm
!          endif
!       endif
!    enddo
    !***(*)_* Temp: define first bpm as 0 phase, etc
    west = 1
    east = NUM_BPMS

    data_struc%loc(west)%a%d_phase_adv = 0.
    data_struc%loc(west)%b%d_phase_adv = 0.
    data_struc%loc(west)%a%gam2_beta_ratio = 1.0
    data_struc%loc(west)%b%gam2_beta_ratio = 1.0
    !dS is the physical distance between the first and last detectors
    !The first term is the total length of the machine (freq is in kHz)
    dS = ( 3e8/FREQ/1000) - data_struc%proc(east)%sPos + &
         data_struc%proc(west)%sPos
    print *, "dS: ", dS


    call phase_adv(data_struc%loc(east)%a, data_struc%loc(west)%a, 1, dS)
    call phase_adv(data_struc%loc(east)%b, data_struc%loc(west)%b, 2, dS) 
print *, "phase A initial: ", data_struc%loc(east)%a%d_phase_adv, "B: ", data_struc%loc(east)%b%d_phase_adv

    data_struc%loc(east)%a%d_phase_adv = &
         mod((data_struc%loc(east)%a%d_phase_adv - &
         data_struc%phi_t(1))+7.0*pi, 2.0*pi) -pi
    data_struc%loc(east)%b%d_phase_adv = &
         mod((data_struc%loc(east)%b%d_phase_adv - &
         data_struc%phi_t(2)+7.0*pi),2.0*pi) - pi
    print *, "phase A: ", data_struc%loc(east)%a%d_phase_adv, "B: ", data_struc%loc(east)%b%d_phase_adv

    do n_bpm = 2, NUM_BPMS

       call phase_adv(data_struc%loc(n_bpm)%a, &
            data_struc%loc(n_bpm-1)%a, 1, &
            (data_struc%proc(n_bpm)%sPos - data_struc%proc(n_bpm-1)%sPos))
       
       call phase_adv(data_struc%loc(n_bpm)%b, &
            data_struc%loc(n_bpm-1)%b, 2, &
            (data_struc%proc(n_bpm)%sPos - data_struc%proc(n_bpm-1)%sPos))
    end do
   
    sum_a = 0.0
    sum_b = 0.0
    !Assigns phi by summing delta phase advances
    do i = 1, NUM_BPMS
       !Reset phi to 0 when switching from west to east

       !Temp disabled--no large gap in detectors anymore
       !Need to have something that detects big gaps and compensates
!       if (data_struc%proc(i)%is_west) then

          sum_a = sum_a + data_struc%loc(i)%a%d_phase_adv
          data_struc%loc(i)%a%phi = sum_a

          sum_b = sum_b + data_struc%loc(i)%b%d_phase_adv
          data_struc%loc(i)%b%phi = sum_b
!       else 
!          if (.not. data_struc%proc(i)%is_west .and. &
!               data_struc%proc(i-1)%is_west) then
!             sum_a = 0.0
!             sum_b = 0.0
!             numWest = i
!          endif
!          sum_a = sum_a + data_struc%loc(NUM_BPMS-(i-numWest))%a%d_phase_adv
!          data_struc%loc(NUM_BPMS-(i-numWest))%a%phi = sum_a

!          sum_b = sum_b + data_struc%loc(NUM_BPMS-(i-numWest))%b%d_phase_adv
!          data_struc%loc(NUM_BPMS-(i-numWest))%b%phi = sum_b
!       endif
    end do

    call intTune(data_struc%loc(NUM_BPMS)%a%phi,1)
    call intTune(data_struc%loc(NUM_BPMS)%b%phi,2)

    do n_bpm = 1, NUM_BPMS

       !
       !finds d_delta (page 15 bottom section)
       !

       call d_delta(data_struc%loc(n_bpm)%a, 2)
       call d_delta(data_struc%loc(n_bpm)%b, 1)

       !
       !Section finds (C bar/gamma of i,j * sqrt(beta-B/beta-A)) (pg 16)
       !

       !A mode - Cbar (1,2) and (2,2)
       do count_x = 1, 2

          count_y = 2
          mag = 2

          call inv_gamma_cbar_sqrt_betas(data_struc%loc(n_bpm)%a, &
               count_x, count_y, mag)

       enddo

       !B mode - Cbar (1,1) and (1,2)
       do count_y = 1, 2 

          count_x = 1
          mag = 1

          call inv_gamma_cbar_sqrt_betas(data_struc%loc(n_bpm)%b, &
               count_x, count_y, mag)

       enddo

    end do

  end subroutine convert_data_from_pi_matrix


  subroutine pi_calc(data_snum, twiss, col_p, col_n, count, n_bpm)
    !
    !Calculates numer, denom, and magnitude^2
    !
    type(data_set) data_snum      !Data from file
    INTEGER :: col_p, col_n, &    !Positive or negative col
         count, n_bpm             !Counter, bpm number
    type(twiss_parameters) twiss  !Location of values in data structure

    !Numer is lambda*pi of positive col, denom same of negative col.
    !Count is 1 for A mode, 2 for B mode
    twiss%numer(count)= &
         data_snum%lambda(col_p) &
         * data_snum%pi_mat(2*n_bpm+(count-2),col_p)

    twiss%denom(count)= &
         data_snum%lambda(col_n) &
         * data_snum%pi_mat(2*n_bpm+(count-2),col_n)

    !Ratio is only used as a graphing option. It may be useless.
    twiss%ratio(count) = twiss%numer(count)/ twiss%denom(count)
    twiss%magnitude2(count) = twiss%numer(count)**2 + twiss%denom(count)**2
  

  end subroutine pi_calc


  subroutine inv_gamma_cbar_sqrt_betas(twiss, ix, iy, mag)
    !
    !Finds (C bar/gamma of i,j * sqrt(beta-B/beta-A))*sin(d_delta) (pg 16)
    !
    type(twiss_parameters) twiss !Location of inv_gam... in data struct
    INTEGER :: ix, iy, &         !Counter for index of Cbar
         mag                     !Either 1 or 2; counter for magnitude^2
    real(rp) :: sdelta, &        !Either sin or cos or d_delta
         temp                    !Temporarily holds +-d_delta (shorter code)

    !?@@
    if (mag ==2) then            !A modes
       temp = -twiss%d_delta
    else                         !B modes
       temp = twiss%d_delta
    endif

    if (ix /= iy) then           !Only applies to Cbar(1,2)--no Cbar(2,1)
       sdelta = -cos(temp)
    else
       sdelta = sin(temp)
    endif

    twiss%inv_gamma_cbar_sqrt_betas(ix, iy) = &
         sqrt(twiss%magnitude2(mag) / &
         twiss%magnitude2(3-mag)) * sdelta          

  end subroutine inv_gamma_cbar_sqrt_betas

  subroutine d_delta(twiss, count)
    !
    !Finds d_delta+
    !
    type(twiss_parameters) twiss !Location of d_delta in dat struct
    INTEGER ::count              !Counter (either 1 or 2)
    !Count is 2 for A mode and 1 for B mode

    !?@@
    !Mod function maps d_delta from -pi to +pi
    twiss%d_delta = mod(atan2(twiss%numer(count), twiss%denom(count)) &
         - atan2(twiss%numer(3-count), twiss%denom(3-count))&
         +7.0*pi+pi/2, 2.0*pi)-pi

  end subroutine d_delta

  subroutine phase_adv(twiss, twiss_old, count, dSPos)
    !
    !Finds phase advance and gamma^2 beta ratio (both pg 15)
    !
    type(twiss_parameters) twiss, & !Location of this numer, denom
         twiss_old                  !Location of previous number, denom
    INTEGER :: count                !Counter (1 for A, 2 for B mode)
    real(rp) :: dSPos, dSMax
    integer :: piFactor


    !Approximate length scale for a phase advance equal to pi for A & B modes
    !Accounts for large gaps between detectors where the phase advance is greater than pi
    !***Disabled for now
!    dSMax = data_struc%intTune(count)*pi

!    piFactor = Floor(dSPos / dSMax)

!    if (piFactor > 0) then
!       Print *, "dS:", dSPos, " factor: ", piFactor
!    end if

    twiss%d_phase_adv =  &
         mod( atan2(twiss%numer(count), twiss%denom(count) ) - &
         atan2(twiss_old%numer(count), twiss_old%denom(count)) + 7.0*pi, 2.0*pi)-pi

!    if (twiss%d_phase_adv < 0) then
!       twiss%d_phase_adv = twiss%d_phase_adv + pi
!       Print *, "Numer 2: ", twiss%numer(count), " 1: ", twiss_old%numer(count)
!       Print *, "Denom 2: ", twiss%denom(count), " 1: ", twiss_old%denom(count)
!       Print *, "---------"
!       Print *, atan2(twiss%numer(count),twiss%denom(count)), " - ", &
!            atan2(twiss_old%numer(count),twiss_old%denom(count))
!       Print *, atan2(twiss%numer(count),twiss%denom(count)) - &
!            atan2(twiss_old%numer(count),twiss_old%denom(count))
!       Print *, twiss%d_phase_adv

!       if (atan2(twiss_old%numer(count),twiss_old%denom(count)) < &
!            atan2(twiss%numer(count),twiss%denom(count))) then
!          Print *,"    ", atan2(twiss%numer(count),twiss%denom(count)) - &
!               (pi + atan2(twiss_old%numer(count),twiss_old%denom(count)))
!       end if
!    end if

!    if (twiss%d_phase_adv < 0) then
!       twiss%d_phase_adv =  &
!            mod( atan2(twiss%denom(count), twiss%numer(count) ) - &
!            atan2(twiss_old%denom(count), twiss_old%numer(count)) + 7.0*pi, 2.0*pi)-pi
!    else 
!       twiss%d_phase_adv = pi - twiss%d_phase_adv 
!    end if

    !Mod function maps d_phase_advance from -pi to +pi
 !   twiss%d_phase_adv =  &
 !        abs(mod( atan2(twiss%numer(count), twiss%denom(count) ) - &
 !        atan2(twiss_old%numer(count), twiss_old%denom(count)) &
 !        +7.0 * pi , 2.0 * pi) - pi) !&
 !        + piFactor * pi


    twiss%gam2_beta_ratio =  &
         twiss_old%magnitude2(count) / twiss%magnitude2(count)

  end subroutine phase_adv


  subroutine intTune(phiMax, count)

    real(rp) :: phiMax
    integer :: count, input
    character(1) :: mode(2)
    mode(1) = "A"
    mode(2) = "B"

    data_struc%intTune(count) = floor(phiMax / (2*pi))

!   Print *, "I think the integer tune for mode ", mode(count), " is ", data_struc%intTune(count)
!    Print *, "Is this correct? Press enter to accept or input the correct value."
!    Read (*,"(i)"), input

    if (input > 0) then
       data_struc%intTune(count) = input
    end if

  end subroutine intTune


  subroutine ring_beta(ring, ring2, i)
    !
    !Calculates beta for a pair of BPMs with known spacing.
    !
    INTEGER :: i                           !Counter for BPM pair
    type(twiss_parameters) ring, ring2     !Ring(i)
    !First BPM:
    ring%beta = bpm_pairs(i)%length / &
         sin(ring2%d_phase_adv) * sqrt(ring2%gam2_beta_ratio)
    !Second BPM:
    ring2%beta = bpm_pairs(i)%length / &
         sin(ring2%d_phase_adv) * sqrt(1/ring2%gam2_beta_ratio)
  end subroutine ring_beta

  subroutine alpha(twiss, d_phase_adv, i, q)
    !
    !Calculates alpha
    !
    type(twiss_parameters) twiss          !A or B mode in ring structure
    real(rp) :: d_phase_adv               !d_phase_adv between BPM pair
    integer :: i, q                       !BPM pair numbers

    twiss%alpha = twiss%beta / bpm_pairs(i)%length - 1.0/tan(d_phase_adv)
    !Second BPM of the pair has the wrong sign.
    if (q == 2) then
       twiss%alpha = -twiss%alpha
    endif
  end subroutine alpha

  subroutine inv_gamma_cbar(ring, index)
    !
    !Calculates cbar/gamma(1,1) and (2,2).
    !
    type(processor_analysis) ring          !Ring(i)
    real(rp) :: beta_ratio, mag2_ratio, &  !Ratios of beta and mag2 A and B
         sin_d_phase_adv                   !sin(+-d_phase_adv)
    integer :: index                       !Indices for Cbar (1 or 2)

    if (index == 1) then
       beta_ratio = ring%b%beta/ring%a%beta
       mag2_ratio = ring%b%magnitude2(1)/ring%b%magnitude2(2) 
       sin_d_phase_adv = sin(ring%b%d_delta)
    else
       beta_ratio = ring%a%beta/ring%b%beta
       mag2_ratio = ring%a%magnitude2(2)/ring%a%magnitude2(1)
       sin_d_phase_adv = sin(-1.0*ring%a%d_delta)
    endif

    ring%inv_gamma_cbar(index, index) = sqrt(beta_ratio*mag2_ratio) * &
         sin_d_phase_adv

  end subroutine inv_gamma_cbar

  subroutine inv_gamma_cbar12(ring)
    !
    !Calculates Cbar/gamma(1,2)
    !
    type(processor_analysis) ring    !Ring(i)
    real(rp) :: beta_ratio, &        !Ratio of betas from A and B modes
         mag2_ratio,&                !Ratio of mag2 from A and B modes
         sin_d_phase_adv,&           !cos(+-d_phase_adv)
         temp(2)                     !Holds two initial values of Cbar/gamma
    integer :: i                     !A (2) or B (1) mode counter

    do i=1, 2
       if (i == 1) then
          beta_ratio = ring%b%beta/ring%a%beta
          mag2_ratio = ring%b%magnitude2(1)/ring%b%magnitude2(2)
          sin_d_phase_adv = cos(ring%b%d_delta)
          temp(i) = ring%b%inv_gamma_cbar_sqrt_betas(1,2) * sqrt(beta_ratio)
       else
          beta_ratio = ring%a%beta/ring%b%beta
          mag2_ratio = ring%a%magnitude2(2)/ring%a%magnitude2(1)
          sin_d_phase_adv = cos(-ring%a%d_delta)
          temp(i) = ring%a%inv_gamma_cbar_sqrt_betas(1,2) * &
               sqrt(beta_ratio)
       endif
    enddo

    ring%inv_gamma_cbar(1,2) = -(temp(1)+temp(2))/2

  end subroutine inv_gamma_cbar12

  subroutine avg(dat_struc, ring, n_ring)
    !
    !Averages twiss values calculated for each file
    !

    real(rp) dat_struc,&          !Average value 
         ring                     !Value for individual file
    integer :: n_ring             !Number of files

    dat_struc = dat_struc + ring/n_ring
    
  end subroutine avg

  subroutine calculate_with_known_spacing (data)

    implicit none
    type (data_set)data(*)        !Data from file
    integer :: i,j, k, ik, jm, &  !Counters
         bpm1, bpm2, &            !Pair of BPMs with known spacing
         n_ring, &                !Number of BPM pairs in use
         r, q,u                   !Counters
    INTEGER :: countx, county, &  !Counters
         count                    !Too many counters
    real(rp) :: ca, cb, &         !Cos of phase advance for A, B
         sa, sb                   !Sin of phase advance for A, B
    Real(rp) :: sqrt_beta_cbar_prelim(2), & !Temporarily holds two values
         inv_gamma_cbar_prelim(2)         !Preliminary values
    real(rp) :: j_amp_a(2), j_amp_b(2)   !Amplitude for both BPMs in a pair
    type(processor_analysis) temp(2)
    type (twiss_parameters), pointer :: ring_a, ring_b,dat_a, dat_b
    type (processor_analysis), pointer :: loc_1, loc_2
    integer :: begin, endpt
    type self_consistent  
       type (processor_analysis), allocatable :: loc(:)
       real(rp) :: j_amp_ave(2) 
    end type self_consistent

    type loc_pointer
       type (processor_analysis), pointer :: loc
    end type loc_pointer

    type twiss_pointer
       type (twiss_parameters), pointer :: twiss
    end type twiss_pointer

    !Pointers to ring%loc for the two bpms in the pair
    type (loc_pointer) :: ring_p(2)
    type (twiss_pointer) :: dat_twiss(2), & 
                                     !Contains data_struc%loc(jm)%a or b
    ring_twiss(2)                    !Same for ring(ik), a or b mode
                           !These pointers are for iterating between a and b
                           !modes with less repeated code.
    !Ring contains values calculated from BPMs with known spacing.
    type (self_consistent), allocatable, target :: ring(:) 

    ring_a => null()
    ring_b => null()
    dat_a => null()
    dat_b => null()
    loc_1 => null()
    loc_2 => null()


    if (.not. bpm_pairs(1)%has_one) return   !Does not continue if there are
    !not BPMs with known spacing
    n_ring = bpm_pairs(1)%number
    allocate (ring(n_ring))

    !Calculate ring parameters from each paired detector separately:
    do i = 1, n_ring
       allocate (ring(i)%loc(NUM_BPMS))
       bpm1 = bpm_pairs(i)%bpm_pntr(1)
       bpm2 = bpm_pairs(i)%bpm_pntr(2)

       if (abs(data_struc%loc(bpm2)%a%d_phase_adv) < 0.01) cycle
       if (abs(data_struc%loc(bpm2)%b%d_phase_adv) < 0.01) cycle

       ring(i)%loc = data_struc%loc

       ring_p(1)%loc => ring(i)%loc(bpm1)
       ring_p(2)%loc => ring(i)%loc(bpm2) 

       !
       !Compute Beta (page 15) +
       !
       call ring_beta(ring_p(1)%loc%a, ring_p(2)%loc%a, i)
       call ring_beta(ring_p(1)%loc%b, ring_p(2)%loc%b, i)

       do q = 1,2

          !
          !Compute Alphas
          !
          call alpha(ring_p(q)%loc%a,ring_p(2)%loc%a%d_phase_adv, i, q)
          call alpha(ring_p(q)%loc%b,ring_p(2)%loc%b%d_phase_adv, i, q)
          !
          !Compute Inv Gamma Cbars (pg 17)+
          !
          call inv_gamma_cbar(ring_p(q)%loc, 1)
          call inv_gamma_cbar12(ring_p(q)%loc)
          call inv_gamma_cbar(ring_p(q)%loc, 2)

          !
          !Compute Inv Gamma Cbar (2,1) (pg 14)
          !
          ca = cos(ring(i)%loc(bpm2)%a%d_phase_adv)
          cb = cos(ring(i)%loc(bpm2)%b%d_phase_adv)
          sa = sin(ring(i)%loc(bpm2)%a%d_phase_adv)
          sb = sin(ring(i)%loc(bpm2)%b%d_phase_adv)

          !Is different for BPMs 1 and 2--can't use ring_p?
          ring_p(1)%loc%inv_gamma_cbar(2,1) =  &
               (cb * ring_p(1)%loc%inv_gamma_cbar(1,1) - &
               sb * ring_p(1)%loc%inv_gamma_cbar(1,2) - &
               ca * ring_p(2)%loc%inv_gamma_cbar(1,1)) / sa

          ring(i)%loc(bpm2)%inv_gamma_cbar(2,1) =  &
               (ca * ring(i)%loc(bpm2)%inv_gamma_cbar(2,2) - &
               sa * ring(i)%loc(bpm2)%inv_gamma_cbar(1,2) - &
               cb * ring(i)%loc(bpm1)%inv_gamma_cbar(2,2)) / sb

          !
          !Checking Inv Gamma Cbar Consistency
          !
          !*Check this
          ring(i)%loc(bpm1)%inv_gamma_cbar_check(1) =  &
               (ca * ring(i)%loc(bpm2)%inv_gamma_cbar(1,2) + &
               sa * ring(i)%loc(bpm2)%inv_gamma_cbar(2,2) - &
               sb * ring(i)%loc(bpm1)%inv_gamma_cbar(1,1) - &
               cb * ring(i)%loc(bpm1)%inv_gamma_cbar(1,2))

          ring(i)%loc(bpm1)%inv_gamma_cbar_check(2) =  &
               (-sa * ring(i)%loc(bpm2)%inv_gamma_cbar(1,1) + &
               ca * ring(i)%loc(bpm2)%inv_gamma_cbar(2,1) - &
               cb * ring(i)%loc(bpm1)%inv_gamma_cbar(2,1) + &
               sb * ring(i)%loc(bpm1)%inv_gamma_cbar(2,2))

          !
          !Gamma (pg 16) +
          !
          ring_p(q)%loc%gamma = sqrt(1.0 / (1.0 + &
               (ring_p(q)%loc%inv_gamma_cbar(1,1) * &
               ring_p(q)%loc%inv_gamma_cbar(2,2) - &
               ring_p(q)%loc%inv_gamma_cbar(1,2) * &
               ring_p(q)%loc%inv_gamma_cbar(2,1))))

          !
          !Compute Cbar
          !
          ring_p(q)%loc%cbar =  &
               ring_p(q)%loc%inv_gamma_cbar * ring_p(q)%loc%gamma

          !
          !Compute J_amp (aka the Action) (pg 16) +
          !**Not correct, seems to be off by some scalar
          j_amp_a(q) = 1. * &
               (ring_p(q)%loc%a%magnitude2(1) / &
               ring_p(q)%loc%gamma**2 / ring_p(q)%loc%a%beta)

          j_amp_b(q) = 1. * &
               ring_p(q)%loc%b%magnitude2(2) / &
               (ring_p(q)%loc%gamma**2 * ring_p(q)%loc%b%beta)
Print *, "j_amp a: ", j_amp_a(q), " b: ", j_amp_b(q)

       enddo    !End calculations for BPMs 1, 2 (using q)

       !
       !Compute J_amp_ave +
       !
       ring(i)%j_amp_ave(1) = (j_amp_a(1) + j_amp_a(2)) / 2.0
       ring(i)%j_amp_ave(2) = (j_amp_b(1) + j_amp_b(2)) / 2.0


       !
       !Compute parameters for rest of ring
       !
       do j = 1, NUM_BPMS

          !
          !Compute Gamma**2 Beta (pg 17) +
          !
          !(*)_*
!          if (i==3) then
             ring(i)%loc(j)%a%gam2_beta = ring(i)%loc(j)%a%magnitude2(1)/ &
                  ring(i)%j_amp_ave(1)
             ring(i)%loc(j)%b%gam2_beta = ring(i)%loc(j)%b%magnitude2(2)/ &
                  ring(i)%j_amp_ave(2)
!          end if
             !
             !Compute Inv Gamma Cbar (pg 17) +
             !
             ring(i)%loc(j)%inv_gamma_cbar(2,2) = &
                  ring(i)%loc(j)%a%inv_gamma_cbar_sqrt_betas(2,2) * &
                  sqrt(ring(i)%loc(j)%a%gam2_beta / &
                  ring(i)%loc(j)%b%gam2_beta)

             ring(i)%loc(j)%inv_gamma_cbar(1,1) = &
                  ring(i)%loc(j)%b%inv_gamma_cbar_sqrt_betas(1,1) * &
                  sqrt(ring(i)%loc(j)%b%gam2_beta / &
                  ring(i)%loc(j)%a%gam2_beta)

             inv_gamma_cbar_prelim(1) = &
                  ring(i)%loc(j)%a%inv_gamma_cbar_sqrt_betas(1,2) * &
                  sqrt(ring(i)%loc(j)%a%gam2_beta / &
                  ring(i)%loc(j)%b%gam2_beta)
             inv_gamma_cbar_prelim(2) = &
                  ring(i)%loc(j)%b%inv_gamma_cbar_sqrt_betas(1,2) * &
                  sqrt(ring(i)%loc(j)%b%gam2_beta / &
                  ring(i)%loc(j)%a%gam2_beta)

             ring(i)%loc(j)%inv_gamma_cbar(1,2) = &
                  (inv_gamma_cbar_prelim(1)+inv_gamma_cbar_prelim(2)) / 2.0

             !
             !Compute sqrt(beta) cbar prelim
             !
             sqrt_beta_cbar_prelim(1) = &
                  - sqrt(data_struc%loc(j)%a%magnitude2(2) / &
                  ring(i)%j_amp_ave(1))*cos(0.-data_struc%loc(j)%a%d_delta)

             sqrt_beta_cbar_prelim(2) = &
                  -sqrt(data_struc%loc(j)%b%magnitude2(1) / &
                  ring(i)%j_amp_ave(2)) * cos(data_struc%loc(j)%b%d_delta)

             !
             !Compute sqrt(beta) cbar
             !
             ring(i)%loc(j)%sqrt_beta_cbar(1,1) = &
                  sqrt(data_struc%loc(j)%b%magnitude2(1) / &
                  ring(i)%j_amp_ave(2)) * sin(data_struc%loc(j)%b%d_delta)

             ring(i)%loc(j)%sqrt_beta_cbar(1,2) = &
                  (sqrt_beta_cbar_prelim(1) + sqrt_beta_cbar_prelim(2)) / 2.0

             ring(i)%loc(j)%sqrt_beta_cbar(2,2) = &
                  sqrt(data_struc%loc(j)%a%magnitude2(2) / &
                  ring(i)%j_amp_ave(1)) * sin(0.-data_struc%loc(j)%a%d_delta)
 !         end if
       end do
    end do

    !
    !Copy averages from rings to data struct
    !

  !this approach cannot work: nothing in this subroutine is file-dependent
!    if (data(1)%noise/data(2)%noise < 2 ) then
!       begin = 2
!       Print *, "Using file 2 only: less noise"
!    else if (data(2)%noise/data(1)%noise < 2 ) then
!       endpt = 1
!       Print *, "Using file 1 only: less noise"
!    end if


!***This line temporarily disables second BPM pair--should not be a permanent feature!
!n_ring = 1
!To ignore data from a bpm pair, change parameters for ik.
    do ik = 1, n_ring
       do jm = 1, NUM_BPMS
          ring_twiss(1)%twiss => ring(ik)%loc(jm)%a
          ring_twiss(2)%twiss => ring(ik)%loc(jm)%b
          dat_twiss(1)%twiss => data_struc%loc(jm)%a
          dat_twiss(2)%twiss => data_struc%loc(jm)%b

          do r = 1, 2
             !This loop calculates some values that depend upon mode.
             !r=1 does A mode, r=2 does B mode
!             data_struc%j_amp_ave(r,ik) = ring(ik)%j_amp_ave(r)
             !
             !Compute Gamma**2 Beta
             !
             !(*)_*
             dat_twiss(r)%twiss%gam2_beta = &
                  dat_twiss(r)%twiss%gam2_beta + &
                  (ring_twiss(r)%twiss%gam2_beta / n_ring )

             dat_twiss(r)%twiss%j_amp(r) = 1000. * &
                  dat_twiss(r)%twiss%magnitude2(r) / &
                  dat_twiss(r)%twiss%gam2_beta
          end do


          !Sum j_amp_ave for A and B modes
          data_struc%j_amp_ave(1) = data_struc%j_amp_ave(1) + &
               (ring(ik)%j_amp_ave(1) / n_ring )
          data_struc%j_amp_ave(2) = data_struc%j_amp_ave(2) + &
               (ring(ik)%j_amp_ave(2) / n_ring )





          !
          !Compute sqrt(beta) cbar
          !
!          call avg(data_struc%loc(jm)%sqrt_beta_cbar(1,1), &
!               ring(ik)%loc(jm)%sqrt_beta_cbar(1,1), n_ring)
!          call avg(data_struc%loc(jm)%sqrt_beta_cbar(1,2), &
!               ring(ik)%loc(jm)%sqrt_beta_cbar(1,2), n_ring)
!          call avg(data_struc%loc(jm)%sqrt_beta_cbar(2,2), &
!               ring(ik)%loc(jm)%sqrt_beta_cbar(2,2), n_ring)

          !
          !Compute Inv Gamma Cbar
          !
          call avg(data_struc%loc(jm)%inv_gamma_cbar(1,1), &
               ring(ik)%loc(jm)%inv_gamma_cbar(1,1), n_ring)
          call avg(data_struc%loc(jm)%inv_gamma_cbar(1,2), &
               ring(ik)%loc(jm)%inv_gamma_cbar(1,2), n_ring)
          call avg(data_struc%loc(jm)%inv_gamma_cbar(2,2), &
               ring(ik)%loc(jm)%inv_gamma_cbar(2,2), n_ring)

          !
          !Compute the special parameters
          !

          !
          !Copies beta into data_struc if it was calculated
          !
          data_struc%loc(jm)%a%beta = ring(ik)%loc(jm)%a%beta + &
               data_struc%loc(jm)%a%beta
          data_struc%loc(jm)%b%beta = ring(ik)%loc(jm)%b%beta + &
               data_struc%loc(jm)%b%beta

          !
          !Copy alpha into data_struc if available
          !
          data_struc%loc(jm)%a%alpha = ring(ik)%loc(jm)%a%alpha + &
               data_struc%loc(jm)%a%alpha
          data_struc%loc(jm)%b%alpha = ring(ik)%loc(jm)%b%alpha + &
               data_struc%loc(jm)%b%alpha

          !
          !Compute Inv Gamma Cbar (2,1)
          !
          data_struc%loc(jm)%inv_gamma_cbar(2,1)  = &
               min( data_struc%loc(jm)%inv_gamma_cbar(2,1), &
               ring(ik)%loc(jm)%inv_gamma_cbar(2,1)) + &
               max( data_struc%loc(jm)%inv_gamma_cbar(2,1), &
               ring(ik)%loc(jm)%inv_gamma_cbar(2,1)) 

          !
          !Checking Inv Gamma Cbar Consistency
          !*Check this
          !
!          DO count = 1, 2                               
!             data_struc%loc(jm)%inv_gamma_cbar_check(count)  = &
!                  min( data_struc%loc(jm)%inv_gamma_cbar_check(count), &
!                  ring(ik)%loc(jm)%inv_gamma_cbar_check(count)) + &
!                  max( data_struc%loc(jm)%inv_gamma_cbar_check(count), &
!                  ring(ik)%loc(jm)%inv_gamma_cbar_check(count)) 
!          ENDDO

          !
          !Gamma
          !
          data_struc%loc(jm)%gamma  = ring(ik)%loc(jm)%gamma + &
               data_struc%loc(jm)%gamma
          !
          !Cbar
          !
          data_struc%loc(jm)%cbar  = ring(ik)%loc(jm)%cbar + &
               data_struc%loc(jm)%cbar

       end do
    end do
    deallocate (ring)
  end subroutine calculate_with_known_spacing



  subroutine output (data)

    implicit none
    type(data_set)data(*)
    integer :: i,j, openstatus
    integer :: n_pairs                  !Number of BPM pairs
    integer :: bpm
    n_pairs = bpm_pairs(1)%number


!51  if (.not. outfile) then
       open (unit = 27, file = "./data/mia.out", &
            action = "write", position = "rewind",&
            iostat = openstatus)
       if (openstatus > 0) print *, "*** Cannot open output file ***",&
            openstatus
!*For easy analysis of beta/cbar error with a plethora of files:
      open (unit = 29, file = "./data/beta.out", &
            action = "write", position = "rewind",&
            iostat = openstatus)
       if (openstatus > 0) print *, "*** Cannot open output file ***",&
            openstatus
!       open (unit = 31, file = "./data/out.cbar", &
!            action = "write", position = "rewind",&
!            iostat = openstatus)
!       if (openstatus > 0) print *, "*** Cannot open output file ***",&
!*            openstatus
!    else
!       open (unit = 27, file = trim(outname), &
!            action = "write", position = "rewind",&
!            iostat = openstatus)
!       if (openstatus > 0) then
!          print *, "*** Cannot open output file ***", openstatus
!          Print *, "Enter a new output filename:"
!          Read *, outfile
!          goto 51
!       endif
!    endif


9   format(1x,a5,a1,2x, f11.7, a3, f11.7, a7, 5x, f11.7, a3, f11.7)
47  format(5x, f11.7, a2,2x, f11.7)
11  format(2x,a5,a1,2x, f11.7, a2,2x, f11.7)
12  format(3(4x, a8, 7x))
13  format(5x, e14.7, a2,2x, e14.7)
17  format(a10, 3x, f11.7, 1x,",", 1x, f11.7)
18  format(4x, a20, 4x, 20x, a5)
34  format(1x, 3(e14.7, a2, 3x))
56  format(a10,3x,f10.7)

    write (27, *) "A mode filename: ", trim(xfile)
    write (27, *) "B mode filename: ", trim(yfile)
    write (27, *) "Number of BPMs: ", NUM_BPMS
    write (27, *) "Number of Turns: ", NUM_TURNS
    write (27, *) "Frequency of machine: ", FREQ
    write (27, *) "Tune A: ", data_struc%tune(1)
    write (27, *), "Tune B: ", data_struc%tune(2)
    write (27, *), "Phi(t) A = ", data_struc%phi_t(1)
    write (27, *), "Phi(t) B = ", data_struc%phi_t(2)
    write (27, *), " Values of Lambda"
!Might want to start beta file with J_amp instead of tune
    write (29,"(4x, a22)") "<J>"
    write (29, 13) data_struc%j_amp_ave(1),",", &
         data_struc%j_amp_ave(2)
!
!*    write (29, *) "Tune A: ", data_struc%tune(1)
!   write (29, *), "Tune B: ", data_struc%tune(2)
!   write (31, *) "Tune A: ", data_struc%tune(1)
!*   write (31, *), "Tune B: ", data_struc%tune(2)

    do i=1, 2*NUM_BPMS
       write (27, *), data(1)%lambda(i), " , ", data(2)%lambda(i)
    enddo

    do i=1, n_pairs
       do j=1, 2
          bpm = bpm_pairs(i)%bpm_pntr(j)
          write (27, *) data_struc%proc(bpm)%label
          write (27, 17) "Beta: ", data_struc%loc(bpm)%a%beta, &
               data_struc%loc(bpm)%b%beta
          write (27, 17) "Alpha: ", data_struc%loc(bpm)%a%alpha, &
               data_struc%loc(bpm)%b%alpha
          write (27, 56) "Gamma: ", data_struc%loc(bpm)%gamma
          write (27, 56) "Cbar (1,1): ", data_struc%loc(bpm)%cbar(1,1)
          write (27,56) "(1,2):", data_struc%loc(bpm)%cbar(1,2)
          write (27,56) "(2,1):", data_struc%loc(bpm)%cbar(2,1)
          write (27,56) "(2,2):", data_struc%loc(bpm)%cbar(2,2)
       enddo
    enddo

    write (27, 18) "Delta phase advance", "Phi"
    do i = 1, NUM_BPMS
       write (27, 9) data_struc%proc(i)%label,":", &
            data_struc%loc(i)%a%d_phase_adv, &
            ",", data_struc%loc(i)%b%d_phase_adv, "||", &
            data_struc%loc(i)%a%phi, ",", data_struc%loc(i)%b%phi
    enddo

!Haven't used this data...disabled output for now.
!-    write (27, "(8x, a12)") "d_delta"
!    do i = 1, NUM_BPMS
!       write (27, 11) data_struc%loc(i)%a%d_delta, &
!            ",", data_struc%loc(i)%b%d_delta
!-    enddo

    write (27,"(10x, a12)") "Gamma^2 Beta"
    do i = 1, NUM_BPMS
       write (27,11)  data_struc%proc(i)%label,":", &
            data_struc%loc(i)%a%gam2_beta, ",", &
            data_struc%loc(i)%b%gam2_beta
    enddo

    write (29,"(10x, a12)") "Gamma^2 Beta"
    do i = 1, NUM_BPMS
       write (29,47) data_struc%loc(i)%a%gam2_beta, ",", &
            data_struc%loc(i)%b%gam2_beta
    enddo

    write (27,"(20x, a15)") "1/gamma * Cbar"
    write (27,12) "(1,1)", "(1,2)", "(2,2)"

    do i = 1, NUM_BPMS
       write (27,34) data_struc%loc(i)%inv_gamma_cbar(1,1), ",", &
            data_struc%loc(i)%inv_gamma_cbar(1,2), ",", &
            data_struc%loc(i)%inv_gamma_cbar(2,2)
    enddo
!*    write (31,"(20x, a15)") "1/gamma * Cbar"
!    write (31,12) "(1,1)", "(1,2)", "(2,2)"
!
!    do i = 1, NUM_BPMS
!       write (31,34) data_struc%loc(i)%inv_gamma_cbar(1,1), ",", &
!            data_struc%loc(i)%inv_gamma_cbar(1,2), ",", &
!            data_struc%loc(i)%inv_gamma_cbar(2,2)
!*    enddo

    write (27,"(4x, a22)") "<J>"
!    do i=1, 2
!       write (27, "(7x, a5, i2)") "File ", i
       write (27, 13) data_struc%j_amp_ave(1),",", &
            data_struc%j_amp_ave(2)
!    enddo


    close(27)
    close(29)
!*    close(31)

  end subroutine output

  subroutine fftOut(data)
    !
    !Outputs fft information to a file
    !
    type(data_set) :: data(*)
    integer :: i,j


    open (unit = 9, file = "./data/fft-a.dat", &
         action = "write", position = "rewind")

    !*Only outputs half; to get all columns change limit to 2*NUM_BPMS
    do i=1, NUM_BPMS
       write(9,*), "Column ", i
       do j=1, NUM_TURNS

          write (9,*) j,"  ,  ",  data(1)%spectrum(j,i)

       end do
       write (9,*) ""
       write (9,*) ""
       write (9,*) ""

    end do
    close(9)

    open (unit = 19, file = "./data/fft-b.dat", &
         action = "write", position = "rewind")

    !*Only outputs half; to get all columns change limit to 2*NUM_BPMS
    do i=1, NUM_BPMS
       write(19,*), "Column ", i
       do j=1, NUM_TURNS

          write (19,*) j,"  ,  ", data(2)%spectrum(j,i)

       end do
       write (19,*) ""
       write (19,*) ""
       write (19,*) ""

    end do
 
    close(19)

  end subroutine fftOut

  subroutine cesrv_out()
    real(rp), allocatable :: cbar_11(:), cbar_12(:),cbar_22(:)
    character(1), allocatable :: valid(:) !True if number was calculated,
                                          !false if not (ex. for cbar)
    real(rp), allocatable :: writeMe(:,:)
    logical :: known_spacing  !temporary
    integer :: i, openstatus
    character(1),allocatable :: allTrueIsm(:) !Contains all true; for variables
                                          !MIA can always calculate
    integer :: cesrv                    !Unit number for cesrv file
    integer, allocatable :: eleNum(:)
    integer :: bpm1, bpm2
    integer, allocatable :: driftSpaces(:) !ele # of detectors separated by
                                           !only a drift space
    cesrv = 51 

    allocate(cbar_11(2*bpm_pairs(1)%number))
    allocate(cbar_12(2*bpm_pairs(1)%number))
    allocate(cbar_22(2*bpm_pairs(1)%number))
    allocate(eleNum(NUM_BPMS))
    allocate(writeMe(2,NUM_BPMS))
    allocate (driftSpaces(2*bpm_pairs(1)%number))

    known_spacing = .false.
    open (unit = cesrv, file = "./data/cesrv.dat", &
         action = "write", position = "rewind",&
         iostat = openstatus)
    if (openstatus > 0) print *, "*** Cannot open output file ***",&
         openstatus
   
    write (cesrv,*) "&DATA_PARAMETERS"
    write (cesrv,*) "  file_type = 'ALL DATA'"
!    write (cesrv,*) "  lattice = 'CTA_2085MEV_20090516'"
    write (cesrv,*) "  lattice = 'CHESS_20090225'"
    write (cesrv,*) "  comment = ", "'MIA'"
    write (cesrv,*) "/"
    Write (cesrv,*) ""

    write (cesrv,*) "&PHASE_PARAMETERS"
    write (cesrv,*) "  species = 1"
    write (cesrv,*) "  horiz_freq = ", data_struc%intTune(1)+data_struc%tune(1)/FREQ
    write (cesrv,*) "  vert_freq = ", data_struc%intTune(2)+data_struc%tune(2)/FREQ
    write (cesrv,*) "/"
    Write (cesrv,*) ""

    Write (cesrv,*) "&all_data"
    Write (cesrv, *) "!    det       X            Y        Valid?"

    close(cesrv)

    do i=1, num_bpms
       eleNum(i) = data_struc%proc(i)%eleNum
    end do

!Added arbitrary phase offset for comparison with CESRV--remove these later
    writeMe(1,:) = data_struc%loc(:)%a%phi +0.88226148
    writeMe(2,:) = data_struc%loc(:)%b%phi +0.14484130

    call writeArray(eleNum, writeMe, "phase")
    writeMe(1,:) = data_struc%loc(:)%a%gam2_beta
    writeMe(2,:) = data_struc%loc(:)%b%gam2_beta
    call writeArray(eleNum, writeMe, "beta")

    writeMe(1,:) = data_struc%loc(:)%inv_gamma_cbar(1,1)
    call writeVector(eleNum, writeMe(1,:), "cbar11")
    writeMe(1,:) = data_struc%loc(:)%inv_gamma_cbar(1,2)
    call writeVector(eleNum, writeMe(1,:), "cbar12")
    writeMe(1,:) = data_struc%loc(:)%inv_gamma_cbar(2,2)
    call writeVector(eleNum, writeMe(1,:), "cbar22")

    open (unit = cesrv, file = "./data/cesrv.dat", &
         action = "write", position = "append",&
         iostat = openstatus)
    if (.not. openstatus == 0) print *, "*** Cannot open output file ***",&
         openstatus


    write (cesrv,*) "/"
    close(cesrv)
  end subroutine cesrv_out

  subroutine writeArray(eleNum, vector, name)

    real(rp) :: vector(:,:)
    integer :: fileNum, i, eleNum(:)
    character(*) :: name
    integer :: openstatus

    fileNum = 55
98  format (a,a1,i4.1,a4,f12.7, 2x, f12.7, 2x, a1)

    open (unit = fileNum, file = "./data/cesrv.dat", &
         action = "write", position = "append",&
         iostat = openstatus)
    if (.not. openstatus == 0) print *, "*** Cannot open output file ***",&
         openstatus

    do i=1, size(eleNum)
       write (fileNum, 98) name,"(", eleNum(i), ") = ", vector(1,i),&
            vector(2,i), "T"
    end do



  end subroutine writeArray

  subroutine writeVector(eleNum, vector, name)

    real(rp) :: vector(:)
    integer :: fileNum, i, eleNum(:)
    character(*) :: name
    integer :: openstatus

    fileNum = 55
98  format (a,a1,i3.1,a4,e12.4, 2x, a1)

    open (unit = fileNum, file = "./data/cesrv.dat", &
         action = "write", position = "append",&
         iostat = openstatus)
    if (.not. openstatus == 0) print *, "*** Cannot open output file ***",&
         openstatus

    do i=1, size(eleNum)
       write (fileNum, 98) name,"(", eleNum(i), ") = ", vector(i), "T"
    end do



  end subroutine writeVector

end module orbit_mia
