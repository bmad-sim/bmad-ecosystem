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

  subroutine bpm_ops(data,iset, nset, first_run)
    !
    !Calls functions that locate and match BPMs
    !Remove if bpm numbers, pairs, etc. are not obtained
    !from input files (ex. from CESRV or another program)
    !
    type(data_set) data(*)
    integer :: iset, nset
    logical :: first_run
 
    call locate_bpm (iset, data(iset))

    if (first_run .and. iset == 1) then
       call find_L(nset,data)
    else
       call match_processors (data)
    end if

  end subroutine bpm_ops


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

    call tune(data) 
    Print *, "set num a ", data_struc%set_num_a
    Print *, "set num b ", data_struc%set_num_b

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
    end do

    !
    !Calculates phase adv around ring in clockwise direction(pg15 Delta Phi)
    !

    !*****************************************************************
    !Phase begins is 0 at BPM 0W.
    !*****************************************************************

    do bpm = 1, data(1)%bpmproc !Find BPMs 0E and 0W
       if (data_struc%proc(bpm)%number == 0) then
          if (data_struc%proc(bpm)%is_west) then
             west = bpm
            ! print *, "Found BPM 0W at ", bpm
          else
             east = bpm
            ! print *, "Found BPM 0E at ", bpm
          endif
       endif
    enddo

    !Initially set east bpm's phase advance to 0
    data_struc%loc(west)%a%d_phase_adv = 0.
    data_struc%loc(west)%b%d_phase_adv = 0.
    data_struc%loc(west)%a%gam2_beta_ratio = 1.0
    data_struc%loc(west)%b%gam2_beta_ratio = 1.0

    call phase_adv(data_struc%loc(east)%a, data_struc%loc(west)%a, 1)
    call phase_adv(data_struc%loc(east)%b, data_struc%loc(west)%b, 2) 

    data_struc%loc(east)%a%d_phase_adv = &
         mod((data_struc%loc(east)%a%d_phase_adv - &
         data_struc%phi_t(1))+7.0*pi, 2.0*pi) -pi
    data_struc%loc(east)%b%d_phase_adv = &
         mod((data_struc%loc(east)%b%d_phase_adv - &
         data_struc%phi_t(2)+7.0*pi),2.0*pi) - pi

    !BPMs must be in clockwise order (first BPM is west, last BPM is east)
    !Gives an error if they are not.
    if (data_struc%proc(1)%is_west) then
       do n_bpm = 2, data(1)%bpmproc
          !Calculate phase advance clockwise around the ring (skipping 0W):
          if (data_struc%proc(n_bpm)%is_west) then 

             call phase_adv(data_struc%loc(n_bpm)%a, &
                  data_struc%loc(n_bpm-1)%a, 1)

             call phase_adv(data_struc%loc(n_bpm)%b, &
                  data_struc%loc(n_bpm-1)%b, 2)
          else
             !Calculate phase counterclockwise for east BPMs (skipping 0E):
             if (.not. data_struc%proc(data(1)%bpmproc)%is_west) then
               ! Print *, "INVERT"
                do n_bpm_e = data(1)%bpmproc-1, n_bpm, -1
                   call phase_adv(data_struc%loc(n_bpm_e)%a, &
                        data_struc%loc(n_bpm_e+1)%a, 1)
                   call phase_adv(data_struc%loc(n_bpm_e)%b, &
                        data_struc%loc(n_bpm_e+1)%b, 2)
                enddo

                Exit
             else
                print *, "ERROR: Found BPMs that are not east, but last BPM is west."
             endif
          end if
       end do
    else
       Print *, "ERROR: BPMs do not start in the west."
    endif
    
    sum_a = 0.0
    sum_b = 0.0
    !Assigns phi by summing delta phase advances
    do i = 1, NUM_BPMS
       !Reset phi to 0 when switching from west to east

       if (data_struc%proc(i)%is_west) then
          sum_a = sum_a + data_struc%loc(i)%a%d_phase_adv
          data_struc%loc(i)%a%phi = sum_a

          sum_b = sum_b + data_struc%loc(i)%b%d_phase_adv
          data_struc%loc(i)%b%phi = sum_b
       else 
          if (.not. data_struc%proc(i)%is_west .and. &
               data_struc%proc(i-1)%is_west) then
             sum_a = 0.0
             sum_b = 0.0
             numWest = i
          endif
          sum_a = sum_a + data_struc%loc(NUM_BPMS-(i-numWest))%a%d_phase_adv
          data_struc%loc(NUM_BPMS-(i-numWest))%a%phi = sum_a

          sum_b = sum_b + data_struc%loc(NUM_BPMS-(i-numWest))%b%d_phase_adv
          data_struc%loc(NUM_BPMS-(i-numWest))%b%phi = sum_b
       endif
    end do

!    do i = 1, data(1)%bpmproc
!    end do

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


  !**********************************************************************
  !The next four subroutines are for use in convert_data_from_pi_matrix
  !**********************************************************************

  subroutine pi_calc(data_snum, twiss, col_p, col_n, count, n_bpm)

    !
    !Calculates numer, denom, and magnitude^2
    !

    type(data_set) data_snum      !Data from file
    INTEGER :: col_p, col_n, &    !Positive or negative col
         count, n_bpm             !Counter, bpm number
    type(twiss_parameters) twiss  !Location of values in data structure

    !Numer is lambda*pi of positive col, denom same of negative col.
    !Numer is x, denom is y
    !Count is 1 or 2.
    twiss%numer(count)= &
         data_snum%lambda(col_p) &
         * data_snum%pi_mat(2*n_bpm+(count-2),col_p)

    twiss%denom(count)= &
         data_snum%lambda(col_n) &
         * data_snum%pi_mat(2*n_bpm+(count-2),col_n)

    !Ratio is a graph option, but is never used elsewhere. What is its purpose?
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
    !Count should be 2 for A mode and 1         for B mode

    !?@@
    twiss%d_delta = mod(atan2(twiss%numer(count), twiss%denom(count)) &
         - atan2(twiss%numer(3-count), twiss%denom(3-count))&
         +7.0*pi+pi/2, 2.0*pi)-pi
    !Mod function maps d_delta from -pi to +pi

  end subroutine d_delta

  subroutine phase_adv(twiss, twiss_old, count)
    !
    !Finds phase advance and gamma^2 beta ratio (both pg 15)
    !
    type(twiss_parameters) twiss, & !Location of this numer, denom
         twiss_old                  !Location of previous number, denom
    INTEGER :: count                !Counter (1 for A, 2 for B mode)

    twiss%d_phase_adv =  &
         mod( atan2(twiss%numer(count), twiss%denom(count) ) - &
         atan2(twiss_old%numer(count), &
         twiss_old%denom(count)) +7.0 * pi , 2.0 * pi) - pi
    !Mod function maps d_phase_advance from -pi to +pi

    !@@ Inverted gam2_beta_ratio. It can be done either here
    !or in ring_beta.

    twiss%gam2_beta_ratio =  &
         twiss%magnitude2(count) / twiss_old%magnitude2(count)

  end subroutine phase_adv

  subroutine ring_beta(ring, ring2, i)
    !
    !Calculates beta for a pair of BPMs with known spacing.
    !Beta is accurate to 0.1 as compared with CESRV.
    !
    INTEGER :: i                           !Counter for BPM pair
    type(twiss_parameters) ring, ring2     !Ring(i)
    !First BPM:
    ring%beta = bpm_pairs(i)%length / &
         sin(ring2%d_phase_adv) * sqrt(1/ring2%gam2_beta_ratio)
    !Second BPM:
    ring2%beta = bpm_pairs(i)%length / &
         sin(ring2%d_phase_adv) * sqrt(ring2%gam2_beta_ratio)

  end subroutine ring_beta

  subroutine alpha(twiss_ring, bpm_l, twiss, i)
    type(twiss_parameters) twiss_ring, &  !Location of ring in data struct
         twiss                            !Location for phase advance
    real(rp) :: bpm_l                     !Distance between BPMs
    integer :: i                          !+-1
    !i changes the sign of the second term

    twiss_ring%alpha = twiss_ring%beta / bpm_l +i*1.0/tan(twiss%d_phase_adv)
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
         bpm1a, bpm2a, &          !Pair of BPMs from set_num_a
         bpm1b, bpm2b, &          !Pair of BPMs from set_num_b
!         set_num_a, set_num_b, &  !Takes set_num_* from data_struc
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

!    set_num_a = data_struc%set_num_a
!    set_num_b = data_struc%set_num_b

    if (.not. bpm_pairs(1)%has_one) return   !Does not continue if there are
    !not BPMs with known spacing

    n_ring = bpm_pairs(1)%number
    allocate (ring(n_ring))

    do i = 1, n_ring

       allocate (ring(i)%loc(NUM_BPMS))

       !?@@           
       !ring(i)%loc(i) = data_struc%loc(i)  !loc # is not = i.
       !Moved below          
       !
       ! Determine info about the known pair of BPMs first.
       !

       bpm1a = bpm_pairs(i)%bpm_pntr(1)                 !A mode BPM pair  
       bpm2a = bpm_pairs(i)%bpm_pntr(2)

       bpm1b = bpm_pairs(i)%bpm_pntr(1) !B mode BPM pair    
       bpm2b = bpm_pairs(i)%bpm_pntr(2)
       !BPM1A = BPM1B, BPM2A = BPM2B (why have A and B?)

       if (abs(data_struc%loc(bpm2a)%a%d_phase_adv) < 0.01) cycle
       if (abs(data_struc%loc(bpm2b)%b%d_phase_adv) < 0.01) cycle

       ring(i)%loc = data_struc%loc

       ring_p(1)%loc => ring(i)%loc(bpm1a)
       ring_p(2)%loc => ring(i)%loc(bpm2a) 

       !
       !Compute Beta (page 15) +
       !
       !Betas are accurate
       call ring_beta(ring_p(1)%loc%a, ring_p(2)%loc%a, i)
       call ring_beta(ring_p(1)%loc%b, ring_p(2)%loc%b, i)
       print *, "Beta B: ", ring_p(1)%loc%b%beta
       print *, "Beta B: ", ring_p(2)%loc%b%beta
       do q = 1,2


          !
          !Compute Alphas
          !
          !**Alphas are wrong
          call alpha(ring_p(q)%loc%a, bpm_pairs(i)%length, &
               ring(i)%loc(2)%a, -1)
          call alpha(ring_p(q)%loc%b, bpm_pairs(i)%length, &
               ring_p(2)%loc%b, -1)

          !
          !Compute Inv Gamma Cbars (pg 17)+
          !

          call inv_gamma_cbar(ring_p(q)%loc, 1)
    !      print *,"Cbar(1,1):",  ring_p(q)%loc%inv_gamma_cbar(1,1)
          call inv_gamma_cbar12(ring_p(q)%loc)
    !      print *,"Cbar(1,2):",  ring_p(q)%loc%inv_gamma_cbar(1,2)     
          call inv_gamma_cbar(ring_p(q)%loc, 2)
    !      print *,"Cbar(2,2):",  ring_p(q)%loc%inv_gamma_cbar(2,2)     

          !
          !Compute Inv Gamma Cbar (2,1) (pg 14)
          !

          ca = cos(ring(i)%loc(bpm2a)%a%d_phase_adv)
          cb = cos(ring(i)%loc(bpm2b)%b%d_phase_adv)
          sa = sin(ring(i)%loc(bpm2a)%a%d_phase_adv)
          sb = sin(ring(i)%loc(bpm2b)%b%d_phase_adv)

          !Is different for BPMs 1 and 2--can't use ring_p?
          ring_p(1)%loc%inv_gamma_cbar(2,1) =  &
               (cb * ring_p(1)%loc%inv_gamma_cbar(1,1) - &
               sb * ring_p(1)%loc%inv_gamma_cbar(1,2) - &
               ca * ring_p(2)%loc%inv_gamma_cbar(1,1)) / sa

          ring(i)%loc(bpm2a)%inv_gamma_cbar(2,1) =  &
               (ca * ring(i)%loc(bpm2a)%inv_gamma_cbar(2,2) - &
               sa * ring(i)%loc(bpm2a)%inv_gamma_cbar(1,2) - &
               cb * ring(i)%loc(bpm1a)%inv_gamma_cbar(2,2)) / sb


          !
          !Checking Inv Gamma Cbar Consistency
          !

          !*Check this
          ring(i)%loc(bpm1a)%inv_gamma_cbar_check(1) =  &
               (ca * ring(i)%loc(bpm2a)%inv_gamma_cbar(1,2) + &
               sa * ring(i)%loc(bpm2a)%inv_gamma_cbar(2,2) - &
               sb * ring(i)%loc(bpm1a)%inv_gamma_cbar(1,1) - &
               cb * ring(i)%loc(bpm1a)%inv_gamma_cbar(1,2))
          !          Print *, "Inv gamma cbar check: ", &
          !               ring(i)%loc(bpm1a)%inv_gamma_cbar_check(1)
          ring(i)%loc(bpm1a)%inv_gamma_cbar_check(2) =  &
               (-sa * ring(i)%loc(bpm2a)%inv_gamma_cbar(1,1) + &
               ca * ring(i)%loc(bpm2a)%inv_gamma_cbar(2,1) - &
               cb * ring(i)%loc(bpm1a)%inv_gamma_cbar(2,1) + &
               sb * ring(i)%loc(bpm1a)%inv_gamma_cbar(2,2))
          !          Print *, "Inv gamma cbar check: ", &
          !               ring(i)%loc(bpm1a)%inv_gamma_cbar_check(2)

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


          !?@@
          ring_p(q)%loc%cbar =  &
               ring_p(q)%loc%inv_gamma_cbar * ring_p(q)%loc%gamma

          !
          !Compute J_amp (aka the Action) (pg 16) +
          !In mm
          !

          j_amp_a(q) = 1000. * 1000. * &
               ring_p(q)%loc%a%magnitude2(1) / &
               ring_p(q)%loc%gamma**2 / ring_p(q)%loc%a%beta

          j_amp_b(q) = 1000. * 1000.* &
               ring_p(q)%loc%b%magnitude2(2) / &
               ring_p(q)%loc%gamma**2 / ring_p(q)%loc%b%beta

       enddo    !End calculations for BPMs 1, 2 (using q)
       !
       !Compute J_amp_ave +
       !

       ring(i)%j_amp_ave(1) = (j_amp_a(1) + j_amp_a(2)) / 2.0

       ring(i)%j_amp_ave(2) = (j_amp_b(1) + j_amp_b(2)) / 2.0

       !
       !Compute parameters for rest of ring
       !

       do j = 1, data(1)%bpmproc
       !***This is always true. Does it have a purpose or should it be .and.?
          if (j /= bpm1a .or. j /= bpm2a) then 

             !
             !Compute Gamma**2 Beta (pg 17) +
             !

             ring(i)%loc(j)%a%gam2_beta = ring(i)%loc(j)%a%magnitude2(1)/ &
                  ring(i)%j_amp_ave(1)*1000.*1000.
             ring(i)%loc(j)%b%gam2_beta = ring(i)%loc(j)%b%magnitude2(2)/ &
                  ring(i)%j_amp_ave(2)*1000.*1000.

             !
             !Compute Inv Gamma Cbar (pg 17) +
             !
             !Does not calculate (2,1)
             !These Cbars are accurate.
             ring(i)%loc(j)%inv_gamma_cbar(2,2) = &
                  ring(i)%loc(j)%a%inv_gamma_cbar_sqrt_betas(2,2) * &
                  sqrt(ring(i)%loc(j)%a%gam2_beta / ring(i)%loc(j)%b%gam2_beta)
             !@@ Was sqrt(1): needs to be beta(a)/beta(b)

             ring(i)%loc(j)%inv_gamma_cbar(1,1) = &
                  ring(i)%loc(j)%b%inv_gamma_cbar_sqrt_betas(1,1) * &
                  sqrt(ring(i)%loc(j)%b%gam2_beta / ring(i)%loc(j)%a%gam2_beta)

             inv_gamma_cbar_prelim(1) = &
                  ring(i)%loc(j)%a%inv_gamma_cbar_sqrt_betas(1,2) * &
                  sqrt(ring(i)%loc(j)%a%gam2_beta / ring(i)%loc(j)%b%gam2_beta)
             inv_gamma_cbar_prelim(2) = &
                  ring(i)%loc(j)%b%inv_gamma_cbar_sqrt_betas(1,2) * &
                  sqrt(ring(i)%loc(j)%b%gam2_beta / ring(i)%loc(j)%a%gam2_beta)

             ring(i)%loc(j)%inv_gamma_cbar(1,2) = &
                  (inv_gamma_cbar_prelim(1) + inv_gamma_cbar_prelim(2)) / 2.0

             !
             !Compute sqrt(beta) cbar prelim
             !

             sqrt_beta_cbar_prelim(1) = &
                  - sqrt(data_struc%loc(j)%a%magnitude2(2) / &
                  ring(i)%j_amp_ave(1)) * cos(0.-data_struc%loc(j)%a%d_delta)

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
          end if
       end do
    end do
    !
    !Copy averages from rings to data struct
    !
    !**Some of these don't actually do anything...

    do ik = 1, n_ring
       do jm = 1, data(1)%bpmproc
          ring_twiss(1)%twiss => ring(ik)%loc(jm)%a
          ring_twiss(2)%twiss => ring(ik)%loc(jm)%b
          dat_twiss(1)%twiss => data_struc%loc(jm)%a
          dat_twiss(2)%twiss => data_struc%loc(jm)%b
    
!          dat_twiss(1) => data_struc%loc(jm)%a
!          dat_twiss(2) => data_struc%loc(jm)%b
!          ring_twiss(1) => ring(ik)%loc(jm)%a
!          ring_twiss(2) => ring(ik)%loc(jm)%b

          do r = 1, 2
             !This loop calculates some values that depend upon mode.
             !r=1 does A mode, r=2 does B mode
             data_struc%j_amp_ave(r,ik) = ring(ik)%j_amp_ave(r)
             !
             !Compute Gamma**2 Beta
             !
             !This takes half of ring_twiss(r)%gam2_beta for each r 
             !(averages gam2_beta from the two files)
             !Final value is accurate to about 0.01
             dat_twiss(r)%twiss%gam2_beta = &
                  dat_twiss(r)%twiss%gam2_beta + &
                  ring_twiss(r)%twiss%gam2_beta / n_ring

          end do

          !
          !Compute sqrt(beta) cbar
          !
          call avg(data_struc%loc(jm)%sqrt_beta_cbar(1,1), &
               ring(ik)%loc(jm)%sqrt_beta_cbar(1,1), n_ring)
          call avg(data_struc%loc(jm)%sqrt_beta_cbar(1,2), &
               ring(ik)%loc(jm)%sqrt_beta_cbar(1,2), n_ring)
          call avg(data_struc%loc(jm)%sqrt_beta_cbar(2,2), &
               ring(ik)%loc(jm)%sqrt_beta_cbar(2,2), n_ring)

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
          !Compute Beta
          !
          !***Does not compute beta.
!          data_struc%loc(jm)%a%beta  = &
!               min( data_struc%loc(jm)%a%beta, &
!               ring(ik)%loc(jm)%a%beta) + &
!               max( data_struc%loc(jm)%a%beta, &
!               ring(ik)%loc(jm)%a%beta) 

!          data_struc%loc(jm)%b%beta  = &
!               min( data_struc%loc(jm)%b%beta, &
!               ring(ik)%loc(jm)%b%beta) + &
!               max( data_struc%loc(jm)%b%beta, &
!               ring(ik)%loc(jm)%b%beta) 

          data_struc%loc(jm)%a%beta = ring(ik)%loc(jm)%a%beta + &
               data_struc%loc(jm)%a%beta
          data_struc%loc(jm)%b%beta = ring(ik)%loc(jm)%b%beta + &
               data_struc%loc(jm)%b%beta

          !
          !Compute Alphas
          !
          !***What is the purpose of min and max?
!          data_struc%loc(jm)%a%alpha  = &
!               min( data_struc%loc(jm)%a%alpha, &
!               ring(ik)%loc(jm)%a%alpha) + &
!               max( data_struc%loc(jm)%a%alpha, &
!               ring(ik)%loc(jm)%a%alpha) 

!          data_struc%loc(jm)%b%alpha  = &
!               min( data_struc%loc(jm)%b%alpha, &
!               ring(ik)%loc(jm)%b%alpha) + &
!               max( data_struc%loc(jm)%b%alpha, &
!               ring(ik)%loc(jm)%b%alpha) 

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
          !


          !*Check this
          DO count = 1, 2                               
             data_struc%loc(jm)%inv_gamma_cbar_check(count)  = &
                  min( data_struc%loc(jm)%inv_gamma_cbar_check(count), &
                  ring(ik)%loc(jm)%inv_gamma_cbar_check(count)) + &
                  max( data_struc%loc(jm)%inv_gamma_cbar_check(count), &
                  ring(ik)%loc(jm)%inv_gamma_cbar_check(count)) 
          ENDDO

          !
          !Gamma
          !

!          data_struc%loc(jm)%gamma  = &
!               min( data_struc%loc(jm)%gamma, &
!               ring(ik)%loc(jm)%gamma) + &
!               max( data_struc%loc(jm)%gamma, &
!               ring(ik)%loc(jm)%gamma) 
          data_struc%loc(jm)%gamma  = ring(ik)%loc(jm)%gamma + &
               data_struc%loc(jm)%gamma
          !
          !Cbar
          !

!          data_struc%loc(jm)%cbar  = &
!               min( data_struc%loc(jm)%cbar, &
!               ring(ik)%loc(jm)%cbar) + &
!               max( data_struc%loc(jm)%cbar, &
!               ring(ik)%loc(jm)%cbar) 

          data_struc%loc(jm)%cbar  = ring(ik)%loc(jm)%cbar + &
               data_struc%loc(jm)%cbar

       end do
    end do
    deallocate (ring)
  end subroutine calculate_with_known_spacing



  subroutine output (data)

    implicit none
    type(data_set)data
!    type (data_file)file
    integer :: i,j, openstatus
    integer :: n_pairs                  !Number of BPM pairs
   ! type (processor_analysis) :: proc_data
    integer :: bpm


    n_pairs = bpm_pairs(1)%number

    open (unit = 27, file = "./data/mia.out", &
         action = "write", position = "rewind",&
         iostat = openstatus)
    if (openstatus > 0) print *, "*** Cannot open output file ***",&
         openstatus
9   format(1x, f10.7, a3, f10.7, a7, 5x, f10.7, a3, f10.7)
11  format(5x, f11.7, a2,2x, f11.7)
12  format(3(4x, a8, 7x))
17  format(a10, 3x, f11.7, 1x,",", 1x, f11.7)
18  format(4x, a20, 4x, 20x, a5)
34  format(1x, 3(e14.7, a2, 3x))
56  format(a10,3x,f10.7)
    write (27, *) "Number of BPMs: ", NUM_BPMS
    write (27, *) "Number of Turns: ", NUM_TURNS
    write (27, *) "Frequency of machine: ", FREQ
    write (27, *) "Tune A: ", data_struc%tune(1)
    write (27, *), "Tune B: ", data_struc%tune(2)
    write (27, *), "Phi(t) A = ", data_struc%phi_t(1)
    write (27, *), "Phi(t) B = ", data_struc%phi_t(2)
    write (27, *), " Values of Lambda"

    do i=1, 2*NUM_BPMS
       write (27, *), data%lambda(i)
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
    do i = 1, data%bpmproc
       write (27, 9) data_struc%loc(i)%a%d_phase_adv, &
            ",", data_struc%loc(i)%b%d_phase_adv, "||", &
            data_struc%loc(i)%a%phi, ",", data_struc%loc(i)%b%phi
    enddo

    write (27, "(8x, a12)") "d_delta"
    do i = 1, data%bpmproc
       write (27, 11) data_struc%loc(i)%a%d_delta, &
            ",", data_struc%loc(i)%b%d_delta
    enddo

    write (27,"(10x, a12)") "Gamma^2 Beta"
    do i = 1, data%bpmproc
       write (27,11) data_struc%loc(i)%a%gam2_beta, ",", &
            data_struc%loc(i)%b%gam2_beta
    enddo

    write (27,"(20x, a15)") "1/gamma * Cbar"
    write (27,12) "(1,1)", "(1,2)", "(2,2)"

    do i = 1, data%bpmproc
       write (27,34) data_struc%loc(i)%inv_gamma_cbar(1,1), ",", &
            data_struc%loc(i)%inv_gamma_cbar(1,2), ",", &
            data_struc%loc(i)%inv_gamma_cbar(2,2)
    enddo

    write (27,"(4x, a22)") "Average amplitude (J)"
    do i=1, 2
       write (27, "(7x, a5, i2)") "File ", i
       write (27, 11) data_struc%j_amp_ave(1,i),",", &
            data_struc%j_amp_ave(2,i)
    enddo

    close(27)

  end subroutine output

end module orbit_mia
