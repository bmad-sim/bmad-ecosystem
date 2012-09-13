module beam_utils

use beam_def_struct
use spin_mod
use eigen_mod
use wake_mod
use bookkeeper_mod

private init_random_distribution, init_grid_distribution
private init_ellipse_distribution, init_kv_distribution
private recenter_bunch, combine_bunch_distributions, calc_this_emit

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track1_bunch_hom (bunch_start, ele, param, bunch_end)
!
! Subroutine to track a bunch of particles through an element.
!
! Note: This routine is overloaded by the routine track1_bunch. See this
! routine for more details.
!
! Each particle experiences a different longitudinal short-range wakefield.
! ele%value(grad_loss_sr_wake$) is used to tell track1_bmad the appropriate loss
! for each particle.
!
! Modules needed:
!   use beam_mod
!
! Input:
!   bunch_start -- bunch_struct: Starting bunch position.
!   ele         -- Ele_struct: The element to track through.
!   param       -- lat_param_struct: General parameters.
!
! Output:
!   bunch_end -- Bunch_struct: Ending bunch position.
!-

subroutine track1_bunch_hom (bunch_start, ele, param, bunch_end)

implicit none

type (bunch_struct) bunch_start, bunch_end
type (ele_struct) ele, half_ele

type (lat_param_struct) param

real(rp) charge
integer i, j, n, ix_z
logical err_flag

character(20) :: r_name = 'track1_bunch_hom'

! Charge and center

bunch_end = bunch_start

!------------------------------------------------
! Without wakefields just track through.

if (ele%key /= lcavity$ .or. .not. associated(ele%rf_wake) .or. &
            (.not. bmad_com%sr_wakes_on .and. .not. bmad_com%lr_wakes_on)) then

  do j = 1, size(bunch_start%particle)
    call track1 (bunch_start%particle(j), ele, param, bunch_end%particle(j))
  enddo

  bunch_end%charge = sum (bunch_end%particle(:)%charge, &
                      mask = (bunch_end%particle(:)%state == alive$))
  return
endif

!------------------------------------------------
! This calculation is for an cavity with wakefields.
! Put the sr wakefield transverse kicks at the half way point.
! First track half way through. This includes the sr longitudinal wakes 

call transfer_ele (ele, half_ele, .true.)
call create_element_slice (half_ele, ele, ele%value(l$)/2, 0.0_rp, param, .true., .false., err_flag)
if (err_flag) then
  if (bmad_status%exit_on_error) call err_exit
  return
endif

call order_particles_in_z (bunch_end)
do j = 1, size(bunch_end%particle)
  ix_z = bunch_end%ix_z(j) ! z-ordered index of the particles
  if (bunch_end%particle(ix_z)%state /= alive$) cycle
  call add_sr_long_wake (ele, param, bunch_end, j-1, ix_z)
  call track1 (bunch_end%particle(ix_z), half_ele, param, bunch_end%particle(ix_z))
enddo

ele%value(grad_loss_sr_wake$) = 0.0

! Put in the transverse wakefields

call track1_sr_wake (bunch_end, ele)
call track1_lr_wake (bunch_end, ele)

! Track the last half of the cavity. This includes the sr longitudinal wakes 

call create_element_slice (half_ele, ele, ele%value(l$)/2, ele%value(l$)/2, param, .false., .true., err_flag, half_ele)

call order_particles_in_z (bunch_end)
do j = 1, size(bunch_end%particle)
  ix_z = bunch_end%ix_z(j) ! z-ordered index of the particles
  if (bunch_end%particle(ix_z)%state /= alive$) cycle
  call add_sr_long_wake (ele, param, bunch_end, j-1, ix_z)
  call track1 (bunch_end%particle(ix_z), half_ele, param, bunch_end%particle(ix_z))
enddo

ele%value(grad_loss_sr_wake$) = 0.0

bunch_end%charge = sum (bunch_end%particle(:)%charge, &
                         mask = (bunch_end%particle(:)%state == alive$))

end subroutine track1_bunch_hom

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine add_sr_long_wake (ele, param, bunch, num_in_front, follower)
!
! Adds the longitudinal wake for all particles in front of the follower.
!
! Input:
!  ele      -- Ele_struct: Element with wakefields.
!  param    -- lat_param_struct: For param%particle.
!  bunch    -- Bunch_struct: Bunch of particles
!  num_in_front -- Integer: number of particles in front of this one
!                   This will be the bunch%particle index number right before
!                   the follower
!  follower -- Integer: index of particle wakes being applied to.
!
! Output:
!  ele%value(grad_loss_sr_wake$) -- Real(rp): net gradient loss due to the leaders.
!-

subroutine add_sr_long_wake (ele, param, bunch, num_in_front, ix_follower)

implicit none

type (ele_struct) ele
type (lat_param_struct) param
type (bunch_struct) bunch
type (coord_struct), pointer :: leader

integer ix_follower, i, num_in_front
integer n_sr_table, n_sr_mode_long, n_sr_mode_trans, k_start

!-----------------------------------
! If there is no wake for this element, or the sr wakes are turned off, then just 
! use the e_loss attribute (as set in track1_bmad).

ele%value(grad_loss_sr_wake$) = 0.0

n_sr_table = size(ele%rf_wake%sr_table) 
n_sr_mode_long = size(ele%rf_wake%sr_mode_long)
n_sr_mode_trans = size(ele%rf_wake%sr_mode_trans)

if ((n_sr_table == 0 .and. n_sr_mode_long == 0 .and. n_sr_mode_trans == 0) .or. &
                                          .not. bmad_com%sr_wakes_on) then 
  ele%value(grad_loss_sr_wake$) = 0.0
  return 
endif

! the self wake only sees the charge of each real particle, not the "macro"
! charge of the simulated particle

if (n_sr_table > 0) then

  ele%value(grad_loss_sr_wake$) = ele%value(grad_loss_sr_wake$) + &
                   ele%rf_wake%sr_table(0)%long * e_charge * abs(charge_of(param%particle)) / 2.0

  !-----------------------------------
  ! add up all wakes from front of bunch to follower

  do i = 1, num_in_front
    if (bunch%particle(bunch%ix_z(i))%state == alive$) &
      call sr_table_add_long_kick (ele, bunch%particle(bunch%ix_z(i)), &
               bunch%particle(bunch%ix_z(i))%charge, &
               bunch%particle(ix_follower))
  enddo

endif

end subroutine add_sr_long_wake

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track1_sr_wake (bunch, ele)
!
! Subroutine to apply the short range wake fields to a bunch. 
!
! Modules needed:
!   use beam_mod
!
! Input:
!   bunch -- Bunch_struct: Bunch of particles.
!   ele   -- Ele_struct: Element with wakefields.
!
! Output:
!   bunch -- Bunch_struct: Bunch with wakefields applied to the particles.
!-

subroutine track1_sr_wake (bunch, ele)

implicit none

type (bunch_struct), target :: bunch
type (ele_struct) ele
type (coord_struct), pointer :: particle, leader
type (coord_struct), pointer :: p(:)

real(rp) dz_sr_table, sr02, z_sr_table_max
integer i, j, k, i1, i2, i_sr_mode, n_sr_table, n_sr_mode_long, n_sr_mode_trans, k_start

logical wake_here
character(16) :: r_name = 'track1_sr_wake'

!-----------------------------------

if (.not. associated(ele%rf_wake)) return
  
p => bunch%particle
  
! error check and zero wake sums and order particles in z

call order_particles_in_z (bunch)  
if (size(ele%rf_wake%sr_mode_long) /= 0) then
  i1 = bunch%ix_z(1) 
  i2 = bunch%ix_z(size(p))
  if (p(i1)%vec(5) - p(i2)%vec(5) > ele%rf_wake%z_sr_mode_max) then
    call out_io (s_abort$, r_name, &
        'Bunch longer than sr_mode wake can handle for element: ' // ele%name)
    if (bmad_status%exit_on_error) call err_exit
  endif
endif

do i = 1, size(ele%rf_wake%sr_mode_long)
  ele%rf_wake%sr_mode_long%b_sin = 0
  ele%rf_wake%sr_mode_long%b_cos = 0
  ele%rf_wake%sr_mode_long%a_sin = 0
  ele%rf_wake%sr_mode_long%a_cos = 0
enddo

!

n_sr_table = size(ele%rf_wake%sr_table) 
z_sr_table_max = 0
if (n_sr_table > 0) then
  z_sr_table_max = ele%rf_wake%sr_table(n_sr_table-1)%z
  dz_sr_table = z_sr_table_max / (n_sr_table - 1)
endif

! loop over all particles in the bunch and apply the wake

i_sr_mode = 1  ! index of next particle to be added to the sr_mode wake sums.

do j = 1, size(p)
  particle => p(bunch%ix_z(j))

  ! Particle_j is kicked by particles k = 1, ..., j-1.
  ! The particles 1, ... i_sr_mode-1 have already had their wakes added to the 
  ! sr_mode wake sums so the loop is from i_sr_mode, ..., j-1.

  k_start = i_sr_mode
  do k = k_start, j-1
    leader => p(bunch%ix_z(k))
    if ((particle%vec(5) - leader%vec(5)) > z_sr_table_max) then
      ! use sr_table table to add to particle j the wake of particle k
      call sr_table_apply_trans_kick (ele, leader, leader%charge, particle)
    else
      ! add contribution of particle(k) to wake sums
      i_sr_mode = k  ! update i_sr_mode
      call sr_mode_long_wake_add_to (ele, leader, leader%charge)
      call sr_mode_trans_wake_add_to (ele, leader, leader%charge)
    endif
  enddo

  ! apply wake to particle(j)

  ! apply longitudinal self wake

  call sr_mode_long_wake_apply_kick (ele, particle%charge, particle)
  call sr_mode_trans_wake_apply_kick(ele, particle)

enddo

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track1_lr_wake (bunch, ele)
!
! Subroutine to put in the long-range wakes for particle tracking.
!
! Note: It is the responsibility of the calling routine to zero the wakefield
! components before the first bunch is sent through. The wakefield components 
! are:
!     ele%rf_wake%lr%b_sin
!     ele%rf_wake%lr%b_cos
!     ele%rf_wake%lr%a_sin
!     ele%rf_wake%lr%a_cos
!
! Modules needed:
!   use beam_mod
!
! Input:
!   bunch -- Bunch_struct: Bunch of particles.
!   ele   -- Ele_struct: Element with wakefields.
!
! Output:
!   bunch -- Bunch_struct: Bunch with wakefields applied to the particles.
!   ele   -- Ele_struct: Element with updated wakefields.
!-

subroutine track1_lr_wake (bunch, ele)

implicit none

type (bunch_struct), target :: bunch
type (ele_struct) ele
type (coord_struct), pointer :: particle

integer n_mode, j, k

if (.not. bmad_com%lr_wakes_on) return
if (.not. associated(ele%rf_wake)) return
  
! Check to see if we need to do any calc

if (.not. associated(ele%rf_wake)) return
n_mode = size(ele%rf_wake%lr)
if (n_mode == 0) return  

call order_particles_in_z (bunch)  ! needed for wakefield calc.

! Give the particles a kick

do k = 1, size(bunch%particle)
  j = bunch%ix_z(k)
  particle => bunch%particle(j)
  if (particle%state /= alive$) cycle
  call lr_wake_apply_kick (ele, bunch%t_center, particle, particle%charge)
enddo

! Add the wakes left by this bunch to the existing wakes.

do k = 1, size(bunch%particle)
  j = bunch%ix_z(k)
  particle => bunch%particle(j)
  if (particle%state /= alive$) cycle
  call lr_wake_add_to (ele, bunch%t_center, particle, particle%charge)
enddo

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine order_particles_in_z (bunch)
!
! Routine to order the particles longitudinally in terms of decreasing %vec(5).
! That is from large z (head of bunch) to small z.
!
! Modules needed:
!   use beam_mod
!
! Input:
!   bunch     -- Bunch_struct: collection of particles.
!     %particle(j)%vec(5) -- Longitudinal position of j^th particle.
!
! Output:
!   bunch     -- bunch_struct: collection of particles.
!     %ix_z(:)     -- Index for the ordering. 
!                     Order is from large z (head of bunch) to small z.
!                     That is: %bunch%ix_z(1) is the particle at the bunch head. 
!-

Subroutine order_particles_in_z (bunch)

implicit none

type (bunch_struct), target :: bunch
type (coord_struct), pointer :: particle(:)
type (coord_struct) temp
integer i, k, nm, i0, i1
real(rp) z1, z2
logical ordered

! Init if needed

particle => bunch%particle
nm = size(particle)

if (bunch%ix_z(1) == 0) then
  forall (i = 1:nm) bunch%ix_z(i) = i
endif

! Order is from large z (head of bunch) to small z.

do
  ordered = .true.
  do i = 1, nm-1
    i0 = bunch%ix_z(i); i1 = bunch%ix_z(i+1)
    if (particle(i0)%vec(5) < particle(i1)%vec(5)) then
      bunch%ix_z(i:i+1) = bunch%ix_z(i+1:i:-1)
      ordered = .false.
    endif
  enddo
  if (ordered) exit
enddo

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine reallocate_beam (beam, n_bunch, n_particle)
! 
! Subroutine to reallocate memory within a beam_struct.
!
! If n_bunch = 0 then all macro beam pointers will be deallocated.
! Rule: If beam%bunch(:) is allocated, beam%bunch(i)%particle(:) will be allocated.
!
! Modules needed:
!   use beam_mod
!
! Input:
!   n_bunch    -- Integer: Number of bunches.
!   n_particle -- Integer: Number of particles. Must be non-negative.
!
! Output:
!   beam -- beam_struct: Allocated beam_struct structure.
!-

subroutine reallocate_beam (beam, n_bunch, n_particle)

implicit none

type (beam_struct) beam

integer i, n_bunch, n_particle

! Deallocate if needed

if (allocated(beam%bunch)) then
  if (n_bunch == 0 .or. size(beam%bunch) /= n_bunch) deallocate (beam%bunch)
endif

if (n_bunch == 0) return
  
! Allocate

if (.not. allocated (beam%bunch)) allocate (beam%bunch(n_bunch))

do i = 1, n_bunch
  call reallocate_bunch (beam%bunch(i), n_particle)
enddo

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine reallocate_bunch (bunch, n_particle)
! 
! Subroutine to reallocate particles within a bunch_struct.
!
! Modules needed:
!   use bunch_mod
!
! Input:
!   n_particle -- Integer: Number of particles. Must be non-negative.
!
! Output:
!   bunch -- bunch_struct: Allocated bunch_struct structure.
!-

subroutine reallocate_bunch (bunch, n_particle)

implicit none

type (bunch_struct) bunch

integer i, n_particle

! Deallocate if needed

if (allocated(bunch%particle)) then
  if (size(bunch%particle) /= n_particle) deallocate (bunch%particle, bunch%ix_z)
endif

if (.not. allocated(bunch%particle)) then
  allocate (bunch%particle(n_particle), bunch%ix_z(n_particle))
  bunch%ix_z = 0
endif

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine init_beam_distribution (ele, param, beam_init, beam)
!
! Subroutine to initialize a beam of particles. 
! 
! Note: This routine sets the random number generator according to the settings
! in beam_int and at the end resets things to their initial state.
!
! For more information on individual bunch initialization, see the 
! init_bunch_distribution routine.
! 
! Modules needed:
!   use beam_mod
!
! Input:
!   ele         -- Ele_struct: element to initialize distribution at
!   param       -- Lat_param_struct: Lattice parameters
!     %particle      -- Type of particle.
!   beam_init   -- beam_init_struct: Use "getf beam_init_struct" for more details.
!
! Output:
!   beam        -- Beam_struct: Structure with initialized particles.
!
!-

subroutine init_beam_distribution (ele, param, beam_init, beam)
 
use random_mod

implicit none

type (ele_struct) ele
type (lat_param_struct) param
type (beam_init_struct), target :: beam_init
type (beam_struct), target :: beam
type (bunch_struct), pointer :: bunch

integer i_bunch, i, n, n_kv
real(rp) old_cutoff

character(16) old_engine, old_converter  
character(22) :: r_name = "init_beam_distribution"

call reallocate_beam (beam, beam_init%n_bunch, 0)

! Save and set the random number generator parameters.

call ran_engine (beam_init%random_engine, old_engine)
call ran_gauss_converter (beam_init%random_gauss_converter, &
                  beam_init%random_sigma_cutoff, old_converter, old_cutoff)

! Loop over all bunches
! Note z_center is negative and t_center is posive for trailing bunches.

do i_bunch = 1, size(beam%bunch)
  bunch => beam%bunch(i_bunch)
  call init_bunch_distribution (ele, param, beam_init, bunch)

  bunch%t_center = (i_bunch-1) * beam_init%dt_bunch
  bunch%z_center = -bunch%t_center * c_light * ele%value(e_tot$) / ele%value(p0c$)
  bunch%ix_bunch = i_bunch

  bunch%particle(:)%t = bunch%particle(:)%t + bunch%t_center

enddo
  
! Reset the random number generator parameters.

call ran_engine (old_engine)
call ran_gauss_converter (old_converter, old_cutoff)

end subroutine init_beam_distribution

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine init_bunch_distribution (ele, param, beam_init, bunch)
!
! Subroutine to initialize a distribution of particles of a bunch.
!
! There are four distributions available: 
!   '', or 'ran_gauss' -- Random gaussian distribution.
!   'ellipse'  -- concentric ellipses representing a Gaussian distribution
!   'grid'     -- uniform rectangular grid
!   'KV'       -- Kapchinsky-Vladimirsky distribution
! See the Bmad manual for more information.
!
! The distribution is matched to the Twiss parameters, centroid position,
! and Energy - z correlation as specified. Coupling in the element ele is 
! incorporated into the distribution.
!
! Note: Except for the random number seed, the random number generator 
! parameters used for this routine are set from the beam_init argument.
! That is, these parameters are independent of what is used everywhere else.
!
! Note: Make sure: |beam_init%dpz_dz| < mode%sigE_E / mode%sig_z
!
! Note: To get good results, It is important to make sure that for 
! circular rings that beam_init%center is the correct closed orbit. 
! The closed orbit will shift if, for example, radiation damping is
! turned on.
!
! Modules needed:
!   use beam_mod
!
! Input:
!   ele         -- Ele_struct: element to initialize distribution at
!   param       -- Lat_param_struct: Lattice parameters
!   beam_init   -- beam_init_struct: Use "getf beam_init_struct" for more details.
!
! Output:
!   bunch        -- bunch_struct: Structure with initialized particles.
!-

subroutine init_bunch_distribution (ele, param, beam_init, bunch)

use mode3_mod

implicit none

type (ele_struct) ele
type (ele_struct) twiss_ele
type (lat_param_struct) param
type (beam_init_struct), target :: beam_init
type (bunch_struct), target :: bunch
type (coord_struct) p_temp
type (coord_struct), pointer :: p
type (kv_beam_init_struct), pointer :: kv

real(rp) beta(3), alpha(3), emit(3), covar
real(rp) v_mat(4,4), v_inv(4,4), beta_vel

real(rp) tunes(1:3)
REAL(rp) G6mat(6,6)
REAL(rp) G6inv(6,6)
REAL(rp) V6mat(6,6)
REAL(rp) t6(6,6)

integer i, j, k, n, species
integer :: n_kv     ! counts how many phase planes are of KV type
integer :: ix_kv(3) ! indices (1,2,3) of the two KV planes or 0 if uninitialized

character(22) :: r_name = "init_bunch_distribution"

logical ok, correct_for_coupling(6)
logical ran_gauss_here

! Checking that |beam_init%dpz_dz| < mode%sigE_E / mode%sig_z

if (abs(beam_init%dPz_dz * beam_init%sig_z) > beam_init%sig_e) then
  call out_io (s_abort$, r_name, "|dpz_dz| MUST be < mode%sigE_E / mode%sig_z")
  if (bmad_status%exit_on_error) call err_exit
endif

! Compute the Twiss parameters beta and alpha, and the emittance for each plane
! 1 = (x,px), 2 = (y,py), 3 = (z,pz)

if (beam_init%full_6D_coupling_calc) then
  twiss_ele%a%beta = 1.0d0
  twiss_ele%a%alpha = 0.0d0
  twiss_ele%b%beta = 1.0d0
  twiss_ele%b%alpha = 0.0d0
  twiss_ele%value(e_tot$) = ele%value(e_tot$)
else
  call transfer_ele (ele, twiss_ele, .true.)
endif

call calc_this_emit (beam_init, twiss_ele, param)

covar = beam_init%dPz_dz * beam_init%sig_z**2
twiss_ele%z%emit = sqrt((beam_init%sig_z*beam_init%sig_e)**2 - covar**2)
if (twiss_ele%z%emit == 0 .or. beam_init%full_6D_coupling_calc) then
  twiss_ele%z%beta = 1
  twiss_ele%z%alpha = 0
else
  twiss_ele%z%beta = beam_init%sig_z**2 / twiss_ele%z%emit
  twiss_ele%z%alpha = - covar / twiss_ele%z%emit
endif

! Init

n_kv = 0
ix_kv = 0
ran_gauss_here = .false.
correct_for_coupling = .true.

! Fill the corresponding struct and generate the distribution for each phase plane.
! init_random_distribution must be called last since the other distributions "multiply"
! the number of particles (see the combine_bunch_distribution routine).

call reallocate_bunch (bunch, 0)

do i = 1, 3
  call str_upcase (beam_init%distribution_type(i), beam_init%distribution_type(i))
  select case (beam_init%distribution_type(i))
  case ('', 'RAN_GAUSS')
    ran_gauss_here = .true.
  case ('ELLIPSE')
    call init_ellipse_distribution (i, beam_init%ellipse(i), twiss_ele, bunch)
  case ('GRID')
    call init_grid_distribution (i, beam_init%grid(i), bunch)
    ! The grid distribution ignores the local twiss parameters.
    correct_for_coupling(2*i-1:2*i) = .false.
  case ('KV') 
    n_kv = n_kv + 1
    ix_kv(n_kv) = i
  case default
    call out_io (s_abort$, r_name, 'PHASE SPACE DISTRIBUTION TYPE NOT RECOGNIZED')
    if (bmad_status%exit_on_error) call err_exit
    return
  end select
enddo

if (n_kv == 2) call init_KV_distribution (ix_kv(1), ix_kv(2), beam_init%kv, twiss_ele, bunch)

if (ran_gauss_here) then
  call init_random_distribution (twiss_ele, param, beam_init, bunch)
endif

if (beam_init%full_6D_coupling_calc) then
  call transfer_matrix_calc(ele%branch%lat, .true., t6, ix1=ele%ix_ele)
  call normal_mode3_calc(t6, tunes, G6mat, V6mat, .true.)
  call mat_inverse(G6mat, G6inv, ok)
  do i = 1, size(bunch%particle)
    p => bunch%particle(i)
    p%vec = matmul(G6inv, p%vec)
  enddo
else
  call make_v_mats(ele, v_mat, v_inv)
endif

bunch%charge = beam_init%bunch_charge
bunch%ix_ele = ele%ix_ele
bunch%species = param%particle

do i = 1, size(bunch%particle)
  p => bunch%particle(i)
  p%charge = bunch%charge * p%charge
  p%state = alive$

  ! Include Dispersion and coupling
  if (beam_init%full_6D_coupling_calc) then
    p%vec(1:6) = matmul(V6mat, p%vec(1:6))

  else
    p_temp%vec(1:4) = p%vec(1:4) + p%vec(6) * [ele%a%eta, ele%a%etap, ele%b%eta, ele%b%etap]
    p_temp%vec(1:4) = matmul(v_mat, p_temp%vec(1:4))
    where (correct_for_coupling(1:4)) p%vec(1:4) = p_temp%vec(1:4)
  endif

enddo

! recenter the bunch and include beam jitter

call recenter_bunch (beam_init, bunch)

bunch%z_center = 0  ! Default
bunch%t_center = 0  ! Default
bunch%ix_bunch = 1  ! Default

! particle spin

call init_spin_distribution (beam_init, bunch)

! Photons:
! For now just give one half e_field_x = 1 and one half e_field_y = 1

species = param%particle

if (species == photon$) then
  n = size(bunch%particle)
  bunch%particle(1:n:2)%e_field_x = 1
  bunch%particle(2:n:2)%e_field_y = 1
endif

! Fill in %species, %t and %s

do i = 1, size(bunch%particle)
  p => bunch%particle(i)
  p%s = ele%s
  call convert_pc_to (ele%value(p0c$) * (1 + p%vec(6)), species, beta = beta_vel)
  p%t = ele%ref_time - p%vec(5) / (beta_vel * c_light)
  call init_coord (p, p, ele, .true.)
enddo

end subroutine init_bunch_distribution

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine calc_this_emit (beam_init, ele, param)
!
! Private routine to calculate the emittances
!
! Input:
!   beam_init -- beam_init_struct: 
!   ele       -- ele_struct:
!   param     -- lat_param_struct:
!
! Ouput:
!   ele%a  -- Real(rp): a emittance
!   ele%b  -- Real(rp): b emittance
!-

subroutine calc_this_emit (beam_init, ele, param)

implicit none

type (beam_init_struct) beam_init
type (ele_struct) ele
type (lat_param_struct) param

real(rp) ran_g(2)

character(16) :: r_name = 'calc_this_emit'

! Check

if ((beam_init%a_norm_emit /= 0 .and. beam_init%a_emit /= 0) .or. &
    (beam_init%b_norm_emit /= 0 .and. beam_init%b_emit /= 0)) then
  call out_io (s_fatal$, r_name, 'SETTING BOTH NORM_EMIT AND EMIT IN BEAM_INIT STRUCTURE IS NOT ALLOWED.')
  if (bmad_status%exit_on_error) call err_exit
endif

!

if (beam_init%a_norm_emit /= 0) then
  ele%a%emit = beam_init%a_norm_emit * mass_of(param%particle) / ele%value(e_tot$)
else
  ele%a%emit = beam_init%a_emit
endif

if (beam_init%b_norm_emit /= 0) then
  ele%b%emit = beam_init%b_norm_emit * mass_of(param%particle) / ele%value(e_tot$)
else
  ele%b%emit = beam_init%b_emit 
endif

! Add jitter if needed

if (any(beam_init%emit_jitter /= 0)) then
  call ran_gauss(ran_g) ! ran(3:4) for z and e jitter used below
  ele%a%emit = ele%a%emit * (1 + beam_init%emit_jitter(1) * ran_g(1))
  ele%b%emit = ele%b%emit * (1 + beam_init%emit_jitter(2) * ran_g(2))
endif

end subroutine calc_this_emit 

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine init_random_distribution (ele, param, beam_init, bunch)
!
! Subroutine to initialize a random bunch of particles matched to
! the Twiss parameters, centroid position, and Energy - z correlation
! as specified. Coupling in the element ele is incorporated into the
! distribution.
!
! Note: This routine is private. Use init_bunch_distribution instead.
!-

subroutine init_random_distribution (ele, param, beam_init, bunch)
 
use random_mod

implicit none

type (ele_struct) ele
type (lat_param_struct) param
type (beam_init_struct) beam_init
type (bunch_struct), target :: bunch
type (coord_struct), allocatable :: p(:)
  
real(rp) dpz_dz, denom
real(rp) a_emit, b_emit, y, a, b
real(rp) ave(6), sigma(6), alpha(6), sig_mat(6,6), r(6)
real(rp) center(6), ran_g(2), old_cutoff

integer i, j, j2, n, i_bunch, n_particle

logical is_ran_plane(3)

character(16) old_engine, old_converter  
character(28) :: r_name = "init_random_distribution"

! If random is to be combined with other distributions, the number
! of particles is set by the other distributions.

is_ran_plane = (beam_init%distribution_type == '' .or. beam_init%distribution_type == 'RAN_GAUSS')

n_particle = beam_init%n_particle
if (any(.not. is_ran_plane)) n_particle = size(bunch%particle)

if (n_particle < 1) then
  call out_io (s_abort$, r_name, 'NUMBER OF PARTICLES TO INIT FOR BEAM IS ZERO!')
  if (bmad_status%exit_on_error) call err_exit
endif

allocate(p(n_particle))

sig_mat = 0
ave = 0
do n = 1, n_particle
  call ran_gauss(r)
  p(n)%vec = r
  ave = ave + r
  forall (i=1:6, j=1:6) sig_mat(i,j) = sig_mat(i,j) + r(i) * r(j)
enddo  

ave = ave / n_particle
sig_mat = sig_mat / n_particle

! Since we are dealing with a finite number of particles, 
! the sigmas of the distributions in each dimension
! will not be exactly 1 and there will 
! be some correlation between different dimensions.
! If beam_init%renorm_sigma = True then take this out.
! That is, we want to make sig_mat = the unit matrix.
! Exception: Ignore if n_particle = 1.
! Also: If n_particle < 7 then cannot remove correlations between 
! dimensions so, in this case, only normalize sig_mat(i,i) to 1.

if (beam_init%renorm_sigma) then

  ! Make sure average is zero.

  do n = 1, n_particle
    p(n)%vec = p(n)%vec - ave
  enddo

  forall (i = 1:6, j = 1:6) sig_mat(i,j) = sig_mat(i,j) - ave(i) * ave(j)

  ! The first step is to zero the off-diagonal elements.
  ! We have to do this in the correct order otherwise zeroing one element
  ! might unzero others that have already been zeroed.

  if (n_particle > 6) then
    do i = 5, 1, -1
      do j = i+1, 6

        b = -sig_mat(i,j) / sig_mat(j,j)
        ! Transform the distribution
        do n = 1, n_particle
          p(n)%vec(i) = p(n)%vec(i) + b * p(n)%vec(j)
        enddo
        ! Since we have transformed the distribution we need to transform
        ! sig_mat to keep things consistant.
        sig_mat(i,i) = sig_mat(i,i) + 2 * b * sig_mat(i,j) + b**2 * sig_mat(j,j)
        do j2 = 1, 6
          if (j2 == i) cycle
          sig_mat(i,j2) = sig_mat(i,j2) + b * sig_mat(j ,j2)
          sig_mat(j2,i) = sig_mat(i,j2)
        enddo

      enddo
    enddo
  endif

  ! Now we make the diagonal elements unity

  if (n_particle > 1) then
    do i = 1, 6
      alpha(i) = sqrt(1/sig_mat(i,i))
    enddo

    do n = 1, n_particle
      p(n)%vec = p(n)%vec * alpha
    enddo
  endif

endif

! Compute sigmas

call calc_this_emit(beam_init, ele, param)

dpz_dz = beam_init%dpz_dz
  
call ran_gauss(ran_g) 
sigma(1) = sqrt(ele%a%emit * ele%a%beta)
sigma(2) = sqrt(ele%a%emit / ele%a%beta)
sigma(3) = sqrt(ele%b%emit * ele%b%beta)
sigma(4) = sqrt(ele%b%emit / ele%b%beta)
if (beam_init%full_6D_coupling_calc) then
  sigma(5) = sqrt(beam_init%sig_z*beam_init%sig_e)
  sigma(6) = sqrt(beam_init%sig_z*beam_init%sig_e)
else
  sigma(5) = beam_init%sig_z * (1 + beam_init%sig_z_jitter*ran_g(1))
  sigma(6) = beam_init%sig_e * (1 + beam_init%sig_e_jitter*ran_g(2))
endif

if (sigma(6) == 0 .or. dpz_dz == 0) then
  a = 0
else if (abs(dpz_dz * sigma(5)) > sigma(6)) then
  call out_io (s_abort$, r_name, "|dpz_dz| MUST be < mode%sigE_E / mode%sig_z")
  if (bmad_status%exit_on_error) call err_exit
else
  a = dpz_dz * sigma(5) / sigma(6)
endif

b = sqrt(1-a**2)
     
! Put everything together to distribute the particles.

do i = 1, n_particle
  r = p(i)%vec
  p(i)%vec(1) =  sigma(1) *  r(1)
  p(i)%vec(2) = -sigma(2) * (r(2) + r(1) * ele%a%alpha)
  p(i)%vec(3) =  sigma(3) *  r(3)
  p(i)%vec(4) = -sigma(4) * (r(4) + r(3) * ele%b%alpha)
  p(i)%vec(5) =  sigma(5) *  r(5)
  p(i)%vec(6) =  sigma(6) * (r(6) * b + r(5) * a)
end do

! Set particle charge and transfer info the the bunch

p(:)%charge = 1.0_rp / n_particle
call combine_bunch_distributions (bunch, p, is_ran_plane, .false.)

end subroutine init_random_distribution

!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
!+
! Subroutine init_grid_distribution (ix_plane, grid, bunch)
!
! Subroutine to initialize a uniform rectangular grid as the phase space
! distribution of a bunch.
!
! Input:
!    ix_plane       -- Integer: Index of plane of this distribution: 1, 2, or 3
!   grid            -- grid_beam_init_struct: Grid info.
!     %n_x              -- number of columns
!     %n_px             -- number of rows
!     %x_min, %x_max    -- upper and lower limits in beam size
!     %px_min, %px_max  -- upper and lower limits in divergence
!
! Output:
!   bunch     -- Bunch_struct: Bunch structure
!-

subroutine init_grid_distribution (ix_plane, grid, bunch)

implicit none

type (grid_beam_init_struct) grid
type (coord_struct), allocatable :: p(:)
type (bunch_struct) bunch

integer i, j, k, ix_plane, n_particle

real(rp) x, px

logical where(3)

character(28) :: r_name = 'init_grid_distribution'

!

n_particle = grid%n_x * grid%n_px       ! total number of particles
if (n_particle < 1) then
  call out_io (s_abort$, r_name, 'NUMBER OF PARTICLES TO INIT FOR BEAM IS ZERO!')
  if (bmad_status%exit_on_error) call err_exit
endif

allocate (p(n_particle))

k = 1

do i = 1, grid%n_x
   if (grid%n_x == 1) then
      x = grid%x_min
   else
      x = grid%x_min + real(i - 1)/(grid%n_x - 1) * (grid%x_max - grid%x_min)
   endif

   do j = 1, grid%n_px
      if (grid%n_px == 1) then
         px = grid%px_min
      else
         px = grid%px_min + real(j - 1)/(grid%n_px - 1) * (grid%px_max - grid%px_min)
      endif

      p(k)%vec(2*ix_plane-1) = x
      p(k)%vec(2*ix_plane)   = px
      p(k)%charge = 1.0_rp / n_particle     ! total charge = 1

      k = k + 1
   enddo
enddo

! Combine with bunch distribution

where = .false.;  where(ix_plane) = .true.
call combine_bunch_distributions (bunch, p, where, .true.)

end subroutine init_grid_distribution


!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
!+
! Subroutine init_ellipse_distribution (ix_plane, ellipse, ix_plane, twiss_ele, bunch)
!
! Subroutine to initalize a phase space distribution as a set of concentric
! ellipses of macroparticles representing a Gaussian distribution.
!
! Input:
!   ellipse           -- ellipse_distribution_struct: Init info.
!     %n_ellipse         -- number of ellipses (>= 1)
!     %part_per_ellipse  -- number of particles per ellipse
!     %sigma_cutoff      -- sigma cutoff of the representation
!   ix_plane          -- Integer: Plane of distribution. 1, 2, or 3.
!   twiss_ele%a       -- Twiss parameters
!   twiss_ele%b       -- Twiss parameters
!   twiss_ele%z       -- Twiss parameters
!   twiss_ele%a%emit  -- emittance
!   twiss_ele%b%emit  -- emittance
!   twiss_ele%z%emit  -- emittance
!
! Output:
!   bunch     -- Bunch_struct: Bunch structure
!
! See manual for more details.
!-

subroutine init_ellipse_distribution (ix_plane, ellipse, twiss_ele, bunch)

implicit none

type (bunch_struct) bunch
type (coord_struct), allocatable :: p(:)
type (ellipse_beam_init_struct), target :: ellipse
type (ele_struct) twiss_ele
type (ellipse_beam_init_struct), pointer :: e

real(rp) beta, alpha, emit

integer ix_plane, n_particle
integer n, m, k

real(rp) b_inner, b_outer                  ! B_{n-1}/epsilon  and B_{n}/epsilon in the bmad manual

real(rp) J, phi
real(rp) x, px, charge

logical where(3)

character(28) :: r_name = 'init_ellipse_distribution'

!

if (ix_plane == 1) then
  beta = twiss_ele%a%beta
  alpha = twiss_ele%a%alpha
  emit = twiss_ele%a%emit
elseif (ix_plane == 2) then
  beta = twiss_ele%b%beta
  alpha = twiss_ele%b%alpha
  emit = twiss_ele%b%emit
elseif (ix_plane == 3) then
  beta = twiss_ele%z%beta
  alpha = twiss_ele%z%alpha
  emit = twiss_ele%z%emit
endif

e => ellipse
n_particle = e%n_ellipse * e%part_per_ellipse
if (n_particle < 1) then
  call out_io (s_abort$, r_name, 'NUMBER OF PARTICLES TO INIT FOR BEAM IS ZERO!')
  if (bmad_status%exit_on_error) call err_exit
endif

allocate (p(n_particle))

k = 0
b_outer = 0

do n = 1, e%n_ellipse

  b_inner = b_outer
  b_outer = e%sigma_cutoff**2/2.0 * (real(n)/e%n_ellipse)**2

  if (n == e%n_ellipse) then
    ! This is the ellipse that represents the distribution out to infinity
    charge = exp(-b_inner)       ! q_n
    J = emit * (b_inner + 1.0) * exp(-b_inner) / charge   ! J_n
  else
    charge = exp(-b_inner) - exp(-b_outer)   ! q_n
    J = emit * ((b_inner + 1.0) * exp(-b_inner) - (b_outer + 1.0) * exp(-b_outer)) / charge   ! J_n
  endif

  do m = 1, e%part_per_ellipse
    phi = (twopi * m) / e%part_per_ellipse
    k = k + 1
    p(k)%vec(2*ix_plane-1) =  sqrt(2 * J * beta) * cos(phi)
    p(k)%vec(2*ix_plane)   = -sqrt(2 * J / beta) * (alpha * cos(phi) + sin(phi))
    p(k)%charge = charge / e%part_per_ellipse
  enddo

enddo

! Combine with bunch distribution

where = .false.;  where(ix_plane) = .true.
call combine_bunch_distributions (bunch, p, where, .true.)

end subroutine init_ellipse_distribution


!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
!+
! Subroutine init_KV_distribution (ix1_plane, ix2_plane, kv, ele, bunch)
!
! Subroutine to initalize a phase space distribution as a set of concentric
! ellipses of macroparticles representing a Kapchinsky-Vladimirsky distribution.
!
! See manual for more details.
!
! Input:
!   ix1_plane        -- Integer: Index of first plane.
!   ix2_plane        -- Integer: Index of second plane.
!   kv               -- kv_beam_init_struct: KV info.
!   twiss_ele%a,b,z  -- Twiss parameters of each phase plane
!   twiss_ele%a%emit -- emittance of a phase plane
!   twiss_ele%b%emit -- emittance of b phase plane
!   twiss_ele%z%emit -- emittance of z phase plane
!
! Output:
!   bunch     -- Bunch_struct: Bunch structure
!-

subroutine init_KV_distribution (ix1_plane, ix2_plane, kv, ele, bunch)

implicit none

type (bunch_struct) bunch
type (kv_beam_init_struct) kv
type (ele_struct) ele
type (coord_struct), allocatable :: p(:)

real(rp) beta1, beta2, alpha1, alpha2, emit1, emit2

integer i_I2, i_phi1, i_phi2, k, n_particle, ix1_plane, ix2_plane, n_p1, n_p2

real(rp) emit_tot
real(rp) I1, I2
real(rp) J1, J2, phi1, phi2
real(rp) x1, x2, px1, px2, charge

logical where(3)

character(28) :: r_name = 'init_kv_distribution'

!

if (ix1_plane == 1) then
  beta1 = ele%a%beta
  alpha1 = ele%a%alpha
  emit1 = ele%a%emit
elseif (ix1_plane == 2) then
  beta1 = ele%b%beta
  alpha1 = ele%b%alpha
  emit1 = ele%b%emit
elseif (ix1_plane == 3) then
  beta1 = ele%z%beta
  alpha1 = ele%z%alpha
  emit1 = ele%z%emit
endif

if (ix2_plane == 1) then
  beta2 = ele%a%beta
  alpha2 = ele%a%alpha
  emit2 = ele%a%emit
elseif (ix2_plane == 2) then
  beta2 = ele%b%beta
  alpha2 = ele%b%alpha
  emit2 = ele%b%emit
elseif (ix2_plane == 3) then
  beta2 = ele%z%beta
  alpha2 = ele%z%alpha
  emit2 = ele%z%emit
endif

n_p1 = kv%part_per_phi(1)
n_p2 = kv%part_per_phi(2)

n_particle = kv%n_i2 * n_p1 * n_p2
if (n_particle < 1) then
  call out_io (s_abort$, r_name, 'NUMBER OF PARTICLES TO INIT FOR BEAM IS ZERO!')
  if (bmad_status%exit_on_error) call err_exit
endif

allocate (p(n_particle))

emit_tot = 1.0 / sqrt(1.0 / emit1**2 + 1.0 / emit2**2)
I1 = kv%A * emit_tot

k = 1

do i_I2 = 1, kv%n_i2
  I2 = -emit1/emit2 * I1 + emit1*emit2/emit_tot**2 * I1 * real(i_I2 - 0.5)/kv%n_i2

  J1 = (I1/emit1 - I2/emit2) * emit_tot
  J2 = (I1/emit2 + I2/emit1) * emit_tot
   
  do i_phi1 = 1, n_p1
    phi1 = 2.0 * pi * real(i_phi1 - 1)/n_p1
 
    do i_phi2 = 1, n_p2
      phi2 = 2.0 * pi * real(i_phi2 - 1)/n_p2

      x1 = sqrt(2.0 * J1 * beta1) * cos(phi1)
      px1 = -sqrt(2.0 * J1 / beta1) * (alpha1 * cos(phi1) + sin(phi1))
      x2 = sqrt(2.0 * J2 * beta2) * cos(phi2)
      px2 = -sqrt(2.0 * J2 / beta2) * (alpha2 * cos(phi2) + sin(phi2))
     
      p(k)%vec(2*ix1_plane-1) = x1
      p(k)%vec(2*ix1_plane)   = px1
      p(k)%vec(2*ix2_plane-1) = x2
      p(k)%vec(2*ix2_plane)   = px2
      p(k)%charge = 1.0_rp / n_particle

      k = k + 1
    enddo
  enddo
enddo

! Combine with bunch distribution

where = .false.;  where(ix1_plane) = .true.;  where(ix2_plane) = .true.
call combine_bunch_distributions (bunch, p, where, .true.)

end subroutine init_KV_distribution

!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
!+
! Subroutine combine_bunch_distributions (bunch, particle, where, do_multiply)
!
! This subroutine combines two bunch distributions together.
!
! If do_multiply is True:
!   size(combined distribution) = size(bunch%particle) * size(particle)
! If do_multiply is False:
!   size(combined distribution) = size(bunch%particle) = size(particle)
!
! Input:
!   bunch       -- bunch_struct: Structure holding the old distribution
!   particle(:) -- coord_struct, allocatable: A new distribution.
!                   This array will be deallocated.
!   where(3)    -- logical: Which planes of particle have the new distribution.
!   do_multiply -- logical: Determines type of combination.
!
! Output:
!   bunch       -- bunch_struct: Structure holding the combined distribution.
!-

subroutine combine_bunch_distributions (bunch, particle, where, do_multiply)

implicit none

type (bunch_struct) bunch
type (coord_struct), allocatable :: particle(:), p(:)

integer i, j, k, m

logical where(:), do_multiply

! If bunch%particle do not contain a distribution, just transfer particle to it.

if (size(bunch%particle) == 0) then
  call reallocate_bunch (bunch, size(particle))
  bunch%particle = particle
  deallocate(particle)
  return
endif

!------------------------
! Multiply combination

if (do_multiply) then
  allocate (p(size(bunch%particle) * size(particle)))
  m = 0
  do i = 1, size(bunch%particle)
    do j = 1, size(particle)
      m = m + 1
      p(m)%charge = bunch%particle(i)%charge * particle(j)%charge
      p(m)%vec = bunch%particle(i)%vec
      do k = 1, 3
        if (.not. where(k)) cycle
        p(m)%vec(2*k-1:2*k) = particle(j)%vec(2*k-1:2*k)
      enddo
    enddo
  enddo

  ! Transfer to bunch

  call reallocate_bunch (bunch, size(p))
  bunch%particle = p
  deallocate (p, particle)

! Overlap combination

else
  do i = 1, size(bunch%particle)
    do k = 1, 3
      if (.not. where(k)) cycle
      bunch%particle(i)%vec(2*k-1:2*k) = particle(i)%vec(2*k-1:2*k)
    enddo
  enddo
  deallocate(particle)
endif

end subroutine combine_bunch_distributions

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine recenter_bunch (beam_init, bunch)
!
! Recenters the bunch and puts in the beam jitter
!
! Input:
!   beam_init  -- beam_init_struct: Use "getf beam_init_struct" for more details.
!     %center(6)          -- Bunch center offset relative to reference
!     %center_jitter(6)   -- Bunch center rms jitter
!   bunch      -- bunch centered at the reference

! Output:
!   bunch      -- bunch recentered
!-

subroutine recenter_bunch (beam_init, bunch)

use random_mod

implicit none

type (beam_init_struct) beam_init
type (bunch_struct), target :: bunch
type (coord_struct), pointer :: p

real(rp) ran(6), center(6)
integer i

call ran_gauss(ran)
center(1) = beam_init%center(1) + beam_init%center_jitter(1)*ran(1)
center(2) = beam_init%center(2) + beam_init%center_jitter(2)*ran(2) 
center(3) = beam_init%center(3) + beam_init%center_jitter(3)*ran(3)
center(4) = beam_init%center(4) + beam_init%center_jitter(4)*ran(4)
center(5) = beam_init%center(5) + beam_init%center_jitter(5)*ran(5)
center(6) = beam_init%center(6) + beam_init%center_jitter(6)*ran(6)

do i = 1, beam_init%n_particle
   p => bunch%particle(i)
   p%vec = p%vec + center
enddo

end subroutine recenter_bunch


!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine init_spin_distribution (beam_init, bunch)
!
! Initializes a spin distribution according to beam_init%spin.
! Note: For non 100% polarizations this routine uses a crude model where
! (1+p)/2 of the spins are pointing in the (theta, phi) direction and the rest
! are pointing in the opposite direction.
!
! Input:
!   beam_init -- beam_init_struct: Initialization parameters 
!     %spin%theta         -- Theta orientation angle.
!     %spin%phi           -- Phi orientation angle.
!     %spin%polarization  -- Polarization percentage in (theta, phi) direction.
!
! Output:
!  bunch    -- bunch_struct: Bunch of particles.
!   %particle(:)%spin
!-

subroutine init_spin_distribution (beam_init, bunch)

implicit none

type (beam_init_struct) beam_init
type (bunch_struct) bunch
type (spin_polar_struct) :: polar

real(rp) pol
integer i, n_diff

! This is a crude way to get the desired polarization.

polar%xi = 0.0 ! spinor phase is zero
n_diff = 0
pol = beam_init%spin%polarization 

do i = 1, size(bunch%particle)
  ! Create spin up if creating one will get us nearer to the desired polarization.
  if (abs(n_diff + 1 - pol * i) < abs(n_diff - 1 - pol * i)) then
    polar%theta = beam_init%spin%theta
    polar%phi = beam_init%spin%phi
    n_diff = n_diff + 1
  ! Else create spin down.
  else
    polar%theta = beam_init%spin%theta + pi
    polar%phi = beam_init%spin%phi
    n_diff = n_diff - 1
  endif

  call polar_to_spinor (polar, bunch%particle(i))
enddo

end subroutine init_spin_distribution

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! subroutine calc_bunch_params_slice (bunch, bunch_params, 
!                           plane, slice_center, slice_spread, err, print_err)
!
! Finds all bunch parameters for a slice through the beam distribution.
!
! Modules needed:
!  use beam_mod
!
! Input:
!   bunch        -- Bunch_struct
!   plane        -- Integer: plane to slice through (x$, px$, & etc...)
!   slice_center -- Real(rp): Center to take slice about
!   slice_spread -- Real(rp): hard-wall spread in slice about center
!   print_err -- Logical, optional: If present and False then suppress 
!                  "no eigen-system found" messages.
!
! Output     
!   params -- bunch_params_struct:
!   err    -- Logical: Set True if there is an error in mat_eigen routine.
! -

subroutine calc_bunch_params_slice (bunch, bunch_params, &
                          plane, slice_center, slice_spread, err, print_err)

implicit none

type (bunch_struct), intent(in) :: bunch
type (beam_struct) :: beam
type (bunch_params_struct) bunch_params

real(rp) slice_center, slice_spread

integer plane
integer i, n_part

logical, optional :: print_err
logical err

!

n_part = 0
do i = 1, size(bunch%particle)
  if (bunch%particle(i)%vec(plane) .le. slice_center + abs(slice_spread) .and. &
      bunch%particle(i)%vec(plane) .ge. slice_center - abs(slice_spread)) &
            n_part = n_part + 1
enddo

call reallocate_beam (beam, 1, n_part)

beam%bunch(1)%charge = bunch%charge
beam%bunch(1)%z_center = bunch%z_center
beam%bunch(1)%t_center = bunch%t_center

n_part = 1
do i = 1, size(bunch%particle)
  if (bunch%particle(i)%vec(plane) .le. slice_center + abs(slice_spread) .and. &
      bunch%particle(i)%vec(plane) .ge. slice_center - abs(slice_spread)) then
            beam%bunch(1)%particle(n_part) = bunch%particle(i)
            n_part = n_part + 1
  endif
enddo

call calc_bunch_params (beam%bunch(1), bunch_params, err, print_err)

end subroutine calc_bunch_params_slice

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine calc_bunch_params (bunch, bunch_params, err, print_err)
!
! Finds all bunch parameters defined in bunch_params_struct, both normal-mode
! and projected. Projected parameters are found purely from the geometrical
! distribution of the beam. Normal-Mode parameters are found using the method
! developed in:
!   "Alternate approach to general coupled linear optics" 
!    A. Wolski, PRST AB 9, 024001 (2006)
!
! Note: If less than two particle remain then the various parameters will be
! set to zero.
! 
! Modules needed:
!  use beam_mod
!
! Input:
!   bunch     -- Bunch_struct
!   print_err -- Logical, optional: If present and False then suppress 
!                  "no eigen-system found" messages.
!
! Output     
!   bunch_params -- bunch_params_struct:
!     %a,%b,%z       -- Projected parameters
!       %alpha; %beta; %gamma
!       %eta, %etap, %norm_emit
!     %a,%b,%c       -- Normal-Mode parameters
!       %alpha; %beta; %gamma
!       %eta, %etap, %norm_emit
!     %sigma(1:21)   -- Projected Sigma Matrix
!     %centroid
!     %spin      
!       %polarization -- Polarization
!       %theta        -- Polar Angle of polarization vector
!       %phi          -- Polar Angle of polarization vector
!     %n_particle_tot
!     %n_particle_live
!     %charge_live
!   err   -- Logical: Set True if there is an error in mat_eigen routine.
!-

subroutine calc_bunch_params (bunch, bunch_params, err, print_err)

implicit none

type (bunch_struct), intent(in) :: bunch
type (bunch_params_struct) bunch_params

real(rp) exp_x2, exp_px2, exp_x_px, exp_x_d, exp_px_d
real(rp) avg_energy, temp6(6), eta, etap
real(rp) :: sigma_s(6,6), s(6,6), sigma_s_save(6,6) = 0.0, sigma(6,6) = 0.0
real(rp) :: d_r(6) = 0.0, d_i(6) = 0.0, e_r(6,6) = 0.0, e_i(6,6) = 0.0
real(rp) :: u(6,6), n_real(6,6), beta_66_iii, charge_live
real(rp), allocatable, save :: charge(:)

complex(rp) :: sigma_s_complex(6,6) = 0.0
complex(rp) :: n(6,6), e(6,6), q(6,6)

integer i, j

logical, optional :: print_err
logical err, err1

character(18) :: r_name = "calc_bunch_params"

! Init

s = 0.0

s(1,2) =  1.0 
s(2,1) = -1.0
s(3,4) =  1.0 
s(4,3) = -1.0
s(5,6) =  1.0 
s(6,5) = -1.0

call re_allocate (charge, size(bunch%particle))

! n_particle and centroid

bunch_params%n_particle_tot = size(bunch%particle)
bunch_params%n_particle_live = count(bunch%particle%state == alive$)
bunch_params%charge_live = sum(bunch%particle%charge, mask = (bunch%particle%state == alive$))

bunch_params%centroid%e_field_x = sum(bunch%particle%e_field_x, mask = (bunch%particle%state == alive$))
bunch_params%centroid%e_field_y = sum(bunch%particle%e_field_y, mask = (bunch%particle%state == alive$))

if (bunch%species == photon$) then
  charge = bunch%particle%e_field_x**2 + bunch%particle%e_field_y**2
else
  charge = bunch%particle%charge
endif

charge_live = sum(charge, mask = (bunch%particle%state == alive$))

!

if (charge_live == 0) then
  bunch_params%centroid%vec = 0.0     ! zero everything
  bunch_params%sigma = 0
  call zero_plane (bunch_params%x)
  call zero_plane (bunch_params%y)
  call zero_plane (bunch_params%z)
  call zero_plane (bunch_params%a)
  call zero_plane (bunch_params%b)
  call zero_plane (bunch_params%c)
  err = .false.
  return
endif
  
if (bmad_com%spin_tracking_on) call calc_spin_params (bunch, bunch_params)
  
! average the energy

avg_energy = sum((1+bunch%particle%vec(6)) * charge, mask = (bunch%particle%state == alive$))
avg_energy = avg_energy * bunch%particle(1)%p0c / charge_live

! Convert to geometric coords and find the sigma matrix

call find_bunch_sigma_matrix (bunch%particle, charge, bunch_params%centroid%vec, bunch_params%sigma, sigma_s)

! X, Y, & Z Projected Parameters
call projected_twiss_calc ('X', bunch_params%x, bunch_params%sigma(s11$), bunch_params%sigma(s22$), &
                      bunch_params%sigma(s12$), bunch_params%sigma(s16$), bunch_params%sigma(s26$))

call projected_twiss_calc ('Y', bunch_params%y, bunch_params%sigma(s33$), bunch_params%sigma(s44$), &
                      bunch_params%sigma(s34$), bunch_params%sigma(s36$), bunch_params%sigma(s46$))

call projected_twiss_calc ('Z', bunch_params%z, bunch_params%sigma(s55$), bunch_params%sigma(s66$), &
                      bunch_params%sigma(s56$), 0.0_rp, 0.0_rp)
     
! Normal-Mode Parameters.
! Use Andy Wolski's eigemode method to find normal-mode beam parameters.
! find eigensystem of sigma.S 

sigma_s_save = sigma_s
call mat_eigen (sigma_s, d_r, d_i, e_r, e_i, err, print_err)
if (err) return

! The eigen-values of Sigma.S are the normal-mode emittances (eq. 32)

bunch_params%a%emit = d_i(1)
bunch_params%b%emit = d_i(3)
bunch_params%c%emit = d_i(5)

bunch_params%a%norm_emit = d_i(1) * (avg_energy/mass_of(bunch%species))
bunch_params%b%norm_emit = d_i(3) * (avg_energy/mass_of(bunch%species))
bunch_params%c%norm_emit = d_i(5) * (avg_energy/mass_of(bunch%species))

! Now find normal-mode sigma matrix and twiss parameters
! N = E.Q from eq. 44

e(1,:) = e_r(1,:) + i_imaginary * e_i(1,:)
e(2,:) = e_r(2,:) + i_imaginary * e_i(2,:)
e(3,:) = e_r(3,:) + i_imaginary * e_i(3,:)
e(4,:) = e_r(4,:) + i_imaginary * e_i(4,:)
e(5,:) = e_r(5,:) + i_imaginary * e_i(5,:)
e(6,:) = e_r(6,:) + i_imaginary * e_i(6,:)

! Eq. 14
! mat_eigen finds row vectors, so switch to column vectors

call normalize_e (e, err)
if (err) return
e = transpose(e)

q = 0.0
q(1,1) = 1.0/sqrt(2.0)
q(2,1) = 1.0/sqrt(2.0)
q(3,3) = 1.0/sqrt(2.0)
q(4,3) = 1.0/sqrt(2.0)
q(5,5) = 1.0/sqrt(2.0)
q(6,5) = 1.0/sqrt(2.0)
q(1,2) =  i_imaginary / sqrt(2.0) 
q(2,2) = -i_imaginary / sqrt(2.0)
q(3,4) =  i_imaginary / sqrt(2.0) 
q(4,4) = -i_imaginary / sqrt(2.0)
q(5,6) =  i_imaginary / sqrt(2.0) 
q(6,6) = -i_imaginary / sqrt(2.0)

! compute N in eq. 44
n = matmul(e,q)
! N is now a real matrix
n_real = real(n)

! Twiss parameters come from equations 59, 63 and 64

bunch_params%a%beta = n_real(1,1)**2 + n_real(1,2)**2
bunch_params%b%beta = n_real(3,3)**2 + n_real(3,4)**2
bunch_params%c%beta = n_real(5,5)**2 + n_real(5,6)**2

bunch_params%a%alpha = -(n_real(1,1)*n_real(2,1) + n_real(1,2)*n_real(2,2))
bunch_params%b%alpha = -(n_real(3,3)*n_real(4,3) + n_real(3,4)*n_real(4,4))
bunch_params%c%alpha = -(n_real(5,5)*n_real(6,5) + n_real(5,6)*n_real(6,6))

bunch_params%a%gamma = n_real(2,1)**2 + n_real(2,2)**2
bunch_params%b%gamma = n_real(4,3)**2 + n_real(4,4)**2
bunch_params%c%gamma = n_real(6,5)**2 + n_real(6,6)**2

! Dispersion comes from equations 69 and 70

beta_66_iii   = n_real(6,5)*n_real(6,5) + n_real(6,6)*n_real(6,6)

bunch_params%a%eta  = n_real(1,5)*n_real(6,5) + n_real(1,6)*n_real(6,6)
bunch_params%a%etap = n_real(2,5)*n_real(6,5) + n_real(2,6)*n_real(6,6)

bunch_params%b%eta  = n_real(3,5)*n_real(6,5) + n_real(3,6)*n_real(6,6)
bunch_params%b%etap = n_real(4,5)*n_real(6,5) + n_real(4,6)*n_real(6,6)

bunch_params%c%eta  = n_real(5,5)*n_real(6,5) + n_real(5,6)*n_real(6,6)
bunch_params%c%etap = n_real(6,5)*n_real(6,5) + n_real(6,6)*n_real(6,6)

!----------------------------------------------------------------------
contains
subroutine zero_plane (twiss)

implicit none

type (twiss_struct), intent(out) :: twiss

twiss%beta       = 0
twiss%alpha      = 0
twiss%gamma      = 0
twiss%eta        = 0
twiss%etap       = 0
twiss%norm_emit  = 0
twiss%emit       = 0

end subroutine zero_plane
  
!----------------------------------------------------------------------
! contains

subroutine projected_twiss_calc (plane, twiss, exp_x2, exp_px2, exp_x_px, exp_x_d, exp_px_d)

implicit none

type (twiss_struct) :: twiss

real(rp), intent(in) :: exp_x2, exp_px2, exp_x_px, exp_x_d, exp_px_d
real(rp) emit, x2, x_px, px2

logical err

character(*) plane

!

if (bunch_params%sigma(s66$) /= 0) then
  twiss%eta   = exp_x_d / bunch_params%sigma(s66$)
  twiss%etap  = exp_px_d / bunch_params%sigma(s66$)
endif

if (bunch_params%sigma(s66$) == 0) then
  x2   = exp_x2 
  x_px = exp_x_px 
  px2  = exp_px2  
else
  x2   = exp_x2 - exp_x_d**2 / bunch_params%sigma(s66$)
  x_px = exp_x_px - exp_x_d * exp_px_d / bunch_params%sigma(s66$)
  px2  = exp_px2  - exp_px_d**2 / bunch_params%sigma(s66$)
endif

emit = sqrt(max(0.0_rp, x2*px2 - x_px**2)) ! Roundoff may give negative argument.

twiss%emit      = emit
twiss%norm_emit = (avg_energy/mass_of(bunch%species)) * emit

if (emit /= 0) then
  twiss%alpha = -x_px / emit
  twiss%beta  = x2 / emit
  twiss%gamma = px2 / emit
endif

end subroutine projected_twiss_calc

!----------------------------------------------------------------------
! contains
! Eq. 14 But using row vectors to conform to BMAD's mat_eigen

subroutine normalize_e (e, err)

implicit none

complex(rp) :: s(6,6)

complex(rp) :: e(6,6), temp(6)
complex(rp) :: wronsk, factor

integer i, j, k

logical err

!

err = .true.

s = 0.0
s(1,2) = ( 1.0,0.0) 
s(2,1) = (-1.0,0.0)
s(3,4) = ( 1.0,0.0) 
s(4,3) = (-1.0,0.0)
s(5,6) = ( 1.0,0.0) 
s(6,5) = (-1.0,0.0)

do i = 1, 6, 2
  e(i,:) = conjg(e(i+1,:)) ! Eq. 14b
  ! set up the normaization factor
  temp = matmul(s,e(i+1,:))
  wronsk = 0.0
  do j = 1, 6
    wronsk = e(i,j)*temp(j) + wronsk
  enddo
  if (wronsk == 0) return ! Error   
  factor = sqrt(i_imaginary) / sqrt(wronsk)
  ! this next step is the actual normalization (Eq. 14a)
  e(i+1,:) = e(i+1,:) * factor 
  e(i,:) = conjg(e(i+1,:))
enddo

err = .false.

end subroutine
  
end subroutine calc_bunch_params
  
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine calc_spin_params (bunch, bunch_params)
!
! Rotine to calculate spin averages
!
! Module needed:
!   use beam_utils
!
! Input:
!   bunch -- bunch_struct: Bunch of spins
!
! Output:
!   bunch_params -- bunch_param_struct: Structure holding average
!     %spin%polarization  -- Polarization (0.0 - 1.0)
!     %spin%theta         -- Theta polar angle of average polarization. 
!     %spin%phi           -- Phi polar angle of average polarization.
subroutine calc_spin_params (bunch, bunch_params)

implicit none

type (bunch_params_struct) bunch_params
type (bunch_struct) bunch
type (spin_polar_struct) polar, ave_polar

real(rp) vec(3), ave_vec(3), charge_live

integer i

! polarization vector

bunch_params%spin%theta = 0.0
bunch_params%spin%phi   = 0.0
charge_live = 0

ave_vec = 0.0
do i = 1, size(bunch%particle)
  if (bunch%particle(i)%state /= alive$) cycle
  call spinor_to_vec (bunch%particle(i), vec)
  ave_vec = ave_vec + vec * bunch%particle(i)%charge
  charge_live = charge_live + bunch%particle(i)%charge
enddo

ave_vec = ave_vec / charge_live
call vec_to_polar (ave_vec, ave_polar)
bunch_params%spin%theta = ave_polar%theta
bunch_params%spin%phi   = ave_polar%phi

! polarization

bunch_params%spin%polarization = sqrt(dot_product(ave_vec, ave_vec))

end subroutine calc_spin_params

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine find_bunch_sigma_matrix (particle, charge, ave, sigma, sigma_s)
!
! Routine to find the sigma matrix elements of a particle distribution.
! 
! Modules needed:
!   use beam_mod
!
! Input:
!   particle(:) -- Coord_struct: Array of particles.
!   charge(:)   -- real(rp): Particle charge or photon intensity.
! Output:
!   sigma(21)    -- Real(rp): Sigma matrix elements.
!   ave(6)       -- Real(rp): Bunch Centroid.
!   sigma_S(6,6) -- Sigma x S matrix for Wolski normal-modes
!-

subroutine find_bunch_sigma_matrix (particle, charge, avg, sigma, sigma_s)

implicit none

type (coord_struct) :: particle(:)

real(rp) charge_live, avg(6), sigma(21)
real(rp) sigma_s(6,6), s(6,6), charge(:)

integer i

!

charge_live = sum(charge, mask = (particle%state == alive$))

do i = 1, 6
  avg(i) = sum(particle(:)%vec(i) * charge, mask = (particle(:)%state == alive$)) / charge_live
enddo

sigma(s11$) = exp_calc (particle, charge, 1, 1, avg)
sigma(s12$) = exp_calc (particle, charge, 1, 2, avg)
sigma(s13$) = exp_calc (particle, charge, 1, 3, avg)
sigma(s14$) = exp_calc (particle, charge, 1, 4, avg)
sigma(s15$) = exp_calc (particle, charge, 1, 5, avg)
sigma(s16$) = exp_calc (particle, charge, 1, 6, avg)
sigma(s22$) = exp_calc (particle, charge, 2, 2, avg)
sigma(s23$) = exp_calc (particle, charge, 2, 3, avg)
sigma(s24$) = exp_calc (particle, charge, 2, 4, avg)
sigma(s25$) = exp_calc (particle, charge, 2, 5, avg)
sigma(s26$) = exp_calc (particle, charge, 2, 6, avg)
sigma(s33$) = exp_calc (particle, charge, 3, 3, avg)
sigma(s34$) = exp_calc (particle, charge, 3, 4, avg)
sigma(s35$) = exp_calc (particle, charge, 3, 5, avg)
sigma(s36$) = exp_calc (particle, charge, 3, 6, avg)
sigma(s44$) = exp_calc (particle, charge, 4, 4, avg)
sigma(s45$) = exp_calc (particle, charge, 4, 5, avg)
sigma(s46$) = exp_calc (particle, charge, 4, 6, avg)
sigma(s55$) = exp_calc (particle, charge, 5, 5, avg)
sigma(s56$) = exp_calc (particle, charge, 5, 6, avg)
sigma(s66$) = exp_calc (particle, charge, 6, 6, avg)

! make sigma.S matrix

sigma_s(1,1) = sigma(s11$)
sigma_s(1,2) = sigma(s12$)
sigma_s(1,3) = sigma(s13$)
sigma_s(1,4) = sigma(s14$)
sigma_s(1,5) = sigma(s15$)
sigma_s(1,6) = sigma(s16$)
sigma_s(2,1) = sigma(s12$)
sigma_s(2,2) = sigma(s22$)
sigma_s(2,3) = sigma(s23$)
sigma_s(2,4) = sigma(s24$)
sigma_s(2,5) = sigma(s25$)
sigma_s(2,6) = sigma(s26$)
sigma_s(3,1) = sigma(s13$)
sigma_s(3,2) = sigma(s23$)
sigma_s(3,3) = sigma(s33$)
sigma_s(3,4) = sigma(s34$)
sigma_s(3,5) = sigma(s35$)
sigma_s(3,6) = sigma(s36$)
sigma_s(4,1) = sigma(s14$)
sigma_s(4,2) = sigma(s24$)
sigma_s(4,3) = sigma(s34$)
sigma_s(4,4) = sigma(s44$)
sigma_s(4,5) = sigma(s45$)
sigma_s(4,6) = sigma(s46$)
sigma_s(5,1) = sigma(s15$)
sigma_s(5,2) = sigma(s25$)
sigma_s(5,3) = sigma(s35$)
sigma_s(5,4) = sigma(s45$)
sigma_s(5,5) = sigma(s55$)
sigma_s(5,6) = sigma(s56$)
sigma_s(6,1) = sigma(s16$)
sigma_s(6,2) = sigma(s26$)
sigma_s(6,3) = sigma(s36$)
sigma_s(6,4) = sigma(s46$)
sigma_s(6,5) = sigma(s56$)
sigma_s(6,6) = sigma(s66$)

s = 0.0

s(1,2) =  1.0 
s(2,1) = -1.0
s(3,4) =  1.0 
s(4,3) = -1.0
s(5,6) =  1.0 
s(6,5) = -1.0

sigma_s = matmul(sigma_s, s)

!----------------------------------------------------------------------
contains

function exp_calc (particle, charge, ix1, ix2, avg) result (this_sigma)

implicit none

type (coord_struct) particle(:)
real(rp) charge(:), avg(:)
real(rp) this_sigma

integer ix1, ix2

!
                                    
this_sigma = sum((particle(:)%vec(ix1) - avg(ix1)) * (particle(:)%vec(ix2) - avg(ix2)) * charge(:), &
                               mask = (particle%state == alive$))

this_sigma = this_sigma / charge_live

end function exp_calc

end subroutine find_bunch_sigma_matrix 
                                    
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine bunch_equal_bunch (bunch1, bunch2)
!
! Subroutine to set one particle bunch equal to another taking care of
! pointers so that they don't all point to the same place.
!
! Note: This subroutine is called by the overloaded equal sign:
!    bunch1 = bunch2
!
! Input: 
!   bunch2 -- bunch_struct: Input bunch
!
! Output
!   bunch1 -- bunch_struct: Output bunch
!-

subroutine bunch_equal_bunch (bunch1, bunch2)

implicit none

type (bunch_struct), intent(inout) :: bunch1
type (bunch_struct), intent(in)    :: bunch2

integer i, n_particle

!

n_particle = size(bunch2%particle)

if (size(bunch1%particle) /= size(bunch2%particle)) then
  deallocate (bunch1%particle)
  allocate (bunch1%particle(n_particle))
endif

bunch1%particle  = bunch2%particle
bunch1%charge    = bunch2%charge
bunch1%z_center  = bunch2%z_center
bunch1%t_center  = bunch2%t_center

end subroutine bunch_equal_bunch

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine beam_equal_beam (beam1, beam2)
!
! Subroutine to set one particle beam equal to another taking care of
! pointers so that they don't all point to the same place.
!
! Note: This subroutine is called by the overloaded equal sign:
!    beam1 = beam2
!
! Input: 
!  beam2 -- beam_struct: Input beam
!
! Output
!  beam1 -- beam_struct: Output beam
!
!-

subroutine beam_equal_beam (beam1, beam2)

implicit none

type (beam_struct), intent(inout) :: beam1
type (beam_struct), intent(in)    :: beam2

integer i, j, n_bun, n_particle
logical allocate_this

! The following rule must be observed: If beam%bunch is allocated then
! beam%bunch%particle must be also.

n_bun = size(beam2%bunch)

allocate_this = .true.
if (allocated(beam1%bunch)) then
  if (size(beam1%bunch) /= size(beam2%bunch)) then
    do i = 1, size(beam1%bunch)
      deallocate (beam1%bunch(i)%particle)
    enddo
    deallocate (beam1%bunch)
  else
    allocate_this = .false.
  endif
endif

if (allocate_this) then
  allocate (beam1%bunch(n_bun))
  do i = 1, n_bun
    n_particle = size(beam2%bunch(i)%particle)
    allocate (beam1%bunch(i)%particle(n_particle))
  enddo
endif

do i = 1, n_bun
  beam1%bunch(i) = beam2%bunch(i)
enddo

end subroutine beam_equal_beam

end module
