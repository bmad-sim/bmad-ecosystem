module beam_utils

  use beam_def_struct
  use bmad_interface
  use spin_mod
  use eigen_mod
  use wake_mod

type phase_space_struct    ! this is only used in beam_utils
   integer n_particle
   real(rp), allocatable :: vec(:,:)
end type

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track1_particle (start, ele, param, end)
!
! Subroutine to track a particle through an element.
!
! Modules needed:
!   use beam_mod
!
! Input:
!   start  -- struct: Starting coords.
!   ele    -- Ele_struct: Element to track through.
!   param  -- lat_param_struct: Global parameters.
!
! Output:
!   end    -- struct: Ending coords.
!-

subroutine track1_particle (start, ele, param, end)

  implicit none

  type (particle_struct) :: start
  type (particle_struct) :: end
  type (ele_struct) :: ele
  type (lat_param_struct), intent(inout) :: param

! transfer z-order index, charge, etc

  end = start
  if (start%ix_lost /= not_lost$) return
  if (ele%key == marker$ .or. ele%key == photon_branch$ .or. ele%key == branch$) return

  call track1 (start%r, ele, param, end%r)
  if (param%lost) end%ix_lost = ele%ix_ele

  if (end%ix_lost /= not_lost$) then
    end%r%vec = 0
    end%charge = 0
    return
  endif

end subroutine

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
! bmad_com%grad_loss_sr_wake is used to tell track1_bmad the appropriate loss
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
type (ele_struct) ele

type (lat_param_struct) param

real(rp), pointer :: a_pole_save(:)
real(rp) charge, value_save(n_attrib_maxx)
integer i, j, n, ix_z

character(20) :: r_name = 'track1_bunch_hom'

! Charge and center

bunch_end = bunch_start

!------------------------------------------------
! Without wakefields just track through.

if (ele%key /= lcavity$ .or. .not. associated(ele%wake) .or. &
            (.not. bmad_com%sr_wakes_on .and. .not. bmad_com%lr_wakes_on)) then

  do j = 1, size(bunch_start%particle)
    call track1_particle (bunch_start%particle(j), &
                                      ele, param, bunch_end%particle(j))
  enddo



  bunch_end%charge = sum (bunch_end%particle(:)%charge, &
                      mask = (bunch_end%particle(:)%ix_lost == not_lost$))
  return
endif

!------------------------------------------------
! This calculation is for an lcavity with wakefields.
! Put the sr wakefield transverse kicks at the half way point.

! first offset the cavity
! wakes applied in cononical coords so don't do canonical coord conversion

do i = 1, size(bunch_end%particle)
  call offset_particle (ele, param, bunch_end%particle(i)%r, set$, &
      set_canonical = .false.)
enddo

! Modify ele temporarily so we can track through half the cavity.

value_save = ele%value
ele%value(l$)      = ele%value(l$) / 2
ele%value(e_tot$)  = (ele%value(e_tot_start$) + ele%value(e_tot$)) / 2
ele%value(p0c$)    = (ele%value(p0c_start$) + ele%value(p0c$)) / 2
ele%value(e_loss$) = 0

! zero all offsets and kicks (offsetting already performed above)

call zero_ele_offsets (ele)
ele%value(hkick$) = 0.0
ele%value(vkick$) = 0.0
if (associated(ele%a_pole)) then
   a_pole_save => ele%a_pole
   nullify(ele%a_pole)
else
   nullify(a_pole_save)
endif
    
! Track half way through. This includes the sr longitudinal wakes 

call order_particles_in_z (bunch_end)
do j = 1, size(bunch_end%particle)
  ix_z = bunch_end%particle(j)%ix_z ! z-ordered index of the particles
  if (bunch_end%particle(ix_z)%ix_lost /= not_lost$) cycle
  call add_sr_long_wake (ele, bunch_end, j-1, ix_z)
  call track1_particle (bunch_end%particle(ix_z), ele, param, bunch_end%particle(ix_z))
enddo

bmad_com%grad_loss_sr_wake = 0.0

! Put in the transverse wakefields

call track1_sr_wake (bunch_end, ele)
call track1_lr_wake (bunch_end, ele)

! Track the last half of the lcavity.  This includes the sr longitudinal wakes 

ele%value(e_tot_start$)  = ele%value(e_tot$)
ele%value(p0c_start$)    = ele%value(p0c$)
ele%value(e_tot$)        = value_save(e_tot$)
ele%value(p0c$)          = value_save(p0c$)

call order_particles_in_z (bunch_end)
do j = 1, size(bunch_end%particle)
  ix_z = bunch_end%particle(j)%ix_z ! z-ordered index of the particles
  if (bunch_end%particle(ix_z)%ix_lost /= not_lost$) cycle
  call add_sr_long_wake (ele, bunch_end, j-1, ix_z)
  call track1_particle (bunch_end%particle(ix_z), ele, param, bunch_end%particle(ix_z))
enddo

bmad_com%grad_loss_sr_wake = 0.0

bunch_end%charge = sum (bunch_end%particle(:)%charge, &
                         mask = (bunch_end%particle(:)%ix_lost == not_lost$))

! Unset the cavity offset.
! Wakes applied in cononical coords so don't do canonical coord conversion

do i = 1, size(bunch_end%particle)
  call offset_particle (ele, param, bunch_end%particle(i)%r, unset$, &
      set_canonical = .false.)
enddo

! Unmodify ele

ele%value = value_save
if (associated(a_pole_save)) ele%a_pole => a_pole_save
      
end subroutine track1_bunch_hom

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine add_sr_long_wake (ele, bunch, num_in_front, follower)
!
! Adds the longitudinal wake for all particles in front of the follower.
!
! Input:
!  ele      -- Ele_struct: Element with wakefields.
!  bunch    -- Bunch_struct: Bunch of particles
!  num_in_front -- Integer: number of particles in front of this one
!                   This will be the bunch%particle index number right before
!                   the follower
!  follower -- Integer: index of particle wakes being applied to.
!
! Output:
!  bmad_com%grad_loss_sr_wake -- Real(rp): net gradient loss due to the leaders.
!-

subroutine add_sr_long_wake (ele, bunch, num_in_front, ix_follower)

implicit none

type (ele_struct) ele
type (bunch_struct) bunch
type (coord_struct), pointer :: leader

integer ix_follower, i, num_in_front
integer n_sr_table, n_sr_mode_long, n_sr_mode_trans, k_start

!-----------------------------------
! If there is no wake for this element, or the sr wakes are turned off, then just 
! use the e_loss attribute (as set in track1_bmad).

bmad_com%grad_loss_sr_wake = 0.0

n_sr_table = size(ele%wake%sr_table) 
n_sr_mode_long = size(ele%wake%sr_mode_long)
n_sr_mode_trans = size(ele%wake%sr_mode_trans)

if ((n_sr_table == 0 .and. n_sr_mode_long == 0 .and. n_sr_mode_trans == 0) .or. &
                                          .not. bmad_com%sr_wakes_on) then 
  bmad_com%grad_loss_sr_wake = 0.0
  return 
endif

! the self wake only sees the charge of each real particle, not the "macro"
! charge of the simulated particle

if (n_sr_table > 0) then

  bmad_com%grad_loss_sr_wake = bmad_com%grad_loss_sr_wake &
                        + ele%wake%sr_table(0)%long * e_charge / 2.0

  !-----------------------------------
  ! add up all wakes from front of bunch to follower

  do i = 1, num_in_front
    if (bunch%particle(bunch%particle(i)%ix_z)%ix_lost == not_lost$) &
      call sr_table_add_long_kick (ele, bunch%particle(bunch%particle(i)%ix_z)%r, &
               bunch%particle(bunch%particle(i)%ix_z)%charge, &
               bunch%particle(ix_follower)%r)
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
type (particle_struct), pointer :: particle, leader
type (particle_struct), pointer :: p(:)

real(rp) dz_sr_table, sr02, z_sr_table_max
integer i, j, k, i1, i2, i_sr_mode, n_sr_table, n_sr_mode_long, n_sr_mode_trans, k_start

logical wake_here
character(16) :: r_name = 'track1_sr_wake'

!-----------------------------------

if (.not. associated(ele%wake)) return
  
p => bunch%particle
  
! error check and zero wake sums and order particles in z

call order_particles_in_z (bunch)  
if (size(ele%wake%sr_mode_long) /= 0) then
  i1 = p(1)%ix_z 
  i2 = p(size(p))%ix_z
  if (p(i1)%r%vec(5) - p(i2)%r%vec(5) > ele%wake%z_sr_mode_max) then
    call out_io (s_abort$, r_name, &
        'Bunch longer than sr_mode wake can handle for element: ' // ele%name)
    call err_exit
  endif
endif

do i = 1, size(ele%wake%sr_mode_long)
  ele%wake%sr_mode_long%b_sin = 0
  ele%wake%sr_mode_long%b_cos = 0
  ele%wake%sr_mode_long%a_sin = 0
  ele%wake%sr_mode_long%a_cos = 0
enddo

!

n_sr_table = size(ele%wake%sr_table) 
z_sr_table_max = 0
if (n_sr_table > 0) then
  z_sr_table_max = ele%wake%sr_table(n_sr_table-1)%z
  dz_sr_table = z_sr_table_max / (n_sr_table - 1)
endif

! loop over all particles in the bunch and apply the wake

i_sr_mode = 1  ! index of next particle to be added to the sr_mode wake sums.

do j = 1, size(p)
  particle => p(p(j)%ix_z)
  ! apply longitudinal self wake

  if (z_sr_table_max .ge. 0) then
    call sr_mode_long_self_wake_apply_kick (ele, particle%charge, particle%r)
  endif

  ! Particle_j is kicked by particles k = 1, ..., j-1.
  ! The particles 1, ... i_sr_mode-1 have already had their wakes added to the 
  ! sr_mode wake sums so the loop is from i_sr_mode, ..., j-1.

  k_start = i_sr_mode
  do k = k_start, j-1
    leader => p(p(k)%ix_z)
    if ((particle%r%vec(5) - leader%r%vec(5)) > z_sr_table_max) then
      ! use sr_table table to add to particle j the wake of particle k
      call sr_table_apply_trans_kick (ele, leader%r, leader%charge, particle%r)
    else
      ! add contribution of particle(k) to wake sums
      i_sr_mode = k  ! update i_sr_mode
      call sr_mode_long_wake_add_to (ele, leader%r, leader%charge)
      call sr_mode_trans_wake_add_to(ele, leader%r, leader%charge)
    endif
  enddo

  ! apply wake to particle(j)
  call sr_mode_long_wake_apply_kick (ele, particle%r)
  call sr_mode_trans_wake_apply_kick(ele, particle%r)

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
!     ele%wake%lr%b_sin
!     ele%wake%lr%b_cos
!     ele%wake%lr%a_sin
!     ele%wake%lr%a_cos
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
type (particle_struct), pointer :: particle

integer n_mode, j, k

if (.not. bmad_com%lr_wakes_on) return
if (.not. associated(ele%wake)) return
  
! Check to see if we need to do any calc

if (.not. associated(ele%wake)) return
n_mode = size(ele%wake%lr)
if (n_mode == 0 .or. .not. bmad_com%lr_wakes_on) return  

call order_particles_in_z (bunch)  ! needed for wakefield calc.

! Give the particles a kick

do k = 1, size(bunch%particle)
  j = bunch%particle(k)%ix_z
  particle => bunch%particle(j)
  if (particle%ix_lost /= not_lost$) cycle
  call lr_wake_apply_kick (ele, bunch%t_center, particle%r)
enddo

! Add the wakes left by this bunch to the existing wakes.

do k = 1, size(bunch%particle)
  j = bunch%particle(k)%ix_z
  particle => bunch%particle(j)
  if (particle%ix_lost /= not_lost$) cycle
  call lr_wake_add_to (ele, bunch%t_center, particle%r, particle%charge)
enddo

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine order_particles_in_z (bunch)
!
! Subroutine to order the particles longitudinally 
! The ordering uses the centroid of the particles:
!   %vec(5) 
!
! Modules needed:
!   use beam_mod
!
! Input:
!   bunch     -- Bunch_struct: collection of particles.
!     %particle(j)%r%vec(5) -- Longitudinal position of j^th particle.
!
! Output:
!   bunch     -- bunch_struct: collection of particles.
!     %particle(j) -- particle ordered using %vec(5).
!                     Order is from large z (head of bunch) to small z.
!                     That is: %particle(1)%ix_z is the particle at the bunch head. 
!       %ix_z        -- Index for the ordering
!-

Subroutine order_particles_in_z (bunch)

implicit none

type (bunch_struct), target :: bunch
type (particle_struct), pointer :: particle(:)
type (particle_struct) temp
integer i, k, nm, i0, i1
real(rp) z1, z2
logical ordered

! Init if needed

particle => bunch%particle
nm = size(particle)

if (particle(1)%ix_z == 0) then
  forall (i = 1:nm) particle(i)%ix_z = i
endif

! Order is from large z (head of bunch) to small z.

do
  ordered = .true.
  do i = 1, nm-1
    i0 = particle(i)%ix_z; i1 = particle(i+1)%ix_z
    if (particle(i0)%r%vec(5) < particle(i1)%r%vec(5)) then
      particle(i:i+1)%ix_z = particle(i+1:i:-1)%ix_z
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
! Subroutine angle_to_canonical_coords (particle, energy0)
!
! Subroutine to convert particle coords from 
!     (x, x', y, y', z, E)
! to
!     (x, px, y, py, z, pz)
!
! Note: the reverse routine is called:
!   canonical_to_angle_coords (particle, energy0)
!
! Modules needed:
!   use beam_mod
!
! Input:
!   particle -- struct: particleparticle with angular coords.
!   energy0  -- real(rp): Reference energy.
!
! Output:
!   particle -- struct: particle-particle with momentum coords.
!-

subroutine angle_to_canonical_coords (particle, energy0)

implicit none

type (particle_struct), target :: particle

real(rp), pointer :: s(:)
real(rp), intent(in) :: energy0
real(rp) f, f2, e, xp0, yp0

!

f = particle%r%vec(6) / energy0
f2 = f * f
e = energy0

xp0 = particle%r%vec(2)
yp0 = particle%r%vec(4)

particle%r%vec(2) = particle%r%vec(2) * f
particle%r%vec(4) = particle%r%vec(4) * f
particle%r%vec(6) = f - 1

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine canonical_to_angle_coords (particle, energy0)
!
! Subroutine to convert particleparticle coords from 
!     (x, px, y, py, z, pz)
! to
!     (x, x', y, y', z, E)
!
! Note: the reverse routine is called:
!   angle_to_canonical_coords (particle, energy0)
!
! Modules needed:
!   use beam_mod
!
! Input:
!   particle -- struct: particle with momentum coords.
!   energy0  -- real(rp): Reference energy.
!
! Output:
!   particle -- struct: particle with angular coords.
!-

subroutine canonical_to_angle_coords (particle, energy0)

implicit none

type (particle_struct), target :: particle

real(rp), pointer :: s(:)
real(rp), intent(in) :: energy0
real(rp) f, f2, e, xp0, yp0

!

f = 1 + particle%r%vec(6)
f2 = f * f
e = energy0

particle%r%vec(2) = particle%r%vec(2) / f
particle%r%vec(4) = particle%r%vec(4) / f
particle%r%vec(6) = energy0 * f 

xp0 = particle%r%vec(2) / f2
yp0 = particle%r%vec(4) / f2

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
!
! Modules needed:
!   use beam_mod
!
! Input:
!   n_bunch -- Integer: Number of bunches.
!   n_particle -- Integer: Number of particles.
!
! Output:
!   beam -- beam_struct: Allocated beam_struct structure.
!-

subroutine reallocate_beam (beam, n_bunch, n_particle)

implicit none

type (beam_struct) beam

integer i, j
integer n_bunch, n_particle

logical de_bunch, de_particle

! Deallocate

de_bunch = .false.
de_particle = .false.

if (allocated(beam%bunch)) then
  if (n_bunch == 0) then
    de_bunch = .true.
    de_particle = .true.
  else
    if (size(beam%bunch) /= n_bunch) then
      de_bunch = .true.
      de_particle = .true.
    endif
    if (size(beam%bunch(1)%particle) /= n_particle) then
      de_particle= .true.
    endif
  endif

  do i = 1, size(beam%bunch)
    if (de_particle) deallocate (beam%bunch(i)%particle)
  enddo
  if (de_bunch) deallocate (beam%bunch)

endif

if (n_bunch == 0) return
  
! Allocate

if (.not. allocated (beam%bunch)) allocate (beam%bunch(n_bunch))
do i = 1, n_bunch
  if (.not. allocated (beam%bunch(i)%particle)) &
                    allocate (beam%bunch(i)%particle(n_particle))
enddo

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine init_beam_distribution (ele, param, beam_init, beam)
!
! Subroutine to initialize either a random or tail-weighted distribution of particles.  
!
! There are three tail-weighted distributions available: rectangular grid, 
! Gaussian represented by concentric ellipses, and Kapchinsky-Vladimirsky (KV).
! See manual for more information.
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
!     %particle      -- Type of particle.
!   beam_init   -- beam_init_struct: Use "getf beam_init_struct" for more details.
!     For Gaussian random distributions:
!       %renorm_center -- Logical: If True then distribution is rescaled to
!                      the desired centroid (to take care of
!                      possible statistical errors in distribution).
!                      (default is True)
!       %renorm_sigma  -- Logical: If True then rescale the distribution to the desired sigmas.
!                      (default is True)
!       %init_spin     -- Logical: If True then the particle spinors will be
!                      initialized with the parameters in beam_init%spin
!                      (default is False)
!       %random_engine          -- Character: "pseudo" (default) or "quasi".
!       %random_gauss_converter -- Character: "exact" (default) or "limited".
!       %random_sigma_cutoff    -- Real: Used with "limited" converter. Default is 4.0.
!     For tail-weighted distributions:
!       %tw_beam_init  -- tail_weighted_beam_init_struct: 
!                         Use "getf tail_weighted_beam_init_struct" for more details.
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
type (beam_init_struct) beam_init
type (beam_struct), target :: beam
type (bunch_struct), pointer :: bunch
  
integer i_bunch
real(rp) old_cutoff

character(16) old_engine, old_converter  
character(22) :: r_name = "init_beam_distribution"

! resize the beam to the number of particles and the number of bunches.

call reallocate_beam (beam, beam_init%n_bunch, beam_init%n_particle)

! Set and save the random number generator parameters.

if (all(beam_init%distribution_type == '')) then
  call ran_engine (beam_init%random_engine, old_engine)
  call ran_gauss_converter (beam_init%random_gauss_converter, &
                  beam_init%random_sigma_cutoff, old_converter, old_cutoff)
endif

! Loop over all bunches
! Note z_center is negative and t_center is posive for trailing bunches.

do i_bunch = 1, size(beam%bunch)

  bunch => beam%bunch(i_bunch)
  call init_bunch_distribution (ele, param, beam_init, bunch)

  bunch%t_center = (i_bunch-1) * beam_init%dt_bunch
  bunch%z_center = -bunch%t_center * c_light * ele%value(e_tot$) / ele%value(p0c$)
  bunch%ix_bunch = i_bunch

enddo
  
! Reset the random number generator parameters.

if (all(beam_init%distribution_type == '')) then
  call ran_engine (old_engine)
  call ran_gauss_converter (old_converter, old_cutoff)
endif

end subroutine init_beam_distribution

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine init_bunch_distribution (ele, param, beam_init, bunch)
!
! Subroutine to initialize either a random or tail-weighted distribution of particles.  
!
! There are three tail-weighted distributions available: rectangular grid, 
! Gaussian represented by concentric ellipses, and Kapchinsky-Vladimirsky (KV).
! See manual for more information.
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
!     %particle      -- Type of particle.
!   beam_init   -- beam_init_struct: Use "getf beam_init_struct" for more details.
!     For gaussian random distributions:
!       %renorm_center -- Logical: If True then distribution is rescaled to
!                      the desired centroid (to take care of
!                      possible statistical errors in distribution).
!                      (default is True)
!       %renorm_sigma  -- Logical: If True then rescale the distribution to the desired sigmas.
!                      (default is True)
!       %init_spin     -- Logical: If True then the particle spinors will be
!                      initialized with the parameters in beam_init%spin
!                      (default is False)
!       %random_engine          -- Character: "pseudo" (default) or "quasi".
!       %random_gauss_converter -- Character: "exact" (default) or "limited".
!       %random_sigma_cutoff    -- Real: Used with "limited" converter. Default is 4.0.
!     For tail-weighted distributions:
!       %tw_beam_init  -- tail_weighted_beam_init_struct: 
!                         Use "getf tail_weighted_beam_init_struct" for more details.
!
! Output:
!   bunch        -- bunch_struct: Structure with initialized particles.
!
!-

subroutine init_bunch_distribution (ele, param, beam_init, bunch)
 
use random_mod

implicit none

type (ele_struct) ele
type (lat_param_struct) param
type (beam_init_struct), target :: beam_init
type (tail_weighted_beam_init_struct), pointer :: tw_init
type (bunch_struct) :: bunch
  
integer i_bunch, i
logical :: is_KV = .false.
real(rp) old_cutoff

character(16) old_engine, old_converter  
character(22) :: r_name = "init_bunch_distribution"

! If we are making a tail-weighted distribution, we must compute the number of particles

if (any(beam_init%distribution_type /= '')) then
   beam_init%n_particle = 1
   tw_init => beam_init%tw_beam_init
   do i = 1, 3
      call str_upcase(beam_init%distribution_type(i), beam_init%distribution_type(i))
      select case (beam_init%distribution_type(i))
      case ('ELLIPSE')
         beam_init%n_particle = beam_init%n_particle * tw_init%n_ellipse(i) * tw_init%part_per_ellipse(i)
      case ('GRID')
         beam_init%n_particle = beam_init%n_particle * tw_init%n_x(i) * tw_init%n_px(i)
      case ('KV')
         is_KV = .true.
         beam_init%n_particle = beam_init%n_particle * tw_init%part_per_ellipse(i)
      case default
         call out_io (s_abort$, r_name, 'PHASE SPACE DISTRIBUTION TYPE NOT RECOGNIZED')
         call err_exit
      end select
   enddo
endif

! number of particles for a 4D KV distribution = n_I2 * n_phi1 * n_phi2

if (is_KV) then
   beam_init%n_particle = beam_init%n_particle * tw_init%n_I2
endif

! Initialize the bunch.

if (all(beam_init%distribution_type == '')) then
   call init_random_bunch (ele, param, beam_init, bunch)
else
   call init_tail_weighted_bunch (ele, param, beam_init, bunch)
endif

end subroutine init_bunch_distribution

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine init_random_bunch (ele, param, beam_init, bunch)
!
! Subroutine to initialize a random bunch of particles matched to
! the Twiss parameters, centroid position, and Energy - z correlation
! as specified. Coupling in the element ele is incorporated into the
! distribution.
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
!     %particle      -- Type of particle.
!   beam_init   -- beam_init_struct: Use "getf beam_init_struct" for more details.
!     %renorm_center -- Logical: If True then distribution is rescaled to
!                    the desired centroid (to take care of
!                    possible statistical errors in distribution).
!                    (default is True)
!     %renorm_sigma  -- Logical: If True then rescale the distribution to the desired sigmas.
!                    (default is True)
!     %init_spin     -- Logical: If True then the particle spinors will be
!                    initialized with the parameters in beam_init%spin
!                    (default is False)
!     %random_engine          -- Character: "pseudo" (default) or "quasi".
!     %random_gauss_converter -- Character: "exact" (default) or "limited".
!     %random_sigma_cutoff    -- Real: Used with "limited" converter. Default is 4.0.
!
! Output:
!   bunch        -- Bunch_struct: Structure with initialized particles.
!
!-

subroutine init_random_bunch (ele, param, beam_init, bunch)
 
use random_mod

implicit none

type (ele_struct) ele
type (lat_param_struct) param
type (beam_init_struct) beam_init
type (bunch_struct), target :: bunch
type (particle_struct), pointer :: p
  
real(rp) dpz_dz, denom
real(rp) a_emitt, b_emitt
real(rp) ave(6), sigma(6), a, b
real(rp) r(6), v_mat(4,4), v_inv(4,4)
real(rp) y, alpha(6), sig_mat(6,6)
real(rp) center(6) ! includes jitter
real(rp) ran(6), old_cutoff

integer i, j, j2, n, i_bunch

character(16) old_engine, old_converter  
character(22) :: r_name = "init_random_bunch"

! make sure the beam%particle array is the correct size.

if (allocated(bunch%particle)) then
  if (size(bunch%particle) /= beam_init%n_particle) deallocate (bunch%particle)
endif
if (.not. allocated(bunch%particle)) allocate (bunch%particle(beam_init%n_particle))

!

sig_mat = 0
ave = 0
do n = 1, beam_init%n_particle
  p => bunch%particle(n)
  call ran_gauss(r)
  p%r%vec = r
  ave = ave + r
  forall (i=1:6, j=1:6) sig_mat(i,j) = sig_mat(i,j) + r(i) * r(j)
enddo  

ave = ave / beam_init%n_particle
sig_mat = sig_mat / beam_init%n_particle

! Now the distribution of bunch%particle(:)%r%vec(n) for fixed n has
! on average, unit sigma and the distribution for n = n1 is uncorrelated
! with the distribution for n = n2, n1 /= n2.

! However, since we are dealing with a finite number of particles, 
! the sigmas of the distributions will not be exactly 1, and there will 
! be some correlation between distributions.
! If beam_init%renorm_sigma = True then take this out.

! Zero the average for now

do n = 1, beam_init%n_particle
  bunch%particle(n)%r%vec = bunch%particle(n)%r%vec - ave
enddo

! renormalize the beam sigmas. Ignore if n_particle = 1.

if (beam_init%renorm_sigma .and. beam_init%n_particle > 1) then

  if (beam_init%n_particle < 7) then
    call out_io (s_abort$, r_name, &
        'INITIALIZATION WITH RENORM_SIGMA MUST USE AT LEAST 7 PARTICLES!')
    call err_exit
  endif

  ! This accounts for subtracting off the average
  forall (i = 1:6, j = 1:6) sig_mat(i,j) = sig_mat(i,j) - ave(i) * ave(j)

  ! To renormalize we want to make sig_mat = the unit matrix.
  ! The first step is to zero the off-diagonal elements.
  ! We have to do this in the correct order otherwise zeroing one element
  ! might unzero others that have already been zeroed.

  do i = 5, 1, -1
    do j = i+1, 6
      b = -sig_mat(i,j) / sig_mat(j,j)
      ! Transform the distribution
      do n = 1, beam_init%n_particle
        p => bunch%particle(n)
        p%r%vec(i) = p%r%vec(i) + b * p%r%vec(j)
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

  ! Now we make the diagonal elements unity

  forall (i = 1:6) alpha(i) = sqrt(1/sig_mat(i,i))
  do n = 1, beam_init%n_particle
    p => bunch%particle(n)
    p%r%vec = p%r%vec * alpha
  enddo

endif

! In general, since we are dealing with a finite number of particles, 
! the averages will not be zero.
! Put back the non-zero center if beam_init%renorm_center = False.

if (.not. beam_init%renorm_center) then
  do n = 1, beam_init%n_particle
    bunch%particle(n)%r%vec = bunch%particle(n)%r%vec + ave
  enddo
endif

! scale by the emittances, etc. and put in jitter

call ran_gauss(ran(1:4)) ! ran(3:4) for z and e jitter used below
a_emitt = beam_init%a_norm_emitt*(1+beam_init%emitt_jitter(1)*ran(1)) &
                                        * mass_of(param%particle) / ele%value(e_tot$)
b_emitt = beam_init%b_norm_emitt*(1+beam_init%emitt_jitter(2)*ran(2)) &
                                        * mass_of(param%particle) / ele%value(e_tot$)
  
dpz_dz = beam_init%dpz_dz
  
call make_v_mats(ele, v_mat, v_inv)

sigma(1) = sqrt(a_emitt * ele%a%beta)
sigma(2) = sqrt(a_emitt / ele%a%beta)
sigma(3) = sqrt(b_emitt * ele%b%beta)
sigma(4) = sqrt(b_emitt / ele%b%beta)
sigma(5) = beam_init%sig_z * (1 + beam_init%sig_z_jitter*ran(3))
sigma(6) = beam_init%sig_e * (1 + beam_init%sig_e_jitter*ran(4))

if (sigma(6) == 0 .or. dpz_dz == 0) then
  a = 0
else if (abs(dpz_dz * sigma(5)) > sigma(6)) then
  call out_io (s_abort$, r_name, &
              "|dpz_dz| MUST be < mode%sigE_E / mode%sig_z")
  call err_exit
else
  a = dpz_dz * sigma(5) / sigma(6)
endif

b = sqrt(1-a**2)
     
! Put everything together to distribute the particles.

do i = 1, beam_init%n_particle

  p => bunch%particle(i)
  r = p%r%vec

  p%r%vec(1) = sigma(1) *  r(1)
  p%r%vec(2) = - sigma(2) * (r(2) + r(1) * ele%a%alpha)
  p%r%vec(3) = sigma(3) *  r(3)
  p%r%vec(4) = - sigma(4) * (r(4) + r(3) * ele%b%alpha)
  p%r%vec(5) = sigma(5) *  r(5)
  p%r%vec(6) = sigma(6) * (r(6) * b + r(5) * a)
      
  ! Include Dispersion
  p%r%vec(1:4) =  p%r%vec(1:4) + &
              p%r%vec(6) * (/ ele%a%eta, ele%a%etap, ele%b%eta, ele%b%etap /)
      
  ! Include Coupling
  p%r%vec(1:4) = matmul(v_mat, p%r%vec(1:4))

end do

! recenter bunch (includes beam jitter)
call recenter_bunch (beam_init, bunch)

! set particle charge

bunch%charge = beam_init%bunch_charge
bunch%particle(:)%charge = beam_init%bunch_charge / beam_init%n_particle
bunch%particle(:)%ix_lost = not_lost$
bunch%ix_ele = ele%ix_ele  ! Current position of the bunch

bunch%z_center = 0  ! Default
bunch%t_center = 0  ! Default
bunch%ix_bunch = 1  ! Default

! particle spin

call init_spin_distribution (beam_init, bunch)

end subroutine init_random_bunch

!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
!+
! Subroutine init_tail_weighted_bunch (ele, param, beam_init, bunch)
!
! Subroutine to initialize a tail_weighted bunch of particles.
!
! Note: Make sure: |beam_init%dpz_dz| < mode%sigE_E / mode%sig_z
!
! Input:
!   ele         -- Ele_struct: element to initialize distribution at
!   param       -- Lat_param_struct: Lattice parameters
!     %particle      -- Type of particle.
!   beam_init   -- beam_init_struct: Use "getf beam_init_struct" for more details.
!     %tw_beam_init -- tail_weighted_beam_init_struct: 
!                      Use "getf tail_weighted_beam_init_struct" for more details.
!
! Output:
!   bunch       -- Bunch_struct: Structure with initialized particles.
!
!
! Phase space distribution types:
!   ellipse$  -- concentric ellipses representing a Gaussian distribution
!   grid$     -- uniform rectangular grid
!   KV$       -- Kapchinsky-Vladimirsky distribution
!
! See manual for more details.   
!
!-

subroutine init_tail_weighted_bunch (ele, param, beam_init, bunch)

implicit none

type (ele_struct) ele
type (lat_param_struct) param
type (beam_init_struct), target :: beam_init
type (tail_weighted_beam_init_struct), pointer :: tw

type (bunch_struct), target :: bunch
type (particle_struct), pointer :: p

type (phase_space_struct) space2D(3), space4D, space6D

real(rp) beta(3), alpha(3), emitt(3), covar, e1, e2

integer i,j,k
integer :: KV_counter = 0     ! counts how many phase planes are of KV type
integer :: KV_indices(2) = 0  ! indices (1,2,3) of the two KV planes or 0 if uninitialized

character(22) :: r_name = "init_tail_weighted_bunch"

tw => beam_init%tw_beam_init

! Checking that |beam_init%dpz_dz| < mode%sigE_E / mode%sig_z
if (abs(beam_init%dPz_dz * beam_init%sig_z) > beam_init%sig_e) then
  call out_io (s_abort$, r_name, &
              "|dpz_dz| MUST be < mode%sigE_E / mode%sig_z")
  call err_exit
endif

! Compute the Twiss parameters beta and alpha, and the emittance for each plane
! 1 = (x,px), 2 = (y,py), 3 = (z,pz)
beta(1) = ele%a%beta
beta(2) = ele%b%beta
alpha(1) = ele%a%alpha
alpha(2) = ele%b%alpha
emitt(1) = beam_init%a_norm_emitt * mass_of(param%particle) / ele%value(E_TOT$)
emitt(2) = beam_init%b_norm_emitt * mass_of(param%particle) / ele%value(E_TOT$)

covar = beam_init%dPz_dz * beam_init%sig_z**2
emitt(3) = sqrt((beam_init%sig_z*beam_init%sig_e)**2 - covar**2)
beta(3) = beam_init%sig_z**2 / emitt(3)
alpha(3) = - covar / emitt(3)


! Fill the corresponding struct and generate the distribution for each phase plane
do i = 1, 3
  call str_upcase (tw%type(i), tw%type(i))
   select case (tw%type(i))
   case ('ELLIPSE')
      call init_ellipse_distribution (tw%n_ellipse(i), tw%part_per_ellipse(i), &
                             tw%sigma_cutoff(i), beta(i), alpha(i), emitt(i), space2D(i))
   case ('GRID')
      call init_grid_distribution (tw%n_x(i), tw%n_px(i), tw%minima(2*i-1), &
                        tw%maxima(2*i-1), tw%minima(2*i), tw%maxima(2*i), space2D(i))
   case ('KV') 
      KV_counter = KV_counter + 1
      ! If we are doing KV, we keep track of which phase planes are KV.
      ! For the first plane, we only store the plane's index and compute the emittance,
      ! but for the second plane, we also call the subroutine to make the distribution.
      if (KV_counter .eq. 1) then
         KV_indices(1) = i
      else if (KV_counter .eq. 2) then
         KV_indices(2) = i
         call init_KV_distribution (tw%n_I2, tw%part_per_ellipse(KV_indices(1)), &
                  tw%part_per_ellipse(i), tw%A, beta(KV_indices(1)), beta(i), &
                  alpha(KV_indices(1)), alpha(i), emitt(KV_indices(1)), emitt(i), space4D)
      endif
   case default
      call out_io (s_abort$, r_name, &
           'PHASE SPACE DISTRIBUTION TYPE NOT RECOGNIZED')
      call err_exit
   end select
enddo

! Check that there are either 0 or 2 planes involved in KV
if (KV_counter .ne. 0 .and. KV_counter .ne. 2) then
   call out_io (s_abort$, r_name, &
        'KV DISTRIBUTION TYPE REQUIRES 2 PHASE PLANES')
   call err_exit
endif

! Combine the planes into one collection of particle positions and weights
if (KV_counter .eq. 2) then
   call multiply_2planes (KV_indices, space4D, space2D(6-KV_indices(1)-KV_indices(2)), space6D)
else
   call multiply_3planes (space2D, space6D)
endif

bunch%charge = beam_init%bunch_charge
bunch%ix_ele = ele%ix_ele

do i = 1, space6D%n_particle
   p => bunch%particle(i)

   p%r%vec = space6D%vec(i,1:6)
   p%charge = space6D%vec(i,7) * bunch%charge
   p%ix_lost = not_lost$
enddo

! recenter the bunch and include beam jitter
call recenter_bunch (beam_init, bunch)

bunch%z_center = 0  ! Default
bunch%t_center = 0  ! Default
bunch%ix_bunch = 1  ! Default

! particle spin

call init_spin_distribution (beam_init, bunch)

end subroutine init_tail_weighted_bunch


!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
!+
! Subroutine init_grid_distribution (n_x, n_px, x_min, x_max, px_min, px_max, space2D)
!
! Subroutine to initialize a uniform rectangular grid as the phase space
! distribution of a bunch.
!
! Input:
!   n_x             -- number of columns
!   n_px            -- number of rows
!   x_min, x_max    -- upper and lower limits in beam size
!   px_min, px_max  -- upper and lower limits in divergence
!
! Output:
!   space2D -- 2D phase_space_struct of particle positions and charges
!
!-

subroutine init_grid_distribution (n_x, n_px, x_min, x_max, px_min, px_max, space2D)

implicit none

integer n_x, n_px
real(rp) x_min, x_max, px_min, px_max

type (phase_space_struct) space2D

integer i, j     ! counters for the x and px directions
integer k        ! particle counter

real(rp) x, px, charge

space2D%n_particle = n_x * n_px       ! total number of particles
charge = 1.0 / space2D%n_particle     ! total charge = 1

if (allocated(space2D%vec)) deallocate(space2D%vec)
allocate (space2D%vec(space2D%n_particle, 3))

k = 1

do i = 1, n_x
   if (n_x == 1) then
      x = x_min
   else
      x = x_min + real(i - 1)/(n_x - 1) * (x_max - x_min)
   endif

   do j = 1, n_px
      if (n_px == 1) then
         px = px_min
      else
         px = px_min + real(j - 1)/(n_px - 1) * (px_max - px_min)
      endif

      space2D%vec(k,1) = x
      space2D%vec(k,2) = px
      space2D%vec(k,3) = charge

      k = k + 1
   enddo
enddo

end subroutine init_grid_distribution


!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
!+
! Subroutine init_ellipse_distribution (n_ellipse, part_per_ellipse, sigma_cutoff,
!                                       beta, alpha, emitt, space2D)
!
! Subroutine to initalize a phase space distribution as a set of concentric
! ellipses of macroparticles representing a Gaussian distribution.
!
! Input:
!   n_ellipse        -- number of ellipses (>= 1)
!   part_per_ellipse  -- number of particles per ellipse
!   sigma_cutoff      -- sigma cutoff of the representation
!   beta, alpha       -- Twiss parameters
!   emitt             -- emittance
!
! Output:
!   space2D -- 2D phase_space_struct of particle positions and charges
!
! See manual for more details.
!-

subroutine init_ellipse_distribution (n_ellipse, part_per_ellipse, sigma_cutoff, beta, alpha, emitt, space2D)

implicit none

integer n_ellipse, part_per_ellipse        ! N and M in the bmad manual
real(rp) sigma_cutoff
real(rp) beta, alpha, emitt

type (phase_space_struct) space2D

integer n, m     ! counters for the J and phi directions
integer k        ! particle counter

real(rp) b_inner, b_outer                  ! B_{n-1}/epsilon  and B_{n}/epsilon in the bmad manual

real(rp) J, phi
real(rp) x, px, charge

space2D%n_particle = n_ellipse * part_per_ellipse

if (allocated(space2D%vec)) deallocate(space2D%vec)
allocate (space2D%vec(space2D%n_particle, 3))

k = 1

do n = 1, n_ellipse
   if (n == n_ellipse) then
      J = emitt * (sigma_cutoff**2/2.0 + 1)                     ! J_N    this is the ellipse that represents
      charge = exp(-sigma_cutoff**2/(2.0 * part_per_ellipse))   ! q_N    the distribution out to infinity
   else
      b_inner = sigma_cutoff**2/2.0 * (real(n - 1)/(n_ellipse - 1))**2
      b_outer = sigma_cutoff**2/2.0 * (real(n)/(n_ellipse - 1))**2

      J = emitt * ((b_inner + 1.0) * exp(-b_inner) - (b_outer + 1.0) * exp(-b_outer)) / (exp(-b_inner) - exp(-b_outer))   ! J_n
      charge = exp(-b_inner) - exp(-b_outer)   ! q_n
   endif

   do m = 1, part_per_ellipse
      phi = (2.0 * pi * m) / part_per_ellipse

      x = sqrt(2.0 * J * beta) * cos(phi)
      px = -sqrt(2.0 * J / beta) * (alpha * cos(phi) + sin(phi))
      
      space2D%vec(k,1) = x
      space2D%vec(k,2) = px
      space2D%vec(k,3) = charge / part_per_ellipse

      k = k + 1
   enddo
enddo

end subroutine init_ellipse_distribution


!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
!+
! Subroutine init_KV_distribution (n_I2, n_phi1, n_phi2, A, beta1, beta2, 
!                                  alpha1, alpha2, e1, e2, space4D)
!
! Subroutine to initalize a phase space distribution as a set of concentric
! ellipses of macroparticles representing a Kapchinsky-Vladimirsky distribution.
!
! Input:
!   n_I2, n_phi1, n_phi2         -- number of points to iterate over in each direction
!   A                            -- = I1/e
!   beta1, beta2, alpha1, alpha2 -- Twiss parameters of each phase plane
!   e1, e2                       -- emittance of each phase plane
!
! Output:
!   space4D -- 4D phase_space_struct of particle positions and charges
!
! See manual for more details.
!
!-

subroutine init_KV_distribution (n_I2, n_phi1, n_phi2, A, beta1, beta2, alpha1, alpha2, e1, e2, space4D)

implicit none

integer n_I2, n_phi1, n_phi2
real(rp) A
real(rp) beta1, beta2, alpha1, alpha2, e1, e2

type (phase_space_struct) space4D

integer i_I2, i_phi1, i_phi2        ! do-loop counters
integer k                           ! particle counter

real(rp) e
real(rp) I1, I2
real(rp) J1, J2, phi1, phi2
real(rp) x1, x2, px1, px2, charge

space4D%n_particle = n_I2 * n_phi1 * n_phi2
charge = 1.0 / space4D%n_particle                       ! constant charge density

if (allocated(space4D%vec)) deallocate(space4D%vec)
allocate (space4D%vec(space4D%n_particle, 5))

e = 1.0 / sqrt(1.0 / e1**2 + 1.0 / e2**2)
I1 = A * e

k = 1

do i_I2 = 1, n_I2
   I2 = -e1/e2 * I1 + e1*e2/e**2 * I1 * real(i_I2 - 0.5)/n_I2

   J1 = (I1/e1 - I2/e2) * e
   J2 = (I1/e2 + I2/e1) * e
   
   do i_phi1 = 1, n_phi1
      phi1 = 2.0 * pi * real(i_phi1 - 1)/n_phi1
 
      do i_phi2 = 1, n_phi2
         phi2 = 2.0 * pi * real(i_phi2 - 1)/n_phi2

         x1 = sqrt(2.0 * J1 * beta1) * cos(phi1)
         px1 = -sqrt(2.0 * J1 / beta1) * (alpha1 * cos(phi1) + sin(phi1))
         x2 = sqrt(2.0 * J2 * beta2) * cos(phi2)
         px2 = -sqrt(2.0 * J2 / beta2) * (alpha2 * cos(phi2) + sin(phi2))
         
         space4D%vec(k,1) = x1
         space4D%vec(k,2) = px1
         space4D%vec(k,3) = x2
         space4D%vec(k,4) = px2
         space4D%vec(k,5) = charge

         k = k + 1
      enddo
   enddo
enddo

end subroutine init_KV_distribution

!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
!+
! Subroutine multiply_3planes (space2D, space6D)
!
! This subroutine multiplies three 2D phase planes by outputting a 6D
! phase_space_struct with all possible combinations of 2D phase space positions.
! Particle weights are conferred accordingly.
!
! Input:
!   space2D(3) -- three 2D phase_space_struct objects representing planes x, y, and z
!
! Output:
!   space6D    -- 6D phase_space_struct representing the full phase space
!
!-

subroutine multiply_3planes (space2D, space6D)

implicit none

type (phase_space_struct) space2D(3), space6D
integer i, j, k
integer :: m = 1          ! particle counter
space6D%n_particle = space2D(1)%n_particle * space2D(2)%n_particle * space2D(3)%n_particle

if (allocated(space6D%vec)) deallocate(space6D%vec)
allocate (space6D%vec(space6D%n_particle, 7))

do i = 1, space2D(1)%n_particle
   do j = 1, space2D(2)%n_particle
      do k = 1, space2D(3)%n_particle
         space6D%vec(m,1) = space2D(1)%vec(i,1)
         space6D%vec(m,2) = space2D(1)%vec(i,2)
         space6D%vec(m,3) = space2D(2)%vec(j,1)
         space6D%vec(m,4) = space2D(2)%vec(j,2)
         space6D%vec(m,5) = space2D(3)%vec(k,1)
         space6D%vec(m,6) = space2D(3)%vec(k,2)
         space6D%vec(m,7) = space2D(1)%vec(i,3) * space2D(2)%vec(j,3) * space2D(3)%vec(k,3)

         m = m + 1
      enddo
   enddo
enddo

end subroutine multiply_3planes


!----------------------------------------------------------
!----------------------------------------------------------
!----------------------------------------------------------
!+
! Subroutine multiply_2planes (index_4D, space4D, space2D, space6D)
!
! This subroutine multiplies a 4D phase plane with a 2D phase plane by 
! outputting a 6D phase_space_struct with all possible combinations of 
! positions between the two planes.  Particle weights are conferred accordingly.
!
! Input:
!   index4D(2) -- tells which indices of the final 6D phase space should be from
!              the 4D phase plane (1 = x-px, 2 = y-py, 3 = z-pz)
!   space4D    -- 4D phase_space_struct
!   space2D    -- 2D phase space struct
!
! Output:
!   space6D    -- 6D phase_space_struct representing the full phase space
!
!-

subroutine multiply_2planes (index_4D, space4D, space2D, space6D)

implicit none

type (phase_space_struct) space4D, space2D, space6D
integer index_4D(2), index_2D

integer i, j, m

m = 1          ! particle counter
index_2D = 6 - index_4D(1) - index_4D(2)   ! 1+2+3=6

space6D%n_particle = space4D%n_particle * space2D%n_particle

if (allocated(space6D%vec)) deallocate(space6D%vec)
allocate (space6D%vec(space6D%n_particle, 7))

do i = 1, space4D%n_particle
   do j = 1, space2D%n_particle
      space6D%vec(m,2*index_4D(1)-1) = space4D%vec(i,1)
      space6D%vec(m,2*index_4D(1)  ) = space4D%vec(i,2)
      space6D%vec(m,2*index_4D(2)-1) = space4D%vec(i,3)
      space6D%vec(m,2*index_4D(2)  ) = space4D%vec(i,4)
      space6D%vec(m,2*index_2D-1   ) = space2D%vec(j,1)
      space6D%vec(m,2*index_2D     ) = space2D%vec(j,2)
      space6D%vec(m,7) = space4D%vec(i,5) * space2D%vec(j,3)
      
      m = m + 1
   enddo
enddo

end subroutine multiply_2planes

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
type (particle_struct), pointer :: p

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
   p%r%vec = p%r%vec + center
enddo

end subroutine recenter_bunch


!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine init_spin_distribution (beam_init, bunch)
!
! Initializes a spin distribution according to beam_init%spin
!
! Input:
!  beam_init -- (beam_init_struct): 
!           %spin  -- (spin_init_struct): spin parameters
!
! Output:
!  bunch          -- (bunch_struct)
!-

subroutine init_spin_distribution (beam_init, bunch)

implicit none

type (beam_init_struct) beam_init
type (bunch_struct) bunch
type (spin_polar_struct) :: polar

real(rp) :: rang, ranl, sigma, vec(3), polarizationvec(3)

integer i

!

polar%xi = 0.0 ! spinor phase is zero

sigma = acos(beam_init%spin%polarization)

if (beam_init%spin%polarization /= 1.0) then
  call out_io (s_error$, "init_spin_distribution", &
                         "Right now, will only set 100% polarization")
endif

! This isn't working correctly yet, so just do %100 polarization for now.
! First set up aroun theta = 0
!   call ran_gauss (rang)
!   call ran_uniform (ranl)
!   polar%theta = sigma * rang
!   polar%phi = 2.0 * pi * ranl

do i = 1, size(bunch%particle)
  polar%theta = beam_init%spin%theta
  polar%phi = beam_init%spin%phi
  call polar_to_spinor (polar, bunch%particle(i)%r)
enddo

end subroutine init_spin_distribution

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! subroutine calc_bunch_params_slice (bunch, ele, param, bunch_params, 
!                           plane, slice_center, slice_spread, err, print_err)
!
! Finds all bunch parameters for a slice through the beam distribution.
!
! Modules needed:
!  use beam_mod
!
! Input:
!   bunch        -- Bunch_struct
!   ele          -- ele_struct: element to find parameters at
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

subroutine calc_bunch_params_slice (bunch, ele, param, bunch_params, &
                          plane, slice_center, slice_spread, err, print_err)

implicit none

type (bunch_struct), intent(in) :: bunch
type (beam_struct) :: beam
type (ele_struct) :: ele
type (lat_param_struct) param
type (bunch_params_struct) bunch_params

real(rp) slice_center, slice_spread

integer plane
integer i, n_part

logical, optional :: print_err
logical err

!

n_part = 0
do i = 1, size(bunch%particle)
  if (bunch%particle(i)%r%vec(plane) .le. slice_center + abs(slice_spread) .and. &
      bunch%particle(i)%r%vec(plane) .ge. slice_center - abs(slice_spread)) &
            n_part = n_part + 1
enddo

call reallocate_beam (beam, 1, n_part)

beam%bunch(1)%charge = bunch%charge
beam%bunch(1)%z_center = bunch%z_center
beam%bunch(1)%t_center = bunch%t_center

n_part = 1
do i = 1, size(bunch%particle)
  if (bunch%particle(i)%r%vec(plane) .le. slice_center + abs(slice_spread) .and. &
      bunch%particle(i)%r%vec(plane) .ge. slice_center - abs(slice_spread)) then
            beam%bunch(1)%particle(n_part) = bunch%particle(i)
            n_part = n_part + 1
  endif
enddo

call calc_bunch_params (beam%bunch(1), ele, param, bunch_params, err, print_err)

end subroutine calc_bunch_params_slice

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine calc_bunch_params (bunch, ele, param, bunch_params, err, print_err)
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
!   ele       -- ele_struct: element to find parameters at
!   param -- Param_struct: lattice parameters.
!     %particle -- Particle being tracked.
!   print_err -- Logical, optional: If present and False then suppress 
!                  "no eigen-system found" messages.
!
! Output     
!   bunch_params -- bunch_params_struct:
!     %a,%b,%z       -- Projected parameters
!       %alpha; %beta; %gamma
!       %eta, %etap, %norm_emitt
!     %a,%b,%c       -- Normal-Mode parameters
!       %alpha; %beta; %gamma
!       %eta, %etap, %norm_emitt
!     %sigma         -- Projected Sigma Matrix
!     %sigma_normal  -- Normal-Mode Sigma Matrix
!     %centroid
!     %spin      
!       %polarization -- Polarization
!       %theta        -- Polar Angle of polarization vector
!       %phi          -- Polar Angle of polarization vector
!     %n_particle ! # particles not lost
!   err   -- Logical: Set True if there is an error in mat_eigen routine.
!-

subroutine calc_bunch_params (bunch, ele, param, bunch_params, err, print_err)

implicit none

type (bunch_struct), intent(in) :: bunch
type (ele_struct) :: ele
type (lat_param_struct) param
type (bunch_params_struct) bunch_params

real(rp) exp_x2, exp_px2, exp_x_px, exp_x_d, exp_px_d
real(rp) avg_energy, temp6(6)
real(rp) eta, etap
real(rp) :: sigma_s(6,6) 
real(rp) :: s(6,6)
real(rp) :: sigma_s_save(6,6) = 0.0
real(rp) :: sigma(6,6) = 0.0
real(rp) :: d_r(6) = 0.0, d_i(6) = 0.0, e_r(6,6) = 0.0, e_i(6,6) = 0.0
real(rp) :: u(6,6)
real(rp) :: n_real(6,6)
real(rp) :: beta_66_iii

complex(rp) :: sigma_s_complex(6,6) = 0.0
complex(rp) :: n(6,6), e(6,6), q(6,6)

integer i, j

logical, optional :: print_err
logical err, err1

character(18) :: r_name = "calc_bunch_params"

!
  
s = 0.0

s(1,2) =  1.0 
s(2,1) = -1.0
s(3,4) =  1.0 
s(4,3) = -1.0
s(5,6) =  1.0 
s(6,5) = -1.0

! n_particle and centroid

bunch_params%n_live_particle = count(bunch%particle%ix_lost == not_lost$)
bunch_params%charge_live = sum(bunch%particle%charge, mask = (bunch%particle%ix_lost == not_lost$))

if (bunch_params%charge_live == 0) then
  bunch_params%centroid%vec = 0.0     ! zero everything
  bunch_params%sigma = 0
  call zero_plane (bunch_params%x)
  call zero_plane (bunch_params%y)
  call zero_plane (bunch_params%z)
  call zero_plane (bunch_params%a)
  call zero_plane (bunch_params%b)
  call zero_plane (bunch_params%c)
  return
endif
  
! average the energy

avg_energy = sum((1+bunch%particle%r%vec(6))*bunch%particle%charge, & 
                              mask = (bunch%particle%ix_lost == not_lost$))
avg_energy = avg_energy * ele%value(E_TOT$) / bunch_params%charge_live

! Convert to geometric coords and find the sigma matrix

call find_bunch_sigma_matrix (bunch%particle, bunch_params%centroid%vec, bunch_params%sigma, sigma_s)

! X, Y, & Z Projected Parameters
call projected_twiss_calc ('X', bunch_params%x, bunch_params%sigma(s11$), bunch_params%sigma(s22$), &
                      bunch_params%sigma(s12$), bunch_params%sigma(s16$), bunch_params%sigma(s26$))

call projected_twiss_calc ('Y', bunch_params%y, bunch_params%sigma(s33$), bunch_params%sigma(s44$), &
                      bunch_params%sigma(s34$), bunch_params%sigma(s36$), bunch_params%sigma(s46$))

call projected_twiss_calc ('Z', bunch_params%z, bunch_params%sigma(s55$), bunch_params%sigma(s66$), &
                      bunch_params%sigma(s56$), bunch_params%sigma(s56$), bunch_params%sigma(s66$))
     
! Normal-Mode Parameters.
! Use Andy Wolski's eigemode method to find normal-mode beam parameters.
! find eigensystem of sigma.S 

sigma_s_save = sigma_s
call mat_eigen (sigma_s, d_r, d_i, e_r, e_i, err, print_err)
if (err) goto 999

! The eigen-values of Sigma.S are the normal-mode emittances (eq. 32)

bunch_params%a%norm_emit = d_i(1) * (avg_energy/mass_of(param%particle))
bunch_params%b%norm_emit = d_i(3) * (avg_energy/mass_of(param%particle))
bunch_params%c%norm_emit = d_i(5) * (avg_energy/mass_of(param%particle))

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

call normalize_e (e)
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

999 continue

if (bmad_com%spin_tracking_on) call calc_spin_params ()
  
! convert back to cannonical coords

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

end subroutine zero_plane
  
!----------------------------------------------------------------------
! contains

subroutine projected_twiss_calc (plane, twiss, exp_x2, exp_px2, exp_x_px, exp_x_d, exp_px_d)

implicit none

type (twiss_struct) :: twiss

real(rp), intent(in) :: exp_x2, exp_px2, exp_x_px, exp_x_d, exp_px_d
real(rp) emitt, x2, x_px, px2

logical err

character(*) plane

!

if (bunch_params%sigma(s66$) /= 0) then
  twiss%eta   = exp_x_d / bunch_params%sigma(s66$)
  twiss%etap  = exp_px_d / bunch_params%sigma(s66$)
endif

x2   = exp_x2   
x_px = exp_x_px 
px2  = exp_px2  

emitt = sqrt(x2*px2 - x_px**2)

twiss%norm_emit = (avg_energy/mass_of(param%particle)) * emitt

if (emitt /= 0) then
  twiss%alpha = -x_px / emitt
  twiss%beta  = x2 / emitt
  twiss%gamma = px2 / emitt
endif

end subroutine projected_twiss_calc

!----------------------------------------------------------------------
! contains
! Eq. 14 But using row vectors to conform to BMAD's mat_eigen

subroutine normalize_e (e)

implicit none

complex(rp) :: s(6,6)

complex(rp) :: e(6,6), temp(6)
complex(rp) :: wronsk, factor

integer i, j, k

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
   factor = sqrt(i_imaginary) / sqrt(wronsk)
   ! this next step is the actual normalization (Eq. 14a)
   e(i+1,:) = e(i+1,:) * factor 
   e(i,:) = conjg(e(i+1,:))
 enddo


end subroutine
  
!----------------------------------------------------------------------
! contains

subroutine calc_spin_params ()

implicit none

type (spin_polar_struct) polar, ave_polar

real(rp) vec(3), ave_vec(3)

! polarization vector

bunch_params%spin%theta = 0.0
bunch_params%spin%phi   = 0.0

ave_vec = 0.0
do i = 1, size(bunch%particle)
  if (bunch%particle(i)%ix_lost /= not_lost$) cycle
  call spinor_to_vec (bunch%particle(i)%r, vec)
  ave_vec = ave_vec + vec * bunch%particle(i)%charge
enddo

ave_vec = ave_vec / bunch_params%charge_live
call vec_to_polar (ave_vec, ave_polar)
bunch_params%spin%theta = ave_polar%theta
bunch_params%spin%phi   = ave_polar%phi

! polarization

bunch_params%spin%polarization = 0.0

  
do i = 1, size(bunch%particle)
  if (bunch%particle(i)%ix_lost /= not_lost$) cycle
  call spinor_to_polar (bunch%particle(i)%r, polar)
  bunch_params%spin%polarization = bunch_params%spin%polarization + &
           cos(angle_between_polars (polar, ave_polar)) * bunch%particle(i)%charge
enddo

bunch_params%spin%polarization = bunch_params%spin%polarization / bunch_params%charge_live
    
end subroutine calc_spin_params

end subroutine calc_bunch_params
  
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine find_bunch_sigma_matrix (particle, ave, sigma, sigma_s)
!
! Routine to find the sigma matrix elements of a particle distribution.
! 
! Modules needed:
!   use beam_mod
!
! Input:
!   particle(:) -- Particle_struct: Array of particles.
!
! Output:
!   sigma(21)    -- Real(rp): Sigma matrix elements.
!   ave(6)       -- Real(rp): Bunch Centroid.
!   sigma_S(6,6) -- Sigma x S matrix for Wolski normal-modes
!-

subroutine find_bunch_sigma_matrix (particle, avg, sigma, sigma_s)

implicit none

type (particle_struct) :: particle(:)

real(rp) charge_live, avg(6), sigma(21)
real(rp) sigma_s(6,6), s(6,6)

integer i

!

charge_live = sum(particle%charge, mask = (particle%ix_lost == not_lost$))

do i = 1, 6
  avg(i) = sum(particle(:)%r%vec(i)*particle(:)%charge, &
                      mask = (particle(:)%ix_lost == not_lost$)) / charge_live
enddo

sigma(s11$) = exp_calc (particle, 1, 1, avg)
sigma(s12$) = exp_calc (particle, 1, 2, avg)
sigma(s13$) = exp_calc (particle, 1, 3, avg)
sigma(s14$) = exp_calc (particle, 1, 4, avg)
sigma(s15$) = exp_calc (particle, 1, 5, avg)
sigma(s16$) = exp_calc (particle, 1, 6, avg)
sigma(s22$) = exp_calc (particle, 2, 2, avg)
sigma(s23$) = exp_calc (particle, 2, 3, avg)
sigma(s24$) = exp_calc (particle, 2, 4, avg)
sigma(s25$) = exp_calc (particle, 2, 5, avg)
sigma(s26$) = exp_calc (particle, 2, 6, avg)
sigma(s33$) = exp_calc (particle, 3, 3, avg)
sigma(s34$) = exp_calc (particle, 3, 4, avg)
sigma(s35$) = exp_calc (particle, 3, 5, avg)
sigma(s36$) = exp_calc (particle, 3, 6, avg)
sigma(s44$) = exp_calc (particle, 4, 4, avg)
sigma(s45$) = exp_calc (particle, 4, 5, avg)
sigma(s46$) = exp_calc (particle, 4, 6, avg)
sigma(s55$) = exp_calc (particle, 5, 5, avg)
sigma(s56$) = exp_calc (particle, 5, 6, avg)
sigma(s66$) = exp_calc (particle, 6, 6, avg)

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

function exp_calc (particle, ix1, ix2, avg) result (this_sigma)

implicit none

type (particle_struct) particle(:)
real(rp) avg(:)
real(rp) this_sigma

integer ix1, ix2

!
                                    
this_sigma = sum((particle(:)%r%vec(ix1) - avg(ix1)) * &
                               (particle(:)%r%vec(ix2) - avg(ix2)) * particle(:)%charge, &
                               mask = (particle%ix_lost == not_lost$))

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
!		bunch1 = bunch2
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
!		beam1 = beam2
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
