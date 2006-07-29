module beam_mod

use wake_mod
use csr_mod
use spin_mod

interface assignment (=)
  module procedure bunch_equal_bunch
  module procedure beam_equal_beam
end interface

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track1_beam 
!
! Subroutine to track a beam of particles through a single element.
!
! This routine is an overloaded name for:
!   track1_beam_lat (beam_start, lat, ix_ele, beam_end)
!   track1_beam_ele (beam_start, ele, param, beam_end)
!
! Note: For the purposes of the wake calculation it is assumed that the
! bunches are ordered with %bunch(1) being the head bunch (largest s).
!
! Note: If longitudinal space charge is on then track1_beam_lat must be used.
!
! Modules needed:
!   use beam_mod
!
! Input:
!   beam_start  -- Beam_struct: Starting beam position.
!   lat         -- Ring_struct: Lattice containing element to be tracked through.
!   ix_ele      -- Integer: Index of element to track through.
!   ele         -- Ele_struct: Element to be tracked through.
!   param       -- Param_struct: General parameters.
!
! Output:
!   beam_end    -- beam_struct: Ending beam position.
!-

interface track1_beam
  module procedure track1_beam_ele
  module procedure track1_beam_lat
end interface


!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track1_bunch
!
! Subroutine to track a bunch of particles through an element.
!
! This routine is an overloaded name for:
!   track1_bunch_lat (beam_start, lat, ix_ele, beam_end)
!   track1_bunch_ele (beam_start, ele, param, beam_end)
!
! Modules needed:
!   use beam_mod
!
! Input:
!   bunch_start -- bunch_struct: Starting bunch position.
!   lat         -- Ring_struct: Lattice containing element to be tracked through.
!   ix_ele      -- Integer: Index of element to track through.
!   ele         -- Ele_struct: Element to be tracked through.
!   param       -- Param_struct: General parameters.
!
! Output:
!   bunch_end -- Bunch_struct: Ending bunch position.
!-

interface track1_bunch
  module procedure track1_bunch_ele
  module procedure track1_bunch_lat
end interface

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track_beam (lat, beam, ix1, ix2)
!
! Subroutine to track a beam of particles from the end of
! lat%ele_(ix1) Through to the end of lat%ele_(ix2).
!
! Modules needed:
!   use beam_mod
!
! Input:
!   lat    -- Ring_struct: Lattice to track through.
!   beam   -- Beam_struct: Beam at end of element ix1.
!   ix1    -- Integer, optional: Index of starting element (this element 
!               is NOT tracked through). Default is 0.
!   ix2    -- Integer, optional: Index of ending element.
!               Default is lat%n_ele_use.
!
! Output:
!   beam   -- beam_struct: Beam at end of element ix2.
!-

subroutine track_beam (lat, beam, ix1, ix2)

implicit none

type (ring_struct) :: lat
type (beam_struct) :: beam

integer, optional, intent(in) :: ix1, ix2
integer i, i1, i2, j

! Init

i1 = 0
if (present(ix1)) i1 = ix1
i2 = lat%n_ele_use
if (present(ix2)) i2 = ix2

! Loop over all elements in the lattice

do i = i1+1, i2
  call track1_beam (beam, lat, i, beam)
enddo

end subroutine track_beam

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track1_beam_lat (beam_start, lat, ix_ele, beam_end)
!
! Subroutine to track a beam of particles through a single element.
!
! Note: This routine is overloaded by the routine track1_beam. See this
! routine for more details.
!
! Modules needed:
!   use beam_mod
!
! Input:
!   beam_start  -- Beam_struct: Starting beam position.
!   lat         -- Ring_struct: Lattice containing element to be tracked through.
!   ix_ele      -- Integer: Index of element to track through.
!
! Output:
!   beam_end    -- beam_struct: Ending beam position.
!-

subroutine track1_beam_lat (beam_start, lat, ix_ele, beam_end)

implicit none

type (beam_struct) beam_start
type (beam_struct) :: beam_end
type (ring_struct) :: lat

integer i, ix_ele, n_mode

! zero the long-range wakes if they exist.

if (associated(lat%ele_(ix_ele)%wake)) then
  lat%ele_(ix_ele)%wake%lr%norm_sin = 0; lat%ele_(ix_ele)%wake%lr%norm_cos = 0
  lat%ele_(ix_ele)%wake%lr%skew_sin = 0; lat%ele_(ix_ele)%wake%lr%skew_cos = 0
endif

! loop over all bunches in a beam

do i = 1, size(beam_start%bunch)
  call track1_bunch (beam_start%bunch(i), lat, ix_ele, beam_end%bunch(i))
enddo

end subroutine track1_beam_lat

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track1_beam_ele (beam_start, ele, param, beam_end)
!
! Subroutine to track a beam of particles through a single element.
!
! Note: This routine is overloaded by the routine track1_beam. See this
! routine for more details.
!
! Modules needed:
!   use beam_mod
!
! Input:
!   beam_start  -- beam_struct: starting beam position
!   ele         -- Ele_struct: The element to track through.
!   param       -- Param_struct: General parameters.
!
! Output:
!   beam_end    -- beam_struct: ending beam position.
!-

subroutine track1_beam_ele (beam_start, ele, param, beam_end)

implicit none

type (beam_struct) beam_start
type (beam_struct), target :: beam_end
type (ele_struct) ele
type (param_struct) param

integer i, n_mode

! zero the long-range wakes if they exist.

if (associated(ele%wake)) then
  ele%wake%lr%norm_sin = 0; ele%wake%lr%norm_cos = 0
  ele%wake%lr%skew_sin = 0; ele%wake%lr%skew_cos = 0
endif

! loop over all bunches in a beam

do i = 1, size(beam_start%bunch)
  call track1_bunch (beam_start%bunch(i), ele, param, beam_end%bunch(i))
enddo

end subroutine track1_beam_ele

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track1_bunch_lat (bunch_start, lat, ix_ele, bunch_end)
!
! Subroutine to track a bunch of particles through an element.
!
! Note: This routine is overloaded by the routine track1_bunch. See this
! routine for more details.
!
! Modules needed:
!   use beam_mod
!
! Input:
!   bunch_start -- bunch_struct: Starting bunch position.
!   lat         -- Ring_struct: Lattice containing element to be tracked through.
!   ix_ele      -- Integer: Index of element to track through.
!
! Output:
!   bunch_end -- Bunch_struct: Ending bunch position.
!-

subroutine track1_bunch_lat (bunch_start, lat, ix_ele, bunch_end)

implicit none

type (bunch_struct) bunch_start, bunch_end
type (ring_struct), target :: lat
type (ele_struct), pointer :: ele
type (ele_struct), save :: rf_ele

real(rp) charge
integer i, j, n, ix_ele

!------------------------------------------------
! space charge tracking will also include wakes if they are on too.

if (bmad_com%coherent_synch_rad_on) then
  call track1_bunch_csr (bunch_start, lat, ix_ele, bunch_end)

else
  call track1_bunch_ele (bunch_start, lat%ele_(ix_ele), lat%param, bunch_end)

endif

end subroutine track1_bunch_lat

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track1_bunch_ele (bunch_start, ele, param, bunch_end)
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
!   param       -- Param_struct: General parameters.
!
! Output:
!   bunch_end -- Bunch_struct: Ending bunch position.
!-

subroutine track1_bunch_ele (bunch_start, ele, param, bunch_end)


implicit none

type (bunch_struct) bunch_start, bunch_end
type (ele_struct) ele
type (ele_struct), save :: rf_ele
type (param_struct) param

real(rp) charge
integer i, j, n

character(20) :: r_name = 'track1_bunch_ele'

! It is not possible to calculate space charge from this

if (bmad_com%coherent_synch_rad_on) then
  call out_io (s_abort$, r_name, 'CANNOT COMPUTE CSR WITHOUT ENTIRE LATTICE!')
  call err_exit
endif

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

! rf_ele is half the cavity

rf_ele = ele
rf_ele%value(l$) = ele%value(l$) / 2
rf_ele%value(beam_energy$) = &
        (ele%value(energy_start$) + ele%value(beam_energy$)) / 2
rf_ele%value(p0c$) = &
            (ele%value(p0c_start$) + ele%value(p0c$)) / 2
rf_ele%value(e_loss$) = 0


! zero all offsets and kicks (offsetting already performed above)
call zero_ele_offsets (rf_ele)
rf_ele%value(hkick$) = 0.0
rf_ele%value(vkick$) = 0.0
if (associated(ele%a)) deallocate(ele%a)
    
! Track half way through. This includes the sr longitudinal wakes 

call order_particles_in_z (bunch_end)
do j = 1, size(bunch_end%particle)
  if (bunch_end%particle(j)%ix_lost .ne. not_lost$) cycle
  call add_sr_long_wake (rf_ele, bunch_end, j-1, bunch_end%particle(j)%ix_z)
  call track1_particle (bunch_end%particle(bunch_end%particle(j)%ix_z), &
                        rf_ele, param, &
                        bunch_end%particle(bunch_end%particle(j)%ix_z))
enddo

bmad_com%grad_loss_sr_wake = 0.0

! Put in the transverse wakefields

rf_ele%value(l$) = ele%value(l$)  ! restore the correct length for the moment
call track1_sr_wake (bunch_end, rf_ele)
call track1_lr_wake (bunch_end, rf_ele)

! Track the last half of the lcavity.  This includes the sr longitudinal wakes 

rf_ele%value(l$)            = ele%value(l$) / 2
rf_ele%value(energy_start$) = rf_ele%value(beam_energy$)
rf_ele%value(p0c_start$)    = rf_ele%value(p0c$)
rf_ele%value(beam_energy$)  = ele%value(beam_energy$)
rf_ele%value(p0c$)          = ele%value(p0c$)

call order_particles_in_z (bunch_end)
do j = 1, size(bunch_end%particle)
  if (bunch_end%particle(j)%ix_lost .ne. not_lost$) cycle
  call add_sr_long_wake (rf_ele, bunch_end, j-1, bunch_end%particle(j)%ix_z)
  call track1_particle (bunch_end%particle(bunch_end%particle(j)%ix_z), &
                        rf_ele, param, &
                        bunch_end%particle(bunch_end%particle(j)%ix_z))
enddo

bmad_com%grad_loss_sr_wake = 0.0

bunch_end%charge = sum (bunch_end%particle(:)%charge, &
                         mask = (bunch_end%particle(:)%ix_lost == not_lost$))

bmad_com%grad_loss_sr_wake = 0.0

! unset the cavity offset
! wakes applied in cononical coords so don't do canonical coord conversion
do i = 1, size(bunch_end%particle)
  call offset_particle (ele, param, bunch_end%particle(i)%r, unset$, &
      set_canonical = .false.)
enddo
      
end subroutine track1_bunch_ele

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
integer n_sr1, n_sr2_long, n_sr2_trans, k_start

!-----------------------------------
! If there is no wake for this element or the sr wakes are turned off then just 
! use the e_loss attribute (as set in track1_bmad).

bmad_com%grad_loss_sr_wake = 0.0

n_sr1 = size(ele%wake%sr1) 
n_sr2_long = size(ele%wake%sr2_long)
n_sr2_trans = size(ele%wake%sr2_trans)

if ((n_sr1 == 0 .and. n_sr2_long == 0 .and. n_sr2_trans == 0) .or. &
                                          .not. bmad_com%sr_wakes_on) then 
  bmad_com%grad_loss_sr_wake = 0.0
  return 
endif

! the self wake only sees the charge of each real particle, not the "macro"
! charge of the simulated particle

if (n_sr1 > 0) then
  bmad_com%grad_loss_sr_wake = bmad_com%grad_loss_sr_wake &
                        + ele%wake%sr1(0)%long * e_charge / 2.0

!-----------------------------------
! add up all wakes from front of bunch to follower
  do i = 1, num_in_front
    if (bunch%particle(bunch%particle(i)%ix_z)%ix_lost .eq. not_lost$) &
      call sr1_add_long_kick (ele, bunch%particle(bunch%particle(i)%ix_z)%r, &
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

real(rp) dz_sr1, sr02, z_sr1_max
integer i, j, k, i1, i2, i_sr2, n_sr1, n_sr2_long, n_sr2_trans, k_start

logical wake_here
character(16) :: r_name = 'track1_sr_wake'

!-----------------------------------

if (.not. associated(ele%wake)) return
  
p => bunch%particle
  
! error check and zero wake sums and order particles in z

call order_particles_in_z (bunch)  
if (size(ele%wake%sr2_long) /= 0) then
  i1 = p(1)%ix_z 
  i2 = p(size(p))%ix_z
  if (p(i1)%r%vec(5) - p(i2)%r%vec(5) > ele%wake%z_sr2_max) then
    call out_io (s_abort$, r_name, &
        'Bunch longer than SR2 wake can handle for element: ' // ele%name)
    call err_exit
  endif
endif

do i = 1, size(ele%wake%sr2_long)
  ele%wake%sr2_long%norm_sin = 0
  ele%wake%sr2_long%norm_cos = 0
  ele%wake%sr2_long%skew_sin = 0
  ele%wake%sr2_long%skew_cos = 0
enddo

!

n_sr1 = size(ele%wake%sr1) 
z_sr1_max = 0
if (n_sr1 > 0) then
  z_sr1_max = ele%wake%sr1(n_sr1-1)%z
  dz_sr1 = z_sr1_max / (n_sr1 - 1)
endif

! loop over all particles in the bunch and apply the wake

i_sr2 = 1  ! index of next particle to be added to the sr2 wake sums.

do j = 1, size(p)
  particle => p(p(j)%ix_z)
  ! apply longitudinal self wake

  if (z_sr1_max .ge. 0) then
    call sr2_long_self_wake_apply_kick (ele, particle%charge, particle%r)
  endif

  ! Particle_j is kicked by particles k = 1, ..., j-1.
  ! The particles 1, ... i_sr2-1 have already had their wakes added to the 
  ! sr2 wake sums so the loop is from i_sr2, ..., j-1.

  k_start = i_sr2
  do k = k_start, j-1
    leader => p(p(k)%ix_z)
    if ((particle%r%vec(5) - leader%r%vec(5)) > z_sr1_max) then
      ! use sr1 table to add to particle j the wake of particle k
      call sr1_apply_trans_kick (ele, leader%r, leader%charge, particle%r)
    else
      ! add contribution of particle(k) to wake sums
      i_sr2 = k  ! update i_sr2
      call sr2_long_wake_add_to (ele, leader%r, leader%charge)
      call sr2_trans_wake_add_to(ele, leader%r, leader%charge)
    endif
  enddo

  ! apply wake to particle(j)
  call sr2_long_wake_apply_kick (ele, particle%r)
  call sr2_trans_wake_apply_kick(ele, particle%r)

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
!     ele%wake%lr%norm_sin
!     ele%wake%lr%norm_cos
!     ele%wake%lr%skew_sin
!     ele%wake%lr%skew_cos
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
!   ele   -- Ele_struct: Element with wakefields.
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
  call lr_wake_apply_kick (ele, bunch%s_center, particle%r)
enddo

! Add the wakes left by this bunch to the existing wakes.

do k = 1, size(bunch%particle)
  j = bunch%particle(k)%ix_z
  particle => bunch%particle(j)
  if (particle%ix_lost /= not_lost$) cycle
  call lr_wake_add_to (ele, bunch%s_center, particle%r, particle%charge)
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

if (associated(beam%bunch)) then
  if (n_bunch .eq. 0) then
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

if (n_bunch .eq. 0) return
  
! Allocate

if (.not. associated (beam%bunch)) allocate (beam%bunch(n_bunch))
do i = 1, n_bunch
  if (.not. associated (beam%bunch(i)%particle)) &
                    allocate (beam%bunch(i)%particle(n_particle))
enddo

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine init_beam_distribution (ele, beam_init, beam)
!
! Subroutine to initialize a distribution of particles matched to
! the Twiss parameters, centroid position, and Energy - z correlation
! as specified. Coupling in the element ele is incorporated into the
! distribution
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
!   beam_init   -- beam_init_struct
!     %renorm_center -- Logical: If True then distribution is rescaled to
!                    the desired centroid (to take care of
!                    possible statistical errors in distribution).
!                    (default is True)
!     %renorm_sigma  -- Logical: If True then rescale the distribution to the desired sigmas.
!                    (default is True)
!     %preserve_dist -- Logical: If True then if first call to routine the
!                    distribution is saved. Each next call to routine will just 
!                    load this same distribution
!                    (default is False)
!     %init_spin     -- Logical: If True then the particle spinors will be
!                    initialized with the parameters in beam_init%spin
!                    (default is False)
!
! Output:
!   beam        -- beam_struct
!
!-
 
subroutine init_beam_distribution (ele, beam_init, beam)
 
use random_mod

implicit none

type (ele_struct) ele
type (beam_init_struct) beam_init
type (beam_struct), target :: beam
type (bunch_struct), pointer :: bunch
type (bunch_params_struct) :: params
type (particle_struct), pointer :: p
  
real(rp) dpz_dz, denom
real(rp) a_emitt, b_emitt
real(rp) ave(6), sigma(6), a, b
real(rp) r(6), v_mat(4,4), v_inv(4,4)
real(rp) y, alpha(6), sig_mat(6,6)
real(rp) center(6) ! includes jitter
real(rp) ran(6)
real(rp), save, allocatable :: ran_save(:,:)

integer i, j, j2, n
  
character(22) :: r_name = "init_beam_distribution"

logical, save :: init = .true.

! Generate a set of random numbers.

call reallocate_beam (beam, beam_init%n_bunch, beam_init%n_particle)
bunch => beam%bunch(1)
 
! if preserve_dist and first call to routine then make a distribution

if (beam_init%preserve_dist) then
  if (init) then
    if (allocated(ran_save)) deallocate(ran_save)
    allocate (ran_save(beam_init%n_particle,6))
  endif
endif

  
sig_mat = 0
ave = 0
do n = 1, beam_init%n_particle
  p => bunch%particle(n)
  if (beam_init%preserve_dist .and. init) then
    call ran_gauss(r)
    ran_save(n,:) = r
  elseif (beam_init%preserve_dist .and. .not. init) then
    r = ran_save(n,:)
  else
    call ran_gauss(r)
  endif
  p%r%vec = r
  ave = ave + r
  forall (i=1:6, j=1:6) sig_mat(i,j) = sig_mat(i,j) + r(i) * r(j)
enddo  

if (beam_init%preserve_dist .and. init) init = .false.
      
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

if (beam_init%renorm_sigma) then

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

! Put in beam jitter, include alpha correlations
call ran_gauss(ran)
center(1) = beam_init%center(1) + beam_init%center_jitter(1)*ran(1)
center(2) = beam_init%center(2) + beam_init%center_jitter(2)*ran(2) + &
                 (ele%x%alpha/ele%x%beta) * beam_init%center_jitter(1)*ran(1)
center(3) = beam_init%center(3) + beam_init%center_jitter(3)*ran(3)
center(4) = beam_init%center(4) + beam_init%center_jitter(4)*ran(4) + &
                 (ele%y%alpha/ele%y%beta) * beam_init%center_jitter(3)*ran(3)
center(5) = beam_init%center(5) + beam_init%center_jitter(5)*ran(5)
center(6) = beam_init%center(6) + beam_init%center_jitter(6)*ran(6) + &
                 beam_init%dpz_dz * beam_init%center_jitter(5)*ran(5)
  
! Now scale by the emittances, etc. and put in jitter

call ran_gauss(ran(1:4)) ! ran(3:4) for z and e jitter used below
denom = (1 + center(6)) * ele%value(beam_energy$)
a_emitt = beam_init%a_norm_emitt*(1+beam_init%emitt_jitter(1)*ran(1)) &
                                                      * m_electron / denom
b_emitt = beam_init%b_norm_emitt*(1+beam_init%emitt_jitter(2)*ran(2)) &
                                                      * m_electron / denom
  
dpz_dz = beam_init%dpz_dz
  
call make_v_mats(ele, v_mat, v_inv)

sigma(1) = sqrt(a_emitt * ele%x%beta)
sigma(2) = sqrt(a_emitt / ele%x%beta)
sigma(3) = sqrt(b_emitt * ele%y%beta)
sigma(4) = sqrt(b_emitt / ele%y%beta)
sigma(5) = beam_init%sig_z * (1 + beam_init%sig_z_jitter*ran(3))
sigma(6) = beam_init%sig_e * (1 + beam_init%sig_e_jitter*ran(4))

a = dpz_dz * sigma(5) / sigma(6)

if (abs(a) > 1)  then
  call out_io (s_abort$, r_name, &
              "|dpz_dz| MUST be < mode%sigE_E / mode%sig_z")
  call err_exit
endif
     
b = sqrt(1-a**2)
     
!

do i = 1, beam_init%n_particle

  p => bunch%particle(i)
  r = p%r%vec

  p%r%vec(1) = sigma(1) *  r(1)
  p%r%vec(2) = - sigma(2) * (r(2) + r(1) * ele%x%alpha)
  p%r%vec(3) = sigma(3) *  r(3)
  p%r%vec(4) = - sigma(4) * (r(4) + r(3) * ele%y%alpha)
  p%r%vec(5) = sigma(5) *  r(5)
  p%r%vec(6) = sigma(6) * (r(6) * b + r(5) * a)
      
  ! Include Dispersion
  p%r%vec(1:4) =  p%r%vec(1:4) + &
              p%r%vec(6) * (/ ele%x%eta, ele%x%etap, ele%y%eta, ele%y%etap /)
      
  ! Include Coupling
  p%r%vec(1:4) = matmul(v_mat, p%r%vec(1:4))

  p%r%vec = p%r%vec + center
      
end do
     
! set particle charge

bunch%particle(:)%charge = beam_init%bunch_charge / beam_init%n_particle
bunch%particle(:)%ix_lost = not_lost$
    
! particle spin
call init_spin_distribution (beam_init, bunch)
  
! init all bunches
  
bunch%s_center = 0.0

do i = 2, size(beam%bunch)
  call bunch_equal_bunch (beam%bunch(i), beam%bunch(1))
  beam%bunch(i)%s_center = (1-i) * beam_init%ds_bunch
enddo
  
end subroutine init_beam_distribution

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine init_spin_distribution (beam_init, bunch)
!
! Initializes a spin distribution according to init_beam%spin
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

if (beam_init%spin%polarization .ne. 1.0) then
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
! Subroutine calc_bunch_params (bunch, ele, params)
!
! Finds all bunch parameters defined in bunch_params_struct, both normal-mode
! and projected
!
! Modules needed:
!  use beam_mod
!
! Input:
!   bunch     -- Bunch_struct
!   ele       -- ele_struct: element to find parameters at
!
! Output     
!   params -- bunch_params_struct:
!     %x%alpha; %x%beta; %x%gamma
!     %x%sigma
!     %x%emitt; %x%dpx_dx
!     %y%alpha; %y%beta; %y%gamma
!     %y%sigma
!     %y%emitt; %y%dpx_dx
!     %z%alpha; %z%beta; %z%gamma
!     %z%sigma
!     %z%emitt; %z%dpx_dx
!     %a%alpha; %a%beta; %a%gamma
!     %a%sigma
!     %a%emitt; %a%dpx_dx
!     %b%alpha; %b%beta; %b%gamma
!     %b%sigma
!     %b%emitt; %b%dpx_dx
!     %centroid
!     %n_particle ! # particle not lost
!     %spin
!-

subroutine calc_bunch_params (bunch, ele, params)

implicit none

type (bunch_struct), intent(in) :: bunch
type (ele_struct), intent(in) :: ele
type (bunch_params_struct) params
type (particle_struct), allocatable, save :: a_mode(:)

real(rp) exp_x2, exp_p_x2, exp_x_p_x, exp_x_d, exp_px_d
real(rp) avg_energy, temp6(6)
real(rp) v_mat(4,4), v_inv_mat(4,4)
real(rp) eta, etap

integer i
  
! n_particle and centroid

params%n_particle = count(bunch%particle%ix_lost == not_lost$)
if (params%n_particle == 0) then
  params%centroid%vec = 0.0     ! zero everything
  call zero_plane (params%x)
  call zero_plane (params%y)
  call zero_plane (params%z)
  call zero_plane (params%a)
  call zero_plane (params%b)
endif
  
! average the energy

avg_energy = sum((1+bunch%particle%r%vec(6)), & 
                              mask = (bunch%particle%ix_lost == not_lost$))
avg_energy = avg_energy * ele%value(beam_energy$) / params%n_particle

! Convert to geometric coords and find the sigma matrix

call switch_to_geometric (bunch%particle, set$)
    
call find_bunch_sigma_matrix (bunch%particle, params%centroid%vec, params%sigma)

!--------
! Projected Parameters
! X
call param_stuffit (params%x, params%sigma(s11$), params%sigma(s22$), &
                      params%sigma(s12$), params%sigma(s16$), params%sigma(s26$))
! Y
call param_stuffit (params%y, params%sigma(s33$), params%sigma(s44$), &
                      params%sigma(s34$), params%sigma(s36$), params%sigma(s46$))
! Z
call param_stuffit (params%z, params%sigma(s55$), params%sigma(s66$), &
                      params%sigma(s56$), params%sigma(s56$), params%sigma(s66$))
     
!***
! Normal-Mode Parameters
  
! take out coupling
call make_v_mats (ele, v_mat, v_inv_mat)

if (.not. allocated(a_mode)) allocate (a_mode(size(bunch%particle)))
if (size(a_mode) .ne. size(bunch%particle)) then
  deallocate(a_mode)
  allocate(a_mode(size(bunch%particle)))
endif
do i = 1, size(a_mode)
  a_mode(i)%r%vec(1:4) = matmul(v_inv_mat, bunch%particle(i)%r%vec(1:4))
  ! longitudinal plane doesn't change
  a_mode(i)%r%vec(5:6) = bunch%particle(i)%r%vec(5:6)
enddo 
  
call find_bunch_sigma_matrix (a_mode, temp6, params%sigma_normal)
  
! Take out linear dispersion for normal mode sigma
! x
eta  = params%sigma_normal(s16$) / params%sigma_normal(s66$)
etap = params%sigma_normal(s26$) / params%sigma_normal(s66$)
params%sigma_normal(s11$) = params%sigma_normal(s11$) &
                             - 2*eta*params%sigma_normal(s16$) &
                             + (eta**2)*params%sigma_normal(s66$)
params%sigma_normal(s22$) = params%sigma_normal(s22$) &
                             - 2*etap*params%sigma_normal(s26$) &
                             + (etap**2)*params%sigma_normal(s66$)
params%sigma_normal(s12$) = params%sigma_normal(s12$) &
                             - etap*params%sigma_normal(s16$) &
                             - eta*params%sigma_normal(s26$) &
                             + eta*etap*params%sigma_normal(s66$)
! y
eta  = params%sigma_normal(s36$) / params%sigma_normal(s66$)
etap = params%sigma_normal(s46$) / params%sigma_normal(s66$)
params%sigma_normal(s33$) = params%sigma_normal(s33$) &
                             - 2*eta*params%sigma_normal(s36$) &
                             + (eta**2)*params%sigma_normal(s66$)
params%sigma_normal(s44$) = params%sigma_normal(s44$) &
                             - 2*etap*params%sigma_normal(s46$) &
                             + (etap**2)*params%sigma_normal(s66$)
params%sigma_normal(s34$) = params%sigma_normal(s34$) &
                             - etap*params%sigma_normal(s36$) &
                             - eta*params%sigma_normal(s46$) &
                             + eta*params%y%etap*params%sigma_normal(s66$)
                             
! A
call param_stuffit (params%a, params%sigma_normal(s11$), params%sigma_normal(s22$), &
                      params%sigma_normal(s12$), params%sigma_normal(s16$), &
                      params%sigma_normal(s26$))
! B
call param_stuffit (params%b, params%sigma_normal(s33$), params%sigma_normal(s44$), &
                      params%sigma_normal(s34$), params%sigma_normal(s36$), &
                      params%sigma_normal(s46$))

if (bmad_com%spin_tracking_on) call calc_spin_params ()
  
! convert back to cannonical coords
call switch_to_geometric (bunch%particle, unset$)
  
contains
!----------------------------------------------------------------------
subroutine zero_plane (param)

implicit none

type (bunch_param_struct), intent(out) :: param

param%beta       = 0.0
param%alpha      = 0.0
param%gamma      = 0.0
param%eta        = 0.0
param%etap       = 0.0
param%norm_emitt = 0.0

end subroutine zero_plane
  
!----------------------------------------------------------------------
! contains
! This assumes no vector potential or curvature function
subroutine switch_to_geometric (particle, set)

implicit none

type (particle_struct) particle(:)

integer i

logical set

if (set) then
    particle(:)%r%vec(2) = particle(:)%r%vec(2) / (1 + particle(:)%r%vec(6))
    particle(:)%r%vec(4) = particle(:)%r%vec(4) / (1 + particle(:)%r%vec(6))
else
    particle(:)%r%vec(2) = particle(:)%r%vec(2) * (1 + particle(:)%r%vec(6))
    particle(:)%r%vec(4) = particle(:)%r%vec(4) * (1 + particle(:)%r%vec(6))
endif

end subroutine switch_to_geometric  

!----------------------------------------------------------------------
! contains

subroutine param_stuffit (param, exp_x2, exp_p_x2, exp_x_p_x, exp_x_d, exp_px_d)

implicit none

type (bunch_param_struct), intent(out) :: param
real(rp), intent(in) :: exp_x2, exp_p_x2, exp_x_p_x, exp_x_d, exp_px_d
real(rp) emitt

emitt = SQRT(exp_x2*exp_p_x2 - exp_x_p_x**2)

param%alpha = -exp_x_p_x / emitt
param%beta  = exp_x2 / emitt
param%gamma = exp_p_x2 / emitt
  
param%eta   = exp_x_d / params%sigma(s66$)
param%etap  = exp_px_d / params%sigma(s66$)

param%norm_emitt = (avg_energy/m_electron) * emitt

end subroutine param_stuffit

!----------------------------------------------------------------------
! contains

subroutine calc_spin_params ()

implicit none

type (spin_polar_struct) polar, ave_polar

real(rp) angle

! polarization vector

params%spin%theta = 0.0
params%spin%phi   = 0.0

do i = 1, size(bunch%particle)
  if (bunch%particle(i)%ix_lost .ne. not_lost$) cycle
  call spinor_to_polar (bunch%particle(i)%r, polar)
  params%spin%theta = params%spin%theta + polar%theta 
  params%spin%phi = params%spin%phi + polar%phi
enddo

params%spin%theta = params%spin%theta / params%n_particle
params%spin%phi = params%spin%phi / params%n_particle

ave_polar%xi = 0.0
ave_polar%theta = params%spin%theta
ave_polar%phi = params%spin%phi
  
! polarization

params%spin%polarization = 0.0

  
do i = 1, size(bunch%particle)
  if (bunch%particle(i)%ix_lost .ne. not_lost$) cycle
  call spinor_to_polar (bunch%particle(i)%r, polar)
  params%spin%polarization = params%spin%polarization + &
               cos(angle_between_polars (polar, ave_polar))
enddo

params%spin%polarization = params%spin%polarization / params%n_particle
    
end subroutine calc_spin_params

end subroutine calc_bunch_params
  
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine find_bunch_sigma_matrix (particle, ave, sigma)
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
!   sigma(21)   -- Real(rp): Sigma matrix elements.
!   ave(6)      -- Real(rp): Bunch Centroid.
!-

subroutine find_bunch_sigma_matrix (particle, avg, sigma)

implicit none

type (particle_struct) :: particle(:)
real(rp) sigma(21)
real(rp) avg(6)

integer i, n

!

n = count(particle(:)%ix_lost == not_lost$)

do i = 1, 6
  avg(i) = sum(particle(:)%r%vec(i), mask = (particle(:)%ix_lost == not_lost$)) / n
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
                               (particle(:)%r%vec(ix2) - avg(ix2)), &
                               mask = (particle%ix_lost .eq. not_lost$))

this_sigma = this_sigma / n

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
!  bunch2 -- bunch_struct: Input bunch
!
! Output
!  bunch1 -- bunch_struct: Output bunch
!
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
bunch1%s_center  = bunch2%s_center

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

! The following rule must be observed: If beam%bunch is associated then
! beam%bunch%particle must be also.

n_bun = size(beam2%bunch)

allocate_this = .true.
if (associated(beam1%bunch)) then
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
