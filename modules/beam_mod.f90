module beam_mod

use bmad_struct
use bmad_interface
use wake_mod

integer, parameter :: not_lost$ = -1

type particle_struct
  type (coord_struct) r   ! Center of the particle
  real(rp) charge         ! charge in a particle (Coul).
  integer :: ix_lost = not_lost$  ! Has the particle been lost in tracking?
end type

type bunch_struct
  type (particle_struct), pointer :: particle(:) => null()
  real(rp) charge   ! total charge in a bunch (Coul).
  real(rp) s_center ! longitudinal center of bunch (m).
end type

type beam_struct
  type (bunch_struct), pointer :: bunch(:) => null()
end type

type beam_init_struct
  real(rp) a_norm_emitt ! a-mode emittance
  real(rp) b_norm_emitt ! b-mode emittance
  real(rp) :: dPz_dz = 0  ! Correlation of Pz with long position.
  real(rp) :: center(6) = 0 ! Bench center offset relative to reference.
  real(rp) ds_bunch    ! Distance between bunches.
  real(rp) sig_z       ! Z sigma in m.
  real(rp) sig_e       ! e_sigma in dE/E.
  real(rp) sig_e_cut   ! Energy cut in sigmas.
  real(rp) sig_z_cut   ! Z cut in sigmas.
  real(rp) sig_trans_cut ! Transverse sigma cutoff
  real(rp) bunch_charge ! charge in a bunch.
  integer n_particle   ! Number of simulated particles per bunch.
  integer n_bunch      ! Number of bunches.
end type

type bunch_param_struct
  real(rp) beta, alpha, gamma
  real(rp) eta, etap
  real(rp) sigma, p_sigma
  real(rp) dpx_dx ! x x' correlation
  real(rp) norm_emitt ! normalized emittance
end type

type bunch_params_struct
  type (bunch_param_struct) :: x, y, z, a, b
  type (coord_struct) :: centroid  ! Lab frame
  integer n_particle               ! all non-lost particles
end type

interface assignment (=)
  module procedure bunch_equal_bunch
  module procedure beam_equal_beam
end interface

real(rp) x_max ! used by renorm_probability_funct

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track_beam (ring, beam, ix1, ix2)
!
! Subroutine to track a beam of particles from the end of
! ring%ele_(ix1) Through to the end of ring%ele_(ix2).
!
! Modules needed:
!   use beam_mod
!
! Input:
!   ring   -- Ring_struct: Lattice to track through.
!   beam   -- Beam_struct: Beam at end of element ix1.
!   ix1    -- Integer, optional: Index of starting element (this element 
!               is NOT tracked through). Default is 0.
!   ix2    -- Integer, optional: Index of ending element.
!               Default is ring%n_ele_use.
!
! Output:
!   beam   -- beam_struct: Beam at end of element ix2.
!-

subroutine track_beam (ring, beam, ix1, ix2)

  implicit none

  type (ring_struct) :: ring
  type (beam_struct) :: beam

  integer, optional, intent(in) :: ix1, ix2
  integer i, i1, i2, j

! Init

  i1 = 0
  if (present(ix1)) i1 = ix1
  i2 = ring%n_ele_use
  if (present(ix2)) i2 = ix2

! Loop over all elements in the lattice

  do i = i1+1, i2
    call track1_beam (beam, ring%ele_(i), ring%param, beam)
  enddo

end subroutine track_beam

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track1_beam (beam_start, ele, param, beam_end)
!
! Subroutine to track a beam of particles through a single element.
!
! Note: For the purposes of the wake calculation it is assumed that the
! bunches are ordered with %bunch(1) being the head bunch (largest s).
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

subroutine track1_beam (beam_start, ele, param, beam_end)

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
    ele%wake%lr%s_ref = 0
  endif

! loop over all bunches in a beam

  do i = 1, size(beam_start%bunch)
    call track1_bunch (beam_start%bunch(i), ele, param, beam_end%bunch(i))
  enddo

end subroutine track1_beam

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track1_bunch (bunch_start, ele, param, bunch_end)
!
! Subroutine to track a bunch of particles through an element.
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

Subroutine track1_bunch (bunch_start, ele, param, bunch_end)

  implicit none

  type (bunch_struct) bunch_start, bunch_end
  type (ele_struct) ele
  type (ele_struct), save :: rf_ele
  type (param_struct) param

  real(rp) charge
  integer i, j

! Charge and center

  bunch_end%s_center = bunch_start%s_center
  bunch_end%charge   = bunch_start%charge

  call order_particles_in_z (bunch_start)  ! needed for wakefield calc.

!------------------------------------------------
! Without wakefields just track through

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
! Put the wakefield kicks at the half way point.

! rf_ele is half the cavity

  rf_ele = ele
  rf_ele%value(l$) = ele%value(l$) / 2
  rf_ele%value(beam_energy$) = &
            (ele%value(energy_start$) + ele%value(beam_energy$)) / 2
  rf_ele%value(p0c$) = &
            (ele%value(p0c$) + ele%value(p0c$)) / 2
  rf_ele%value(e_loss$) = 0

! Track half way through. 
! This includes the short-range longitudinal wakefields

  do j = 1, size(bunch_start%particle)
    call track1_particle (bunch_start%particle(j), &
                                    rf_ele, param, bunch_end%particle(j))
  enddo

! Put in the short-range transverse wakefields

  rf_ele%value(l$) = ele%value(l$)  ! restore the correct length for the moment
  call track1_sr_wake (bunch_end, rf_ele)
  call track1_lr_wake (bunch_end, rf_ele)

! Track the last half of the lcavity. 
! This includes the short-range longitudinal wakes.

  rf_ele%value(l$)            = ele%value(l$) / 2
  rf_ele%value(energy_start$) = rf_ele%value(beam_energy$)
  rf_ele%value(beam_energy$)  = ele%value(beam_energy$)
  rf_ele%value(p0c$)          = ele%value(p0c$)

  do j = 1, size(bunch_start%particle)
    call track1_particle (bunch_end%particle(j), &
                                    rf_ele, param, bunch_end%particle(j))
  enddo

  bunch_end%charge = sum (bunch_end%particle(:)%charge, &
                         mask = (bunch_end%particle(:)%ix_lost == not_lost$))

end subroutine track1_bunch

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track1_sr_wake (bunch, ele)
!
! Subroutine to calculate the longitudinal component of the
! short range wake fields. 
! This routine is not really meant for general use.
!-

subroutine track1_sr_wake (bunch, ele)

  implicit none

  type (bunch_struct), target :: bunch
  type (ele_struct) ele
  type (particle_struct), pointer :: particle, leader

  real(rp) dz_sr1, sr02, z_cut
  integer i, j, k, i_sr2, n_sr1, n_sr2_long, n_sr2_trans, k_start

  logical wake_here

!-----------------------------------
! If there is no wake for this element then just use the e_loss attribute.

  n_sr1 = size(ele%wake%sr1) 
  n_sr2_long = size(ele%wake%sr2_long)
  n_sr2_trans = size(ele%wake%sr2_trans)

  if (n_sr1 == 0 .and. n_sr2_long == 0 .and. n_sr2_trans == 0) then 
    bunch%particle(:)%r%vec(6) = bunch%particle(:)%r%vec(6) - &
                       ele%value(e_loss$) * bunch%charge / ele%value(p0c$) 
    return 
  endif

!-----------------------------------
! zero wake sums

  do i = 1, size(ele%wake%sr2_long)
    ele%wake%sr2_long%norm_sin = 0
    ele%wake%sr2_long%norm_cos = 0
    ele%wake%sr2_long%skew_sin = 0
    ele%wake%sr2_long%skew_cos = 0
  enddo

!

  if (n_sr1 > 0) then
    dz_sr1 = ele%wake%sr1(n_sr1-1)%z / (n_sr1 - 1)
    sr02 = ele%wake%sr1(0)%long / 2
  endif

! loop over all particles in the bunch and apply the wake

  i_sr2 = 1  ! index of next particle to be added to the sr2 wake sums.
  z_cut = ele%wake%z_cut_sr ! remember z_cut is negative

  do j = 1, size(bunch%particle)
    particle => bunch%particle(j)
    ! apply longitudinal self wake

    if (z_cut < 0) then
      particle%r%vec(6) = particle%r%vec(6) - sr02 * particle%charge 
    else
      call sr2_long_self_wake_apply_kick (ele, particle%charge, particle%r)
    endif

    ! Particle_j is kicked by particles k = 1, ..., j-1.
    ! The particles 1, ... i_sr2-1 have already had their wakes added to the 
    ! sr2 wake sums so the loop is from i_sr2, ..., j-1.

    k_start = i_sr2
    do k = k_start, j-1
      leader => bunch%particle(k)
      if (leader%r%vec(5) - particle%r%vec(5) < z_cut) then
        ! use sr1 table to add to particle j the wake of particle k
        call sr1_apply_kick (ele, leader%r, leader%charge, particle%r)
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
! This routine is not really meant  for general use.
!-

subroutine track1_lr_wake (bunch, ele)

  implicit none

  type (bunch_struct), target :: bunch
  type (ele_struct) ele
  type (particle_struct), pointer :: particle

  integer n_mode, j, k

! Check to see if we need to do any calc

  if (.not. associated(ele%wake)) return
  n_mode = size(ele%wake%lr)
  if (n_mode == 0 .or. .not. bmad_com%lr_wakes_on) return  

! Give the particles a kick

  do k = 1, size(bunch%particle)
    particle => bunch%particle(k)
    if (particle%ix_lost /= not_lost$) cycle
    call lr_wake_apply_kick (ele, bunch%s_center, particle%r)
  enddo

! Add the wakes left by this bunch to the existing wakes.

  do k = 1, size(bunch%particle)
    particle => bunch%particle(k)
    if (particle%ix_lost /= not_lost$) cycle
    call lr_wake_add_to (ele, bunch%s_center, particle%r, particle%charge)
  enddo


end subroutine

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
!   param  -- Param_struct: Global parameters.
!
! Output:
!   end    -- struct: Ending coords.
!-

subroutine track1_particle (start, ele, param, end)

  implicit none

  type (particle_struct) :: start
  type (particle_struct) :: end
  type (ele_struct) :: ele
  type (param_struct), intent(inout) :: param

! transfer z-order index, charge, etc

  end = start
  if (start%ix_lost /= not_lost$) return
  if (ele%key == marker$) return

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
!                   Order is from large z (head of bunch) to small z.
!                   That is: %particle(1) is the particle at the bunch head. 
!-

Subroutine order_particles_in_z (bunch)

  implicit none

  type (bunch_struct), target :: bunch
  type (particle_struct), pointer :: particle(:)
  type (particle_struct) temp
  integer i, k, nm
  real(rp) z1, z2
  logical ordered

! Order is from large z (head of bunch) to small z.

  particle => bunch%particle
  nm = size(particle)

  do
    ordered = .true.
    do i = 1, nm-1
      if (particle(i)%r%vec(5) < particle(i+1)%r%vec(5)) then
        particle(i:i+1) = particle(i+1:i:-1)
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

subroutine particle_to_angle_coords (particle, energy0)

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
    if (size(beam%bunch) /= n_bunch) then
      de_bunch = .true.
      de_particle = .true.
    endif
    if (size(beam%bunch(1)%particle) /= n_particle) then
      de_particle= .true.
    endif

    do i = 1, size(beam%bunch)
      if (de_particle) deallocate (beam%bunch(i)%particle)
    enddo
    if (de_bunch) deallocate (beam%bunch)

  endif

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
! Subroutine init_beam_distribution (ele, beam_init, beam, renormalize, random_dist)
!
! Subroutine to initialize a distribution of particles matched to
! the Twiss parameters, centroid position, and Energy - z correlation
! as specified. Coupling in the element ele is incorporated into the
! distribution
!
! Make sure: dpz_dz < mode%sigE_E / mode%sig_z
!
! Modules needed:
!   use random_mod
!   use bmad
!
! Input:
!   ele         -- Ele_struct: element to initialize distribution at
!   beam_init   -- beam_init_struct
!   beam        -- beam_struct
!   renormalize -- Logical: If True then distribution is rescaled to
!                   the desired centroid (to take care of
!                   possible statistical errors in distribution).
!   rendom_dist -- Logical: If True then a random gaussian distribution will be
!                  used. Otherwise the inverse error function will be used to get 
!                  an ordered distribution
!
!
! Output:
!   beam        -- beam_struct
!
!-
 
subroutine init_beam_distribution (ele, beam_init, beam, renormalize, random_dist)
 
  use random_mod
  use bmad
  
  implicit none

  type (ele_struct) ele
  type (beam_init_struct) beam_init
  type (beam_struct), target :: beam
  type (bunch_struct), pointer :: bunch
  type (bunch_params_struct) :: params
  
  real(rp) dpz_dz, denom
  real(rp) a_emitt, b_emitt
  real(rp) ave(6), sigma(6), sigma_cut(6), dpz_dz2, a, b
  real(rp) center(6), r(6), v_mat(4,4), v_inv(4,4)
  real(rp) y, phi(6)
  real(rp) :: tol = 1.0e-6_rp
  real(rp) :: emitt_tol = 0.01 ! 1% of design emittance
  real(rp) :: x_domain = 20.0 ! something really big
  real(rp) a_factor, b_factor

  integer i,j,i_p
  integer n_ellipse, n_ellipse_parts, i_part, n_more, n_too_many, n_factor
  
  logical renormalize, random_dist

  character(22) :: r_name = "init_beam_distribution"

!

  call reallocate_beam (beam, beam_init%n_bunch, beam_init%n_particle)
  
  ave = 0.
 
  denom = (1 + beam_init%center(6)) * ele%value(beam_energy$)
  a_emitt = beam_init%a_norm_emitt * m_electron / denom
  b_emitt = beam_init%b_norm_emitt * m_electron / denom
  
  sigma_cut(1:4) = beam_init%sig_trans_cut
  sigma_cut(5)   = beam_init%sig_z_cut
  sigma_cut(6)   = beam_init%sig_e_cut
                                                 
  center = beam_init%center
  dpz_dz = beam_init%dpz_dz
  
  call make_v_mats(ele,v_mat,v_inv)

  bunch => beam%bunch(1)

  if (random_dist) then
    sigma(1) = sqrt(a_emitt * ele%x%beta)
    sigma(2) = sqrt(a_emitt / ele%x%beta)
    sigma(3) = sqrt(b_emitt * ele%y%beta)
    sigma(4) = sqrt(b_emitt / ele%y%beta)
    sigma(5) = beam_init%sig_z
    sigma(6) = beam_init%sig_e
 
    a = dpz_dz * sigma(5) / sigma(6)
    if (a > 1)  then
      call out_io (s_abort$, r_name, &
                   "dpz_dz MUST be < mode%sigE_E / mode%sig_z")
      call err_exit
    endif
 
    b = sqrt(1-a**2)
 
    call ran_seed(0)

    do i = 1, beam_init%n_particle
      do j = 1, 6
        do
          call ran_gauss(r(j))
          if (abs(r(j)) .le. sigma_cut(j)) exit
        enddo
      enddo
       
      !Initialize the Distribution
      bunch%particle(i)%r%vec(1) = sigma(1) *  r(1)
      bunch%particle(i)%r%vec(2) = - sigma(2) * (r(2) + r(1) * ele%x%alpha)
      bunch%particle(i)%r%vec(3) = sigma(3) *  r(3)
      bunch%particle(i)%r%vec(4) = - sigma(4) * (r(4) + r(3) * ele%y%alpha)
      bunch%particle(i)%r%vec(5) = sigma(5) *  r(5)
      bunch%particle(i)%r%vec(6) = sigma(6) * (r(6) * b + r(5) * a)
    
      !Include Dispersion
      bunch%particle(i)%r%vec(1:4) =  bunch%particle(i)%r%vec(1:4) + bunch%particle(i)%r%vec(6) * &
           (/ ele%x%eta, ele%x%etap, ele%y%eta, ele%y%etap /)
    
      !Include Coupling
      bunch%particle(i)%r%vec(1:4) = matmul(v_mat, bunch%particle(i)%r%vec(1:4))
    
      !Calculate the Distribution Centroid
      ave = ave + bunch%particle(i)%r%vec
    end do
  else ! ordered distribution
    n_ellipse = int(sqrt(real(beam_init%n_particle)))
    n_ellipse_parts = n_ellipse
    
    ! fudge particle distribution so total number of particles is correct
    n_too_many = 0
    do 
      n_more = beam_init%n_particle - (n_ellipse * n_ellipse_parts - n_too_many)
      if (n_more .eq. 0) exit
      n_factor = int(n_more / n_ellipse_parts)
      if (n_factor .eq. 0) then
        ! only one more ellipse
        n_ellipse = n_ellipse + 1
        n_too_many = n_ellipse_parts - n_more
      else
        ! several more ellipses
        n_ellipse = n_ellipse + n_factor
      endif
    enddo
      
    sigma(1:2) = a_emitt
    sigma(3:4) = b_emitt
    sigma(5) = beam_init%sig_z
    sigma(6) = beam_init%sig_e
    
    a = dpz_dz * sigma(5) / sigma(6)
    if (a > 1)  then
      call out_io (s_abort$, r_name, &
                   "dpz_dz MUST be < mode%sigE_E / mode%sig_z")
      call err_exit
    endif
    
    b = sqrt(1-a**2)
 
    x_max = (2.0*real(n_ellipse) - 1.0) / (4.0 * real(n_ellipse))
    i_part = 0
    do i = 1, n_ellipse
      y = (2.0*real(i) - 1.0) / (4.0 * real(n_ellipse))
      do j = 1, 6
        r(j) = inverse(renorm_probability_funct, y, -x_domain, x_domain, tol)
      enddo

      do j = 1, n_ellipse_parts
        i_part = i_part + 1
 
        phi(:) = (2.0*pi / real(n_ellipse_parts)) * real(j-1)
         
        if (i .eq. n_ellipse .and. n_too_many .ne. 0) then
          if (j .eq. (n_ellipse_parts - n_too_many + 1) ) exit
          phi(:) = (2.0*pi / real(n_ellipse_parts - n_too_many)) * (real(j-1))
        endif
        
        !Initialize the Distribution
        bunch%particle(i_part)%r%vec(1) = sqrt(sigma(1) * r(1)*ele%x%beta) * cos(phi(1))
        bunch%particle(i_part)%r%vec(2) = &
            - sqrt(sigma(1) * r(1) / ele%x%beta) * (ele%x%alpha * cos(phi(1)) + sin(phi(1)))

        bunch%particle(i_part)%r%vec(3) = sqrt(sigma(3) * r(3)*ele%y%beta) * cos(phi(3))
        bunch%particle(i_part)%r%vec(4) = &
            - sqrt(sigma(3) * r(3) / ele%y%beta) * (ele%y%alpha * cos(phi(3)) + sin(phi(3)))

        bunch%particle(i_part)%r%vec(5) = sigma(5) * r(5) * cos(phi(5))
        bunch%particle(i_part)%r%vec(6) = &
            sigma(6) * r(5) * (a * cos(phi(5)) + b * sin(phi(5)))
     
        !Include Dispersion
        bunch%particle(i_part)%r%vec(1:4) =  bunch%particle(i_part)%r%vec(1:4) + &
            bunch%particle(i_part)%r%vec(6) * (/ ele%x%eta, ele%x%etap, ele%y%eta, ele%y%etap /)
     
        !Include Coupling
        bunch%particle(i_part)%r%vec(1:4) = matmul(v_mat, bunch%particle(i_part)%r%vec(1:4))
     
        !Calculate the Distribution Centroid
        ave = ave + bunch%particle(i_part)%r%vec
      enddo
    end do

    ! Now scale to get correct emittances
!   do 
!     call calc_bunch_params (bunch, ele, params)
!     a_factor = 
!     if (abs(params%a%norm_emitt - beam_init%a_norm_emitt) .lt.
!                                     emitt_tol*beam_init%a_norm_emitt) exit
!     a_factor = beam_init%a_norm_emitt / params%a%norm_emitt
!     bunch%particle(:)%r%vec(1) = bunch%particle(:)%r%vec(1) + a_factor
      
    
  endif

  ave = ave / beam_init%n_particle
     
  ! Recenter the Distribution
  if (renormalize) then
    do i = 1, beam_init%n_particle
      bunch%particle(i)%r%vec = bunch%particle(i)%r%vec - ave
    enddo
  endif

  bunch%particle(:)%ix_lost = not_lost$
  
  ! move to specified center
  center(1:4) = matmul(v_mat, center(1:4))
  do i = 1, beam_init%n_particle
    bunch%particle(i)%r%vec = bunch%particle(i)%r%vec + center
  end do

  ! set particle charge
  bunch%particle(:)%charge = beam_init%bunch_charge / beam_init%n_particle
  
  ! init all bunches
  
  bunch%s_center = 0.0

  do i = 2, size(beam%bunch)
    call bunch_equal_bunch (beam%bunch(i), beam%bunch(1))
    beam%bunch(i)%s_center = (1-i) * beam_init%ds_bunch
  enddo
  
end subroutine init_beam_distribution

!--------------------------------------------------------------------------
! This is the probability function renormalized for integration only out to
! n_max

function renorm_probability_funct(x)

  use nr


  use precision_def

  implicit none

  real(rp) renorm_probability_funct, x

  renorm_probability_funct = probability_funct(x) / (2.0 * probability_funct(x_max))

end function
 
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
!  bunch     -- Bunch_struct
!  ele       -- ele_struct: element to find parameters at
!
! Output     -- bunch_params_struct
!                 %x%alpha; %x%beta; %x%gamma
!                 %x%sigma; %x%p_sigma
!                 %x%emitt; %x%dpx_dx
!                 %y%alpha; %y%beta; %y%gamma
!                 %y%sigma; %y%p_sigma
!                 %y%emitt; %y%dpx_dx
!                 %z%alpha; %z%beta; %z%gamma
!                 %z%sigma; %z%p_sigma
!                 %z%emitt; %z%dpx_dx
!                 %a%alpha; %a%beta; %a%gamma
!                 %a%sigma; %a%p_sigma
!                 %a%emitt; %a%dpx_dx
!                 %b%alpha; %b%beta; %b%gamma
!                 %b%sigma; %b%p_sigma
!                 %b%emitt; %b%dpx_dx
!                 %centroid
!-

subroutine calc_bunch_params (bunch, ele, params)

  implicit none

  type (bunch_struct), intent(in) :: bunch
  type (ele_struct), intent(in) :: ele
  type (bunch_params_struct) params
  type (coord_struct), automatic :: a_mode(size(bunch%particle))

  real(rp) exp_x2, exp_p_x2, exp_x_p_x, exp_x_d, exp_px_d
  real(rp) avg_energy 
  real(rp) avg_delta, exp_delta2
  real(rp) v_mat(4,4), v_inv_mat(4,4)

  integer i
  
! centroid and n_particle

  params%n_particle = count(bunch%particle%ix_lost == not_lost$)
  if (params%n_particle == 0) then
    ! zero everything
    params%centroid%vec = 0.0
    call zero_plane (params%x)
    call zero_plane (params%y)
    call zero_plane (params%z)
    call zero_plane (params%a)
    call zero_plane (params%b)
  endif
  
  do i = 1, 6
    params%centroid%vec(i) = sum(bunch%particle%r%vec(i), &
                              mask = (bunch%particle%ix_lost == not_lost$))
  enddo
  
  params%centroid%vec = params%centroid%vec / params%n_particle
  
  ! average energy
  avg_energy = sum((1+bunch%particle%r%vec(6)), & 
                              mask = (bunch%particle%ix_lost == not_lost$))
  avg_energy = avg_energy * ele%value(beam_energy$) / params%n_particle

  ! delta spread and center
  avg_delta = sum(bunch%particle%r%vec(6), & 
                              mask = (bunch%particle%ix_lost == not_lost$))
  avg_delta = avg_delta  / params%n_particle
  
  exp_delta2 = sum((bunch%particle%r%vec(6) - avg_delta)**2, &
                              mask = (bunch%particle%ix_lost == not_lost$))
  exp_delta2 = exp_delta2 / params%n_particle
  
  ! Projected Parameters
  ! X
  call find_expectations (bunch, bunch%particle(:)%r%vec(1), bunch%particle(:)%r%vec(2), &
                          exp_x2, exp_p_x2, exp_x_p_x, exp_x_d, exp_px_d, .false.)

  call param_stuffit (params%x, exp_x2, exp_p_x2, exp_x_p_x, exp_x_d, exp_px_d)
     
  ! Y
  call find_expectations (bunch, bunch%particle(:)%r%vec(3), bunch%particle(:)%r%vec(4), &
                          exp_x2, exp_p_x2, exp_x_p_x, exp_x_d, exp_px_d, .false.)

  call param_stuffit (params%y, exp_x2, exp_p_x2, exp_x_p_x, exp_x_d, exp_px_d)
  
  ! Z
  call find_expectations (bunch, bunch%particle(:)%r%vec(5), bunch%particle(:)%r%vec(6), &
                          exp_x2, exp_p_x2, exp_x_p_x, exp_x_d, exp_px_d, .false.)

  call param_stuffit (params%z, exp_x2, exp_p_x2, exp_x_p_x, exp_x_d, exp_px_d)
  
  !***
  ! Normal-Mode Parameters
  
  ! take out coupling
  call make_v_mats (ele, v_mat, v_inv_mat)
  do i = 1, size(a_mode)
    a_mode(i)%vec(1:4) = matmul(v_inv_mat, bunch%particle(i)%r%vec(1:4))
  enddo 
  
  ! A
  call find_expectations (bunch, a_mode(:)%vec(1), a_mode(:)%vec(2), &
                          exp_x2, exp_p_x2, exp_x_p_x, exp_x_d, exp_px_d,  .true.)

  call param_stuffit (params%a, exp_x2, exp_p_x2, exp_x_p_x, exp_x_d, exp_px_d)
     
  ! B
  call find_expectations (bunch, a_mode(:)%vec(3), a_mode(:)%vec(4), &
                          exp_x2, exp_p_x2, exp_x_p_x, exp_x_d, exp_px_d, .true.)

  call param_stuffit (params%b, exp_x2, exp_p_x2, exp_x_p_x, exp_x_d, exp_px_d)
     
  
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
  param%sigma      = 0.0
  param%p_sigma    = 0.0
  param%dpx_dx     = 0.0
  param%norm_emitt = 0.0

end subroutine zero_plane
  
!----------------------------------------------------------------------
subroutine find_expectations (bunch, x, p_x, exp_x2, exp_p_x2, exp_x_p_x, &
                              exp_x_d, exp_px_d, normal_mode_flag)

  implicit none

  type (bunch_struct), intent(in) :: bunch
  real(rp), intent(in)  :: x(:)
  real(rp), intent(in)  :: p_x(:)
  real(rp), intent(out) ::  exp_x2, exp_p_x2, exp_x_p_x, exp_x_d, exp_px_d
  real(rp) avg_x, avg_p_x, eta, etap
  
  integer i

  logical normal_mode_flag

  if (size(x) .ne. size(p_x)) then
    exp_x2     = 0.0
    exp_p_x2   = 0.0
    exp_x_p_x = 0.0
    return
  endif

  avg_x = sum(x, mask = (bunch%particle%ix_lost == not_lost$))/params%n_particle
  avg_p_x = sum(p_x, mask = (bunch%particle%ix_lost == not_lost$))/params%n_particle
 
  ! take out dispersion
  exp_x_d   = sum((x - avg_x)*(bunch%particle(:)%r%vec(6) - avg_delta),&
                                    mask = (bunch%particle%ix_lost .eq. not_lost$))
  exp_px_d  = sum((p_x - avg_p_x)*(bunch%particle(:)%r%vec(6) - avg_delta),&
                                    mask = (bunch%particle%ix_lost .eq. not_lost$))
  exp_x2    = sum((x - avg_x)**2, mask = (bunch%particle%ix_lost .eq. not_lost$))
  exp_p_x2  = sum((p_x - avg_p_x)**2, mask = (bunch%particle%ix_lost .eq. not_lost$))
  exp_x_p_x = sum((x - avg_x)*(p_x - avg_p_x), mask = (bunch%particle%ix_lost .eq. not_lost$))
   
  exp_x2    = exp_x2    / params%n_particle 
  exp_p_x2  = exp_p_x2  / params%n_particle
  exp_x_p_x = exp_x_p_x / params%n_particle
  exp_x_d   = exp_x_d   / params%n_particle
  exp_px_d  = exp_px_d  / params%n_particle
  
  
  if (normal_mode_flag) then
    eta   = exp_x_d / exp_delta2
    etap  = exp_px_d / exp_delta2

    exp_x2    = exp_x2 - 2*eta*exp_x_d + (eta**2)*exp_delta2
    exp_p_x2  = exp_p_x2 - 2*etap*exp_px_d + (etap**2)*exp_delta2
    exp_x_p_x = exp_x_p_x - etap*exp_x_d - eta*exp_px_d + eta*etap*exp_delta2
  endif

end subroutine find_expectations

!----------------------------------------------------------------------
! contains

subroutine param_stuffit (param, exp_x2, exp_p_x2, exp_x_p_x, exp_x_d, exp_px_d)

  implicit none

  type (bunch_param_struct), intent(out) :: param
  real(rp), intent(in) :: exp_x2, exp_p_x2, exp_x_p_x, exp_x_d, exp_px_d
  real(rp) emitt

  emitt = SQRT(exp_x2*exp_p_x2 - exp_x_p_x**2)

  param%alpha = exp_x_p_x / emitt
  param%beta  = exp_x2 / emitt
  param%gamma = exp_p_x2 / emitt
  
  param%eta   = exp_x_d / exp_delta2
  param%etap  = exp_px_d / exp_delta2

  param%norm_emitt = (avg_energy/m_electron) * emitt

  param%sigma = SQRT(exp_x2)
  param%p_sigma = SQRT(exp_p_x2)

  param%dpx_dx = exp_x_p_x / exp_x2

end subroutine param_stuffit

end subroutine calc_bunch_params
  
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

!

  n_bun = size(beam2%bunch)

  if (size(beam1%bunch) /= size(beam2%bunch)) then
    do i = 1, size(beam1%bunch)
      deallocate (beam1%bunch(i)%particle)
    enddo
    deallocate (beam1%bunch)
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
