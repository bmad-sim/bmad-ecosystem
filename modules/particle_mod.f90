module particle_mod

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine create_distribution (ele, mode, bunch, renormalize, center, dpz_dz)
!
! Subroutine to initialize a distribution of particles matched to 
! the Twiss parameters, centroid position, and Energy - z correlation
! as specified.
!
! Note:  This routine requires calculation of the Twiss parameter.
! Use twiss_propagate_all.
!
! Note:  This routine requires calculation of emittances and 
! longitudinal beam parameters.  Use radiation_integrals.
!
! Modules needed:
!   particle_mod
!
! Input:
!   ele         -- Ele_struct: element to initialize distribution at.
!     %x          -- Twiss parameters.
!     %y          -- Twiss parameters.
!   mode        -- Modes_struct: longitudinal beam parameters and 
!                   transverse emittances to match distrubition to.
!   bunch(:)    -- Coord_struct, allocatable:  Array of coord_struct
!                   allocated to the size of how many particles you want.
!   renormalize -- Logical: If True then distribution is rescaled to 
!                   the desired centroid (to take care of 
!                   possible statistical errors in distribution).
!   center(:)   -- Real(6), optional: desired 6-D centroid of beam.
!                   If renormalize is True then beam will be centered here.
!                   Distribution will be centered on zero if not specified
!                   and renormalize is True.
!   dpz_dz      -- Real, optional: desired Energy-z correlation of beam.  
!                   Default is no correlation.
!
!                   Make sure: dpz_dz < mode%sigE_E / mode%sig_z
!
! Output:
!   bunch(:) -- Coord_struct, allocatable: Array of coord_structs
!                holding the resulting particle distribution.
!
!-

subroutine create_distribution(ele, mode, bunch, renormalize, center, dpz_dz)

  use random_mod
  use bmad

  implicit none

  type (ele_struct) ele
  type (modes_struct) mode
  type (coord_struct), allocatable :: bunch(:)

  real(rp) center2(6), ave(6),ave_sq(6)
  real(rp) sigma(6), dpz_dz2, a, b
  real(rp), optional :: center(6)
  real(rp), optional :: dpz_dz
  real(rp) r(7), v_mat(4,4), v_inv(4,4)

  integer number, i, j
  logical renormalize

!

  number = size(bunch)
  ave = 0.
  ave_sq = 0.

  sigma(1) = sqrt(mode%a%emittance*ele%x%beta)
  sigma(2) = sigma(1) / ele%x%beta
  sigma(3) = sqrt(mode%b%emittance*ele%y%beta)
  sigma(4) = sigma(3) / ele%y%beta
  sigma(5) = mode%sig_z
  sigma(6) = mode%sigE_E

  if (present(center)) then
     center2 = center
  else
     center2 = 0.
  end if

  if (present(dpz_dz)) then
     dpz_dz2 = dpz_dz
  else
     dpz_dz2 = 0.
  end if

  a = dpz_dz2 * sigma(5) / sigma(6)
  if (a > 1) print*, 'Program Crashing'
  b = sqrt(1-a**2)

  call make_v_mats(ele,v_mat,v_inv)
  call ran_seed(2)

  do i = 1, number
     do j = 1, 7
        call ran_gauss(r(j))
     end do
     
     !Initialize the Distribution
     bunch(i)%vec(1) = sigma(1) *  r(1)
     bunch(i)%vec(2) = sigma(2) * (r(2) - r(1) * ele%x%alpha)
     bunch(i)%vec(3) = sigma(3) *  r(3)
     bunch(i)%vec(4) = sigma(4) * (r(4) - r(3) * ele%y%alpha)
     bunch(i)%vec(5) = sigma(5) *  r(5)
     bunch(i)%vec(6) = sigma(6) * (r(6) * b + r(5) * a)

     !Include Dispersion
     bunch(i)%vec(1:4) =  bunch(i)%vec(1:4) + bunch(i)%vec(6) * &
          (/ ele%x%eta, ele%x%etap, ele%y%eta, ele%y%etap /)

     !Include Coupling
     bunch(i)%vec(1:4) = matmul(v_mat, bunch(i)%vec(1:4))

     !Calculate the Distribution Centroid
     ave = ave + bunch(i)%vec
     ave_sq = ave_sq + bunch(i)%vec**2
  end do
  
  ave = ave / number
  ave_sq = ave_sq / number

  if (renormalize) then
     !Recenter the Distribution
     do i = 1, number
        bunch(i)%vec = bunch(i)%vec - ave + center2
     end do
  endif

end subroutine create_distribution





end module
