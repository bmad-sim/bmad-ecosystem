module macro_particle_mod

  use bmad_struct
  use bmad_interface

  type macro_particle_struct
    type (coord_struct) r
    real(rp) sigma(21)
    real(rp) k_loss
    real(rp) charge
    real(rp) effective_charge ! effective charge for long sr wake calc
    integer :: iz = 0         ! index to ordering of particles in z.
  end type

  type bunch_struct
    type (macro_particle_struct), allocatable :: macro(:)
    real(rp) charge
  end type

  type beam_struct
    type (bunch_struct), allocatable :: bunch(:)
  end type

  integer, parameter :: s11$ = 1, s12$ = 2, s13$ = 3, s14$ =  4, s15$ =  5
  integer, parameter :: s16$ = 6, s22$ = 7, s23$ = 8, s24$ = 9
  integer, parameter :: s25$ = 10, s26$ = 11, s33$ = 12, s34$ = 13, s35$ = 14
  integer, parameter :: s36$ = 15, s44$ = 16, s45$ = 17, s46$ = 18
  integer, parameter :: s55$ = 19, s56$ = 20, s66$ = 21

  type macro_particle_com_struct
    real(rp) ::  dz_bin = 0.5 ! in units of the spacing between wakefield pts
  end type

  type (macro_particle_com_struct), save :: mp_com

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track_all_beam (ring, beam, ix1, ix2)
!
! Subroutine to track a beam of macroparticles from the end of
! ring%ele_(ix1) Through to the end of ring%ele_(ix2).
!
! Modules needed:
!   use macro_particle_mod
!
! Input:
!   ring -- Ring_struct: Lattice to track through.
!   beam -- Beam_struct: Collection of macroparticles.
!   ix1  -- Index of starting element (this element is NOT tracked through).
!   ix2  -- Index of ending element.
!
! Output:
!   beam -- Beam_struct: macroparticles at end of the tracking.
!-

subroutine track_all_beam (ring, beam, ix1, ix2)

  implicit none

  type (ring_struct), intent(in) :: ring
  type (beam_struct) :: beam

  integer, optional, intent(in) :: ix1, ix2
  integer i, j, i1, i2

! Init

  i1 = 0
  if (present(ix1)) i1 = ix1
  i2 = ring%n_ele_ring
  if (present(ix2)) i2 = ix2

!

  do i = i1+1, i2
    do j = 1, size(beam%bunch)
      call track1_bunch (beam%bunch(j), ring%ele_(i), ring%param)
    enddo
  enddo

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
!-

subroutine track1_bunch (bunch, ele, param)

  implicit none

  type (bunch_struct) bunch
  type (ele_struct) ele, rf_ele
  type (param_struct) param

  real(rp) charge
  integer i, ix

  logical sr_wake_on

!------------------------------------------------
! Without wakefields just track through

  if (ele%key /= lcavity$ .or. .not. bmad_com%sr_wakes_on) then
    do i = 1, size(bunch%macro)
      call track1_macro_particle (bunch%macro(i), ele, param, bunch%macro(i))
    enddo
    return
  endif

!------------------------------------------------
! For an lcavity without a wakefield file use the e_loss attribute 
! to calculate energy loss. The effective charge of a bunch is:
!   2*sum[charge_of_previous_bunches] + charge_of_this_bunch

  if (.not. associated(ele%wake%sr)) then
    call order_macro_particles_in_z (bunch)
    charge = param%charge
    param%charge = 0
    do ix = 1, size(bunch%macro)
      i = bunch%macro(ix)%iz
      param%charge = param%charge + bunch%macro(i)%charge
      call track1_macro_particle (bunch%macro(i), ele, param, bunch%macro(i))
      param%charge = param%charge + bunch%macro(i)%charge
    enddo
    param%charge = charge
    return
  endif

!------------------------------------------------
! This calculation is for an lcavity with full wakefields.
! With wakefields put them in at the half way point.
! First track half way through.

  rf_ele = ele
  rf_ele%value(l$) = ele%value(l$) / 2
  rf_ele%value(beam_energy$) = &
            (ele%value(energy_start$) + ele%value(beam_energy$)) / 2

  charge = param%charge
  param%charge = 0

  call sr_long_wake_calc (bunch, ele)

  do i = 1, size(bunch%macro)
    bmad_com%k_loss = bunch%macro(i)%k_loss
    call track1_macro_particle (bunch%macro(i), rf_ele, param, bunch%macro(i))
  enddo

! Put in the wakefields

  rf_ele%value(l$) = ele%value(l$)
  call track1_sr_trans_wake (bunch, rf_ele)

! Track the last half of the lcavity

  rf_ele%value(l$) = ele%value(l$) / 2
  rf_ele%value(energy_start$) = rf_ele%value(beam_energy$)
  rf_ele%value(beam_energy$) = ele%value(beam_energy$)

  param%charge = 0
  do i = 1, size(bunch%macro)
    bmad_com%k_loss = bunch%macro(i)%k_loss
    call track1_macro_particle (bunch%macro(i), rf_ele, param, bunch%macro(i))
  enddo

  param%charge = charge
  bmad_com%k_loss = 0

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine sr_long_wake_calc (bunch, ele)
!
! Subroutine to put in the longitudinal component of the
! short range wake fields. 
! This routine is not really meant  for general use.
!-

subroutine sr_long_wake_calc (bunch, ele)

  implicit none

  type (bunch_struct) bunch
  type (ele_struct) ele

  real(rp) z_ave, charge, z0
  real(rp) f1, f2, z, dE, E_rel, fact, dz_wake

  integer i, ix, i0, i1, ix0, ix1, iw, n_wake, n_macro

! Init

  call order_macro_particles_in_z (bunch)

  n_wake = ubound(ele%wake%sr, 1)
  dz_wake = ele%wake%sr(n_wake)%z / n_wake
  n_macro = size(bunch%macro)

  ix0 = 1
  bunch%macro(:)%k_loss = 0

! loop over all bins
! ix0, ix1 give the range of the current bin

  do

    i0 = bunch%macro(ix0)%iz

    z0 = bunch%macro(i0)%r%vec(5)

    charge = 0
    do ix = ix0, n_macro
      i = bunch%macro(ix)%iz
      if (bunch%macro(i)%r%vec(5) < z0 - mp_com%dz_bin * dz_wake) exit
      charge = charge + bunch%macro(i)%charge 
      z_ave = z_ave + bunch%macro(i)%charge * bunch%macro(i)%r%vec(5)
      ix1 = ix 
    enddo

    do ix = ix0, ix1
      i = bunch%macro(ix)%iz
      bunch%macro(i)%k_loss = bunch%macro(i)%k_loss + &
                                      abs(charge) * ele%wake%sr(0)%long / 2
    enddo

! now apply the wakefields to the other macro-particles.

    do ix = ix1+1, n_macro
      i = bunch%macro(ix)%iz
      z = z_ave - bunch%macro(i)%r%vec(5)
      iw = z / dz_wake
      f2 = z/dz_wake - iw
      f1 = 1 - f2
      bunch%macro(i)%k_loss = bunch%macro(i)%k_loss + abs(charge) * &
                  (ele%wake%sr(iw)%long * f1 + ele%wake%sr(iw+1)%long * f2)
    enddo

    if (ix1 == n_macro) return
    ix0 = ix1 + 1

  enddo

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track1_sr_trans_wake (bunch, ele)
!
! Subroutine to put in the transverse component of the
! short range wake fields. 
! This routine is not really meant  for general use.
!-

subroutine track1_sr_trans_wake (bunch, ele)

  implicit none

  type (bunch_struct) bunch
  type (ele_struct) ele

  real(rp) x_ave, y_ave, z_ave, charge, z0
  real(rp) f1, f2, z, dE, E_rel, fact, dz_wake

  integer i, ix, i0, i1, ix0, ix1, iw, n_wake, n_macro

! Init

  call order_macro_particles_in_z (bunch)

  n_wake = ubound(ele%wake%sr, 1)
  dz_wake = ele%wake%sr(n_wake)%z / n_wake
  n_macro = size(bunch%macro)

! loop over all bins
! ix0, ix1 give the range of the current bin

  ix0 = 1

  do

    i0 = bunch%macro(ix0)%iz

    z0 = bunch%macro(i0)%r%vec(5)

    charge = 0
    x_ave = 0
    y_ave = 0
    z_ave = 0
    do ix = ix0, n_macro
      i = bunch%macro(ix)%iz
      if (bunch%macro(i)%r%vec(5) < z0 - mp_com%dz_bin * dz_wake) exit
      charge = charge + bunch%macro(i)%charge / 2
!      bunch%macro(i)%r%vec(6) = bunch%macro(i)%r%vec(6) - &
!          abs(charge) * ele%wake%sr(0)%long * ele%value(l$) / &
!          ele%value(beam_energy$)
      charge = charge + bunch%macro(i)%charge / 2
      x_ave = x_ave + bunch%macro(i)%charge * bunch%macro(i)%r%vec(1)
      y_ave = y_ave + bunch%macro(i)%charge * bunch%macro(i)%r%vec(3)
      z_ave = z_ave + bunch%macro(i)%charge * bunch%macro(i)%r%vec(5)
      ix1 = ix 
    enddo

    if (charge /= 0) then
      x_ave = x_ave / charge
      y_ave = y_ave / charge
      z_ave = z_ave / charge
    endif

! now apply the wakefields to the other macro-particles.

    do ix = ix1+1, n_macro
      i = bunch%macro(ix)%iz
      z = z_ave - bunch%macro(i)%r%vec(5)
      iw = z / dz_wake
      f2 = z/dz_wake - iw
      f1 = 1 - f2
      fact = abs(charge) * ele%value(l$) 
!      dE = fact * (ele%wake%sr(iw)%long * f1 + ele%wake%sr(iw+1)%long * f2)
!      bunch%macro(i)%r%vec(6) = bunch%macro(i)%r%vec(6) - &
!                                         dE / ele%value(beam_energy$)
      fact = fact * (ele%wake%sr(iw)%trans*f1 + ele%wake%sr(iw+1)%trans*f2) / &
                                                  ele%value(beam_energy$)
      bunch%macro(i)%r%vec(2) = bunch%macro(i)%r%vec(2) + fact * x_ave 
      bunch%macro(i)%r%vec(4) = bunch%macro(i)%r%vec(4) + fact * y_ave
    enddo

    if (ix1 == n_macro) return
    ix0 = ix1 + 1

  enddo

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track1_macro_particle (start, ele, param, end)
!
! Subroutine to track a macro-particle through an element.
! Note: the transfer ele%mat6 matrix is changed and not restored
! during traking.
!
! Modules needed:
!   use macro_particle_mod
!
! Input:
!   start  -- Macro_particle_struct: Starting coords.
!   ele    -- Ele_struct: Element to track through.
!   param  -- Param_struct: Global parameters.
!
! Output:
!   end    -- macro_particle_struct: Ending coords.
!-

subroutine track1_macro_particle (start, ele, param, end)

  implicit none

  type (macro_particle_struct) start
  type (macro_particle_struct) end
  type (ele_struct) ele
  type (param_struct) param

  real(rp) l, l2, s(21), m4(4,4), s_mat4(4,4), s_mat6(6,6)

! Very simple cases

  select case (ele%key)
  case (drift$, ecollimator$, elseparator$, hkicker$, instrument$, &
                    kicker$, marker$, monitor$, rcollimator$)

    call track1 (start%r, ele, param, end%r)

    s = start%sigma
    end%sigma = s
    l = ele%value(l$) / (1 + start%r%vec(6))
    l2 = l*l

    end%sigma(s11$) = s(s11$) + 2 * l * s(s12$) + l2 * s(s22$)
    end%sigma(s12$) = s(s12$) +     l * s(s22$)  
    end%sigma(s13$) = s(s13$) +     l * s(s14$) + l * s(s23$) + l2 * s(s24$)
    end%sigma(s14$) = s(s14$) +     l * s(s24$)  
    end%sigma(s15$) = s(s15$) +     l * s(s25$)
    end%sigma(s16$) = s(s16$) +     l * s(s26$)
    end%sigma(s23$) = s(s23$) +     l * s(s24$)
    end%sigma(s33$) = s(s33$) + 2 * l * s(s34$) + l2 * s(s44$)
    end%sigma(s34$) = s(s34$) +     l * s(s44$)  
    end%sigma(s35$) = s(s35$) +     l * s(s45$)
    end%sigma(s36$) = s(s36$) +     l * s(s46$)
    return

  end select

! Simple case where longitudinal motion can be ignored.

  if (start%sigma(s55$) == 0 .and. start%sigma(s66$) == 0 .and. &
                                                  ele%key /= sbend$) then

    call make_mat6 (ele, param, start%r, end%r)  

    s = start%sigma
    end%sigma = s
    m4 = ele%mat6(1:4,1:4)

    if (all(ele%mat6(1:2,3:4) == 0) .and. all(ele%mat6(3:4,1:2) == 0)) then
      end%sigma(s11$) = m4(1,1)*m4(1,1)*s(s11$) + m4(1,2)*m4(1,2)*s(s22$) + &
                                              2*m4(1,1)*m4(1,2)*s(s12$) 
      end%sigma(s12$) = m4(1,2)*m4(2,1)*s(s12$) + m4(1,2)*m4(2,2)*s(s22$) + &
                        m4(1,1)*m4(2,1)*s(s11$) + m4(1,1)*m4(2,2)*s(s12$)
      end%sigma(s13$) = m4(1,1)*m4(3,4)*s(s14$) + m4(1,2)*m4(3,3)*s(s23$) + &
                        m4(1,2)*m4(3,4)*s(s24$) + m4(1,1)*m4(3,3)*s(s13$)
      end%sigma(s14$) = m4(1,1)*m4(4,3)*s(s13$) + m4(1,2)*m4(4,4)*s(s24$) + &
                        m4(1,2)*m4(4,3)*s(s23$) + m4(1,1)*m4(4,4)*s(s14$) 
      end%sigma(s22$) = m4(2,2)*m4(2,2)*s(s22$) + m4(2,1)*m4(2,1)*s(s11$) + &
                                              2*m4(2,2)*m4(2,1)*s(s12$) 
      end%sigma(s23$) = m4(2,2)*m4(3,4)*s(s24$) + m4(2,1)*m4(3,3)*s(s13$) + &
                        m4(2,1)*m4(3,4)*s(s14$) + m4(2,2)*m4(3,3)*s(s23$)
      end%sigma(s24$) = m4(2,2)*m4(4,3)*s(s23$) + m4(2,1)*m4(4,4)*s(s14$) + &
                        m4(2,1)*m4(4,3)*s(s13$) + m4(2,2)*m4(4,4)*s(s24$) 
      end%sigma(s33$) = m4(3,3)*m4(3,3)*s(s33$) + m4(3,4)*m4(3,4)*s(s44$) + &
                                              2*m4(3,3)*m4(3,4)*s(s34$) 
      end%sigma(s34$) = m4(3,4)*m4(4,3)*s(s34$) + m4(3,4)*m4(4,4)*s(s44$) + &
                        m4(3,3)*m4(4,3)*s(s33$) + m4(3,3)*m4(4,4)*s(s34$)
      end%sigma(s44$) = m4(4,4)*m4(4,4)*s(s44$) + m4(4,3)*m4(4,3)*s(s33$) + &
                                              2*m4(4,4)*m4(4,3)*s(s34$) 
    else
      call mp_sigma_to_mat (s, s_mat4)
      s_mat4 = matmul(m4, matmul(s_mat4, transpose(m4)))
      call mat_to_mp_sigma (s_mat4, end%sigma)
    endif

    return
  endif

! Full tracking. 

    call make_mat6 (ele, param, start%r, end%r)  
    call mp_sigma_to_mat (start%sigma, s_mat6)
    s_mat6 = matmul(ele%mat6, matmul(s_mat6, transpose(ele%mat6)))
    call mat_to_mp_sigma (s_mat6, end%sigma)

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine order_macro_particles_in_z (bunch)
!
! Subroutine to order the macro particles longitudinally.
! This routine is not really meant for general use.
!-

subroutine order_macro_particles_in_z (bunch)

  implicit none

  type (bunch_struct) bunch
  integer i, j1, j2
  logical ordered

! Init if needed

  if (any (bunch%macro(:)%iz == 0)) then
    do i = 1, size(bunch%macro)
      bunch%macro(i)%iz = i
    enddo
  endif

! Order is from large z (head of bunch) to small z.

  do 

    ordered = .true.

    do i = 1, size(bunch%macro)-1
      j1 = bunch%macro(i)%iz
      j2 = bunch%macro(i+1)%iz
      if (bunch%macro(j1)%r%vec(5) < bunch%macro(j2)%r%vec(5)) then
        bunch%macro(i:i+1)%iz = bunch%macro(i+1:i:-1)%iz
        ordered = .false.
      endif
    enddo

    if (ordered) return

  enddo 

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine mp_sigma_to_mat (sigma, mat)
!
! Subroutine to convert a linear array of macro-particle sigmas to a 
! sigma matrix. 
!
! Modules needed:
!   use macro_particle_mod
!
! Input:
!   sigma(21) -- Real(rp): array of sigmas.
!
! Output:
!   mat(:,:)  -- Real(rp): Sigma matrix. 
!                 If size(mat, 1) = 4 then only the transverse sigma matrix
!                                          is constructed.
!                 If size(mat, 1) = 6 then the full sigma matrix
!                                          is constructed.
!-

subroutine mp_sigma_to_mat (s, mat)

  implicit none

  real(rp), intent(in) :: s(:)
  real(rp), intent(out) :: mat(:,:)

!  

  if (size(mat, 1) == 4 .and. size(mat, 2) == 4) then
    mat(1,:) = (/ s(s11$), s(s12$), s(s13$), s(s14$) /)
    mat(2,:) = (/ s(s12$), s(s22$), s(s23$), s(s24$) /)
    mat(3,:) = (/ s(s13$), s(s23$), s(s33$), s(s34$) /)
    mat(4,:) = (/ s(s14$), s(s24$), s(s34$), s(s44$) /)
  elseif (size(mat, 1) == 6 .and. size(mat,2) == 6) then
    mat(1,:) = (/ s(s11$), s(s12$), s(s13$), s(s14$), s(s15$), s(s16$)/)
    mat(2,:) = (/ s(s12$), s(s22$), s(s23$), s(s24$), s(s25$), s(s26$)/)
    mat(3,:) = (/ s(s13$), s(s23$), s(s33$), s(s34$), s(s35$), s(s36$)/)
    mat(4,:) = (/ s(s14$), s(s24$), s(s34$), s(s44$), s(s45$), s(s46$)/)
    mat(5,:) = (/ s(s15$), s(s25$), s(s35$), s(s45$), s(s55$), s(s56$)/)
    mat(6,:) = (/ s(s16$), s(s26$), s(s36$), s(s46$), s(s56$), s(s66$)/)
  else
    print *, 'ERROR IN MP_SIGMA_TO_MAT: MATRIX SIZE NOT 4 OR 6!'
    call err_exit
  endif

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine mat_to_mp_sigma (sigma, mat)
!
! Subroutine to convert a sigma matrix. to a linear array of 
! macro-particle sigmas 
!
! Modules needed:
!   use macro_particle_mod
!
! Input:
!   mat(:,:)  -- Real(rp): Sigma matrix. 
!                 If size(mat, 1) = 4 then only the transverse sigmas are
!                                          converted.
!                 If size(mat, 1) = 6 then all the sigmas are
!                                          converted.
!
! Output:
!   sigma(21) -- Real(rp): array of sigmas.
!-

subroutine mat_to_mp_sigma (mat, sig)

  implicit none

  real(rp), intent(in) :: mat(:,:)
  real(rp), intent(out) :: sig(:)

  sig(s11$) = mat(1,1)
  sig(s12$) = mat(1,2)
  sig(s13$) = mat(1,3)
  sig(s14$) = mat(1,4)
  sig(s22$) = mat(2,2)
  sig(s23$) = mat(2,3)
  sig(s24$) = mat(2,4)
  sig(s33$) = mat(3,3)
  sig(s34$) = mat(3,4)
  sig(s44$) = mat(4,4)

  if (size(mat, 1) == 4 .and. size(mat, 2) == 4) then
    return
  elseif (size(mat, 1) == 6 .and. size(mat,2) == 6) then
    sig(s15$) = mat(1,5)
    sig(s16$) = mat(1,6)
    sig(s25$) = mat(2,5)
    sig(s26$) = mat(2,6)
    sig(s35$) = mat(3,5)
    sig(s36$) = mat(3,6)
    sig(s45$) = mat(4,5)
    sig(s46$) = mat(4,6)
    sig(s55$) = mat(5,5)
    sig(s56$) = mat(5,6)
    sig(s66$) = mat(6,6)
  else
    print *, 'ERROR IN NAT_TO_MP_SIGMA: MATRIX SIZE IS NOT 4 OR 6!'
    call err_exit
  endif

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine mp_to_momentum_coords (mp, energy0)
!
! Subroutine to convert macro-particle coords from 
!     (x, x', y, y', z, E)
! to
!     (x, px, y, py, z, pz)
!
! Modules needed:
!   use macro_particle_mod
!
! Input:
!   mp      -- macro_particle_struct: Macro-particle with angular coords.
!   energy0 -- real(rp): Reference energy.
!
! Output:
!   mp      -- macro_particle_struct: Macro-particle with momentum coords.
!-

subroutine mp_to_momentum_coords (mp, energy0)

  implicit none

  type (macro_particle_struct), target :: mp

  real(rp), pointer :: s(:)
  real(rp), intent(in) :: energy0
  real(rp) f, f2, e, xp0, yp0

!

  f = mp%r%vec(6) / energy0
  f2 = f * f
  e = energy0

  xp0 = mp%r%vec(2)
  yp0 = mp%r%vec(4)

  mp%r%vec(2) = mp%r%vec(2) * f
  mp%r%vec(4) = mp%r%vec(4) * f
  mp%r%vec(6) = f - 1

  s => mp%sigma

  s(s16$) = s(s16$) / e
  s(s26$) = s(s26$) / e
  s(s36$) = s(s36$) / e
  s(s46$) = s(s46$) / e
  s(s56$) = s(s56$) / e
  s(s66$) = s(s66$) / (e * e)

  s(s12$) = s(s12$) * f  +   s(s16$) * xp0
  s(s14$) = s(s14$) * f  +   s(s16$) * yp0
  s(s22$) = s(s22$) * f2 + 2*s(s26$) * xp0 + s(s66$) * xp0 * xp0
  s(s23$) = s(s23$) * f  +   s(s36$) * xp0
  s(s24$) = s(s24$) * f2 +   s(s46$) * xp0 + s(s26$) * yp0 + s(s66$) * xp0 * yp0
  s(s25$) = s(s25$) * f  +   s(s56$) * xp0
  s(s26$) = s(s26$) * f  +   s(s66$) * xp0
  s(s34$) = s(s34$) * f  +   s(s36$) * yp0
  s(s44$) = s(s44$) * f2 + 2*s(s46$) * yp0 + s(s66$) * yp0 * yp0
  s(s45$) = s(s45$) * f  +   s(s56$) * yp0
  s(s46$) = s(s46$) * f  +   s(s66$) * yp0

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine mp_to_angle_coords (mp, energy0)
!
! Subroutine to convert macro-particle coords from 
!     (x, px, y, py, z, pz)
! to
!     (x, x', y, y', z, E)
!
! Modules needed:
!   use macro_particle_mod
!
! Input:
!   mp      -- macro_particle_struct: Macro-particle with momentum coords.
!   energy0 -- real(rp): Reference energy.
!
! Output:
!   mp      -- macro_particle_struct: Macro-particle with angular coords.
!-

subroutine mp_to_angle_coords (mp, energy0)

  implicit none

  type (macro_particle_struct), target :: mp

  real(rp), pointer :: s(:)
  real(rp), intent(in) :: energy0
  real(rp) f, f2, e, xp0, yp0

!

  f = 1 + mp%r%vec(6)
  f2 = f * f
  e = energy0

  mp%r%vec(2) = mp%r%vec(2) / f
  mp%r%vec(4) = mp%r%vec(4) / f
  mp%r%vec(6) = energy0 * f 

  xp0 = mp%r%vec(2) / f2
  yp0 = mp%r%vec(4) / f2

  s => mp%sigma

  s(s12$) = s(s12$) / f  -   s(s16$) * xp0
  s(s14$) = s(s14$) / f  -   s(s16$) * yp0
  s(s22$) = s(s22$) / f2 - 2*s(s26$) * xp0 + s(s66$) * xp0 * xp0
  s(s23$) = s(s23$) / f  -   s(s36$) * xp0
  s(s24$) = s(s24$) / f2 -   s(s46$) * xp0 - s(s26$) * yp0 + s(s66$) * xp0 * yp0
  s(s25$) = s(s25$) / f  -   s(s56$) * xp0
  s(s26$) = s(s26$) / f  -   s(s66$) * xp0
  s(s34$) = s(s34$) / f  -   s(s36$) * yp0
  s(s44$) = s(s44$) / f2 - 2*s(s46$) * yp0 + s(s66$) * yp0 * yp0
  s(s45$) = s(s45$) / f  -   s(s56$) * yp0
  s(s46$) = s(s46$) / f  -   s(s66$) * yp0

  s(s16$) = s(s16$) * e
  s(s26$) = s(s26$) * e
  s(s36$) = s(s36$) * e
  s(s46$) = s(s46$) * e
  s(s56$) = s(s56$) * e
  s(s66$) = s(s66$) * e * e

end subroutine

end module
