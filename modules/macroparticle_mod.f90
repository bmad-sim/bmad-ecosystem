module macroparticle_mod

  use bmad_struct
  use bmad_interface

  type macro_init_twiss_struct
    real(rp) beta, alpha, emit
  end type

  type macro_init_struct
    type (macro_init_twiss_struct) x, y
    real(rp) E_0            ! Nominal Energy (eV).
    real(rp) center(6) ! Bench center.
    real(rp) n_part    ! Number of particles per bunch.
    real(rp) ds_bunch  ! Distance between bunches.
    real(rp) sig_z     ! Z sigma in m.
    real(rp) sig_e     ! e_sigma in eV.
    real(rp) sig_e_cut ! Energy cut in sigmas.
    real(rp) sig_z_cut ! Z cut in sigmas.
    integer n_bunch    ! Number of bunches.
    integer n_slice    ! Number of slices per bunch.
    integer n_macro    ! Number of macroparticles per slice.
  endtype

  type macro_struct
    type (coord_struct) r   ! Center of the macroparticle
    real(rp) sigma(21)      ! Sigma matrix.
    real(rp) :: sig_z = 0   ! longitudinal macroparticle length (m).
    real(rp) k_loss         ! loss factor (V/m). scratch variable for tracking.
    real(rp) charge         ! charge in a macroparticle (Coul).
    integer :: iz = 0       ! index to ordering of particles in z.
  end type

  type slice_struct
    type (macro_struct), allocatable :: macro(:)
    real(rp) charge   ! total charge in a slice (Coul).
  end type

  type bunch_struct
    type (slice_struct), allocatable :: slice(:)
    real(rp) charge   ! total charge in a bunch (Coul).
    real(rp) s_center ! longitudinal center of bunch (m).
  end type

  type beam_struct
    type (bunch_struct), allocatable :: bunch(:)
  end type

  integer, parameter :: s11$ = 1, s12$ = 2, s13$ = 3, s14$ =  4, s15$ =  5
  integer, parameter :: s16$ = 6, s22$ = 7, s23$ = 8, s24$ = 9
  integer, parameter :: s25$ = 10, s26$ = 11, s33$ = 12, s34$ = 13, s35$ = 14
  integer, parameter :: s36$ = 15, s44$ = 16, s45$ = 17, s46$ = 18
  integer, parameter :: s55$ = 19, s56$ = 20, s66$ = 21

  type macroparticle_com_struct
    real(rp) ::  sig_z_min = 5e-6 ! min effective macroparticle length.
  end type

  type (macroparticle_com_struct), save :: mp_com

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
!   use macroparticle_mod
!
! Input:
!   ring     -- Ring_struct: Lattice to track through.
!   beam(0:) -- Beam_struct: Beam at end of i^th element.
!   ix1      -- Integer, optional: Index of starting element (this element 
!                 is NOT tracked through). Default is 0.
!   ix2      -- Integer, optional: Index of ending element.
!                 Default is ring%n_ele_ring.
!
! Output:
!   beam(ix1+1:ix2) -- Beam_struct: Beam at end of elements.
!-

subroutine track_all_beam (ring, beam, ix1, ix2)

  implicit none

  type (ring_struct), intent(in) :: ring
  type (beam_struct) :: beam(0:)

  integer, optional, intent(in) :: ix1, ix2
  integer i, j, i1, i2

! Init

  i1 = 0
  if (present(ix1)) i1 = ix1
  i2 = ring%n_ele_ring
  if (present(ix2)) i2 = ix2

!

  do i = i1+1, i2
    call reallocate_beam (beam(i), size(beam(i1)%bunch), &
                size(beam(i1)%bunch(1)%slice), size(beam(i1)%bunch(1)%slice(1)%macro)) 
    do j = 1, size(beam(i)%bunch)
      call track1_bunch (beam(i-1)%bunch(j), ring%ele_(i), ring%param, beam(i)%bunch(j))
    enddo
  enddo

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track1_bunch (bunch_start, ele, param, bunch_end)
!
! Subroutine to track a bunch of macroparticles through an element.
!
! Modules needed:
!   use macroparticle_mod
!
! Input:
!   bunch_start -- Bunch_struct: Starting bunch position.
!   ele         -- Ele_struct: The element to track through.
!   param       -- Param_struct: General parameters.
!
! Output:
!   bunch_end -- Bunch_struct: Ending bunch position.
!-

Subroutine track1_bunch (bunch_start, ele, param, bunch_end)

  implicit none

  type (bunch_struct) bunch_start, bunch_end
  type (ele_struct) ele, rf_ele
  type (param_struct) param

  real(rp) charge
  integer i, j, ix

  logical sr_wake_on

! Charge and center

  bunch_end%charge   = bunch_start%charge
  bunch_end%s_center = bunch_start%s_center
  do i = 1, size(bunch_start%slice)
    bunch_end%slice(i)%charge = bunch_start%slice(i)%charge
  enddo

!------------------------------------------------
! Without wakefields just track through

  if (ele%key /= lcavity$ .or. .not. bmad_com%sr_wakes_on) then
    do i = 1, size(bunch_start%slice)
      do j = 1, size(bunch_start%slice(i)%macro)
        call track1_macroparticle (bunch_start%slice(i)%macro(j), &
                                      ele, param, bunch_end%slice(i)%macro(j))
      enddo
    enddo
    return
  endif

!------------------------------------------------
! For an lcavity without a wakefield file use the e_loss attribute 

  if (.not. associated(ele%wake%sr)) then
    call order_macroparticles_in_z (bunch_start)
    charge = param%charge
    param%charge = 0
    call sr_long_wake_calc (bunch_start, ele) ! calc %k_loss factors
    do i = 1, size(bunch_start%slice)
      do j = 1, size(bunch_start%slice(i)%macro)
        bmad_com%k_loss = bunch_start%slice(i)%macro(j)%k_loss
        call track1_macroparticle (bunch_start%slice(i)%macro(j), &
                                      ele, param, bunch_end%slice(i)%macro(j))
      enddo
    enddo
    param%charge = charge
    bmad_com%k_loss = 0
    return
  endif

!------------------------------------------------
! This calculation is for an lcavity with full wakefields.
! With wakefields put them in at the half way point.

! rf_ele is half the cavity

  rf_ele = ele
  rf_ele%value(l$) = ele%value(l$) / 2
  rf_ele%value(beam_energy$) = &
            (ele%value(energy_start$) + ele%value(beam_energy$)) / 2

  charge = param%charge
  param%charge = 0

  call sr_long_wake_calc (bunch_start, ele) ! calc %k_loss factors

! Track half way through. This includes the longitudinal wakefields

  do i = 1, size(bunch_start%slice)
    do j = 1, size(bunch_start%slice(i)%macro)
      bmad_com%k_loss = bunch_start%slice(i)%macro(j)%k_loss
      call track1_macroparticle (bunch_start%slice(i)%macro(j), &
                                    rf_ele, param, bunch_end%slice(i)%macro(j))
    enddo
  enddo

! Put in the transverse wakefields

  rf_ele%value(l$) = ele%value(l$)
  call track1_sr_trans_wake (bunch_end, rf_ele)

! Track the last half of the lcavity. This includes the longitudinal wakes.

  rf_ele%value(l$) = ele%value(l$) / 2
  rf_ele%value(energy_start$) = rf_ele%value(beam_energy$)
  rf_ele%value(beam_energy$) = ele%value(beam_energy$)

  do i = 1, size(bunch_start%slice)
    do j = 1, size(bunch_start%slice(i)%macro)
      bmad_com%k_loss = bunch_start%slice(i)%macro(j)%k_loss
      call track1_macroparticle (bunch_end%slice(i)%macro(j), &
                                    rf_ele, param, bunch_end%slice(i)%macro(j))
    enddo
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

  type (bunch_struct), target :: bunch
  type (ele_struct) ele
  type (macro_struct), pointer :: macro(:), macro2(:)

  real(rp) z_ave, charge, z0, sr02, ch_j, sig0, sig, dz
  real(rp) f, f1, f2, z, dE, E_rel, fact, dz_wake

  integer i, j, k, ix, i0, i1, ij, ik, iw, nm, n_wake

! Init

  call order_macroparticles_in_z (bunch)

  n_wake = ubound(ele%wake%sr, 1)
  dz_wake = ele%wake%sr(n_wake)%z / n_wake
  sr02 = ele%wake%sr(0)%long / 2

  do i = 1, size(bunch%slice)
    bunch%slice(i)%macro(:)%k_loss = 0
  enddo

! Loop over all slices

  do i = 1, size(bunch%slice)

    macro => bunch%slice(i)%macro
    nm = size(macro)
    charge = bunch%slice(i)%charge
    if (charge == 0) cycle

! Calculate wakefields within a slice.
! Easy calc is when all the particles have the same position

    i0 = macro(1)%iz
    i1 = macro(nm)%iz
  
    if (macro(i0)%r%vec(5) - macro(i1)%r%vec(5) < mp_com%sig_z_min/100) then
      z_ave = (macro(i0)%r%vec(5) + macro(i1)%r%vec(5)) / 2
      do j = 1, nm
        macro(j)%k_loss = macro(j)%k_loss + charge * sr02
      enddo

! If not all the same position then we need to look at all pairs

    else

      z_ave = 0

      do j = 1, nm
        ij = macro(j)%iz
        ch_j = abs(macro(ij)%charge)
        macro(ij)%k_loss = macro(ij)%k_loss + ch_j * sr02
        sig0 = macro(ij)%sig_z + mp_com%sig_z_min
        z0 = macro(ij)%r%vec(5)
        z_ave = z_ave + (z0 * ch_j) / charge

        do k = j+1, nm
          ik = macro(k)%iz
          sig = macro(ik)%sig_z + sig0
          dz = z0 - macro(ik)%r%vec(5)
          if (dz > sig) then
            f = 1
          else
            f = (1 + dz / sig) / 2
          endif
          macro(ik)%k_loss = macro(ik)%k_loss + f * ch_j * sr02
          macro(ij)%k_loss = macro(ij)%k_loss + &
                                   (1-f) * abs(macro(ik)%charge) * sr02
        enddo

      enddo

    endif

! now apply the wakefields to the other slices.

    do j = i+1, size(bunch%slice)
      macro2 => bunch%slice(j)%macro
      do k = 1, size(macro2)
        z = z_ave - macro2(k)%r%vec(5)
        iw = z / dz_wake
        f2 = z/dz_wake - iw
        f1 = 1 - f2
        macro2(k)%k_loss = macro2(k)%k_loss + charge * &
                  (ele%wake%sr(iw)%long * f1 + ele%wake%sr(iw+1)%long * f2)
      enddo
    enddo

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

  type (bunch_struct), target :: bunch
  type (ele_struct) ele
  type (macro_struct), pointer :: macro(:), macro2(:)

  real(rp) x_ave, y_ave, z_ave, charge, z0, ch_j
  real(rp) f1, f2, z, dE, E_rel, fact, dz_wake

  integer i, j, k, ix, i0, i1, ix0, ix1, iw, n_wake

! Init

  call order_macroparticles_in_z (bunch)

  n_wake = ubound(ele%wake%sr, 1)
  dz_wake = ele%wake%sr(n_wake)%z / n_wake

! Loop over all slices

  do i = 1, size(bunch%slice)

    macro => bunch%slice(i)%macro
    charge = bunch%slice(i)%charge
    if (charge == 0) cycle

! Calculate average position within a slice.
! It is assumed that the affect between particles in a slice is negligible.

    x_ave = 0
    y_ave = 0
    z_ave = 0

    do j = 1, size(macro)
      ch_j = abs(macro(j)%charge)
      z0 = macro(j)%r%vec(5)
      x_ave = x_ave + (macro(j)%charge * macro(j)%r%vec(1)) / charge
      y_ave = y_ave + (macro(j)%charge * macro(j)%r%vec(3)) / charge
      z_ave = z_ave + (macro(j)%charge * macro(j)%r%vec(5)) / charge
    enddo

! now apply the wakefields to the other slices.

    do j = i+1, size(bunch%slice)
      macro2 => bunch%slice(j)%macro
      do k = 1, size(macro2)
        z = z_ave - macro2(k)%r%vec(5)
        iw = z / dz_wake
        f2 = z/dz_wake - iw
        f1 = 1 - f2
        fact = (ele%wake%sr(iw)%trans*f1 + ele%wake%sr(iw+1)%trans*f2) * &
                              charge * ele%value(l$) / ele%value(beam_energy$)
        macro2(k)%r%vec(2) = macro2(k)%r%vec(2) + fact * x_ave 
        macro2(k)%r%vec(4) = macro2(k)%r%vec(4) + fact * y_ave
      enddo
    enddo

  enddo


end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track1_macroparticle (start, ele, param, end)
!
! Subroutine to track a macro-particle through an element.
! Note: the transfer ele%mat6 matrix is changed and not restored
! during traking.
!
! Modules needed:
!   use macroparticle_mod
!
! Input:
!   start  -- macro_struct: Starting coords.
!   ele    -- Ele_struct: Element to track through.
!   param  -- Param_struct: Global parameters.
!
! Output:
!   end    -- macro_struct: Ending coords.
!-

subroutine track1_macroparticle (start, ele, param, end)

  implicit none

  type (macro_struct) start
  type (macro_struct) end
  type (ele_struct) ele
  type (param_struct) param

  real(rp) l, l2, s(21), m4(4,4), s_mat4(4,4), s_mat6(6,6)

! transfer z-order index, charge, etc

  end = start

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

    if (end%sigma(s55$) /= 0) end%sig_z = sqrt(end%sigma(s55$))    

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine order_macroparticles_in_z (bunch)
!
! Subroutine to order in each slice the macro particles longitudinally 
! using the leading edge of the bunch: %vec(5) + %sig_z
!
! Modules needed:
!   use macroparticle_mod
!
! Input:
!   bunch     -- Bunch_struct: collection of macroparticles.
!     %slice(i)  -- i^th slice.
!       %macro(j)%r%vec(5) -- Longitudinal position of j^th macroparticle.
!       %macro(j)%sig_z    -- Longitudinal sigma.
!
! Output:
!   bunch     -- Bunch_struct: collection of macroparticles.
!     %slice(i)  -- i^th slice.
!       %macro(j)%iz -- Index of macroparticle ordered using %vec(5)+%sig_z.
!                   Order is from large z (head of slice) to small z.
!                   That is: slice(i)%macro(1)%iz is the index of the macro 
!                   particle at the head of the slice.
!-

Subroutine order_macroparticles_in_z (bunch)

  implicit none

  type (bunch_struct), target :: bunch
  type (macro_struct), pointer :: macro(:)
  integer i, k, j1, j2
  real(rp) z1, z2
  logical ordered

! Loop over all slices

  do k = 1, size(bunch%slice)

    macro => bunch%slice(k)%macro

! Init if needed

    if (any (macro(:)%iz == 0)) then
      do i = 1, size(macro)
        macro(i)%iz = i
      enddo
    endif

! Order is from large z (head of bunch) to small z.

    do 

      ordered = .true.

      do i = 1, size(macro)-1
        j1 = macro(i)%iz
        j2 = macro(i+1)%iz
        z1 = macro(j1)%r%vec(5) + macro(j1)%sig_z
        z2 = macro(j2)%r%vec(5) + macro(j2)%sig_z
        if (z1 < z2) then
          macro(i:i+1)%iz = macro(i+1:i:-1)%iz
          ordered = .false.
        endif
      enddo

      if (ordered) exit  ! goto next slice.

    enddo 

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
!   use macroparticle_mod
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
!   use macroparticle_mod
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
! Subroutine mp_to_canonical_coords (mp, energy0)
!
! Subroutine to convert macro-particle coords from 
!     (x, x', y, y', z, E)
! to
!     (x, px, y, py, z, pz)
!
! Modules needed:
!   use macroparticle_mod
!
! Input:
!   mp      -- macro_struct: Macro-particle with angular coords.
!   energy0 -- real(rp): Reference energy.
!
! Output:
!   mp      -- macro_struct: Macro-particle with momentum coords.
!-

subroutine mp_to_canonical_coords (mp, energy0)

  implicit none

  type (macro_struct), target :: mp

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
!   use macroparticle_mod
!
! Input:
!   mp      -- macro_struct: Macro-particle with momentum coords.
!   energy0 -- real(rp): Reference energy.
!
! Output:
!   mp      -- macro_struct: Macro-particle with angular coords.
!-

subroutine mp_to_angle_coords (mp, energy0)

  implicit none

  type (macro_struct), target :: mp

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

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine reallocate_beam (beam, n_bunch, n_slice, n_macro)
! 
! Subroutine to reallocate memory within a beam_struct.
!
! Modules needed:
!   use macroparticle_mod
!
! Input:
!   n_bunch -- Integer: Number of bunches.
!   n_slice -- Integer: Number of slices per bunch.
!   n_macro -- Integer: Number of macroparticles per slice.
!
! Output:
!   beam -- Beam_struct: Allocated beam_struct structure.
!-

subroutine reallocate_beam (beam, n_bunch, n_slice, n_macro)

  implicit none

  type (beam_struct) beam

  integer i, j
  integer n_bunch, n_slice, n_macro

  logical de_bunch, de_slice, de_macro

! Deallocate

  de_bunch = .false.
  de_slice = .false.
  de_macro = .false.

  if (allocated(beam%bunch)) then
    if (size(beam%bunch) /= n_bunch) then
      de_bunch = .true.
      de_slice = .true.
      de_macro = .true.
    endif
    if (size(beam%bunch(1)%slice) /= n_slice) then
      de_slice = .true.
      de_macro = .true.
    endif
    if (size(beam%bunch(1)%slice(1)%macro) /= n_macro) then
      de_macro= .true.
    endif

    do i = 1, size(beam%bunch)
      do j = 1, size(beam%bunch(i)%slice)
        if (de_macro) deallocate (beam%bunch(i)%slice(j)%macro)
      enddo
      if (de_slice) deallocate (beam%bunch(i)%slice)
    enddo
    if (de_bunch) deallocate (beam%bunch)

  endif

! Allocate

  if (.not. allocated (beam%bunch)) allocate (beam%bunch(n_bunch))
  do i = 1, n_bunch
    if (.not. allocated (beam%bunch(i)%slice)) allocate (beam%bunch(i)%slice(n_slice))
    do j = 1, n_slice
      if (.not. allocated (beam%bunch(i)%slice(j)%macro)) &
                    allocate (beam%bunch(i)%slice(j)%macro(n_macro))
    enddo
  enddo

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine init_macro_distribution (beam, init, canonical_out)
!
! Subroutine to initialize a macroparticle distribution.
! This routine uses the LIAR algorithm. See the Bmad manual for more details.
!
! Modules needed:
!   use macroparticle_mod
!
! Input:
!   init          -- macroparticle_init_struct: Structure holding the 
!                     parameters of the initial distribution. See the 
!                     definition of the structure for more details.
!   canonical_out -- Logical: If True then convert to canonical coords.
!
! Output:
!   beam -- beam_struct: Initialized Structure
!     %bunch(:) -- %bunch(1) is leading bunch.
!       %slice(:) -- %slice(1) is leading slice in a bunch. 
!-

subroutine init_macro_distribution (beam, init, canonical_out)

  implicit none

  type (beam_struct), target ::  beam
  type (macro_init_struct) init
  type (bunch_struct), pointer :: bunch

  real(rp) z_fudge, e_fudge, z0, dz, e0, de, ex, ey
  integer i, j, k
  logical canonical_out

! Reallocate if needed 

  call reallocate_beam (beam, init%n_bunch, init%n_slice, init%n_macro)

! Initalize distribution of 1st bunch
! z_fudge and e_fudge are used so that the total charge adds up to the 
!  input charge.

  z_fudge = 2 * gauss_int(init%sig_z_cut)
  e_fudge = 2 * gauss_int(init%sig_e_cut)

  bunch => beam%bunch(1)
  bunch%charge = init%n_part * e_charge

  ex = init%x%emit * m_electron / init%E_0
  ey = init%y%emit * m_electron / init%E_0

  do j = 1, init%n_slice
    z0 = (init%n_slice - 2*j + 1) * init%sig_z_cut / init%n_slice 
    dz = init%sig_z_cut / init%n_slice
    bunch%slice(j)%charge = bunch%charge * &
                            (gauss_int(z0+dz) - gauss_int(z0-dz)) / z_fudge

    do k = 1, init%n_macro
      e0 = (2*k - 1 - init%n_macro) * init%sig_e_cut / init%n_macro
      de = init%sig_e_cut / init%n_macro
      bunch%slice(j)%macro(k)%charge = bunch%slice(j)%charge * &
                           (gauss_int(e0+de) - gauss_int(e0-de)) / e_fudge
      bunch%slice(j)%macro(k)%r%vec = init%center
      bunch%slice(j)%macro(k)%r%vec(5) = init%sig_z * z0 + init%center(5)
      bunch%slice(j)%macro(k)%r%vec(6) = init%sig_e * e0 + init%center(6)
      bunch%slice(j)%macro(k)%sigma = 0
      bunch%slice(j)%macro(k)%sigma(s11$) =  ex * init%x%beta
      bunch%slice(j)%macro(k)%sigma(s12$) = -ex * init%x%alpha
      bunch%slice(j)%macro(k)%sigma(s22$) =  ex * (1+init%x%alpha**2) / &
                                                                 init%x%beta
      bunch%slice(j)%macro(k)%sigma(s33$) =  ey * init%y%beta
      bunch%slice(j)%macro(k)%sigma(s34$) = -ey * init%y%alpha
      bunch%slice(j)%macro(k)%sigma(s44$) =  ey * (1+init%y%alpha**2) / &
                                                                 init%y%beta
      bunch%slice(j)%macro(k)%iz = k
      if (canonical_out) &
              call mp_to_canonical_coords (bunch%slice(j)%macro(k), init%E_0)
    enddo
  enddo

! Initialize all bunches

  do i = 1, size(beam%bunch)
    beam%bunch(i)%s_center = (1-i) * init%ds_bunch
    beam%bunch(i) = beam%bunch(1)
  enddo

!-----------------------------------------------------
contains

!+
! Function gauss_int (x)
!
! Subroutine to return the integral of the normalized Gaussian function.
!   gauss_int[x] = Integral [0, x]: e^(-x^2/2) / sqrt(2pi)
! Thus
!   gauss_int[0]        = 0
!   gauss_int{Infinity] = 0.5
!   gauss_int[-x]       = -gauss_int[x]
!
! Input:
!   x    -- Real(rp): Upper limit of the integral.
!
! Output:
!   gauss_int -- Real(rp): Integral of the Normal Curve from 0 to x.
!-

function gauss_int(x) result (int)

  use nr

  real(rp) x, int
  int = erf(x/sqrt_2) / 2  ! Error function from Numerical Recipes.

end function

end subroutine

end module
