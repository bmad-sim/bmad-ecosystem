module macroparticle_mod

use bmad_struct
use bmad_interface
use wake_mod
use beam_def_struct

! Remember: If any of the macroparticle structures change you will need to modify:
!   mp_beam_equal_mp_beam, mp_slice_equal_mp_slice, and mp_bunch_equal_mp_bunch

type macro_struct
  type (coord_struct) r   ! Center of the macroparticle
  real(rp) sigma(21)      ! Sigma matrix.
  real(rp) :: sig_z = 0   ! longitudinal macroparticle length (m).
  real(rp) grad_loss_sr_wake         ! loss factor (V/m). scratch variable for tracking.
  real(rp) charge         ! charge in a macroparticle (Coul).
  logical :: lost = .false.  ! Has the particle been lost in tracking?
end type

type macro_slice_struct
  type (macro_struct), pointer :: macro(:) => null()
  real(rp) charge   ! total charge in a slice (Coul).
end type

type macro_bunch_struct
  type (macro_slice_struct), pointer :: slice(:) => null()
  real(rp) charge   ! total charge in a bunch (Coul).
  real(rp) z_center ! longitudinal center of bunch (m).
end type

type macro_beam_struct
  type (macro_bunch_struct), pointer :: bunch(:) => null()
end type

type macro_init_twiss_struct
  real(rp) norm_emit
end type

type macro_init_struct
  type (macro_init_twiss_struct) a, b
  real(rp) :: dPz_dz = 0  ! Correlation of Pz with long position.
  real(rp) :: center(6) = 0 ! Bench center offset relative to reference.
  real(rp) n_part      ! Number of particles per bunch.
  real(rp) ds_bunch    ! Distance between bunches.
  real(rp) sig_z       ! Z sigma in m.
  real(rp) sig_e       ! e_sigma in dE/E.
  real(rp) sig_e_cut   ! Energy cut in sigmas.
  real(rp) sig_z_cut   ! Z cut in sigmas.
  integer n_bunch      ! Number of bunches.
  integer n_slice      ! Number of slices per bunch.
  integer n_macro      ! Number of macroparticles per slice.
end type

type macroparticle_common_struct
  real(rp) ::  sig_z_min = 5e-6 ! min effective macroparticle length.
end type

type (macroparticle_common_struct), save :: mp_com

interface assignment (=)
  module procedure mp_slice_equal_mp_slice
  module procedure mp_bunch_equal_mp_bunch
  module procedure mp_beam_equal_mp_beam
end interface


contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track_macro_beam (lat, beam, ix1, ix2)
!
! Subroutine to track a beam of macroparticles from the end of
! lat%ele(ix1) Through to the end of lat%ele(ix2).
!
! Modules needed:
!   use macroparticle_mod
!
! Input:
!   lat   -- lat_struct: Lattice to track through.
!   beam   -- Macro_beam_struct: Beam at end of element ix1.
!   ix1    -- Integer, optional: Index of starting element (this element 
!               is NOT tracked through). Default is 0.
!   ix2    -- Integer, optional: Index of ending element.
!               Default is lat%n_ele_track.
!
! Output:
!   beam   -- Macro_beam_struct: Beam at end of element ix2.
!-

subroutine track_macro_beam (lat, beam, ix1, ix2)

  implicit none

  type (lat_struct) :: lat
  type (macro_beam_struct) :: beam

  integer, optional, intent(in) :: ix1, ix2
  integer i, i1, i2

! Init

  i1 = 0
  if (present(ix1)) i1 = ix1
  i2 = lat%n_ele_track
  if (present(ix2)) i2 = ix2

! Loop over all elements in the lattice

  do i = i1+1, i2
    call track1_macro_beam (beam, lat%ele(i), lat%param, beam)
  enddo

end subroutine track_macro_beam

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track1_macro_beam (beam_start, ele, param, beam_end)
!
! Subroutine to track a beam of macroparticles through a single element.
!
! Note: For the purposes of the wake calculation it is assumed that the
! bunches are ordered with %bunch(1) being the head bunch (largest s).
!
! Modules needed:
!   use macroparticle_mod
!
! Input:
!   beam_start  -- Macro_beam_struct: starting beam position
!   ele         -- Ele_struct: The element to track through.
!   param       -- lat_param_struct: General parameters.
!
! Output:
!   beam_end    -- Macro_beam_struct: ending beam position.
!-

subroutine track1_macro_beam (beam_start, ele, param, beam_end)

  implicit none

  type (macro_beam_struct) beam_start
  type (macro_beam_struct), target :: beam_end
  type (ele_struct) ele
  type (lat_param_struct) param

  integer i, n_mode

! zero the long-range wakes if they exist.

  if (associated(ele%rf%wake)) then
    ele%rf%wake%lr%norm_sin = 0; ele%rf%wake%lr%norm_cos = 0
    ele%rf%wake%lr%skew_sin = 0; ele%rf%wake%lr%skew_cos = 0
  endif

! loop over all bunches in a beam

  do i = 1, size(beam_start%bunch)
    call track1_macro_bunch (beam_start%bunch(i), ele, param, beam_end%bunch(i))
  enddo

end subroutine track1_macro_beam

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track1_macro_bunch (bunch_start, ele, param, bunch_end)
!
! Subroutine to track a bunch of macroparticles through an element.
!
! Modules needed:
!   use macroparticle_mod
!
! Input:
!   bunch_start -- Macro_bunch_struct: Starting bunch position.
!   ele         -- Ele_struct: The element to track through.
!   param       -- lat_param_struct: General parameters.
!
! Output:
!   bunch_end -- Bunch_struct: Ending bunch position.
!-

Subroutine track1_macro_bunch (bunch_start, ele, param, bunch_end)

  implicit none

  type (macro_bunch_struct) bunch_start, bunch_end
  type (ele_struct) ele
  type (ele_struct), save :: rf_ele
  type (lat_param_struct) param

  real(rp) charge
  integer i, j

! Charge and center

  bunch_end%z_center = bunch_start%z_center
  bunch_end%charge   = bunch_start%charge
  do i = 1, size(bunch_start%slice)
    bunch_end%slice(i)%charge = bunch_start%slice(i)%charge
  enddo

!------------------------------------------------
! Without wakefields just track through

  if (ele%key /= lcavity$ .or. &
            (.not. bmad_com%sr_wakes_on .and. .not. bmad_com%lr_wakes_on)) then
    do i = 1, size(bunch_start%slice)
      do j = 1, size(bunch_start%slice(i)%macro)
        call track1_macroparticle (bunch_start%slice(i)%macro(j), &
                                      ele, param, bunch_end%slice(i)%macro(j))
      enddo
    enddo
    call recalc_charge
    return
  endif

!------------------------------------------------
! This calculation is for an lcavity with wakefields.
! Put the wakefield kicks at the half way point.

! rf_ele is half the cavity

  rf_ele = ele
  rf_ele%value(l$) = ele%value(l$) / 2
  rf_ele%value(E_TOT$) = &
            (ele%value(E_TOT_START$) + ele%value(E_TOT$)) / 2
  rf_ele%value(p0c$) = &
            (ele%value(p0c_start$) + ele%value(p0c$)) / 2

! param%charge is used with e_loss$ when there is

  call order_macroparticles_in_z (bunch_start)
  call grad_loss_macro_sr_wake_calc (bunch_start, ele) 

! Track half way through. 
! This includes the short-range longitudinal wakefields

  do i = 1, size(bunch_start%slice)
    do j = 1, size(bunch_start%slice(i)%macro)
      call track1_macroparticle (bunch_start%slice(i)%macro(j), &
                                    rf_ele, param, bunch_end%slice(i)%macro(j))
    enddo
  enddo

! Put in the short-range transverse wakefields

  rf_ele%value(l$) = ele%value(l$)  ! restore the correct length for the moment
  call track1_macro_sr_trans_wake (bunch_end, rf_ele)

! Put in the long-range wakes

  call track1_macro_lr_wake (bunch_end, rf_ele)

! Track the last half of the lcavity. 
! This includes the short-range longitudinal wakes.

  rf_ele%value(l$)            = ele%value(l$) / 2
  rf_ele%value(E_TOT_START$) = rf_ele%value(E_TOT$)
  rf_ele%value(p0c_start$)    = rf_ele%value(p0c$)
  rf_ele%value(E_TOT$)  = ele%value(E_TOT$)
  rf_ele%value(p0c$)          = ele%value(p0c$)

  do i = 1, size(bunch_start%slice)
    do j = 1, size(bunch_start%slice(i)%macro)
      call track1_macroparticle (bunch_end%slice(i)%macro(j), &
                                    rf_ele, param, bunch_end%slice(i)%macro(j))
    enddo
  enddo

  call recalc_charge

!--------------------------------------
! This subroutine is needed when particles get lost.

contains

subroutine recalc_charge

  do i = 1, size(bunch_end%slice)
    bunch_end%slice(i)%charge = sum (bunch_end%slice(i)%macro(:)%charge, &
                            mask = .not. bunch_end%slice(i)%macro(:)%lost)
  enddo
  bunch_end%charge = sum (bunch_end%slice(:)%charge)

end subroutine

end subroutine track1_macro_bunch

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine grad_loss_macro_sr_wake_calc (bunch, ele)
!
! Subroutine to put in the longitudinal component of the
! short range wake fields. 
! This routine is not really meant  for general use.
!-

subroutine grad_loss_macro_sr_wake_calc (bunch, ele)

  implicit none

  type (macro_bunch_struct), target :: bunch
  type (ele_struct) ele
  type (macro_struct), pointer :: macro(:), macro2(:)

  real(rp) z_ave, charge, z0, sr02, ch_j, sig0, sig, dz
  real(rp) f, f1, f2, z, dz_wake

  integer i, j, k, iw, nm, n_wake
  character(24) :: r_name = 'grad_loss_sr_wake_calc'
  logical wake_here

! no wake file: just use e_loss$ factor

  wake_here = .true.
  if (.not. associated(ele%rf%wake)) wake_here = .false.
  if (wake_here) n_wake = size (ele%rf%wake%sr_table) - 1
    
  if (.not. wake_here .or. n_wake == -1) then 
    do i = 1, size(bunch%slice)
      bunch%slice(i)%macro(:)%grad_loss_sr_wake = &
                           ele%value(e_loss$) * bunch%charge / ele%value(l$)
    enddo
    return 
  endif

!

  do i = 1, size(bunch%slice)
    bunch%slice(i)%macro(:)%grad_loss_sr_wake = 0
  enddo

  dz_wake = ele%rf%wake%sr_table(n_wake)%z / n_wake
  sr02 = ele%rf%wake%sr_table(0)%long / 2

! Loop over all slices

  do i = 1, size(bunch%slice)

    macro => bunch%slice(i)%macro
    nm = size(macro)
    charge = bunch%slice(i)%charge
    if (charge == 0) cycle

! Calculate wakefields within a slice.
! Easy calc is when all the particles have the same position

    if (macro(1)%r%vec(5) - macro(nm)%r%vec(5) < mp_com%sig_z_min/100) then
      z_ave = (macro(1)%r%vec(5) + macro(nm)%r%vec(5)) / 2
      do j = 1, nm
        macro(j)%grad_loss_sr_wake = macro(j)%grad_loss_sr_wake + charge * sr02
      enddo

! If not all the same position then we need to look at all pairs

    else

      z_ave = 0

      do j = 1, nm
        if (macro(j)%lost) cycle
        ch_j = abs(macro(j)%charge)
        macro(j)%grad_loss_sr_wake = macro(j)%grad_loss_sr_wake + ch_j * sr02
        sig0 = macro(j)%sig_z + mp_com%sig_z_min
        z0 = macro(j)%r%vec(5)
        z_ave = z_ave + (z0 * ch_j) / charge

        do k = j+1, nm
          if (macro(k)%lost) cycle
          sig = macro(k)%sig_z + sig0
          dz = z0 - macro(k)%r%vec(5)
          if (dz > sig) then
            f = 1
          else
            f = (1 + dz / sig) / 2
          endif
          macro(k)%grad_loss_sr_wake = macro(k)%grad_loss_sr_wake + f * ch_j * sr02
          macro(j)%grad_loss_sr_wake = macro(j)%grad_loss_sr_wake + &
                                   (1-f) * abs(macro(k)%charge) * sr02
        enddo

      enddo

    endif

! Now apply the wakefields to the other slices.
! Use linear interpolation.
! If z is larger than the array size then use a linear extrapolation.

    do j = i+1, size(bunch%slice)
      macro2 => bunch%slice(j)%macro
      do k = 1, size(macro2)
        if (macro(k)%lost) cycle
        z = macro2(k)%r%vec(5) - z_ave
        iw = z / dz_wake
        iw = min (iw, n_wake-1)
        if (z > 0) then
          call out_io (s_abort$, r_name, &
             'MACROPARTICLE Z POSITION HAS SHIFTED ACROSS A SLICE BOUNDARY.')
          if (global_com%exit_on_error) call err_exit
        endif
        if (iw+1 > n_wake) then
          call out_io (s_error$, r_name, &
             'SHORT RANGE WAKE ARRAY FROM FILE: ' // ele%rf%wake%sr_file, &
             'DOES NOT COVER MACROPARTICLE LENGTH')
          cycle
        endif
        f2 = z/dz_wake - iw
        f1 = 1 - f2
        macro2(k)%grad_loss_sr_wake = macro2(k)%grad_loss_sr_wake + charge * &
                  (ele%rf%wake%sr_table(iw)%long * f1 + ele%rf%wake%sr_table(iw+1)%long * f2)
      enddo
    enddo

  enddo

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track1_macro_lr_wake (bunch, ele)
!
! Subroutine to put in the long-range wakes for macroparticle tracking.
! This routine is not really meant  for general use.
!-

subroutine track1_macro_lr_wake (bunch, ele)

  implicit none

  type (macro_bunch_struct), target :: bunch
  type (ele_struct) ele
  type (macro_struct), pointer :: macro

  integer n_mode, j, k

! Check to see if we need to do any calc

  if (.not. associated(ele%rf%wake)) return
  n_mode = size(ele%rf%wake%lr)
  if (n_mode == 0 .or. .not. bmad_com%lr_wakes_on) return  

! Give the macroparticles a kick

  do j = 1, size(bunch%slice)
    do k = 1, size(bunch%slice(j)%macro)
      macro => bunch%slice(j)%macro(k)
      if (macro%lost) cycle
      call lr_wake_apply_kick (ele, bunch%z_center, macro%r)
    enddo
  enddo

! Add the wakes left by this bunch to the existing wakes.

  do j = 1, size(bunch%slice)
    do k = 1, size(bunch%slice(j)%macro)
      macro => bunch%slice(j)%macro(k)
      if (macro%lost) cycle
      call lr_wake_add_to (ele, bunch%z_center, macro%r, macro%charge)
    enddo
  enddo


end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine track1_macro_sr_trans_wake (bunch, ele)
!
! Subroutine to put in the transverse component of the
! short range wake fields. 
! This routine is not really meant  for general use.
!-

subroutine track1_macro_sr_trans_wake (bunch, ele)

  implicit none

  type (macro_bunch_struct), target :: bunch
  type (ele_struct) ele
  type (macro_struct), pointer :: macro(:), macro2(:)

  real(rp) x_ave, y_ave, z_ave, charge
  real(rp) f1, f2, z, fact, dz_wake

  integer i, j, k, iw, n_wake

! Init

  if (.not. associated(ele%rf%wake)) return
  n_wake = size(ele%rf%wake%sr_table) - 1
  if (n_wake == -1 .or. .not. bmad_com%sr_wakes_on) return

  call order_macroparticles_in_z (bunch)

  dz_wake = ele%rf%wake%sr_table(n_wake)%z / n_wake

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
      if (macro(j)%lost) cycle
      x_ave = x_ave + (macro(j)%charge * macro(j)%r%vec(1)) / charge
      y_ave = y_ave + (macro(j)%charge * macro(j)%r%vec(3)) / charge
      z_ave = z_ave + (macro(j)%charge * macro(j)%r%vec(5)) / charge
    enddo

    x_ave = x_ave - ele%value(x_offset_tot$)
    y_ave = y_ave - ele%value(y_offset_tot$)

! Now apply the wakefields to the other slices.
! Use linear interpolation.
! If z is larger than the array size then use a linear extrapolation.

    do j = i+1, size(bunch%slice)
      macro2 => bunch%slice(j)%macro
      do k = 1, size(macro2)
        if (macro2(k)%lost) cycle
        z = macro2(k)%r%vec(5) - z_ave  ! distance from slice to particle
        iw = z / dz_wake                ! Index of wake array
        iw = min(n_wake-1, iw)    ! effectively do an extrapolation.
        f2 = z/dz_wake - iw
        f1 = 1 - f2
        fact = (ele%rf%wake%sr_table(iw)%trans*f1 + ele%rf%wake%sr_table(iw+1)%trans*f2) * &
                              charge * ele%value(l$) / ele%value(p0c$)
        macro2(k)%r%vec(2) = macro2(k)%r%vec(2) - fact * x_ave 
        macro2(k)%r%vec(4) = macro2(k)%r%vec(4) - fact * y_ave
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
!
! Modules needed:
!   use macroparticle_mod
!
! Input:
!   start  -- macro_struct: Starting coords.
!   ele    -- Ele_struct: Element to track through.
!   param  -- lat_param_struct: Global parameters.
!
! Output:
!   end    -- macro_struct: Ending coords.
!-

subroutine track1_macroparticle (start, ele, param, end)

  implicit none

  type (macro_struct) :: start
  type (macro_struct) :: end
  type (ele_struct) :: ele
  type (lat_param_struct), intent(inout) :: param

  real(rp) l, l2, s(21), m4(4,4), s_mat4(4,4), s_mat6(6,6)
  real(rp) mat6(6,6), vec0(6)

! transfer z-order index, charge, etc

  end = start

  if (start%lost) return
  if (ele%key == marker$ .or. ele%key == photon_branch$ .or. ele%key == branch$) return

! Init

  bmad_com%grad_loss_sr_wake = start%grad_loss_sr_wake

! Very simple cases

  select case (ele%key)
  case (drift$, ecollimator$, elseparator$, hkicker$, instrument$, &
                    kicker$, monitor$, rcollimator$)

    call track1 (start%r, ele, param, end%r)
    call check_lost

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

    bmad_com%grad_loss_sr_wake = 0
    return

  end select

! Simple case where longitudinal motion can be ignored.

  if (start%sigma(s55$) == 0 .and. start%sigma(s66$) == 0 .and. &
                                                  ele%key /= sbend$) then

    mat6 = ele%mat6; vec0 = ele%vec0  ! save
    call make_mat6 (ele, param, start%r, end%r)  
    call check_lost

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

    bmad_com%grad_loss_sr_wake = 0
    ele%mat6 = mat6; ele%vec0 = vec0
    return

  endif

! Full tracking. 

    mat6 = ele%mat6; vec0 = ele%vec0  ! save
    call make_mat6 (ele, param, start%r, end%r)  
    call check_lost
    call mp_sigma_to_mat (start%sigma, s_mat6)
    s_mat6 = matmul(ele%mat6, matmul(s_mat6, transpose(ele%mat6)))
    call mat_to_mp_sigma (s_mat6, end%sigma)

    bmad_com%grad_loss_sr_wake = 0
    ele%mat6 = mat6; ele%vec0 = vec0

! Sig_z calc. Because of roundoff sigma(s55$) can be negative.

    if (end%sigma(s55$) < 0) end%sigma(s55$) = 0
    if (end%sigma(s55$) /= 0) end%sig_z = sqrt(end%sigma(s55$))    

!--------------------------------------------
contains

subroutine check_lost

  end%lost = param%lost

  if (end%lost) then
    end%r%vec = 0
    end%sigma = 0 
    end%charge = 0
    return
  endif

end subroutine

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine order_macroparticles_in_z (bunch)
!
! Subroutine to order in each slice the macro particles longitudinally 
! The ordering uses the centroid of the macroparticles:
!   %vec(5) 
! If the position of a macroparticle in one slice starts overlaping 
! macroparticles in another slice, this routine will swap pairs of 
! macroparticles to prevent this.
!
! Modules needed:
!   use macroparticle_mod
!
! Input:
!   bunch     -- Macro_Bunch_struct: collection of macroparticles.
!     %slice(i)  -- i^th slice.
!       %macro(j)%r%vec(5) -- Longitudinal position of j^th macroparticle.
!       %macro(j)%sig_z    -- Longitudinal sigma.
!
! Output:
!   bunch     -- Macro_bunch_struct: collection of macroparticles.
!     %slice(i)  -- i^th slice.
!       %macro(j) -- Macroparticle ordered using %vec(5).
!                   Order is from large z (head of slice) to small z.
!                   That is: slice(i)%macro(1) is the macro  
!                   particle at the head of the slice.
!-

Subroutine order_macroparticles_in_z (bunch)

  implicit none

  type (macro_bunch_struct), target :: bunch
  type (macro_struct), pointer :: macro(:)
  type (macro_struct) temp
  integer i, k, nm
  real(rp) z1, z2
  logical ordered

! Loop over all slices
! Order is from large z (head of bunch) to small z.

  do

    ordered = .true.

    do k = 1, size(bunch%slice)

! order macroparticles within a slice

      macro => bunch%slice(k)%macro
      nm = size(macro)

      do i = 1, nm-1
        if (macro(i)%r%vec(5) < macro(i+1)%r%vec(5)) then
          macro(i:i+1) = macro(i+1:i:-1)
          ordered = .false.
        endif
      enddo

! Now make sure that slices do not overlap

      if (k < size(bunch%slice)) then
        z1 = macro(nm)%r%vec(5)
        z2 = bunch%slice(k+1)%macro(1)%r%vec(5)
        if (z1 < z2) then
          temp = macro(nm)
          macro(nm) = bunch%slice(k+1)%macro(1)
          bunch%slice(k+1)%macro(1) = temp
          bunch%slice(k)%charge = sum (bunch%slice(k)%macro(:)%charge, &
                            mask = .not. bunch%slice(k)%macro(:)%lost)
          bunch%slice(k+1)%charge = sum (bunch%slice(k+1)%macro(:)%charge, &
                            mask = .not. bunch%slice(k+1)%macro(:)%lost)
          ordered = .false.
        endif
      endif

    enddo

    if (ordered) exit  ! goto next slice.

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
    mat(1,:) = [s(s11$), s(s12$), s(s13$), s(s14$) ]
    mat(2,:) = [s(s12$), s(s22$), s(s23$), s(s24$) ]
    mat(3,:) = [s(s13$), s(s23$), s(s33$), s(s34$) ]
    mat(4,:) = [s(s14$), s(s24$), s(s34$), s(s44$) ]
  elseif (size(mat, 1) == 6 .and. size(mat,2) == 6) then
    mat(1,:) = [s(s11$), s(s12$), s(s13$), s(s14$), s(s15$), s(s16$)]
    mat(2,:) = [s(s12$), s(s22$), s(s23$), s(s24$), s(s25$), s(s26$)]
    mat(3,:) = [s(s13$), s(s23$), s(s33$), s(s34$), s(s35$), s(s36$)]
    mat(4,:) = [s(s14$), s(s24$), s(s34$), s(s44$), s(s45$), s(s46$)]
    mat(5,:) = [s(s15$), s(s25$), s(s35$), s(s45$), s(s55$), s(s56$)]
    mat(6,:) = [s(s16$), s(s26$), s(s36$), s(s46$), s(s56$), s(s66$)]
  else
    print *, 'ERROR IN MP_SIGMA_TO_MAT: MATRIX SIZE NOT 4 OR 6!'
    if (global_com%exit_on_error) call err_exit
  endif

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine mat_to_mp_sigma (mat, sigma)
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
    if (global_com%exit_on_error) call err_exit
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
! Note: the reverse routine is called:
!   mp_to_angle_coords (mp, energy0)
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
! Note: the reverse routine is called:
!   mp_to_canonical_coords (mp, energy0)
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
! Subroutine reallocate_macro_beam (beam, n_bunch, n_slice, n_macro)
! 
! Subroutine to reallocate memory within a macro_beam_struct.
!
! If n_bunch = 0 then all macro beam pointers will be deallocated.
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
!   beam -- Macro_beam_struct: Allocated beam_struct structure.
!-

subroutine reallocate_macro_beam (beam, n_bunch, n_slice, n_macro)

  implicit none

  type (macro_beam_struct) beam

  integer i, j
  integer n_bunch, n_slice, n_macro

  logical de_bunch, de_slice, de_macro

! Deallocate

  de_bunch = .false.
  de_slice = .false.
  de_macro = .false.

  if (associated(beam%bunch)) then
    if (n_bunch .eq. 0) then
      de_bunch = .true.
      de_slice = .true.
      de_macro = .true.
    else
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
    endif

    do i = 1, size(beam%bunch)
      do j = 1, size(beam%bunch(i)%slice)
        if (de_macro) deallocate (beam%bunch(i)%slice(j)%macro)
      enddo
      if (de_slice) deallocate (beam%bunch(i)%slice)
    enddo
    if (de_bunch) deallocate (beam%bunch)

  endif

  if (n_bunch .eq. 0) return
  
! Allocate

  if (.not. associated (beam%bunch)) allocate (beam%bunch(n_bunch))
  do i = 1, n_bunch
    if (.not. associated (beam%bunch(i)%slice)) allocate (beam%bunch(i)%slice(n_slice))
    do j = 1, n_slice
      if (.not. associated (beam%bunch(i)%slice(j)%macro)) &
                    allocate (beam%bunch(i)%slice(j)%macro(n_macro))
    enddo
  enddo

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine init_macro_distribution (beam, init, ele,
!                                             canonical_out, liar_gaussian)
!
! Subroutine to initialize a macroparticle distribution.
! This routine uses the LIAR algorithm. See the Bmad manual for more details.
!
! Modules needed:
!   use macroparticle_mod
!
! Input:
!   init          -- Macro_init_struct: Structure holding the 
!                      parameters of the initial distribution. See the 
!                      definition of the structure for more details.
!   ele           -- Ele _struct: Structure with dispersion and Twiss parameters.
!     %value(E_TOT$) -- Reference energy.
!     %a%beta              -- Beta
!     %a%alpha             -- Alpha
!     %a%eta               -- Dispersion
!     %c_mat(2,2)          -- Coupling matrix.
!   canonical_out -- Logical: If True then convert to canonical coords.
!   liar_gaussian -- Logical, optional: If true then use the liar algorithm
!                       for calculating macroparticle charge values.
!                       This is inaccurate at the few percent level.
!                       Use this for testing purposes. Default is False.
!
! Output:
!   beam -- Macro_beam_struct: Initialized Structure
!     %bunch(:) -- %bunch(1) is leading bunch.
!       %slice(:) -- %slice(1) is leading slice in a bunch. 
!-

subroutine init_macro_distribution (beam, init, ele, &
                                                canonical_out, liar_gaussian)

  implicit none

  type (macro_beam_struct), target ::  beam
  type (ele_struct) ele
  type (macro_init_struct), intent(in) :: init
  type (macro_bunch_struct), pointer :: bunch
  type (macro_struct), pointer :: macro

  real(rp) z_fudge, e_fudge, z_rel, dz, e_rel, ex, ey, dE_E, z, E0, del_e
  real(rp) mat4(4,4), v_mat(4,4), v_inv_mat(4,4), r(4), E_center

  integer i, j, k

  logical, intent(in) :: canonical_out
  logical, optional, intent(in) :: liar_gaussian

! Reallocate if needed 

  call reallocate_macro_beam (beam, init%n_bunch, init%n_slice, init%n_macro)

! Initalize distribution of 1st bunch
! z_fudge and e_fudge are used so that the total charge adds up to the 
!  input charge.

  z_fudge = 2 * gauss_int(init%sig_z_cut)
  e_fudge = 2 * gauss_int(init%sig_e_cut)
  E0 = ele%value(E_TOT$)
  E_center =  E0 * (1 + init%center(6))

  bunch => beam%bunch(1)
  bunch%charge = init%n_part * e_charge

  ex = init%a%norm_emit * m_electron / E_center
  ey = init%b%norm_emit * m_electron / E_center

  do j = 1, init%n_slice
    dz = init%sig_z_cut / init%n_slice
    z_rel = (init%n_slice - 2*j + 1) * dz
    bunch%slice(j)%charge = bunch%charge * &
                            (gauss_int(z_rel+dz) - gauss_int(z_rel-dz)) / z_fudge

    do k = 1, init%n_macro
      macro => bunch%slice(j)%macro(k)
      del_e = init%sig_e_cut / init%n_macro
      e_rel = (2*k - 1 - init%n_macro) * del_e
      macro%charge = bunch%slice(j)%charge * &
                           (gauss_int(e_rel+del_e) - gauss_int(e_rel-del_e)) / e_fudge
      z = init%center(5) + init%sig_z * z_rel 
      dE_E = init%center(6) + init%dPz_dz * z + e_rel * init%sig_e
      macro%r%vec = init%center
      macro%r%vec(5) = init%center(5) + init%sig_z * z_rel 
      macro%r%vec(6) = dE_E
      r = [ele%a%eta, ele%a%etap, ele%b%eta, ele%b%etap ]* dE_E
      ele%a%gamma = (1+ele%a%alpha**2) / ele%a%beta
      ele%b%gamma = (1+ele%b%alpha**2) / ele%b%beta
      macro%sigma = 0
      macro%sigma(s11$) =  ex * ele%a%beta  + r(1)**2
      macro%sigma(s12$) = -ex * ele%a%alpha + r(1) * r(2)
      macro%sigma(s22$) =  ex * ele%a%gamma + r(2)**2 
      macro%sigma(s33$) =  ey * ele%b%beta  + r(3)**2
      macro%sigma(s34$) = -ey * ele%b%alpha + r(3) * r(4)
      macro%sigma(s44$) =  ey * ele%b%gamma + r(4)**2
      if (any(ele%c_mat /= 0)) then
        call mp_sigma_to_mat (macro%sigma, mat4)
        call make_v_mats (ele, v_mat, v_inv_mat)
        mat4 = matmul (v_mat, matmul (mat4, transpose(v_mat)))
        call mat_to_mp_sigma (mat4, macro%sigma)
      endif
      macro%lost = .false.
      if (.not. canonical_out) &
              call mp_to_angle_coords (macro, E0)
    enddo
  enddo

! Initialize all bunches

  bunch%z_center = 0.0

  do i = 2, size(beam%bunch)
    beam%bunch(i) = beam%bunch(1)
    beam%bunch(i)%z_center = (1-i) * init%ds_bunch
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

  if (logic_option(.false., liar_gaussian)) then
    int = gauss_int_liar(x)
  else
    int = erf(x/sqrt_2) / 2  ! Error function from Numerical Recipes.
  endif

end function

!-----------------------------------------------------

function gauss_int_liar(x) result (g_int)

  implicit none

  real(rp) x
  real(rp) g_int
  integer i
  integer ipos
  real*4 zlo_gauss(46)
  real*4 int_gauss(46)
  real*4 sign
  real*4 slope

! define table of integrated gaussian distribution. integral goes
! from zero to zlo_gauss(i):

  zlo_gauss(01) = 0.0
  int_gauss(01) = 0.0000
  zlo_gauss(02) = 0.1
  int_gauss(02) = 0.0398
  zlo_gauss(03) = 0.2
  int_gauss(03) = 0.0793
  zlo_gauss(04) = 0.3
  int_gauss(04) = 0.1179
  zlo_gauss(05) = 0.4
  int_gauss(05) = 0.1554
  zlo_gauss(06) = 0.5
  int_gauss(06) = 0.1915
  zlo_gauss(07) = 0.6
  int_gauss(07) = 0.2257
  zlo_gauss(08) = 0.7
  int_gauss(08) = 0.2580
  zlo_gauss(09) = 0.8
  int_gauss(09) = 0.2881
  zlo_gauss(10) = 0.9
  int_gauss(10) = 0.3159
  zlo_gauss(11) = 1.0
  int_gauss(11) = 0.3413
  zlo_gauss(12) = 1.1
  int_gauss(12) = 0.3643
  zlo_gauss(13) = 1.2
  int_gauss(13) = 0.3849
  zlo_gauss(14) = 1.3
  int_gauss(14) = 0.4032
  zlo_gauss(15) = 1.4
  int_gauss(15) = 0.4192
  zlo_gauss(16) = 1.5
  int_gauss(16) = 0.4332
  zlo_gauss(17) = 1.6
  int_gauss(17) = 0.4452
  zlo_gauss(18) = 1.7
  int_gauss(18) = 0.4554
  zlo_gauss(19) = 1.8
  int_gauss(19) = 0.4641
  zlo_gauss(20) = 1.9
  int_gauss(20) = 0.4713
  zlo_gauss(21) = 2.0
  int_gauss(21) = 0.4772
  zlo_gauss(22) = 2.1
  int_gauss(22) = 0.4821
  zlo_gauss(23) = 2.2
  int_gauss(23) = 0.4860966
  zlo_gauss(24) = 2.3
  int_gauss(24) = 0.4892759
  zlo_gauss(25) = 2.4
  int_gauss(25) = 0.4918025
  zlo_gauss(26) = 2.5
  int_gauss(26) = 0.4937903
  zlo_gauss(27) = 2.6
  int_gauss(27) = 0.4953388
  zlo_gauss(28) = 2.7
  int_gauss(28) = 0.4965330
  zlo_gauss(29) = 2.8
  int_gauss(29) = 0.4974449
  zlo_gauss(30) = 2.9
  int_gauss(30) = 0.4981342
  zlo_gauss(31) = 3.0
  int_gauss(31) = 0.4986501
  zlo_gauss(32) = 3.1
  int_gauss(32) = 0.4990324
  zlo_gauss(33) = 3.2
  int_gauss(33) = 0.4993129
  zlo_gauss(34) = 3.3
  int_gauss(34) = 0.4995166
  zlo_gauss(35) = 3.4
  int_gauss(35) = 0.4996631
  zlo_gauss(36) = 3.5
  int_gauss(36) = 0.4997674
  zlo_gauss(37) = 3.6
  int_gauss(37) = 0.4998409
  zlo_gauss(38) = 3.7
  int_gauss(38) = 0.4998922
  zlo_gauss(39) = 3.8
  int_gauss(39) = 0.4999276
  zlo_gauss(40) = 3.9
  int_gauss(40) = 0.4999519
  zlo_gauss(41) = 4.0
  int_gauss(41) = 0.4999683
  zlo_gauss(42) = 4.1
  int_gauss(42) = 0.4999793
  zlo_gauss(43) = 4.2
  int_gauss(43) = 0.4999867
  zlo_gauss(44) = 4.3
  int_gauss(44) = 0.4999915
  zlo_gauss(45) = 4.4
  int_gauss(45) = 0.4999946
  zlo_gauss(46) = 4.5
  int_gauss(46) = 0.4999966


  ipos = 0
  do i = 1, 45
    if ( zlo_gauss(i)   .le. abs(x) .and. &
        zlo_gauss(i+1) .gt. abs(x) ) then
      ipos = i
      exit
    endif
  end do

  i     = ipos
  slope = (int_gauss(i+1) - int_gauss(i)) / &
         (zlo_gauss(i+1) - zlo_gauss(i))

  if (x == 0) then
    g_int = 0
  else
    if (x < 0) sign = -1.
    g_int = (int_gauss(i) + (abs(x) - zlo_gauss(i)) * slope )
    if (x < 0) g_int = -g_int
  endif
  
end function

end subroutine

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine mp_slice_equal_mp_slice (slice1, slice2)
!
! Subroutine to set one macroparticle slice equal to another taking care of
! pointers so that they don't all point to the same place.
!
! Note: This subroutine is called by the overloaded equal sign:
!		slice1 = slice2
!
! Input: 
!  slice2 -- macro_slice_struct: Input slice
!
! Output
!  slice1 -- macro_slice_struct: Output slice
!
!-

subroutine mp_slice_equal_mp_slice (slice1, slice2)

  implicit none

  type (macro_slice_struct), intent(inout) :: slice1
  type (macro_slice_struct), intent(in)    :: slice2

!

  if (size(slice1%macro) /= size(slice2%macro)) then
    deallocate(slice1%macro)
    allocate(slice1%macro(size(slice2%macro)))
  endif

  slice1%macro(:)  = slice2%macro(:)
  slice1%charge    = slice2%charge

end subroutine mp_slice_equal_mp_slice

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine mp_bunch_equal_mp_bunch (bunch1, bunch2)
!
! Subroutine to set one macroparticle bunch equal to another taking care of
! pointers so that they don't all point to the same place.
!
! Note: This subroutine is called by the overloaded equal sign:
!		bunch1 = bunch2
!
! Input: 
!  bunch2 -- macro_bunch_struct: Input bunch
!
! Output
!  bunch1 -- macro_bunch_struct: Output bunch
!
!-

subroutine mp_bunch_equal_mp_bunch (bunch1, bunch2)

  implicit none

  type (macro_bunch_struct), intent(inout) :: bunch1
  type (macro_bunch_struct), intent(in)    :: bunch2

  integer i

!

  if (size(bunch1%slice) /= size(bunch2%slice)) then
    do i = 1, size(bunch1%slice)
      deallocate (bunch1%slice(i)%macro)
    enddo
    deallocate (bunch1%slice)
    allocate (bunch1%slice(size(bunch2%slice)))
    do i = 1, size(bunch1%slice)
      allocate (bunch1%slice(i)%macro(size(bunch2%slice(i)%macro)))
    enddo  
  endif

  do i = 1, size(bunch1%slice)
    bunch1%slice(i) = bunch2%slice(i)
  enddo  

  bunch1%charge    = bunch2%charge
  bunch1%z_center  = bunch2%z_center

end subroutine mp_bunch_equal_mp_bunch

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine mp_beam_equal_mp_beam (beam1, beam2)
!
! Subroutine to set one macroparticle beam equal to another taking care of
! pointers so that they don't all point to the same place.
!
! Note: This subroutine is called by the overloaded equal sign:
!		beam1 = beam2
!
! Input: 
!  beam2 -- macro_beam_struct: Input beam
!
! Output
!  beam1 -- macro_beam_struct: Output beam
!
!-

subroutine mp_beam_equal_mp_beam (beam1, beam2)

  implicit none

  type (macro_beam_struct), intent(inout) :: beam1
  type (macro_beam_struct), intent(in)    :: beam2

  integer i, j, n_bun, n_slice, n_macro

!

  n_bun = size(beam2%bunch)

  if (size(beam1%bunch) /= size(beam2%bunch)) then
    do i = 1, size(beam1%bunch)
      do j = 1, size(beam1%bunch(i)%slice)
        deallocate (beam1%bunch(i)%slice(j)%macro)
       enddo
      deallocate (beam1%bunch(i)%slice)
    enddo
    deallocate (beam1%bunch)
    allocate (beam1%bunch(n_bun))
    do i = 1, n_bun
      n_slice = size(beam2%bunch(i)%slice)
      allocate (beam1%bunch(i)%slice(n_slice))
      do j = 1, n_slice
        n_macro = size(beam2%bunch(i)%slice(j)%macro)
        allocate (beam1%bunch(i)%slice(j)%macro(n_macro))
      enddo
    enddo
  endif

  do i = 1, n_bun
    beam1%bunch(i) = beam2%bunch(i)
  enddo

end subroutine mp_beam_equal_mp_beam

end module
