!+                                
! Subroutine ele_to_fibre (ele, fiber, param, integ_order, steps)
!
! Subroutine to convert a BMAD element to a PTC fibre element.
! This subroutine allocates fresh storage for the fibre so after calling
! this routine you need to deallocate at some point with:
!       call kill (fiber)
!
! Note: You need to call set_ptc before using this routine.
!
! Modules Needed:
!   use accelerator
!
! Input:
!   ele    -- Ele_struct: BMAD element.
!   param       -- param_struct: 
!     %energy     -- Beam energy (for wigglers).
!   integ_order -- Integer, optional: Order for the 
!                    sympletic integrator. Possibilities are: 2, 4, or 6
!                    Overrides ele%integration_order.
!                    default = 2 (if not set with set_ptc).
!   steps       -- Integer, optional: Number of integration steps.
!                    Overrides ele%num_steps.
!                    If ele%num_steps = 0 and steps is not present
!                    then the default is used. The default = 10 if not
!                    set with set_ptc. 
!
! Output:
!   fiber -- Fibre: PTC fibre element.
!+

subroutine ele_to_fibre (ele, fiber, param, integ_order, steps)

  use accelerator

  implicit none
 
  type (ele_struct) ele
  type (fibre) fiber
  type (el_list) el
  type (param_struct) :: param

  real(8) mis_rot(6)

  real(rdef) an0(0:n_pole_maxx), bn0(0:n_pole_maxx)
  real(rdef) cos_t, sin_t, len, hk, vk, x_off, y_off

  integer i, n, metd_temp, nstd_temp, key, n_term
  integer, optional :: integ_order, steps

  character name*16

!

  el = 0  ! init: subroutine el_0

  el%name    = ele%name
  el%vorname = ele%type

  el%l    = ele%value(l$)
  el%ld   = ele%value(l$)
  el%lc   = ele%value(l$)

  el%tilt = ele%value(tilt$)

!

  key = ele%key
  len = ele%value(l$)
  if (.not. ele%is_on) key = drift$

  select case (key)

  case (drift$) 
    el%kind = kind1

  case (quadrupole$) 
    el%kind = matrix_kick_matrix  ! kind7
    el%k(2) = ele%value(k1$)

  case (sbend$) 
    el%kind = matrix_kick_matrix  ! kind7
    el%b0 = ele%value(g_design$)
    el%lc = ele%value(l_chord$)
    el%t1 = ele%value(e1$)
    el%t2 = ele%value(e2$)

  case (sextupole$)
    el%kind = drift_kick_drift  ! kind2
    el%k(3) = ele%value(k2$) / 2

  case (octupole$)
    el%kind = drift_kick_drift ! kind2
    el%k(4) = ele%value(k3$) / 6

  case (solenoid$)
    el%kind = kind17    ! kind5 will be quicker but not exact
    el%bsol = ele%value(ks$)

  case (sol_quad$)
    el%kind = kind17    ! kind5 will be quicker but not exact
    el%bsol = ele%value(ks$)
    el%k(2) = ele%value(k1$)

  case (marker$)
    el%kind = kind0

  case (kicker$)
    el%kind = kind2

  case (rfcavity$)
    el%kind = kind4
    el%volt = ele%value(volt$)
    if (param%total_length == 0) then
      el%freq0 = 1
    else
      el%freq0 = c_light / param%total_length
    endif
    el%lag = ele%value(lag$)
    el%delta_e = 0

  case (elseparator$)
!    el%kind = drift_kick_drift  ! kind2
    el%kind = kind15
    if (ele%value(hkick$) == 0 .and. ele%value(vkick$) == 0) then
      el%tilt = 0
    else
      el%tilt = -atan2 (ele%value(hkick$), ele%value(vkick$))
    endif
    el%volt = (1e3 * param%energy / len) * &
                   sqrt(ele%value(hkick$)**2 + ele%value(vkick$)**2)
    call multipole_ele_to_ab (ele, param%particle, an0, bn0, .false.) 
    if (any(an0 /= 0) .or. any(bn0 /= 0)) then
      type *, 'ERROR IN ELE_TO_FIBRE: MULTIPOLES IN AN ELSEPARATOR NOT SUPPORTED IN A FIBRE'
      call err_exit
    endif

  case (ab_multipole$, multipole$)
    el%kind = kind3

  case (beambeam$)
    print *, 'ERROR IN ELE_TO_FIBRE: BEAMBEAM ELEMENT NOT YET IMPLEMENTED!'
    call err_exit

  case (wiggler$)
    el%kind = kinduser2    
    if (ele%sub_key == periodic_type$) then
      print *, 'ERROR IN ELE_TO_FIBRE: OLD STYLE WIGGLER: ', ele%name
      print *, '       CANNOT BE USED WITH TAYLOR.'
      call err_exit
    endif

  case default
    print *, 'ERROR IN ELE_TO_FIBRE: UNKNOWN ELEMENT KEY: ', &
                                                 key_name(ele%key)
    print *, '      FOR ELEMENT: ', ele%name
    call err_exit

  end select

! multipole components
! bmad an and bn are integrated fields. PTC uses just the field.

  if (ele%key /= elseparator$) then
    if (ele%value(hkick$) /= 0 .or. ele%value(vkick$) /= 0) then
      cos_t = cos(ele%value(tilt$))
      sin_t = sin(ele%value(tilt$))
      hk =  ele%value(hkick$) / len
      vk =  ele%value(hkick$) / len
      el%k(1)  = -hk * cos_t - vk * sin_t
      el%ks(1) = -hk * sin_t + vk * cos_t
    endif

    call multipole_ele_to_ab (ele, param%particle, an0, bn0, .false.)
    if (len /= 0) then
      an0 = an0 / len
      bn0 = bn0 / len
    endif

    n = min(n_pole_maxx+1, size(el%k))
    if (n-1 < n_pole_maxx) then
      if (any(an0(n:n_pole_maxx) /= 0) .or. any(bn0(n:n_pole_maxx) /= 0)) then
        print *, 'WARNING IN ELE_TO_FIBRE: MULTIPOLE NOT TRANSFERED TO FIBRE'
        print *, '        FOR: ', ele%name
      endif
    endif
 
    el%ks(1:n) = el%ks(1:n) + an0(0:n-1)
    el%k(1:n) = el%k(1:n) + bn0(0:n-1)

    if (key == sbend$) el%k(1) = el%k(1) + ele%value(g$)

    do n = nmax, 1, -1
      if (el%ks(n) /= 0 .or. el%k(n) /= 0) exit
    enddo
    el%nmul  = n
  endif

  if (el%kind == kind1 .and. el%nmul > 0) el%kind = matrix_kick_matrix  ! kind7

  if (el%kind == matrix_kick_matrix) el%nmul = max(2, el%nmul)

! METD and STEPS are global variables in PTC.

  if (present (integ_order)) then
    el%method = integ_order
  elseif (ele%integration_order /= 0) then
    el%method = ele%integration_order
  else
    el%method = METD
  endif

  if (present (steps)) then
    el%nst = steps
  elseif (ele%num_steps /= 0) then
    el%nst = ele%num_steps
  else
    el%nst = NSTD
  endif

  fiber = el

! wiggler

  if (key == wiggler$) then

    if (hyper_x$ /= hyperbolic_x$ .or. hyper_y$ /= hyperbolic_y$ .or. &
                                          hyper_xy$ /= hyperbolic_xy$) then
      print *, 'ERROR IN ELE_TO_FIBRE: WIGGLER FORM/TYPE MISMATCH!'
      print *, '     ', hyper_y$, hyper_xy$, hyper_x$
      print *, '     ', hyperbolic_y$, hyperbolic_xy$, hyperbolic_x$
      call err_exit
    endif

    n_term = size(ele%wig_term)
    call init_wig_pointers (fiber%mag%u2%w, n_term)   

    fiber%mag%u2%w%a(1:n_term) = 0.2997925 * &
            ele%value(polarity$) * ele%wig_term%coef / param%energy
    fiber%mag%u2%w%k(1,1:n_term)  = ele%wig_term%kx
    fiber%mag%u2%w%k(2,1:n_term)  = ele%wig_term%ky
    fiber%mag%u2%w%k(3,1:n_term)  = ele%wig_term%kz
    fiber%mag%u2%w%f(1:n_term)    = ele%wig_term%phi_z
    fiber%mag%u2%w%form(1:n_term) = ele%wig_term%type

    call copy (fiber%mag, fiber%magp)
  endif

!  misalignments.
! In PTC the reference point for the offsets is the beginning of the element.
! In BMAD the reference point is the center of the element..

  x_off = ele%value(x_offset$)-ele%value(l$)*ele%value(x_pitch$)
  y_off = ele%value(y_offset$)-ele%value(l$)*ele%value(y_pitch$)

  mis_rot = (/ x_off, y_off, 0.0_rdef, &
              -ele%value(y_pitch$), -ele%value(x_pitch$),  0.0_rdef /)

  fiber = mis_rot  ! call fibre_mis

end subroutine
