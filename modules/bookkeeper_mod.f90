#include "CESR_platform.inc"

module bookkeeper_mod

  use bmad_basic_mod
  use bmad_utils_mod
  use multipole_mod

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine control_bookkeeper (ring, ix_ele)
!
! Subroutine to transfer attibute information from lord to slave elements.
! If you want to do the bookkeeping for the entire ring you only need to
! call control_bookkeeper for the lord elements from ring%n_ele_ring+1 
! through ring%n_ele_max.
!
! Modules needed:
!   use bmad
!
! Input:
!   RING   -- Ring_struct: Ring to be used
!   IX_ELE -- Integer: Index of element whose attribute values have been
!               changed.
!-

subroutine control_bookkeeper (ring, ix_ele)

  implicit none

  type (ring_struct), target :: ring
  type (ele_struct), pointer :: ele

  integer ix_ele, i, j, ix, ix1, ix2
  integer ix_eles(300)

! attribute bookkeeping

  call attribute_bookkeeper (ring%ele_(ix_ele), ring%param)

! Make a list of elements to update.
! we do not need to update free elements of group lords.

  ix1 = 0   ! index for processed elements
  ix_eles(1) = ix_ele
  ix2 = 1   ! index for last element in list

  do
    ix1 = ix1 + 1
    ele => ring%ele_(ix_eles(ix1))
    do i = ele%ix1_slave, ele%ix2_slave
      ix = ring%control_(i)%ix_slave
      if (ele%key == group_lord$ .and. ring%ele_(ix)%key == free$) cycle
      if (ix == ix_eles(ix2)) cycle   ! do not use duplicates
      ix2 = ix2 + 1
      ix_eles(ix2) = ix
    enddo
    if (ix1 == ix2) exit
  enddo

! First: Makup lords and group slaves.
! If a super_lord has lords in turn then these lords must be overlay_lords. 
! therefore treat the super_lord as an overlay_slave.
! the same is true if an overlay lord has lords.

  do j = 1, ix2

    ix = ix_eles(j)
    ele => ring%ele_(ix)

    if (ele%control_type == group_lord$) then
      call makeup_group_slaves (ring, ix)

    elseif (ele%control_type == super_lord$ .and. ele%n_lord > 0) then
      call adjust_super_lord_s_position (ring, ix)
      call makeup_overlay_slave (ring, ix)

    elseif (ele%control_type == overlay_lord$ .and. ele%n_lord > 0) then
      call makeup_overlay_slave (ring, ix)

    endif

  enddo

! Second: Makeup all slaves but group slaves.

  do j = 1, ix2

    ix = ix_eles(j)
    ele => ring%ele_(ix)       

    if (ele%control_type == super_slave$) then
      call makeup_super_slave (ring, ix)

    elseif (ele%control_type == overlay_slave$) then
      call makeup_overlay_slave (ring, ix)

    endif

  enddo

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine adjust_super_lord_s_position (ring, ix_lord)
!
! Subroutine to adjust the positions of the slaves of a super_lord due
! to changes in the lord's s_offset.
!-

Subroutine adjust_super_lord_s_position (ring, ix_lord)

  implicit none

  type (ring_struct), target :: ring
  type (ele_struct), pointer :: lord, slave 

  integer ix_lord, ix

  real(rp) s_start, s_start2, s_end

!

  lord => ring%ele_(ix_lord)

  if (lord%control_type /= super_lord$) then
    print *, 'ERROR IN ADJUST_SUPER_LORD_S_POSITION: ELEMENT IS NOT A LORD!'
    call err_exit
  endif

! If a super lord is moved then we just need to adjust the start and end edges.
! Adjust end position.

  s_end = lord%s + lord%value(s_offset$)
  ix = ring%control_(lord%ix2_slave)%ix_slave
  slave => ring%ele_(ix)
  s_start = slave%s - slave%value(l$)
  slave%value(l$) = s_end - s_start
  slave%s = s_end

! Adjust start position

  s_start = s_end - lord%value(l$)
  if (s_start < 0) s_start = s_start + ring%param%total_length
  ix = ring%control_(lord%ix1_slave)%ix_slave
  slave => ring%ele_(ix)
  s_start2 = slave%s - slave%value(l$)
  if (s_start < 0) s_start = s_start + ring%param%total_length 
  slave%value(l$) = slave%value(l$) + s_start2 - s_start

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine makeup_group_slaves (ring, ix_slave)
!
! Subroutine to calculate the attributes of group slave elements
!-

Subroutine makeup_group_slaves (ring, ix_lord)   

  implicit none

  type (ring_struct), target :: ring
  type (ele_struct), pointer :: lord, slave

  real(rp) delta, coef

  integer ix_lord, ix, iv, ict, i

  logical moved

!

  lord => ring%ele_(ix_lord)

  delta = lord%value(command$) - lord%value(old_command$)    ! change
  lord%value(old_command$) = lord%value(command$) ! save old

  moved = .false.   ! have we longitudinally moved an element?

  do i = lord%ix1_slave, lord%ix2_slave

    ix = ring%control_(i)%ix_slave
    iv = ring%control_(i)%ix_attrib
    ict = ring%ele_(ix)%control_type
    slave => ring%ele_(ix)

    if (iv == l$) then
      moved = .true.
      if (ict /= free$ .and. ict /= super_slave$) then
        print *, 'ERROR IN CONTROL_BOOKKEEPER: A GROUP: ', lord%name
        print *, '      CONTROLS THE LENGTH OF A LORD ELEMENT: ', slave%name
        call err_exit
      endif
    endif
    coef = ring%control_(i)%coef
    slave%value(iv) = slave%value(iv) + delta * coef
  enddo

  if (moved) call s_calc (ring)       ! recompute s distances


end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine makeup_super_slave (ring, ix_slave)
!
! Subroutine to calcualte the attributes of overlay slave elements
!-
         
subroutine makeup_super_slave (ring, ix_slave)

  implicit none

  type (ring_struct), target :: ring
  type (ele_struct), pointer :: lord, slave
  type (ele_struct) :: sol_quad

  integer i, ix_con, j, ix, ix_slave

  real(rp) tilt, k1_x, k1_y, x_kick, y_kick, ks, k1, coef
  real(rp) x_o, y_o, x_p, y_p, s_slave, s_del
  real(rp) sin_2, cos_2, x_off, y_off, a(0:n_pole_maxx), b(0:n_pole_maxx)
  real(rp) knl(0:n_pole_maxx), t(0:n_pole_maxx), value(n_attrib_maxx)
  real(rp) a_tot(0:n_pole_maxx), b_tot(0:n_pole_maxx)
  real(rp) sum_1, sum_2, sum_3, sum_4, ks_sum, ks_xp_sum, ks_xo_sum
  real(rp) ks_yp_sum, ks_yo_sum, l_slave, r_off(4)
  real(rp) t_1(4), t_2(4), T_end(4,4), mat4(4,4), mat4_inv(4,4), beta(4)
  real(rp) T_tot(4,4), x_o_sol, x_p_sol, y_o_sol, y_p_sol

  logical, save :: init_needed = .true.

! init

  if (init_needed) then
    call init_ele (sol_quad)
    sol_quad%key = sol_quad$
    init_needed = .false.
  endif

! Super_slave:

  slave => ring%ele_(ix_slave)
                     
  if (slave%control_type /= super_slave$) then
    print *, 'ERROR IN MAKEUP_SUPER_SLAVE: ELEMENT IS NOT AN SUPER SLAVE:'
    print *, '      ', slave%name
    call err_exit
  endif

! If this slave is the last slave for some lord (so that then longitudinal 
! end of the slave matches the end of the lord) then the limits of the lord
! are transfered to the slave.

!-----------------------------------------------------------------------
! 1 super_lord for this super_slave: just transfer attributes except length

  if (slave%n_lord == 1) then

    ix_con = ring%ic_(slave%ic1_lord)  
    ix = ring%control_(ix_con)%ix_lord
    lord => ring%ele_(ix)
    coef = ring%control_(ix_con)%coef  ! = len_slave / len_lord

    value = lord%value
    value(check_sum$) = slave%value(check_sum$) ! do not change the check_sum
    value(l$) = slave%value(l$)                 ! do not change slave length
    value(hkick$) = lord%value(hkick$) * coef
    value(vkick$) = lord%value(vkick$) * coef
    slave%num_steps = lord%num_steps * coef + 1
    if (slave%key == rfcavity$) value(volt$) = lord%value(volt$) * coef
    if (ix_con /= lord%ix2_slave) then   ! if not at end of the lord domain
      value(x_limit$) = 0
      value(y_limit$) = 0
    endif

! s_del is the distance between lord and slave centers

    if (value(x_pitch$) /= 0 .or. value(y_pitch$) /= 0) then
      s_del = (slave%s - slave%value(l$)/2) - &
                  (lord%s + lord%value(s_offset$) - lord%value(l$)/2)
      s_del = modulo2 (s_del, ring%param%total_length/2)
      value(x_offset$) = value(x_offset$) + s_del * value(x_pitch$)
      value(y_offset$) = value(y_offset$) + s_del * value(y_pitch$)
    endif
      
    slave%value = value
    slave%is_on = lord%is_on

! if a wiggler: 
! must keep track of where we are in terms of the unsplit wiggler.
! This is for anything which does not try to make a homogeneous approximation.
! l_original is the length of the unsplit original wiggler.
! l_start is the starting point with respect to the original wiggler.
! l_end is the ending point with respect to the original wiggler.

    if (slave%key == wiggler$) then
      slave%value(n_pole$) = lord%value(n_pole$) * coef
      slave%value(l_original$) = lord%value(l$)
      slave%value(l_start$)    = (slave%s - slave%value(l$)) - &
                                                   (lord%s - lord%value(l$))
      slave%value(l_end$)      = slave%value(l_start$) + slave%value(l$)

      if (associated(lord%wig_term)) then
        if (.not. associated (slave%wig_term) .or. &
                size(slave%wig_term) /= size(lord%wig_term)) then
          deallocate (slave%wig_term)
          allocate (slave%wig_term(size(lord%wig_term)))
        endif
        do i = 1, size(lord%wig_term)
          slave%wig_term(i)%coef = lord%wig_term(i)%coef
          slave%wig_term(i)%kx = lord%wig_term(i)%kx
          slave%wig_term(i)%ky = lord%wig_term(i)%ky
          slave%wig_term(i)%kz = lord%wig_term(i)%kz
          slave%wig_term(i)%phi_z = lord%wig_term(i)%phi_z + &
                               lord%wig_term(i)%kz * slave%value(l_start$)
        enddo
      endif

    endif

! if an sbend:
!     1) renormalize the angles
!     2) zero the face angles next to the split

    if (slave%key == sbend$) then
      if (ix_con == lord%ix1_slave) then   ! first slave bend
        slave%value(e2$) = 0
      elseif (ix_con == lord%ix2_slave) then 
        slave%value(e1$) = 0
      else
        slave%value(e1$) = 0
        slave%value(e2$) = 0
      endif
    endif                       

    return

  endif

!-----------------------------------------------------------------------
! Multiple super_lords for this super_slave: 
! must be a solenoid/quadrupole combo.
! combine the lord elements.
                                           
  k1_x = 0
  k1_y = 0
  x_kick = 0
  y_kick = 0
  a_tot = 0
  b_tot = 0
  sum_1 = 0
  sum_2 = 0
  sum_3 = 0
  sum_4 = 0
  ks_sum = 0
  ks_xp_sum = 0
  ks_xo_sum = 0
  ks_yp_sum = 0
  ks_yo_sum = 0

  value = 0
  value(l$) = slave%value(l$)
  value(check_sum$) = slave%value(check_sum$) ! do not change the check_sum

  s_slave = slave%s - value(l$)/2  ! center of slave

! sum over all lords...

  do j = slave%ic1_lord, slave%ic2_lord

    ix_con = ring%ic_(j)
    ix = ring%control_(ix_con)%ix_lord
    coef = ring%control_(ix_con)%coef
    lord => ring%ele_(ix)

    if (lord%control_type /= super_lord$) then
      print *, 'ERROR IN MAKEUP_SUPER_SLAVE: SUPER_SLAVE HAS A'
      print *, '      CONTROL ELEMENT THAT IS NOT SUPER_LORD'
      print *, '      SLAVE: ', slave%name, ix_slave
      print *, '      LORD:  ', lord%name, ix
      call err_exit
    endif

    if (ix_con == lord%ix2_slave) then      ! longitudinal ends match
      value(x_limit$) = lord%value(x_limit$)
      value(y_limit$) = lord%value(y_limit$)
    endif

    if (lord%is_on) then

      x_p = lord%value(x_pitch$);  x_o = lord%value(x_offset$)
      y_p = lord%value(y_pitch$);  y_o = lord%value(y_offset$)

      s_del = s_slave - (lord%s + lord%value(s_offset$) - lord%value(l$)/2)
      s_del = modulo2 (s_del, ring%param%total_length/2)

      ks = lord%value(ks$)

      ks_sum = ks_sum + ks

      ks_xp_sum = ks_xp_sum + ks * x_p
      ks_yp_sum = ks_yp_sum + ks * y_p

      ks_xo_sum = ks_xo_sum + ks * (x_o + x_p * s_del)
      ks_yo_sum = ks_yo_sum + ks * (y_o + y_p * s_del)

      x_kick = x_kick + lord%value(hkick$) * coef
      y_kick = y_kick + lord%value(vkick$) * coef

      tilt = lord%value(tilt$)
      cos_2 = lord%value(k1$) * cos(2 * tilt)
      sin_2 = lord%value(k1$) * sin(2 * tilt)

      k1_x = k1_x + cos_2
      k1_y = k1_y + sin_2

      sum_1 = sum_1 + cos_2 * x_p + sin_2 * y_p
      sum_2 = sum_2 + sin_2 * x_p - cos_2 * y_p

      sum_3 = sum_3 + cos_2 * (x_o + x_p * s_del) + sin_2 * (y_o + y_p * s_del)
      sum_4 = sum_4 + sin_2 * (x_o + x_p * s_del) - cos_2 * (y_o + y_p * s_del)

      if (associated(lord%a)) then
        call multipole_ele_to_kt (lord, +1, knl, t, .true.)
        call multipole_kt_to_ab (knl/lord%value(l$), t, a, b)
        a_tot = a_tot + a
        b_tot = b_tot + b
      endif

    endif

  enddo

! stuff sums into slave element

  ks = ks_sum
  value(ks$) = ks
  value(hkick$) = x_kick
  value(vkick$) = y_kick

  slave%value = value

  if (k1_x == 0 .and. k1_y == 0 .and. ks == 0) return

  if (ks /= 0) then
    x_o_sol = ks_xo_sum / ks
    x_p_sol = ks_xp_sum / ks
    y_o_sol = ks_yo_sum / ks
    y_p_sol = ks_yp_sum / ks
  endif

  if (k1_x == 0 .and. k1_y == 0) then  ! pure solenoid
    slave%value(k1$) = 0
    slave%value(tilt$) = 0
    deallocate (slave%a, slave%b, stat = ix)
    slave%value(x_offset$) = x_o_sol
    slave%value(y_offset$) = y_o_sol
    slave%value(x_pitch$)  = x_p_sol
    slave%value(y_pitch$)  = y_p_sol
    return
  endif   

! here if have quadrupole component

  k1 = sqrt(k1_x**2 + k1_y**2)
  tilt = atan2(k1_y, k1_x) / 2

  if (tilt > pi/4) then
    k1 = -k1
    tilt = tilt - pi/2
  elseif (tilt < -pi/4) then
    k1 = -k1
    tilt = tilt + pi/2
  endif

  slave%value(k1$) = k1
  slave%value(tilt$) = tilt

  cos_2 = k1_x / (k1_x**2 + k1_y**2)
  sin_2 = k1_y / (k1_x**2 + k1_y**2)

  slave%value(x_pitch$)  = cos_2 * sum_1 + sin_2 * sum_2
  slave%value(y_pitch$)  = sin_2 * sum_1 - cos_2 * sum_2
  slave%value(x_offset$) = cos_2 * sum_3 + sin_2 * sum_4
  slave%value(y_offset$) = sin_2 * sum_3 - cos_2 * sum_4

  deallocate (slave%a, slave%b, stat = ix)
  if (any(a_tot /= 0) .or. any(b_tot /= 0)) then
    allocate (slave%a(0:n_pole_maxx), slave%b(0:n_pole_maxx))
    call multipole_ab_to_kt(a_tot, b_tot, knl, t)
    call multipole_kt_to_ab(knl/k1, t-tilt, a, b)
    slave%a = a
    slave%b = b
    slave%value(radius$) = 1
  endif

! if ks /= 0 then we have to recalculate the offsets and pitches.

  if (ks == 0) return

  x_p = slave%value(x_pitch$) - x_p_sol; x_o = slave%value(x_offset$) - x_o_sol
  y_p = slave%value(y_pitch$) - y_p_sol; y_o = slave%value(y_offset$) - y_o_sol

  if (x_p == 0 .and. x_o == 0 .and. y_p == 0 .and. y_o == 0) return

  t_2 = (/ x_o, x_p, y_o, y_p /)
  call tilt_coords (tilt, t_2, .true.)  ! set

  l_slave = slave%value(l$)

  t_1 = (/ t_2(2), 0.0_rp, t_2(4), 0.0_rp /)
  t_2(1) = t_2(1) + ks * t_2(4) / k1 
  t_2(3) = t_2(3) + ks * t_2(2) / k1
             
  call mat_make_unit (T_end)
  T_end(4,1) =  ks / 2
  T_end(2,3) = -ks / 2

  sol_quad%value(ks$) = ks
  sol_quad%value(k1$) = k1
  sol_quad%value(l$)  = l_slave
  call make_mat6 (sol_quad, ring%param)
  T_tot = sol_quad%mat6(1:4,1:4)

  r_off = matmul (T_end, l_slave * t_1 / 2 - t_2) 
  r_off = matmul (T_tot, r_off) + matmul (T_end, l_slave * t_1 / 2 + t_2)

  call mat_make_unit (mat4)
  mat4(:,2) = mat4(:,2) + l_slave * T_tot(:,1) / 2
  mat4(:,4) = mat4(:,4) + l_slave * T_tot(:,3) / 2
  mat4(1,2) = mat4(1,2) + l_slave / 2
  mat4(3,4) = mat4(3,4) + l_slave / 2
  mat4 = mat4 - T_tot

  call mat_inverse (mat4, mat4_inv)
  beta = matmul (mat4_inv, r_off)

  call tilt_coords (tilt, beta, .false.)  ! unset

  slave%value(x_offset$) = beta(1) + x_o_sol
  slave%value(x_pitch$)  = beta(2) + x_p_sol
  slave%value(y_offset$) = beta(3) + y_o_sol
  slave%value(y_pitch$)  = beta(4) + y_p_sol

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

subroutine makeup_overlay_slave (ring, ix_ele)

  implicit none

  type (ring_struct), target :: ring
  type (ele_struct), pointer :: ele

  real(rp) value(n_attrib_maxx), coef
  integer i, j, ix, iv, ix_ele, icom, ct
  logical used(n_attrib_maxx)

!
                               
  ele => ring%ele_(ix_ele)
  ct = ele%control_type

  if (ct /= super_lord$ .and. ct /= overlay_slave$ .and. &
                                               ct /= overlay_lord$) then
    print *, 'ERROR IN MAKEUP_OVERLAY_SLAVE: ELEMENT IS NOT OF PROPER TYPE.'
    print *, '      RING INDEX:', ix_ele
    call type_ele (ele, .true., 0, .false., 0, .true., ring)
    call err_exit
  endif

  value = 0
  used = .false.

  do i = ele%ic1_lord, ele%ic2_lord
    j = ring%ic_(i)
    ix = ring%control_(j)%ix_lord
    if (ring%ele_(ix)%control_type /= overlay_lord$) then
      print *, 'ERROR IN MAKEUP_OVERLAY_SLAVE:',  &
                          ' THE LORD IS NOT AN OVERLAY_LORD', ix_ele
      call type_ele (ele, .true., 0, .false., 0, .true., ring)
      call err_exit
    endif     
    coef = ring%control_(j)%coef
    iv = ring%control_(j)%ix_attrib
    icom = ring%ele_(ix)%ix_value
    value(iv) = value(iv) + ring%ele_(ix)%value(icom)*coef
    used(iv) = .true.
  enddo

  do i = 1, n_attrib_maxx
    if (used(i)) ele%value(i) = value(i)
  enddo
                                            
end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!+
! Subroutine attribute_bookkeeper (ele, param)
!
! Subroutine to recalculate the dependent attributes of an element.
! If the attributes have changed then any Taylor Maps will be killed.
!
!   BEAMBEAM:   bbi_const$ = param%n_part * m_electron * charge$ * r_e /
!                           (2 * pi * param%beam_energy * (sig_x$ + sig_y$)
!
!   RFCAVITY:   rf_wavelength$ = param%total_length / harmon$
!
!   SBEND:      angle$   = L$ * G$
!               l_chord$ = 2 * sin(Angle$/2) / G$
!               rho$     = 1 / G$
!
!   WIGGLER:    k1$  = -0.5 * (c_light * b_max$ / param%beam_energy)**2
!               rho$ = param%beam_energy / (c_light * b_max$)
!
!   LCAVITY:    e_loss$  = ele%wake%sr(0)%long/2 if ele%wake%sr exists
!               delta_e$ = gradient$ * L$ 
!
! Modules needed:
!   use bmad
!
! Input:
!   ele   -- Ele_struct: Element with attributes 
!   param -- Param_struct: 
!
! Output:
!   ele  -- Ele_struct: Element with self-consistant attributes.
!-

subroutine attribute_bookkeeper (ele, param)

  implicit none

  type (ele_struct) ele
  type (param_struct) param

  real(rp) r, factor, check_sum
  real(rp), save :: old_energy = 0, p
  
! field_master

  if (ele%field_master) then

    if (ele%value(beam_energy$) == 0) then
      factor = 0
    else
      if (old_energy /= ele%value(beam_energy$)) &
           call energy_to_kinetic (ele%value(beam_energy$), &
                                                    param%particle, p0c = p)
      factor = c_light / p
      old_energy = ele%value(beam_energy$)
    endif

    select case (ele%key)
    case (quadrupole$)
      ele%value(k1$) = factor * ele%value(B_gradient$)
    case (sextupole$)
      ele%value(k2$) = factor * ele%value(B_gradient$)
    case (octupole$)
      ele%value(k3$) = factor * ele%value(B_gradient$)
    case (solenoid$)
      ele%value(ks$) = factor * ele%value(B_field$)
    case (sbend$)
      ele%value(g$) = factor * ele%value(B_field$)
    case default
      print *, 'ERROR IN ATTRIBUTE_BOOKKEEPER: ', &
                      '"FIELD_MASTER" NOT IMPLEMENTED FOR: ', trim(ele%name)
      call err_exit
    end select

  else

    if (ele%value(beam_energy$) == 0) then
      factor = 0
    else
      if (old_energy /= ele%value(beam_energy$)) &
           call energy_to_kinetic (ele%value(beam_energy$), &
                                                  param%particle, p0c = p)
      factor = p / c_light
      old_energy = ele%value(beam_energy$)
    endif

    select case (ele%key)
    case (quadrupole$)
       ele%value(B_gradient$) = factor * ele%value(k1$)
    case (sextupole$)
       ele%value(B_gradient$) = factor * ele%value(k2$)
    case (octupole$)
       ele%value(B_gradient$) = factor * ele%value(k3$)
    case (solenoid$)
       ele%value(B_field$) = factor * ele%value(ks$)
    case (sbend$)
       ele%value(B_field$) = factor * ele%value(g$)
    end select

  endif

! Dependent attribute bookkeeping.

  select case (ele%key)

! Bends

  case (sbend$)
    ele%value(angle$) = ele%value(l$) * ele%value(g$)
    if (ele%value(l$) == 0 .or. ele%value(g$) == 0) then
      ele%value(l_chord$) = 0
    else
      ele%value(l_chord$) = 2 * sin(ele%value(angle$)/2) / ele%value(g$)
    endif
    if (ele%value(g$) == 0) then
      ele%value(rho$) = 0
    else
      ele%value(rho$) = 1 / ele%value(g$)
    endif

! Lcavity

  case (lcavity$)
    if (associated (ele%wake%sr)) &
                          ele%value(e_loss$) = ele%wake%sr(0)%long / 2
    ele%value(delta_e$) = ele%value(gradient$) * ele%value(L$) 
    

! RFcavity

  case (rfcavity$)
    if (ele%value(harmon$) /= 0) ele%value(rf_wavelength$) =  &
                                   param%total_length / ele%value(harmon$)

! BeamBeam

  case (beambeam$)

    if (ele%value(n_slice$) == 0) ele%value(n_slice$) = 1.0 ! revert to default

    if (ele%value(charge$) == 0 .or. param%n_part == 0) then
      ele%value(bbi_const$) = 0
      return
    endif

    if (ele%value(sig_x$) == 0 .or. ele%value(sig_y$) == 0) then
      print *, 'ERROR IN ATTRIBUTE_BOOKKEEPER: ZERO SIGMA IN BEAMBEAM ELEMENT!'
      call type_ele(ele, .true., 0, .false., 0, .false.)
      stop
    endif

    ele%value(bbi_const$) = &
        -param%n_part * m_electron * ele%value(charge$) * r_e /  &
        (2 * pi * param%beam_energy * (ele%value(sig_x$) + ele%value(sig_y$)))


! Wiggler

  case (wiggler$) 

    if (param%beam_energy == 0) then
      ele%value(k1$) = 0
    else
      ele%value(k1$) = -0.5 * &
                    (c_light * ele%value(b_max$) / param%beam_energy)**2
    endif

    if (ele%value(b_max$) == 0) then
      ele%value(rho$) = 0
    else
      ele%value(rho$) = param%beam_energy / (c_light * ele%value(b_max$))
    endif

  end select

! We need to kill the Taylor Map, etc. if things have changed.
! calculate a check sum to see if things have changed.
! ele%value(check_sum$) == 0 means that the check_sum has never been 
! computed so in this case do not kill the Taylor Map

  select case (ele%key)

  case (wiggler$)
    if (ele%sub_key == periodic_type$) return
    check_sum = ele%value(polarity$)

  case (quadrupole$)
    check_sum = ele%value(k1$) 

  case (sol_quad$)
    check_sum = ele%value(ks$) + ele%value(k1$)

  case (solenoid$)
    check_sum = ele%value(ks$)

  case (sbend$)
    check_sum = ele%value(g$) + ele%value(delta_g$) + ele%value(e1$) + &
        ele%value(e2$)

  case (sextupole$)
    check_sum = ele%value(k2$) + ele%value(tilt$)

  case (octupole$)
    check_sum = ele%value(k3$) + ele%value(tilt$)

  case (rfcavity$)
    check_sum = ele%value(volt$) + ele%value(phi0$)

  case default
    return

  end select

  check_sum = check_sum + ele%value(l$) + ele%value(x_offset$) + &
        ele%value(y_offset$) + ele%value(x_pitch$) + ele%value(y_pitch$) + &
        ele%num_steps + ele%value(s_offset$)

  if (ele%value(check_sum$) /= check_sum) then

    if (ele%value(check_sum$) == 0) then
      ele%value(check_sum$) = check_sum
      return
    endif

    ele%value(check_sum$) = check_sum
    if (associated(ele%taylor(1)%term)) call kill_taylor(ele%taylor)
    if (associated(ele%gen_field)) call kill_gen_field(ele%gen_field)
    if (ele%key == wiggler$) ele%value(z_patch$) = 0

  endif

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine transfer_ring_taylors (ring_in, ring_out, 
!                                              type_out, transfered_all)
!
! Subroutine to transfer the taylor maps from the elements of one ring to
! the elements of another. The elements are matched between the rings so 
! that the appropriate element in ring_out will get the correct Taylor map
! even if the order of the elements is different in the 2 rings.
!
! Note: The transfered Taylor map will be truncated to bmad_com%taylor_order.
! Note: If the taylor_order of an element in ring_in is less than 
!   bmad_com%taylor_order then it will not be used.  
!
! Modules needed:
!   use bmad
!
! Input:
!   ring_in   -- Ring_struct: Input ring with Taylor maps.
!   type_out  -- Logical: If True then print a message for each Taylor map
!                 transfered.
!
! Output:
!   ring_out  -- Ring_struct: Ring to receive the Taylor maps.
!   transfered_all -- Logical, optional: Set True if a Taylor map is found
!                 for all elements in ring_out that need one. False otherwise.
!-

subroutine transfer_ring_taylors (ring_in, ring_out, type_out, transfered_all)

  implicit none

  type (ring_struct), target, intent(in) :: ring_in
  type (ring_struct), target, intent(inout) :: ring_out
  type (ele_struct), pointer :: ele_in, ele_out

  integer i, j, k, it, ix
  integer n_in, ix_in(ring_in%n_ele_maxx)
 
  logical, intent(in)  :: type_out
  logical, optional :: transfered_all

! check global parameters

  if (present(transfered_all)) transfered_all = .true.

  if (ring_in%param%beam_energy /= ring_out%param%beam_energy) then
    if (type_out) then
      print *, 'TRANSFER_RING_TAYLORS: THE RING ENERGIES ARE DIFFERENT.'
      print *, '    TAYLOR MAPS NOT TRANSFERED.'
    endif
    if (present(transfered_all)) transfered_all = .false.
    return
  endif

! Find the taylor series in the first ring.

  n_in = 0
  do i = 1, ring_in%n_ele_max
    if (associated(ring_in%ele_(i)%taylor(1)%term)) then
      if (bmad_com%taylor_order > ring_in%ele_(i)%taylor_order) cycle
      n_in = n_in + 1
      ix_in(n_in) = i
    endif
  enddo

! Go through ring_out and match elements.
! If we have a match transfer the Taylor map.
! Call attribute_bookkeeper before transfering the taylor map to make sure
! the check_sum is correct. 

  out_loop: do i = 1, ring_out%n_ele_max

    ele_out => ring_out%ele_(i)

    do j = 1, n_in

      ele_in => ring_in%ele_(ix_in(j))

      if (equivalent_eles (ele_in, ele_out)) then
        if (type_out) print *, 'TRANSFER_RING_TAYLORS: ', &
             'Reusing Taylor from: ', trim(ele_in%name), '  to: ', ele_out%name
        call attribute_bookkeeper (ele_out, ring_out%param)
        call transfer_ele_taylor (ele_in, ele_out, bmad_com%taylor_order)
        cycle out_loop
      endif

    enddo

    if (ele_out%tracking_method == taylor$ .or. &
                    ele_out%mat6_calc_method == taylor$ .and. type_out) then
      print *, 'TRANSFER_RING_TAYLORS: NO TAYLOR FOR: ', ele_out%name
      if (present(transfered_all)) transfered_all = .false.
    endif

  enddo out_loop

end subroutine

end module
