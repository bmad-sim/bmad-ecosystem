#include "CESR_platform.inc"

module bookkeeper_mod

  use bmad_interface
  use bmad_utils_mod
  use multipole_mod

  integer, parameter :: off$ = 1, on$ = 2
  integer, parameter :: save_state$ = 3, restore_state$ = 4
        
contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine control_bookkeeper (ring, ix_ele)
!
! Subroutine to transfer attibute information from lord to slave elements.
! Note: To do a complete bookkeeping job on a lattice use:
!   lattice_bookkeeper
!
! Modules needed:
!   use bmad
!
! Input:
!   ring   -- Ring_struct: Ring to be used
!   ix_ele -- Integer, optional: Index of element whose attribute values 
!               have been changed. If not present bookkeeping will be done 
!               for all elements.
!-

subroutine control_bookkeeper (ring, ix_ele)

  implicit none

  type (ring_struct), target :: ring
  type (ele_struct), pointer :: lord, slave

  integer, optional :: ix_ele
  integer ie

! If ix_ele is present we only do bookkeeping for this one element

  if (present(ix_ele)) then
    call this_bookkeeper (ix_ele)

! Else we need to make up all the lords. 
! The group lords must be done last since their slaves don't know about them.

  else
    do ie = ring%n_ele_use+1, ring%n_ele_max
      if (ring%ele_(ie)%control_type /= group_lord$) call this_bookkeeper (ie)
    enddo

    do ie = ring%n_ele_use+1, ring%n_ele_max
      if (ring%ele_(ie)%control_type == group_lord$) call this_bookkeeper (ie)
    enddo

  endif

!--------------------------------------------------------------------------
contains

subroutine this_bookkeeper (ix_ele)

  integer ix_ele, j, k, ix, ix1, ix2, ix_lord
  integer ix_slaves(300)

! Attribute bookkeeping for this element

  call attribute_bookkeeper (ring%ele_(ix_ele), ring%param)
  if (ring%ele_(ix_ele)%n_slave == 0 .and. &
            ring%ele_(ix_ele)%n_lord == 0) return  ! nothing more to do

! Make a list of slave elements to update.
! we do not need to update free elements of group lords.

  ix1 = 0   ! index for processed elements
  ix_slaves(1) = ix_ele
  ix2 = 1   ! index for last element in list

  do
    ix1 = ix1 + 1
    ix_lord = ix_slaves(ix1)
    lord => ring%ele_(ix_lord)
    do j = lord%ix1_slave, lord%ix2_slave
      ix = ring%control_(j)%ix_slave
      if (lord%control_type == group_lord$ .and. ring%ele_(ix)%control_type == free$) cycle
      if (ix == ix_slaves(ix2)) cycle   ! do not use duplicates
      ix2 = ix2 + 1
      ix_slaves(ix2) = ix
    enddo
    if (ix1 == ix2) exit
  enddo

! First: Makup lords and group slaves.
! If an overlay_lord has lords above it then these lords must be overlay_lords.
! Therefore treat the overlay_lord as an overlay_slave.
! The same is true if a super_lord has lords except in this case the lord
! may be a multipass_lord.

  do j = 1, ix2

    ix = ix_slaves(j)
    slave => ring%ele_(ix)

    if (slave%control_type == group_lord$) then
      call makeup_group_slaves (ring, ix)
      call attribute_bookkeeper (slave, ring%param)

    elseif (slave%control_type == super_lord$ .and. slave%n_lord > 0) then
      k =  ring%ic_(slave%ic1_lord)
      lord => ring%ele_(ring%control_(k)%ix_lord)
      if (lord%control_type == multipass_lord$) then
        call makeup_multipass_slave (ring, ix)
      else
        call adjust_super_lord_s_position (ring, ix)
        call makeup_overlay_and_i_beam_slave (ring, ix)
      endif
      call attribute_bookkeeper (slave, ring%param)

    elseif (slave%control_type == multipass_lord$ .and. slave%n_lord > 0) then
      call makeup_overlay_and_i_beam_slave (ring, ix)
      call attribute_bookkeeper (slave, ring%param)

    elseif (slave%control_type == overlay_lord$ .and. slave%n_lord > 0) then
      call makeup_overlay_and_i_beam_slave (ring, ix)
      call attribute_bookkeeper (slave, ring%param)

    endif

  enddo

! Second: Makeup all slaves but group slaves.

  do j = 1, ix2

    ix = ix_slaves(j)
    slave => ring%ele_(ix)       

    if (slave%control_type == super_slave$) then
      call makeup_super_slave (ring, ix)
      call attribute_bookkeeper (slave, ring%param)

    elseif (slave%control_type == overlay_slave$) then
      call makeup_overlay_and_i_beam_slave (ring, ix)
      call attribute_bookkeeper (slave, ring%param)

    elseif (slave%control_type == multipass_slave$) then
      call makeup_multipass_slave (ring, ix)
      call attribute_bookkeeper (slave, ring%param)

    endif

  enddo

end subroutine

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine lattice_bookkeeper (ring)
!
! Subroutine to do a complete bookkeeping job on a lattice.
!
! Modules needed:
!   use bmad
!
! Input:
!   ring   -- Ring_struct: Lattice needing bookkeeping.
!
! Output:
!   ring   -- Ring_struct: Lattice with bookkeeping done.
!-

subroutine lattice_bookkeeper (ring)

  implicit none

  type (ring_struct) ring
  integer i

!

  call control_bookkeeper (ring)
  call compute_element_energy (ring)

  do i = 1, ring%n_ele_use
    call attribute_bookkeeper (ring%ele_(i), ring%param)
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

  character(40) :: r_name = 'adjust_super_lord_s_position'

!

  lord => ring%ele_(ix_lord)

  if (lord%control_type /= super_lord$) then
     call out_io (s_abort$, r_name, 'ELEMENT IS NOT A LORD! ' // lord%name)
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

  character(20) :: r_name = 'makeup_group_slaves'

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
        call out_io (s_abort$, r_name, "A GROUP: " // lord%name, &
                    "CONTROLS THE LENGTH OF A LORD ELEMENT: " // slave%name)
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
! Subroutine makeup_multipass_slave (ring, ix_slave)
!
! Subroutine to calcualte the attributes of multipass slave elements.
! This routine is not meant for general use.
!-

subroutine makeup_multipass_slave (ring, ix_slave)

  implicit none

  type (ring_struct), target :: ring
  type (ele_struct), pointer :: lord, slave

  real(rp) s, val(n_attrib_maxx)
  integer j, ix_slave

!

  slave => ring%ele_(ix_slave)
  j =  ring%ic_(slave%ic1_lord)
  lord => ring%ele_(ring%control_(j)%ix_lord)

  val = slave%value

  slave%value = lord%value
  if (lord%key == lcavity$ .or. lord%key == rfcavity$) then
    slave%value(dphi0$)        = val(dphi0$)
    slave%value(energy_start$) = val(energy_start$)
  endif

  slave%value(beam_energy$) = val(beam_energy$)
  slave%value(p0c$)         = val(p0c$)

  if (associated (slave%a)) then
    slave%a = lord%a
    slave%b = lord%b
  endif

  if (associated (slave%r)) slave%r = lord%r
  if (associated (slave%const)) slave%const = lord%const
  if (associated (slave%wake)) then
    slave%wake%sr1       = lord%wake%sr1
    slave%wake%sr2_long  = lord%wake%sr2_long
    slave%wake%sr2_trans = lord%wake%sr2_trans
    slave%wake%lr        = lord%wake%lr
  endif

  slave%mat6_calc_method = lord%mat6_calc_method
  slave%tracking_method  = lord%tracking_method
  slave%num_steps        = lord%num_steps      
  slave%is_on            = lord%is_on
  slave%aperture_at      = lord%aperture_at

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine makeup_super_slave (ring, ix_slave)
!
! Subroutine to calcualte the attributes of superposition slave elements.
! This routine is not meant for general use.
!-
         
subroutine makeup_super_slave (ring, ix_slave)

  implicit none

  type (ring_struct), target :: ring
  type (ele_struct), pointer :: lord, slave
  type (ele_struct), save :: sol_quad

  integer i, ix_con, j, ix, ix_slave

  real(rp) tilt, k_x, k_y, x_kick, y_kick, ks, k1, coef
  real(rp) x_o, y_o, x_p, y_p, s_slave, s_del, k2, k3
  real(rp) sin_n, cos_n, a(0:n_pole_maxx), b(0:n_pole_maxx)
  real(rp) knl(0:n_pole_maxx), t(0:n_pole_maxx), value(n_attrib_maxx)
  real(rp) a_tot(0:n_pole_maxx), b_tot(0:n_pole_maxx)
  real(rp) sum_1, sum_2, sum_3, sum_4, ks_sum, ks_xp_sum, ks_xo_sum
  real(rp) ks_yp_sum, ks_yo_sum, l_slave, r_off(4)
  real(rp) t_1(4), t_2(4), T_end(4,4), mat4(4,4), mat4_inv(4,4), beta(4)
  real(rp) T_tot(4,4), x_o_sol, x_p_sol, y_o_sol, y_p_sol

  logical, save :: init_needed = .true.

  character(20) :: r_name = 'makeup_super_slave'

! init

  if (init_needed) then
    call init_ele (sol_quad)
    sol_quad%key = sol_quad$
    init_needed = .false.
  endif

! Super_slave:

  slave => ring%ele_(ix_slave)

  if (slave%control_type /= super_slave$) then
     call out_io(s_abort$, r_name, "ELEMENT IS NOT AN SUPER SLAVE: " // slave%name)
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
    if (slave%key == rfcavity$) value(voltage$) = lord%value(voltage$) * coef

    slave%aperture_at = no_end$
    call compute_slave_aperture (value, slave, lord, ix_con)

! s_del is the distance between lord and slave centers

    s_del = (slave%s - slave%value(l$)/2) - &
                  (lord%s + lord%value(s_offset$) - lord%value(l$)/2)
    s_del = modulo2 (s_del, ring%param%total_length/2)
    value(x_pitch$) = value(x_pitch_tot$)
    value(y_pitch$) = value(y_pitch_tot$)
    value(x_offset$) = value(x_offset_tot$) + s_del * value(x_pitch_tot$)
    value(y_offset$) = value(y_offset_tot$) + s_del * value(y_pitch_tot$)
    value(tilt$)         = value(tilt_tot$)

    slave%value = value
    slave%is_on = lord%is_on
    slave%mat6_calc_method = lord%mat6_calc_method
    slave%tracking_method  = lord%tracking_method

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
        slave%value(e2$)    = 0
        slave%value(h2$)    = 0
        slave%value(fintx$) = 0
        slave%value(hgapx$) = 0
      elseif (ix_con == lord%ix2_slave) then 
        slave%value(e1$)    = 0
        slave%value(h1$)    = 0
        slave%value(fint$)  = 0
        slave%value(hgap$)  = 0
      else
        slave%value(e1$)    = 0
        slave%value(h1$)    = 0
        slave%value(fint$)  = 0
        slave%value(hgap$)  = 0
        slave%value(e2$)    = 0
        slave%value(h2$)    = 0
        slave%value(fintx$) = 0
        slave%value(hgapx$) = 0
      endif
    endif                       

    return

  endif

!-----------------------------------------------------------------------
! Multiple super_lords for this super_slave: 
! combine the lord elements.
                                           
  k_x = 0
  k_y = 0
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
  value(beam_energy$) = slave%value(beam_energy$)
  value(p0c$) = slave%value(p0c$)
  value(check_sum$) = slave%value(check_sum$) ! do not change the check_sum

  s_slave = slave%s - value(l$)/2  ! center of slave

! sum over all lords...

  do j = slave%ic1_lord, slave%ic2_lord

    ix_con = ring%ic_(j)
    ix = ring%control_(ix_con)%ix_lord
    coef = ring%control_(ix_con)%coef
    lord => ring%ele_(ix)

    if (lord%control_type /= super_lord$) then
      call out_io (s_abort$, r_name, &
            "SUPER_SLAVE HAS A CONTROL ELEMENT THAT IS NOT A SUPER_LORD", &
            'SLAVE: ' //  slave%name // '  \i\ ', &
            'LORD:  ' //  lord%name  // '  \i\ ', i_array = (/ ix_slave, ix /) )
      call err_exit
    endif

    call compute_slave_aperture (value, slave, lord, ix_con)

    if (j == slave%ic1_lord) then
      slave%mat6_calc_method = lord%mat6_calc_method
      slave%tracking_method  = lord%tracking_method
    else
      if (slave%mat6_calc_method /= lord%mat6_calc_method) then
        ix = ring%control_(ring%ic_(slave%ic1_lord))%ix_lord
        call out_io(s_abort$, r_name, 'MAT6_CALC_METHOD DOES NOT AGREE FOR DIFFERENT', &
             'SUPERPOSITION LORDS: ' // trim(lord%name) // ', ' // trim(ring%ele_(ix)%name))
        call err_exit
      endif
      if (slave%tracking_method /= lord%tracking_method) then
        ix = ring%control_(ring%ic_(slave%ic1_lord))%ix_lord
        call out_io(s_abort$, r_name, ' TRACKING_METHOD DOES NOT AGREE FOR DIFFERENT', &
             'SUPERPOSITION LORDS: ' // trim(lord%name) // ', ' // trim(ring%ele_(ix)%name))
        call err_exit
      endif
    endif

    if (.not. lord%is_on) cycle

    x_kick = x_kick + lord%value(hkick$) * coef
    y_kick = y_kick + lord%value(vkick$) * coef
    tilt = lord%value(tilt_tot$)

    if (associated(lord%a)) then
      call multipole_ele_to_kt (lord, +1, knl, t, .true.)
      call multipole_kt_to_ab (knl/lord%value(l$), t, a, b)
      a_tot = a_tot + a
      b_tot = b_tot + b
    endif

!------

    select case (slave%key)


! sextupole

    case (sextupole$) 

      cos_n = lord%value(k2$) * cos(3 * tilt)
      sin_n = lord%value(k2$) * sin(3 * tilt)            
      
      k_x = k_x + cos_n
      k_y = k_y + sin_n
  
! octupole

    case (octupole$)

      cos_n = lord%value(k3$) * cos(4 * tilt)
      sin_n = lord%value(k3$) * sin(4 * tilt)        
      
      k_x = k_x + cos_n
      k_y = k_y + sin_n

! solenoid/quadrupole combo.

    case (solenoid$, sol_quad$, quadrupole$)

      x_p = lord%value(x_pitch_tot$);  x_o = lord%value(x_offset_tot$)
      y_p = lord%value(y_pitch_tot$);  y_o = lord%value(y_offset_tot$)

      s_del = s_slave - (lord%s + lord%value(s_offset_tot$) - lord%value(l$)/2)
      s_del = modulo2 (s_del, ring%param%total_length/2)

      ks = lord%value(ks$)

      ks_sum = ks_sum + ks

      ks_xp_sum = ks_xp_sum + ks * x_p
      ks_yp_sum = ks_yp_sum + ks * y_p

      ks_xo_sum = ks_xo_sum + ks * (x_o + x_p * s_del)
      ks_yo_sum = ks_yo_sum + ks * (y_o + y_p * s_del)

      cos_n = lord%value(k1$) * cos(2 * tilt)
      sin_n = lord%value(k1$) * sin(2 * tilt)

      k_x = k_x + cos_n
      k_y = k_y + sin_n

      sum_1 = sum_1 + cos_n * x_p + sin_n * y_p
      sum_2 = sum_2 + sin_n * x_p - cos_n * y_p

      sum_3 = sum_3 + cos_n * (x_o + x_p * s_del) + sin_n * (y_o + y_p * s_del)
      sum_4 = sum_4 + sin_n * (x_o + x_p * s_del) - cos_n * (y_o + y_p * s_del)

! bend_sol_quad

    case (bend_sol_quad$)
      call out_io (s_abort$, r_name, &
                 'CODING NOT YET IMPLEMENTED FOR A: ' // key_name(slave%key))
      call err_exit


! default

    case default
      call out_io (s_abort$, r_name, &
                 'CODING NOT YET IMPLEMENTED FOR A: ' // key_name(slave%key))
      call err_exit

    end select

  enddo

!------------------------------
! stuff sums into slave element

  value(hkick$) = x_kick
  value(vkick$) = y_kick

  slave%value = value

  if (any(a_tot /= 0) .or. any(b_tot /= 0)) then
    call multipole_init(slave)
    call multipole_ab_to_kt(a_tot, b_tot, knl, t)
    call multipole_kt_to_ab(knl/k1, t-tilt, a, b)
    slave%a = a
    slave%b = b
    slave%value(radius$) = 1
  elseif (associated(slave%a)) then
    deallocate (slave%a, slave%b)
  endif

!-----------------------------

  select case (slave%key)

  case (sextupole$) 

    if (k_x == 0 .and. k_y == 0) return

    k2 = sqrt(k_x**2 + k_y**2)
    tilt = atan2(k_y, k_x) / 3

    if (tilt > pi/6) then
      k2 = -k2
      tilt = tilt - pi/3
    elseif (tilt < -pi/6) then
      k2 = -k2
      tilt = tilt + pi/3
    endif

    slave%value(k2$) = k2
    slave%value(tilt$) = tilt

! octupole

  case (octupole$)

    if (k_x == 0 .and. k_y == 0 .and. ks == 0) return

    k3 = sqrt(k_x**2 + k_y**2)
    tilt = atan2(k_y, k_x) / 4

    if (tilt > pi/8) then
      k3 = -k3
      tilt = tilt - pi/4
    elseif (tilt < -pi/8) then
      k3 = -k3
      tilt = tilt + pi/4
    endif

    slave%value(k3$) = k3
    slave%value(tilt$) = tilt

! sol_quad, etc.

  case (solenoid$, sol_quad$, quadrupole$)

    ks = ks_sum
    slave%value(ks$) = ks

    if (k_x == 0 .and. k_y == 0 .and. ks == 0) return

    if (ks /= 0) then
      x_o_sol = ks_xo_sum / ks
      x_p_sol = ks_xp_sum / ks
      y_o_sol = ks_yo_sum / ks
      y_p_sol = ks_yp_sum / ks
    endif

    if (k_x == 0 .and. k_y == 0) then  ! pure solenoid
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

    k1 = sqrt(k_x**2 + k_y**2)
    tilt = atan2(k_y, k_x) / 2

    if (tilt > pi/4) then
      k1 = -k1
      tilt = tilt - pi/2
    elseif (tilt < -pi/4) then
      k1 = -k1
      tilt = tilt + pi/2
    endif

    slave%value(k1$) = k1
    slave%value(tilt$) = tilt

    cos_n = k_x / (k_x**2 + k_y**2)
    sin_n = k_y / (k_x**2 + k_y**2)

    slave%value(x_pitch$)  = cos_n * sum_1 + sin_n * sum_2
    slave%value(y_pitch$)  = sin_n * sum_1 - cos_n * sum_2
    slave%value(x_offset$) = cos_n * sum_3 + sin_n * sum_4
    slave%value(y_offset$) = sin_n * sum_3 - cos_n * sum_4

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

! bend_sol_quad

  case (bend_sol_quad$)
    call out_io (s_abort$, r_name, &
                 'CODING NOT YET IMPLEMENTED FOR A: ' // key_name(slave%key))
    call err_exit

  end select

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine compute_slave_aperture (value, slave, lord, ix_con)
!
! This routine is not meant for general use.
!-

subroutine compute_slave_aperture (value, slave, lord, ix_con)

  use bmad_struct

  implicit none

  type (ele_struct) slave, lord
  real(rp) value(n_attrib_maxx)
  integer ix_con

!

  select case (lord%aperture_at)
  case (exit_end$) 
    if (ix_con == lord%ix2_slave) slave%aperture_at = exit_end$
  case (entrance_end$)
    if (ix_con == lord%ix1_slave) slave%aperture_at = entrance_end$
  case (both_ends$)
    if (ix_con == lord%ix1_slave .and. ix_con == lord%ix2_slave) then
      slave%aperture_at = both_ends$
    elseif (ix_con == lord%ix1_slave) then
      slave%aperture_at = entrance_end$
    elseif (ix_con == lord%ix2_slave) then 
      slave%aperture_at = exit_end$
    endif
  end select

  if (slave%aperture_at == no_end$) then
    value(x_limit$) = 0
    value(y_limit$) = 0
    value(aperture$) = 0
  endif

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine makeup_overlay_and_i_beam_slave (ring, ix_ele)
!
! This routine is not meant for general use.
!-

subroutine makeup_overlay_and_i_beam_slave (ring, ix_ele)

  implicit none

  type (ring_struct), target :: ring
  type (ele_struct), pointer :: ele, i_beam

  real(rp) value(n_attrib_maxx), coef, ds
  integer i, j, ix, iv, ix_ele, icom, ct
  logical used(n_attrib_maxx)

  character(40) :: r_name = 'makeup_overlay_and_i_beam_slave'

!
                               
  ele => ring%ele_(ix_ele)
  ct = ele%control_type

  if (ct /= super_lord$ .and. ct /= overlay_slave$ .and. &
           ct /= overlay_lord$ .and. ct /= multipass_lord$) then
    call out_io(s_abort$, r_name, 'ELEMENT IS NOT OF PROPER TYPE. RING INDEX: \i\ ', ix_ele)
    call type_ele (ele, .true., 0, .false., 0, .true., ring)
    call err_exit
  endif

  value = 0
  used = .false.
  ele%on_an_i_beam = .false.

  do i = ele%ic1_lord, ele%ic2_lord
    j = ring%ic_(i)
    ix = ring%control_(j)%ix_lord

    if (ring%ele_(ix)%control_type == i_beam_lord$) then
      i_beam => ring%ele_(ix)
      ds = (ele%s - ele%value(l$)/2) - i_beam%value(s_center$) 
      ele%value(x_offset_tot$) = ele%value(x_offset$) + &
                     ds * i_beam%value(x_pitch$) + i_beam%value(x_offset$)
      ele%value(y_offset_tot$) = ele%value(y_offset$) + &
                     ds * i_beam%value(y_pitch$) + i_beam%value(y_offset$)
      ele%value(s_offset_tot$) = ele%value(s_offset$) + i_beam%value(s_offset$)
      ele%value(x_pitch_tot$)  = ele%value(x_pitch$)  + i_beam%value(x_pitch$)
      ele%value(y_pitch_tot$)  = ele%value(y_pitch$)  + i_beam%value(y_pitch$)
      ele%value(tilt_tot$)     = ele%value(tilt$)     + i_beam%value(tilt$)
      ele%on_an_i_beam = .true.
      cycle
    endif

    if (ring%ele_(ix)%control_type /= overlay_lord$) then
      call out_io (s_abort$, r_name, 'THE LORD IS NOT AN OVERLAY_LORD \i\ ', ix_ele)
      call type_ele (ele, .true., 0, .false., 0, .true., ring)
      call err_exit
    endif     

    coef = ring%control_(j)%coef
    iv = ring%control_(j)%ix_attrib
    icom = ring%ele_(ix)%ix_value
    value(iv) = value(iv) + ring%ele_(ix)%value(icom)*coef
    used(iv) = .true.
  enddo

  where (used) ele%value = value

! If no i_beam then simply transfer tilt to tilt_tot, etc.

  if (.not. ele%on_an_i_beam) then
    ele%value(tilt_tot$)     = ele%value(tilt$)
    ele%value(x_offset_tot$) = ele%value(x_offset$)
    ele%value(y_offset_tot$) = ele%value(y_offset$)
    ele%value(s_offset_tot$) = ele%value(s_offset$)
    ele%value(x_pitch_tot$)  = ele%value(x_pitch$)
    ele%value(y_pitch_tot$)  = ele%value(y_pitch$)
  endif

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
! BEAMBEAM:   
!     bbi_const$ = param%n_part * m_electron * charge$ * r_e /
!                           (2 * pi * beam_energy$ * (sig_x$ + sig_y$)
!
! ELSEPARATOR:
!     e_field$ = sqrt(hkick$**2 + vkick$**2) * beam_energy$ / l$
!     voltage$ = e_field$ * gap$ 
!
! LCAVITY:    
!     delta_e$ = gradient$ * L$ 
!
! RFCAVITY:   
!     rf_frequency$ = harmon$ * c_light / param%total_length (only if harmon$ /= 0)
!
! SBEND:      
!     angle$   = L$ * G$
!     l_chord$ = 2 * sin(Angle$/2) / G$
!     rho$     = 1 / G$
!
! WIGGLER:    
!     k1$  = -0.5 * (c_light * b_max$ / beam_energy$)**2
!     rho$ = beam_energy$ / (c_light * b_max$)
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
!
! Programming Note: If the dependent attributes are changed then 
!       attribute_free must be modified.
!-

subroutine attribute_bookkeeper (ele, param)

  implicit none

  type (ele_struct) ele
  type (param_struct) param

  real(rp) factor, check_sum
  
  character(20) ::  r_name = 'attribute_bookkeeper'

! Transfer tilt to tilt_tot, etc.

  if (.not. ele%on_an_i_beam .and. ele%key /= match$) then
    ele%value(tilt_tot$)     = ele%value(tilt$)
    ele%value(x_offset_tot$) = ele%value(x_offset$)
    ele%value(y_offset_tot$) = ele%value(y_offset$)
    ele%value(s_offset_tot$) = ele%value(s_offset$)
    ele%value(x_pitch_tot$)  = ele%value(x_pitch$)
    ele%value(y_pitch_tot$)  = ele%value(y_pitch$)
  endif

! field_master

  if (ele%field_master) then

    if (ele%value(p0c$) == 0) then
      factor = 0
    else
      factor = c_light / ele%value(p0c$)
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
    case (hkicker$)
      ele%value(kick$) = factor * ele%value(B_field$)
    case (vkicker$)
      ele%value(kick$) = factor * ele%value(B_field$)
    case default
       call out_io(s_abort$,r_name,' "FIELD_MASTER" NOT IMPLEMENTED FOR: ' // trim(ele%name))
      call err_exit
    end select

    ele%value(hkick$) = factor * ele%value(hkick_B_field$)
    ele%value(vkick$) = factor * ele%value(vkick_B_field$)

  else

    if (ele%value(p0c$) == 0) then
      factor = 0
    else
      factor = ele%value(p0c$) / c_light
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
    case (hkicker$)
      ele%value(B_field$) = factor * ele%value(kick$)
    case (vkicker$) 
      ele%value(B_field$) = factor * ele%value(kick$)
    end select

    ele%value(hkick_B_field$) = factor * ele%value(hkick$)
    ele%value(vkick_B_field$) = factor * ele%value(vkick$)

  endif

! Dependent attribute bookkeeping.
! Note: If the dependent attributes are changed then attribute_free 
!       must be modified.

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
    ele%value(delta_e$) = ele%value(gradient$) * ele%value(L$) 
    

! RFcavity

  case (rfcavity$)
    if (ele%value(harmon$) /= 0) ele%value(rf_frequency$) =  &
                              ele%value(harmon$) * c_light / param%total_length 

! BeamBeam

  case (beambeam$)

    if (ele%value(n_slice$) == 0) ele%value(n_slice$) = 1.0 ! revert to default

    if (ele%value(charge$) == 0 .or. param%n_part == 0) then
      ele%value(bbi_const$) = 0

    else

      if (ele%value(sig_x$) == 0 .or. ele%value(sig_y$) == 0) then
        call out_io(s_abort$, r_name, 'ZERO SIGMA IN BEAMBEAM ELEMENT!')
        call type_ele(ele, .true., 0, .false., 0, .false.)
        call err_exit
      endif

      ele%value(bbi_const$) = &
        -param%n_part * m_electron * ele%value(charge$) * r_e /  &
        (2 * pi * ele%value(beam_energy$) * (ele%value(sig_x$) + ele%value(sig_y$)))

    endif

! Elseparator

  case (elseparator$)

    if (ele%value(l$) == 0 .or. ele%value(gap$) == 0) then
      ele%value(e_field$) = 0
      ele%value(voltage$) = 0
    else
      ele%value(e_field$) = sqrt(ele%value(hkick$)**2 + ele%value(vkick$)**2) * &
                                               ele%value(beam_energy$) / ele%value(l$)
      ele%value(voltage$) = ele%value(e_field$) * ele%value(gap$) 
    endif


! Wiggler

  case (wiggler$) 

    if (ele%value(beam_energy$) == 0) then
      ele%value(k1$) = 0
    else
      ele%value(k1$) = -0.5 * &
                    (c_light * ele%value(b_max$) / ele%value(beam_energy$))**2
    endif

    if (ele%value(b_max$) == 0) then
      ele%value(rho$) = 0
    else
      ele%value(rho$) = ele%value(beam_energy$) / (c_light * ele%value(b_max$))
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
    check_sum = ele%value(voltage$) + ele%value(phi0$)

  case (elseparator$)
    check_sum = check_sum + ele%value(hkick$) + ele%value(vkick$)

  case default
    return

  end select

  check_sum = check_sum + ele%value(l$) + ele%value(x_offset$) + &
        ele%value(y_offset$) + ele%value(x_pitch$) + ele%value(y_pitch$) + &
        ele%num_steps + ele%value(s_offset$) + ele%value(tilt$)

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

  integer i, j
  integer n_in, ix_in(ubound(ring_in%ele_, 1))
 
  logical, intent(in)  :: type_out
  logical, optional :: transfered_all

  character(25) :: r_name = 'transfer_ring_taylors'

! check global parameters

  if (present(transfered_all)) transfered_all = .true.

  if (ring_in%ele_(0)%value(beam_energy$) /= &
                              ring_out%ele_(0)%value(beam_energy$)) then
    if (type_out) then
       call out_io (s_warn$, r_name, &
              'THE RING ENERGIES ARE DIFFERENT. TAYLOR MAPS NOT TRANSFERED.')
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
        if (type_out) call out_io (s_info$, r_name, &
            ' Reusing Taylor from: ' // trim(ele_in%name) // '  to: ' //  ele_out%name)
        call attribute_bookkeeper (ele_out, ring_out%param)
        call transfer_ele_taylor (ele_in, ele_out, bmad_com%taylor_order)
        cycle out_loop
      endif

    enddo

    if (ele_out%tracking_method == taylor$ .or. &
                    ele_out%mat6_calc_method == taylor$ .and. type_out) then
      call out_io (s_warn$, r_name, ' NO TAYLOR FOR: ' // ele_out%name)
      if (present(transfered_all)) transfered_all = .false.
    endif

  enddo out_loop

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine set_on_off (key, ring, switch, orb_)
!
! Subroutine to turn on or off a set of elements (quadrupoles, rfcavities,
! etc.) in a ring. An element that is turned off acts like a drift.
! RING_MAKE_MAT6 will be called to remake ring%ele_()%mat6.
!
!
! Modules needed:
!   use bmad
!
! Input:
!   key      -- Integer: Key name of elements to be turned on or off.
!                  [Key = quadrupole$, etc.]
!   ring     -- Ring_struct: Ring structure holding the elements
!   switch   -- Integer: 
!                 on$            => Turn elements on.  
!                 off$           => Turn elements off. 
!                 save_state$    => Save present on/off state. 
!                                     No turning on or off is done.
!                 restore_state$ => Restore saved on/off state.
!   orb_(0:) -- Coord_struct, optional: Needed for ring_make_mat6
!
! Output:
!   ring -- Ring_struct: Modified ring.
!-

#include "CESR_platform.inc"
                                    
subroutine set_on_off (key, ring, switch, orb_)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ring_struct) ring
  type (coord_struct), optional :: orb_(0:)

  integer i, key               
  integer, intent(in) :: switch

  character(20) :: r_name = 'set_on_off'

!

  do i = 1, ring%n_ele_max

    if (ring%ele_(i)%key /= key) cycle

    select case (switch)
    case (on$) 
      ring%ele_(i)%is_on = .true.
    case (off$)
      ring%ele_(i)%is_on = .false.
    case (save_state$)
      ring%ele_(i)%internal_logic = ring%ele_(i)%is_on
      cycle
    case (restore_state$)
      ring%ele_(i)%is_on = ring%ele_(i)%internal_logic
    case default
      call out_io (s_abort$, r_name, 'BAD SWITCH: \i\ ', switch)
      call err_exit
    end select

    call ring_make_mat6(ring, i, orb_)

  enddo

end subroutine

end module
