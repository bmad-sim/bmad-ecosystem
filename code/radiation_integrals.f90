!+
! Subroutine RADIATION_INTEGRALS (RING, ORB_, MODE)
!
! Subroutine to calculate the synchrotron radiation integrals along with the
! emittance, and energy spread.
!
! Modules to use:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!   RING      -- Ring_struct: Ring to use. The calculation assumes that 
!                    the Twiss parameters have been calculated.
!   ORB_(0:*) -- Coord_struct: Closed orbit.
!
! Output:
!   MODE -- Modes_struct: Parameters for the ("horizontal like") a-mode,
!                              ("vertical like") b-mode, and the z-mode
!     %SYNCH_INT(1:3) -- Synchrotron integrals.
!     %SIG_E          -- Sigma_E/E energy spread
!     %SIG_Z          -- Bunch Length
!     %ENERGY_LOSS    -- Energy loss in GeV per turn
!     %A, %B, %Z      -- Amode_struct: Substructure
!       %EMITTANCE      -- Emittance
!       %SYNCH_INT(4:5) -- Synchrotron integrals
!       %J_DAMP         -- Damping partition factor
!       %ALPHA_DAMP     -- Exponential damping coefficient per turn
!
! Notes:
!     1) %SYNCH_INT(1) = momentum_compaction * ring_length
!-       

subroutine radiation_integrals (ring, orb_, mode)
                     
  use nr
  use rad_int_common

  implicit none

  type (ring_struct), target :: ring
  type (coord_struct), target :: orb_(0:*)
  type (modes_struct) mode

  real, parameter :: c_gam = 4.425e-5, c_q = 3.84e-13
  real i1, i2, i3, i4a, i4b, i4z, i5a, i5b, m65, G_max, g3_ave
  real theta, energy, gamma2_factor, energy_loss, arg

  integer ir, key

!----------------------------------------------------------------------

  interface
    function eval_i1 (s_vec)
      use nrtype
      real, intent(in) :: s_vec(:)
      real, dimension(size(s_vec)) :: eval_i1
    end function
  end interface

  interface
    function eval_i4a (s_vec)
      use nrtype
      real, intent(in) :: s_vec(:)
      real, dimension(size(s_vec)) :: eval_i4a
    end function
  end interface

  interface
    function eval_i4b (s_vec)
      real, intent(in) :: s_vec(:)
      real, dimension(size(s_vec)) :: eval_i4b
    end function
  end interface

  interface
    function eval_i5a (s_vec)
      real, intent(in) :: s_vec(:)
      real, dimension(size(s_vec)) :: eval_i5a
    end function
  end interface

  interface
    function eval_i5b (s_vec)
      real, intent(in) :: s_vec(:)
      real, dimension(size(s_vec)) :: eval_i5b
    end function
  end interface

!---------------------------------------------------------------------
! init

  m65 = 0

  rad_com%i1_ = 0;  rad_com%i2_ = 0;  rad_com%i3_ = 0
  rad_com%i4a_ = 0;  rad_com%i4b_ = 0
  rad_com%i5a_ = 0;  rad_com%i5b_ = 0

!---------------------------------------------------------------------
! loop over all elements

  pring => ring

  do ir = 1, ring%n_ele_use     

    ele => ring%ele_(ir)         
    orb1 => orb_(ir)

    if (.not. ele%is_on) cycle

    ele0 => ring%ele_(ir-1)
    orb0 => orb_(ir-1)

    rad_com%ll = ele%value(l$) 
    if (rad_com%ll == 0) cycle

    key = ele%key

    rad_com%g_x0 = -ele%value(hkick$) / rad_com%ll
    rad_com%g_y0 = -ele%value(vkick$) / rad_com%ll

    if (key == rfcavity$) m65 = m65 + ele%mat6(6,5) 

! for a wiggler we make the approximation that the variation of G is
! fast compaired to the variation in eta.

    if (key == wiggler$) then
      G_max = sqrt(2*abs(ele%value(k1$)))       ! 1/rho at max B
      rad_com%i2_(ir) = rad_com%ll * G_max**2 / 2
      g3_ave = 4 * G_max**3 / (3 * pi)
      rad_com%i3_(ir) = rad_com%ll * g3_ave
      rad_com%g_x0 = g3_ave**(1.0/3)
      rad_com%g_y0 = 0
      rad_com%k1 = 0
      rad_com%s1 = 0
      rad_com%i5a_(ir) = qromb_rad_int (eval_i5a, 0.0, rad_com%ll, i5a)
      rad_com%i5b_(ir) = qromb_rad_int (eval_i5b, 0.0, rad_com%ll, i5b)
      cycle
    endif
    
    if (key == custom$) then
      call custom_radiation_integrals (ring, ir, orb_)
      cycle
    endif

    if (rad_com%g_x0 == 0 .and. rad_com%g_y0 == 0 .and. &
          key /= quadrupole$ .and. key /= sol_quad$ .and. key /= sbend$) cycle

    if (key == sbend$) then
      theta = ele%value(tilt$) + ele%value(roll$)
      rad_com%g_x0 = rad_com%g_x0 + cos(theta) / ele%value(rho$) 
      rad_com%g_y0 = rad_com%g_y0 - sin(theta) / ele%value(rho$) 
    endif

    rad_com%g2 = rad_com%g_x0**2 + rad_com%g_y0**2
    rad_com%g = sqrt(rad_com%g2)

    rad_com%i2_(ir)  = rad_com%g2 * rad_com%ll
    rad_com%i3_(ir)  = rad_com%g2 * rad_com%g * rad_com%ll

    if (key == quadrupole$ .or. key == sol_quad$) then
      theta = ele%value(tilt$)
      rad_com%k1 = ele%value(k1$) * cos(2*theta)
      rad_com%s1 = ele%value(k1$) * sin(2*theta)
    elseif (key == sbend$) then
      theta = ele%value(tilt$) + ele%value(roll$)
      rad_com%k1 = ele%value(k1$) * cos(2*theta)
      rad_com%s1 = ele%value(k1$) * sin(2*theta)
    else
      rad_com%k1 = 0
      rad_com%s1 = 0
    endif

! edge effects for a bend. In this case we ignore any rolls.

    if (key == sbend$) then
      call propagate_part_way(0.0)
      rad_com%i4a_(ir) = -rad_com%eta_a(1) * rad_com%g2 * tan(ele%value(e1$))
      rad_com%i4b_(ir) = -rad_com%eta_b(1) * rad_com%g2 * tan(ele%value(e1$))
      call propagate_part_way(rad_com%ll)
      rad_com%i4a_(ir) = rad_com%i4a_(ir) - &
                           rad_com%eta_a(1) * rad_com%g2 * tan(ele%value(e2$))
      rad_com%i4b_(ir) = rad_com%i4a_(ir) - &
                           rad_com%eta_b(1) * rad_com%g2 * tan(ele%value(e2$))
    endif

! integrate 

    rad_com%i1_(ir)  =   qromb_rad_int (eval_i1,  0.0, rad_com%ll, i1)
    rad_com%i4a_(ir) = rad_com%i4a_(ir) + &
                         qromb_rad_int (eval_i4a, 0.0, rad_com%ll, i4a)
    rad_com%i4b_(ir) = rad_com%i4b_(ir) + &
                         qromb_rad_int (eval_i4b, 0.0, rad_com%ll, i4b)
    rad_com%i5a_(ir) =   qromb_rad_int (eval_i5a, 0.0, rad_com%ll, i5a)
    rad_com%i5b_(ir) =   qromb_rad_int (eval_i5b, 0.0, rad_com%ll, i5b)

  enddo

!---------------------------------------------------------------------
! now put everything together

  i1   = sum(rad_com%i1_(1:ring%n_ele_use))
  i2   = sum(rad_com%i2_(1:ring%n_ele_use))
  i3   = sum(rad_com%i3_(1:ring%n_ele_use))
  i4a  = sum(rad_com%i4a_(1:ring%n_ele_use))
  i4b  = sum(rad_com%i4b_(1:ring%n_ele_use))
  i5a  = sum(rad_com%i5a_(1:ring%n_ele_use))
  i5b  = sum(rad_com%i5b_(1:ring%n_ele_use))

  i4z = i4a + i4b

  energy = ring%param%energy
  gamma2_factor = (energy * 1956.95)**2
  energy_loss = c_gam * energy**4 * i2 / pi

  mode%synch_int(1) = i1
  mode%synch_int(2) = i2
  mode%synch_int(3) = i3

  mode%a%synch_int(4) = i4a
  mode%b%synch_int(4) = i4b
  mode%z%synch_int(4) = i4z

  mode%a%synch_int(5) = i5a
  mode%b%synch_int(5) = i5b

  if (i2 /= 0) then

    mode%a%emittance = c_q * gamma2_factor * i5a / (i2 - i4a)
    mode%b%emittance = c_q * gamma2_factor * i5b / (i2 - i4b)

    mode%a%j_damp = 1 - i4a / i2
    mode%b%j_damp = 1 - i4b / i2
    mode%z%j_damp = 2 + i4z / i2

    arg = (c_q * i3 * gamma2_factor / (2*i2 + i4z))
    mode%sig_e = sqrt(max(0.0, arg))

  endif

  mode%a%alpha_damp = energy_loss * mode%a%j_damp / energy
  mode%b%alpha_damp = energy_loss * mode%b%j_damp / energy
  mode%z%alpha_damp = energy_loss * mode%z%j_damp / energy

  mode%energy_loss = energy_loss

  if(abs(m65) > 0. ) then
    mode%sig_z = sqrt( mode%synch_int(1)/abs(m65) ) * mode%sig_e
  else
    mode%sig_z = 0.
  endif

end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

subroutine propagate_part_way (s)

  use rad_int_common
  use bmad_interface

  implicit none

  real s, v(4,4), v_inv(4,4)

!

  if (s == 0) then
    runt%x = ele0%x
    runt%y = ele0%y
    runt%c_mat = ele0%c_mat
    runt%gamma_c = ele0%gamma_c
    orb = orb0
  elseif (s == rad_com%ll) then
    runt = ele
    orb = orb1
  else
    runt = ele
    runt%value(l$) = s
    if (ele%key == sbend$) then
      runt%value(angle$) = ele%value(angle$) * s / rad_com%ll
      runt%value(e2$) = 0
    endif
    call track1 (orb0, runt, pring%param, orb)
    call make_mat6 (runt, pring%param, orb0, orb)
    call twiss_propagate1 (ele0, runt)
  endif

  call make_v_mats (runt, v, v_inv)

  rad_com%eta_a = &
          matmul(v, (/ runt%x%eta, runt%x%etap, 0.0,       0.0        /))
  rad_com%eta_b = &
          matmul(v, (/ 0.0,       0.0,        runt%y%eta, runt%y%etap /))
  rad_com%eta = rad_com%eta_a + rad_com%eta_b

  rad_com%g_x = rad_com%g_x0 + orb%x%pos * rad_com%k1 + orb%y%pos * rad_com%s1
  rad_com%g_y = rad_com%g_y0 - orb%y%pos * rad_com%k1 + orb%x%pos * rad_com%s1
                   
  rad_com%dg2_x = 2 * (rad_com%g_x * rad_com%k1 + rad_com%g_y * rad_com%s1)
  rad_com%dg2_y = 2 * (rad_com%g_x * rad_com%s1 - rad_com%g_y * rad_com%k1) 

  rad_com%g2 = rad_com%g_x**2 + rad_com%g_y**2
  rad_com%g = sqrt(rad_com%g2)

end subroutine
  
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

function  eval_i1 (s_vec)
                    
  use rad_int_common

  implicit none

  real s_vec(:)
  real, dimension(size(s_vec)) :: eval_i1

  integer i

!                      
                                         
  do i = 1, size(s_vec)
    call propagate_part_way (s_vec(i))
    eval_i1(i) = rad_com%g_x * rad_com%eta(1) + rad_com%g_y * rad_com%eta(3)
  enddo

end function
  
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

function  eval_i4a (s_vec)
                    
  use rad_int_common
                       
  implicit none

  real s_vec(:)
  real, dimension(size(s_vec)) :: eval_i4a

  integer i

!

  do i = 1, size(s_vec)
    call propagate_part_way (s_vec(i))
    eval_i4a(i) = rad_com%g2 * (rad_com%g_x * rad_com%eta_a(1) + rad_com%g_y * rad_com%eta_a(3)) + &
                  (rad_com%dg2_x * rad_com%eta_a(1) + rad_com%dg2_y * rad_com%eta_a(3))
  enddo

end function
  
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

function  eval_i4b (s_vec)
                    
  use rad_int_common

  implicit none

  real s_vec(:)
  real, dimension(size(s_vec)) :: eval_i4b

  integer i

!

  call propagate_part_way (s_vec)

  do i = 1, size(s_vec)
    call propagate_part_way (s_vec(i))
    eval_i4b(i) = rad_com%g2 * (rad_com%g_x * rad_com%eta_b(1) + rad_com%g_y * rad_com%eta_b(3)) + &
                  (rad_com%dg2_x * rad_com%eta_b(1) + rad_com%dg2_y * rad_com%eta_b(3))
  enddo

end function
  
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

function  eval_i5a (s_vec)
                    
  use rad_int_common

  implicit none

  real s_vec(:)
  real, dimension(size(s_vec)) :: eval_i5a

!

  integer i

  do i = 1, size(s_vec)
    call propagate_part_way (s_vec(i))
    eval_i5a(i) = rad_com%g2 * rad_com%g * (runt%x%gamma * runt%x%eta**2 + &
                      2 * runt%x%alpha * runt%x%eta * runt%x%etap + &
                      runt%x%beta * runt%x%etap**2)
  enddo             

end function

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

function  eval_i5b (s_vec)
                    
  use rad_int_common

  implicit none

  real s_vec(:)
  real, dimension(size(s_vec)) :: eval_i5b

  integer i

!

  do i = 1, size(s_vec)
    call propagate_part_way (s_vec(i))
    eval_i5b(i) = rad_com%g2 * rad_com%g * (runt%y%gamma * runt%y%eta**2 + &
                      2 * runt%y%alpha * runt%y%eta * runt%y%etap + &
                      runt%y%beta * runt%y%etap**2)
  enddo

end function

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
! This is a modified version of QROMB from Num. Rec.


function qromb_rad_int (func, a, b, sum)

  use nrtype; use nrutil, only : nrerror
  use nr, only : polint,trapzd

  implicit none

  interface
    function func(x)
    use nrtype
    real(sp), dimension(:), intent(in) :: x
    real(sp), dimension(size(x)) :: func
    end function func
  end interface

  integer(i4b), parameter :: jmax = 6, jmaxp = jmax+1, k = 5, km = k-1
  integer(i4b) :: j

  real(sp), intent(in) :: a, b, sum
  real(sp) :: qromb_rad_int
  real(sp) :: dqromb
  real(sp), parameter :: eps = 1.0e-4_sp, eps2 = 1.0e-6_sp
  real(sp), dimension(jmaxp) :: h, s

  logical :: debug = .false.

!

  h(1) = 1.0

  do j = 1, jmax
    call trapzd(func, a, b, s(j), j)
    if (j >=  k) then
      call polint(h(j-km:j), s(j-km:j), 0.0_sp, qromb_rad_int, dqromb)
      if (abs(dqromb) <= eps * abs(qromb_rad_int) + eps2 * abs(sum)) return
    elseif (j >= 3) then
      call polint(h(1:j), s(1:j), 0.0_sp, qromb_rad_int, dqromb)
      if (abs(dqromb) <= eps * abs(qromb_rad_int) + eps2 * abs(sum)) return
    end if
    s(j+1) = s(j)
    h(j+1) = 0.25_sp * h(j)
  end do

  if (debug) then
    type *, 'Warning in RADIATION_INTEGRALS: Integral does not converge.'
  endif

end function qromb_rad_int
