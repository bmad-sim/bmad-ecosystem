#include "CESR_platform.inc"

module reverse_mod

  use bookkeeper_mod
  use bmad_struct
  use bmad_interface

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine lat_reverse (lat_in, lat_rev)
!
! Subroutine to construct a lat structure with the elements in reversed order.
! This may be used for backward tracking through the lat. 
!
! The correspondence between elements in the two lattices is as follows:
!     lat_rev%ele(lat%n_ele_track+1-i) = lat_in%ele(i)  
!                                                for 0 < i <= lat%n_ele_track
!     lat_rev%ele(i)                   = lat_in%ele(i)   
!                                                for lat%n_ele_track < i 
!
! All longitudial quantities (for example, the value of ks for a solenoid) 
! are flipped in sign for the reversed lat. 
! This means that the appropriate transformation from particle coordinates
! in one lat the corrisponding coordinates in the other is:
!     (x, P_x, y, P_y, z, P_z) -> (x, -P_x, y, -P_y, -z, P_z)
!
! Note: The Twiss parameters will not be correct for the reversed lat.
! You will need to compute them.
!
! Modules needed:
!   use bmad
!
! Input:
!   lat_in -- lat_struct: Input lat.
!
! Output:
!   lat_rev -- lat_struct: Lat with the elements in reversed order.
!-

subroutine lat_reverse (lat_in, lat_rev)

  implicit none

  type (lat_struct), intent(in) :: lat_in
  type (lat_struct), intent(out), target :: lat_rev
  type (lat_struct), save :: lat
  type (ele_struct), pointer :: lord, ele
  type (control_struct), pointer :: con

  integer i, n, i1, i2, nr, n_con
  integer :: ix_con(size(lat_in%control))

! Transfer info from lat_in to lat_rev.
! the lat lattice is used since the actual arguments of lat_in and lat_rev
! may be the same

  n_con = size(lat_in%control)

  lat = lat_in 
  lat_rev = lat

  nr = lat_rev%n_ele_track
  lat_rev%ele(1:nr) = lat%ele(nr:1:-1)

! Flip longitudinal stuff, maps

  do i = 1, lat_rev%n_ele_max
    ele => lat_rev%ele(i)
    call reverse_ele (ele, lat%param)
    if (i <= nr) ele%s = lat_rev%param%total_length - (ele%s - ele%value(l$))
  enddo

! Correct control information

  do i = 1, lat_rev%n_control_max
    con => lat_rev%control(i)
    if (con%ix_slave <= nr) con%ix_slave = nr+1-con%ix_slave
    if (con%ix_lord <= nr)  con%ix_lord  = nr+1-con%ix_lord
  enddo

! Slaves of a super lord must be in assending sequence.
! ix_con keeps track of the switching.
! Also: adjust s-position of lords.

  forall (i = 1:n_con) ix_con(i) = i 

  do i = nr+1, lat_rev%n_ele_max
    lord => lat_rev%ele(i)
    if (lord%lord_status /= super_lord$) cycle
    i1 = lord%ix1_slave
    i2 = lord%ix2_slave
    lat_rev%control(i1:i2) = lat_rev%control(i2:i1:-1)
    ix_con(i1:i2) = ix_con(i2:i1:-1)
    if (lord%s > lord%value(l$)) then
      lord%s = lat_rev%param%total_length - (lord%s - lord%value(l$))
    else  ! Lord wraps around zero case
      lord%s = lord%value(l$) - lord%s
    endif
  enddo

  n = lat_rev%n_ic_max
  lat_rev%ic(1:n) = ix_con(lat_rev%ic(1:n))

! Cleanup

  lat_rev%param%t1_with_RF = 0  ! Init
  lat_rev%param%t1_no_RF = 0    ! Init

  call check_lat_controls (lat_rev, .true.)
  call lattice_bookkeeper (lat_rev)

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine reverse_ele (ele, param)
!
! Subroutine to "reverse" an element for backward tracking.
!
! All longitudial quantities (for example, the value of ks for a solenoid) 
! are flipped in sign for the reversed lat. 
! This means that the appropriate transformation that corresponds to the 
! reverse transformation is:
!     (x, P_x, y, P_y, z, P_z) -> (x, -P_x, y, -P_y, -z, P_z)
!
! Modules needed:
!   use bmad
!
! Input:
!   ele   -- Ele_struct: Input element.
!   param -- Lat_param_struct: Lattice parameters
!
! Output:
!   ele -- Ele_struct: Reversed element.
!-

subroutine reverse_ele (ele, param)

  use ptc_interface_mod, only: taylor_inverse

  implicit none

  type (ele_struct) ele
  type (lat_param_struct) param
  type (coord_struct) temp

  integer i, j, sum245

  real(rp) tempp

! Flip map_ref coords

  temp = ele%map_ref_orb_out
  ele%map_ref_orb_out = ele%map_ref_orb_in
  ele%map_ref_orb_in = temp

  ele%map_ref_orb_in%vec(2) = -ele%map_ref_orb_in%vec(2)
  ele%map_ref_orb_in%vec(4) = -ele%map_ref_orb_in%vec(4)
  ele%map_ref_orb_in%vec(5) = -ele%map_ref_orb_in%vec(5)

  ele%map_ref_orb_out%vec(2) = -ele%map_ref_orb_out%vec(2)
  ele%map_ref_orb_out%vec(4) = -ele%map_ref_orb_out%vec(4)
  ele%map_ref_orb_out%vec(5) = -ele%map_ref_orb_out%vec(5)

! Flip aperture limit position

  if (ele%aperture_at == entrance_end$) then
    ele%aperture_at = exit_end$
  elseif (ele%aperture_at == exit_end$) then
    ele%aperture_at = entrance_end$
  endif

! Flip coupler limit position

  if (ele%key == lcavity$) then
    if (nint(ele%value(coupler_at$)) == entrance_end$) then
      ele%value(coupler_at$) = exit_end$
    elseif (nint(ele%value(coupler_at$)) == exit_end$) then
      ele%value(coupler_at$) = entrance_end$
    endif
  endif

! Flip longitudinal attributes

  ele%value(x_pitch$) = -ele%value(x_pitch$)
  ele%value(y_pitch$) = -ele%value(y_pitch$)

  ele%value(x_pitch_tot$) = -ele%value(x_pitch_tot$)
  ele%value(y_pitch_tot$) = -ele%value(y_pitch_tot$)

  select case (ele%key)

  case (solenoid$, sol_quad$)
    ele%value(ks$) = -ele%value(ks$)

  case (lcavity$)
    ele%value(gradient$)     = -ele%value(gradient$) 
    ele%value(gradient_err$) = -ele%value(gradient_err$) 

  case (sbend$)
    tempp = ele%value(e1$)
    ele%value(e1$) = ele%value(e2$)
    ele%value(e2$) = tempp
    tempp = ele%value(h1$)
    ele%value(h1$) = ele%value(h2$)
    ele%value(h2$) = tempp
    tempp = ele%value(fint$)
    ele%value(fint$)  = ele%value(fintx$)
    ele%value(fintx$) = tempp 
    tempp = ele%value(hgap$)
    ele%value(hgap$)  = ele%value(hgapx$)
    ele%value(hgapx$) = tempp

! For wigglers:
!       phi_z -> -phi_z -  k_z * Length
! This transforms:
!       (B_x, B_y, B_z) @ (s) -> (B_x, B_y, -B_z) @ (L-s)
! Also: Since the wiggler trajectory starting on the origin may not end
!   on the origin, z_patch may shift. Zero z_patch so that there are
!   no shifts in tracking when remaking the element.

  case (wiggler$)
    if (associated(ele%wig_term)) then
      do i = 1, size(ele%wig_term)
        ele%wig_term(i)%phi_z = -ele%wig_term(i)%phi_z - &
                                        ele%wig_term(i)%kz * ele%value(l$)
      enddo
      ele%value(z_patch$) = 0
      call attribute_bookkeeper (ele, param)
    endif

  end select

! The inverse of the Taylor map is: 
!   M * T^(-1) * M
! where: 
!   T = the input taylor series
!   M = The map (x, P_x, y, P_y, z, P_z) -> (x, -P_x, y, -P_y, -z, P_z)

  if (associated(ele%taylor(1)%term)) then

    call taylor_inverse (ele%taylor, ele%taylor)

! Apply M to the right

    do i = 1, 6
      do j = 1, size(ele%taylor(i)%term)
        sum245 = ele%taylor(i)%term(j)%exp(2) + ele%taylor(i)%term(j)%exp(4) + &
                                                  ele%taylor(i)%term(j)%exp(5)
        if (mod(sum245, 2) == 1) ele%taylor(i)%term(j)%coef = &
                                                    -ele%taylor(i)%term(j)%coef
      end do
    end do

! Apply M to the left

    ele%taylor(2)%term(:)%coef = -ele%taylor(2)%term(:)%coef
    ele%taylor(4)%term(:)%coef = -ele%taylor(4)%term(:)%coef
    ele%taylor(5)%term(:)%coef = -ele%taylor(5)%term(:)%coef

  endif

! kill any gen_field

  if (associated(ele%gen_field)) call kill_gen_field (ele%gen_field)

! reverse mat6

  call mat_symp_conj (ele%mat6, ele%mat6)
  ele%vec0 = -matmul(ele%mat6, ele%vec0)

  ele%mat6(2,:) = -ele%mat6(2,:)
  ele%mat6(4,:) = -ele%mat6(4,:)
  ele%mat6(5,:) = -ele%mat6(5,:)

  ele%vec0(2) = -ele%vec0(2)
  ele%vec0(4) = -ele%vec0(4)
  ele%vec0(5) = -ele%vec0(5)

  ele%mat6(:,2) = -ele%mat6(:,2)
  ele%mat6(:,4) = -ele%mat6(:,4)
  ele%mat6(:,5) = -ele%mat6(:,5)

end subroutine

end module


