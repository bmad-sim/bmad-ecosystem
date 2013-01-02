module reverse_mod

use bookkeeper_mod

!! private ele_reverse

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
!     lat_rev%ele(lat%n_ele_track+1-i) = lat_in%ele(i)  For 0 < i <= lat%n_ele_track
!     lat_rev%ele(i)                   = lat_in%ele(i)  For lat%n_ele_track < i 
!
! The transformation from particle coordinates
! in one lat to the corrisponding coordinates in the other is:
!     (x, P_x, y, P_y, z, P_z) -> (x, -P_x, -y, P_y, -z, P_z)
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
!               The lat_rev actual argument may not be the same as the lat_in actual argument.
!-

subroutine lat_reverse (lat_in, lat_rev)

implicit none

type (lat_struct), target :: lat_in, lat_rev
type (ele_struct), pointer :: lord, ele
type (control_struct), pointer :: con
type (branch_struct), pointer :: branch, branch_in

integer i, n, i1, i2, nr, n_con, ib
integer :: ix_con(size(lat_in%control))

logical err_flag

! Correct control information

lat_rev = lat_in

do i = 1, lat_rev%n_control_max
  con => lat_rev%control(i)
  if (con%ix_slave <= lat_rev%n_ele_track) con%ix_slave = lat_rev%n_ele_track+1-con%ix_slave
  if (con%ix_lord <= lat_rev%n_ele_track)  con%ix_lord  = lat_rev%n_ele_track+1-con%ix_lord
enddo

! Slaves of a super lord must be in assending sequence.
! ix_con keeps track of the switching.
! Also: adjust s-position of lords.

n_con = size(lat_in%control)
forall (i = 1:n_con) ix_con(i) = i 

do i = lat_rev%n_ele_track+1, lat_rev%n_ele_max
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

! Transfer info from lat_in to lat_rev.
! the lat lattice is used since the actual arguments of lat_in and lat_rev
! may be the same

do ib = 0, ubound(lat_in%branch, 1)

  branch => lat_rev%branch(ib)
  branch_in => lat_in%branch(ib)

  nr = branch%n_ele_track
  branch%ele(1:nr) = branch_in%ele(nr:1:-1)

  ! Flip longitudinal stuff, maps

  do i = 1, branch%n_ele_max
    ele => branch%ele(i)
    call ele_reverse (ele, branch%param)
    if (i <= nr) then
      ele%s = branch%param%total_length - (ele%s - ele%value(l$))
      ele%ref_time = branch_in%ele(nr)%ref_time - branch_in%ele(nr-i)%ref_time
    endif
  enddo

  ! Cleanup

  branch%param%t1_with_RF = 0  ! Init
  branch%param%t1_no_RF = 0    ! Init

enddo

! Finish

call lat_sanity_check (lat_rev, err_flag)
call set_ele_status_stale (lat_rev%ele(0), floor_position_group$)
call lattice_bookkeeper (lat_rev)

end subroutine lat_reverse

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine ele_reverse (ele, param)
!
! Subroutine to "reverse" an element for backward tracking.
!
! All longitudial quantities (for example, the value of ks for a solenoid) 
! are flipped in sign for the reversed lat. 
! This means that the appropriate transformation that corresponds to the 
! reverse transformation is:
!     (x, P_x, y, P_y, z, P_z) -> (x, -P_x, y, -P_y, -z, P_z)
!
! Note: Due to complications occuring when you have, for example, super_slave 
! em_field elements, this routine is private and cannot be called directly.
! Use lat_reverse instead.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele   -- Ele_struct: Input element.
!   param -- Lat_param_struct: Lattice parameters.
!
! Output:
!   ele -- Ele_struct: Reversed element.
!-

subroutine ele_reverse (ele, param)

use ptc_interface_mod, only: taylor_inverse

implicit none

type (ele_struct) ele
type (lat_param_struct) param

integer i, j, sum235

real(rp) temp(6)
real(rp) tempp

! Flip map_ref coords

temp = ele%map_ref_orb_out
ele%map_ref_orb_out = ele%map_ref_orb_in
ele%map_ref_orb_in = temp

ele%map_ref_orb_in(2) = -ele%map_ref_orb_in(2)
ele%map_ref_orb_in(3) = -ele%map_ref_orb_in(3)
ele%map_ref_orb_in(5) = -ele%map_ref_orb_in(5)

ele%map_ref_orb_out(2) = -ele%map_ref_orb_out(2)
ele%map_ref_orb_out(3) = -ele%map_ref_orb_out(3)
ele%map_ref_orb_out(5) = -ele%map_ref_orb_out(5)

! Flip aperture limit position

if (ele%aperture_at == upstream_end$) then
  ele%aperture_at = downstream_end$
elseif (ele%aperture_at == downstream_end$) then
  ele%aperture_at = upstream_end$
endif

! Flip coupler limit position

if (ele%key == lcavity$) then
  if (nint(ele%value(coupler_at$)) == upstream_end$) then
    ele%value(coupler_at$) = downstream_end$
  elseif (nint(ele%value(coupler_at$)) == downstream_end$) then
    ele%value(coupler_at$) = upstream_end$
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
!       phi_z -> -phi_z - k_z * Length
!       coef  -> -coef
! This transforms:
!       (B_x, B_y, B_z) @ (x, y, s) -> (B_x, -B_y, -B_z) @ (x, -y, L-s)

case (wiggler$)
  if (associated(ele%wig)) then
    do i = 1, size(ele%wig%term)
      ele%wig%term(i)%phi_z = -ele%wig%term(i)%phi_z - ele%wig%term(i)%kz * ele%value(l$)
      ele%wig%term(i)%coef = -ele%wig%term(i)%coef
    enddo
  endif

end select

! Since, for example, the wiggler trajectory starting on the origin may not end
!   on the origin, the reference transit time may shift and needs to be recalculated.

if (.not. ele_has_constant_ds_dt_ref(ele)) then
  !...
endif

! The inverse of the Taylor map is: 
!   M * T^(-1) * M
! where: 
!   T = the input taylor series
!   M = The map (x, P_x, y, P_y, z, P_z) -> (x, -P_x, -y, P_y, -z, P_z)

if (associated(ele%taylor(1)%term)) then

  call taylor_inverse (ele%taylor, ele%taylor)

  ! Apply M to the right

  do i = 1, 6
    do j = 1, size(ele%taylor(i)%term)
      sum235 = ele%taylor(i)%term(j)%expn(2) + ele%taylor(i)%term(j)%expn(3) + &
                                                ele%taylor(i)%term(j)%expn(5)
      if (mod(sum235, 2) == 1) ele%taylor(i)%term(j)%coef = -ele%taylor(i)%term(j)%coef
    end do
  end do

  ! Apply M to the left

  ele%taylor(2)%term(:)%coef = -ele%taylor(2)%term(:)%coef
  ele%taylor(3)%term(:)%coef = -ele%taylor(3)%term(:)%coef
  ele%taylor(5)%term(:)%coef = -ele%taylor(5)%term(:)%coef

endif

! kill any ptc_genfield

if (associated(ele%ptc_genfield)) call kill_ptc_genfield (ele%ptc_genfield)

! reverse mat6

call mat_symp_conj (ele%mat6, ele%mat6)
ele%vec0 = -matmul(ele%mat6, ele%vec0)

ele%mat6(2,:) = -ele%mat6(2,:)
ele%mat6(3,:) = -ele%mat6(3,:)
ele%mat6(5,:) = -ele%mat6(5,:)

ele%vec0(2) = -ele%vec0(2)
ele%vec0(3) = -ele%vec0(3)
ele%vec0(5) = -ele%vec0(5)

ele%mat6(:,2) = -ele%mat6(:,2)
ele%mat6(:,3) = -ele%mat6(:,3)
ele%mat6(:,5) = -ele%mat6(:,5)

end subroutine ele_reverse

end module


