module linear_aperture_mod

implicit none

contains

subroutine linear_aperture(ring,da_config)
  !currently only works for round apertures
  use bmad
  use custom_dynamic_aperture_mod, only: custom_aperture_scan_struct, local_cosphi, local_sinphi
  use transfer_map_mod

  implicit none

  type(lat_struct) ring
  type(custom_aperture_scan_struct) da_config

  integer i, k, j
  real(rp) tmat(6,6), t1(6,6)
  real(rp) tvec(6), v1(6)
  real(rp) delta_angle
  real(rp) dE, theta, pax, pay
  real(rp) B1, B2, C1, C2
  real(rp) alpha, beta, gamma
  real(rp) sinphi1, cosphi1
  real(rp) sinphi2, cosphi2
  real(rp) pl, pl1, pl2
  real(rp) vec_length
  logical mask_x(1:ring%n_ele_track)
  logical mask_y(1:ring%n_ele_track)

  delta_angle = (da_config%max_angle - da_config%min_angle)/(da_config%n_angle -1)
  dE = da_config%param%closed_orbit%vec(6)
  pax = ring%ele(1)%value(x1_limit$)
  pay = ring%ele(1)%value(y1_limit$)

  mask_x = abs(ring%ele(1:ring%n_ele_track)%value(x1_limit$)) > 0.0001
  mask_y = abs(ring%ele(1:ring%n_ele_track)%value(y1_limit$)) > 0.0001
  da_config%Sx = minval( ring%ele(1:ring%n_ele_track)%value(x1_limit$) / sqrt(ring%ele(1:ring%n_ele_track)%a%beta), mask_x ) * sqrt(ring%ele(1)%a%beta)
  da_config%Sy = minval( ring%ele(1:ring%n_ele_track)%value(y1_limit$) / sqrt(ring%ele(1:ring%n_ele_track)%b%beta), mask_y ) * sqrt(ring%ele(1)%b%beta)

  !initialize to physical aperture at track point
  do i=1,da_config%n_angle
    theta = (i-1)*delta_angle + da_config%min_angle
    vec_length = sqrt(pax**2 * cos(theta)**2 + pay**2 * sin(theta)**2)
    da_config%aperture(i)%x = vec_length * local_cosphi(theta,da_config%Sx,da_config%Sy)
    da_config%aperture(i)%y = vec_length * local_sinphi(theta,da_config%Sx,da_config%Sy)
    ! vec_length = sqrt(pax**2 * cos(theta)**2 + pay**2 * sin(theta)**2)
    ! da_config%aperture(i)%x = vec_length * cos(theta)
    ! da_config%aperture(i)%y = vec_length * sin(theta)
    da_config%aperture(i)%i_turn = 1  !this routine does not calculate the turn number at which the particle was lost.
  enddo

  do i=1,ring%n_ele_track
    pax = ring%ele(i)%value(x1_limit$)
    pay = ring%ele(i)%value(y1_limit$)

    call transfer_matrix_calc (ring, t1, v1, ix1=i, one_turn=.true.)
    call transfer_matrix_calc (ring, tmat, tvec, ix1=0, ix2=i, one_turn=.true.)
    do j=1,100  !simulate 100 turns
      if (j>1) then
        call concat_transfer_mat (t1, v1, tmat, tvec, tmat, tvec)
      endif
      do k=1,da_config%n_angle
        theta = (k-1)*delta_angle + da_config%min_angle

        B1 = tmat(1,1)*local_cosphi(theta,da_config%Sx,da_config%Sy)+tmat(1,3)*local_sinphi(theta,da_config%Sx,da_config%Sy)
        B2 = tmat(3,1)*local_cosphi(theta,da_config%Sx,da_config%Sy)+tmat(3,3)*local_sinphi(theta,da_config%Sx,da_config%Sy)
        ! B1 = tmat(1,1)*cos(theta)+tmat(1,3)*sin(theta)
        ! B2 = tmat(3,1)*cos(theta)+tmat(3,3)*sin(theta)
        C1 = tmat(1,6)*dE+tvec(1)
        C2 = tmat(3,6)*dE+tvec(3)

        ! ****** Calculate downstream angle
        alpha = B1*pay
        gamma = -B2*pax
        beta = B1*C2-B2*C1
        cosphi1 = (beta*gamma-sqrt(alpha*alpha*(alpha*alpha-beta*beta+gamma*gamma)))/(alpha*alpha+gamma*gamma)
        cosphi2 = (beta*gamma+sqrt(alpha*alpha*(alpha*alpha-beta*beta+gamma*gamma)))/(alpha*alpha+gamma*gamma)
        sinphi1 = (alpha*alpha*beta+gamma*sqrt(alpha*alpha*(alpha*alpha-beta*beta+gamma*gamma)))/(alpha**3+alpha*gamma*gamma)
        sinphi2 = (alpha*alpha*beta-gamma*sqrt(alpha*alpha*(alpha*alpha-beta*beta+gamma*gamma)))/(alpha**3+alpha*gamma*gamma)
        ! ******

        pl1 = (pax*cosphi1+pay*sinphi1-C1-C2)/(B1+B2)
        pl2 = (pax*cosphi2+pay*sinphi2-C1-C2)/(B1+B2)

        if( pl1 .gt. 0.0d0 ) then 
          pl = pl1
        else
          pl = pl2
        endif

        if( pl**2 .lt. da_config%aperture(k)%x**2+da_config%aperture(k)%y**2 ) then
          da_config%aperture(k)%x = abs(pl) * local_cosphi(theta,da_config%Sx,da_config%Sy)
          da_config%aperture(k)%y = abs(pl) * local_sinphi(theta,da_config%Sx,da_config%Sy)
          ! da_config%aperture(k)%x = abs(pl) * cos(theta)
          ! da_config%aperture(k)%y = abs(pl) * sin(theta)
        endif
      enddo
    enddo
  enddo
end subroutine

end module
