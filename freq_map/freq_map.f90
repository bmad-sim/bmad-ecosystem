!........................................................................
!+
! program    :  freq_map
!
! Description:  Compute frequency map of a lattice 
!
! Arguments  :
!
! Mod/Commons:
!
! Calls      :
!
! Author     : 
!
! Modified   :
!-
!........................................................................
!
!
! $Log: freq_map.f90,v $
! Revision 1.3  2007/01/30 16:14:31  dcs
! merged with branch_bmad_1.
!
!
!
!........................................................................
!
#include "CESR_platform.inc"

program freq_map

  use bmad
  use bmadz_interface
  use bsim_interface
  use nr
  use eigen_mod
  use nrutil, only: dpc
  use scan_parameters

  implicit none

  type (lat_struct) ring, ring_in
  type (coord_struct), allocatable :: co(:), orb(:), data(:)
  type (coord_struct) start_coord
  type (normal_modes_struct) mode
  type (scan_params_struct) scan_params

  integer i, j, k, isign, r, s, g
  integer ix
  integer n_turn, particle, i_train, j_car, n_trains_tot, n_cars, n_part, slices
  integer ios
  integer i_x, i_y, i_z
  integer n_x, n_y, n_z
  integer ix_pn, ix_dot
  integer integer_qx, integer_qy
  integer version
  integer i_dim
  integer fft_turns

  real(rp) phy_x_set, phy_y_set
  real(rp) delta_e/0.0/, chromx, chromy, qp_x, qp_y
  real(rp), allocatable :: dk1(:)
  real(rp) eps_rel(4), eps_abs(4)
  real(rp), allocatable :: tune(:,:,:)
  real(rp) Qx, Qy, Qz, current, Q, d, b, c, A
  real(rp) Qx_init, Qy_init
  real(rp) sig_in(3) !sigx, sigy, sigz, initial distribution
  real(rp) x0, y0, e0, x1, y1, e1, dx, dy, de
  real(rp), ALLOCATABLE :: fft1(:), fft2(:), fft3(:)
  real(rp), ALLOCATABLE :: fft4(:), fft5(:), fft6(:)
  complex(dpc), ALLOCATABLE ::  n1(:,:), m1(:,:),  n2(:,:), m2(:,:), n3(:,:), m3(:,:), n4(:,:)
  real(rp), ALLOCATABLE ::  rftxa(:), rftxb(:), rftya(:), rftyb(:), rftea(:), rfteb(:)
  real(rp) :: final_pos_in(1:4)
  real(rp) hanning
  type(coord_struct) :: final_pos_out

  real(rp) one_turn_mat(1:6,1:6)
  real(rp) eval_r(1:6), eval_i(1:6), evec_r(1:6,1:6), evec_i(1:6,1:6)
  real(rp) T(1:6,1:6), Tinv(1:6,1:6), P(1:6,1:6), Pinv(1:6,1:6)
  real(rp) psv(1:6), psv_norm(1:6)

  character*40 lattice, out_file_prefix
  character*180 lat_file, out_file
  character*80 line, last_line, file_name, prefix_name
  character*4 type
  character*1 answer
  character*60 in_file

  logical keep_trying/.true./, error
  logical write_orbit/.false./                                        
  logical beambeam_ip, close_pretz, close_vert, lrbbi, rec_taylor, go
  logical ok, aperture_limits

  namelist / parameters /lat_file, type, out_file_prefix, &
                Qx_init, Qy_init, Qx, Qy, Qz, Qp_x, Qp_y, &
                x0, y0, e0, x1, y1, e1, dx, dy, de, &
                n_turn, particle, aperture_limits,  &
                i_train, j_car, n_trains_tot, n_cars, current, &
                lrbbi, rec_taylor, beambeam_ip, close_pretz, close_vert, &
                slices, sig_in, go, final_pos_in

  do   !read file with input parameters
    type '(a, $)', ' Input command file <CR=freq_map.in>: '
    read  '(a)', file_name
    call string_trim (file_name, file_name, ix)
    if (ix .eq. 0) file_name = 'freq_map.in'
    open (unit= 1, file = file_name, status = 'old', iostat = ios)
    if (ios == 0) then
      exit
    else
      print *
      print *, 'ERROR: CANNOT OPEN FILE: ', trim(file_name)
    endif
  enddo

  particle = positron$
  lrbbi = .false.
  beambeam_ip = .false.
  close_pretz = .false.
  close_vert = .false.
  scan_params%final_pos_in%vec(1:4) = 0.0
  read(1, nml = parameters)
  
  scan_params%Q_z                   = Qz 
  scan_params%lat_file              = lat_file 
  scan_params%n_turn                = n_turn 
  scan_params%particle              = particle
  scan_params%i_train               = i_train 
  scan_params%j_car                 = j_car 
  scan_params%n_trains_tot          = n_trains_tot 
  scan_params%n_cars                = n_cars 
  scan_params%current               = current
  scan_params%lrbbi                 = lrbbi
  scan_params%beambeam_ip           = beambeam_ip 
  scan_params%close_pretz           = close_pretz
  scan_params%close_vert            = close_vert 
  scan_params%slices                = slices 
  scan_params%rec_taylor            = rec_taylor
  scan_params%sig_in                = sig_in
  scan_params%final_pos_in%vec(1:4) = final_pos_in(1:4)

  print *, ' lat_file = ', lat_file

  IF(type .eq. 'xsif') THEN
    CALL xsif_parser(lat_file, ring)
  ELSEIF(type .eq. 'bmad') THEN
    if(lat_file(1:8) == 'digested') then
      call read_digested_bmad_file(lat_file, ring, version)
    else
      call fullfilename(lat_file, lat_file)
      call bmad_parser(lat_file, ring)
    endif
  ELSE
    WRITE(*,*) "Error: unknown lattice type: ", type
    STOP
  ENDIF

  call reallocate_coord(orb, ring%n_ele_max)
  call reallocate_coord(co, ring%n_ele_max)

  if(go)keep_trying=.false.
  do while (keep_trying)
    type '(a, $)', ' FREQ_MAP: element change or GO> '
    read  '(a)', line
    ix = index(line, '!')
    if (ix /= 0) line = line(:ix-1)        ! strip off comments
    call str_upcase(line, line)
    call string_trim(line, line, ix)
    if (ix == 0) then       ! nothing typed. do the same thing
      line = last_line
    endif
    last_line = line
    if(line(1:1) .eq. 'G')exit
    call find_change( line, ring)
  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ring%param%aperture_limit_on = aperture_limits
  CALL reallocate_coord(orb,ring%n_ele_max)
  CALL set_on_off(rfcavity$, ring, on$)
  CALL twiss_and_track(ring,orb)
  CALL calc_z_tune(ring)
  ring%param%particle = particle
  
  WRITE(*,'(A,F6.3,A,F6.3,A,F6.3)') &
    "Original Tune: Qx=", ring%a%tune/twopi, " Qy=", ring%b%tune/twopi, " Qz=", ring%z%tune/twopi

  IF(Qx_init /= 0. .and. Qy_init /= 0.) THEN
    ALLOCATE(dk1(ring%n_ele_max))
    integer_qx = int(ring%ele(ring%n_ele_track)%a%phi / twopi)
    integer_qy = int(ring%ele(ring%n_ele_track)%b%phi / twopi)
    phy_x_set = (integer_qx + qx_init)*twopi
    phy_y_set = (integer_qy + qy_init)*twopi
    CALL choose_quads(ring, dk1)
    DO i=0, ring%n_ele_max
      orb(i)%vec(1:6) = 0.0
    ENDDO
    CALL set_on_off(elseparator$,ring,off$)
    CALL custom_set_tune(phy_x_set,phy_y_set,dk1,ring,orb,ok)
    IF(.not. ok) WRITE(*,*) "Qtune 1 failed"
    CALL set_on_off(elseparator$,ring,on$)
    CALL twiss_and_track(ring,orb)
    WRITE(*,*) "After QTune 1: Qx = ", ring%a%tune/twopi, "  Qy = ", ring%b%tune/twopi
    deallocate(dk1)
  ENDIF

  if(lrbbi)then
    call lrbbi_setup(ring, ring, particle, i_train, j_car, n_trains_tot, n_cars, current, rec_taylor)
    CALL twiss_and_track(ring,orb)
    WRITE(*,*) "Long range bb interaction on!"
  endif

  if(beambeam_ip) then
    call beambeam_setup(ring, particle, current, scan_params, slices)
    call twiss_and_track(ring,orb)
    WRITE(*,*) "Beambeam interation on!"
  endif

  IF(qz /= 0.) THEN
    ring%z%tune = qz * twopi
    CALL set_z_tune(ring)
    CALL twiss_and_track(ring,orb)
    WRITE(*,*) "After set_z_tune: Qz = ", ring%z%tune/twopi
  ENDIF

  i_dim = 6
  if(close_pretz) then
    call close_pretzel(ring,i_dim,scan_params%final_pos_in,final_pos_out)
    call twiss_and_track(ring,orb)
    WRITE(*,*) "Close pretzel on!"
  endif
  
  if(close_vert) then
    call close_vertical(ring,i_dim,scan_params%final_pos_in,final_pos_out)
    call twiss_and_track(ring,orb)
    WRITE(*,*) "Close vertical on!"
  endif

  IF(qx /= 0. .and. qy /= 0.) THEN
    integer_qx = int(ring%ele(ring%n_ele_track)%a%phi / twopi)
    integer_qy = int(ring%ele(ring%n_ele_track)%b%phi / twopi)
    phy_x_set = (integer_qx + qx)*twopi
    phy_y_set = (integer_qy + qy)*twopi
    allocate(dk1(ring%n_ele_max))
    CALL choose_quads(ring, dk1)
    CALL custom_set_tune(phy_x_set,phy_y_set,dk1,ring,orb,ok)
    IF(.not. ok) WRITE(*,*) "Qtune 2 failed"
    CALL twiss_and_track(ring,orb)
    WRITE(*,*) "After QTune 2: Qx = ", ring%a%tune/twopi, "  Qy = ", ring%b%tune/twopi
    deallocate(dk1)
  ENDIF

  CALL chrom_calc(ring,delta_e,chromx,chromy)
  WRITE(*,'(A,F6.3,A,F6.3)') "Initial Chromaticity: qp_x = ", chromx, "   qp_y = ", chromy
  IF(qp_x /= 0. .and. qp_y /= 0.) THEN
    CALL qp_tune(ring,qp_x,qp_y,ok)
    IF(.not. ok) WRITE(*,*) "qp_tune failed"
    CALL chrom_calc(ring,0.0_rp,chromx,chromy)
    CALL twiss_at_start(ring)
    WRITE(*,*) "After QPtune: qp_x = ", chromx, "   qp_y = ", chromy
    WRITE(*,*) "After QPtune: Qx = ", ring%a%tune/twopi, "  Qy = ", ring%b%tune/twopi
  ENDIF
  WRITE(*,'(A,F6.3,A,F6.3,A,F6.3)') &
    "Final Tune: Qx=", ring%a%tune/twopi, " Qy=", ring%b%tune/twopi, " Qz=", ring%z%tune/twopi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !This loop can be changed to set the default beampipe size
  do i = 1, ring%n_ele_max
    if(ring%ele(i)%value(x1_limit$) == 0) ring%ele(i)%value(x1_limit$) = 0
    if(ring%ele(i)%value(x2_limit$) == 0) ring%ele(i)%value(x2_limit$) = 0
    if(ring%ele(i)%value(y1_limit$) == 0) ring%ele(i)%value(y1_limit$) = 0
    if(ring%ele(i)%value(y2_limit$) == 0) ring%ele(i)%value(y2_limit$) = 0
  enddo

  n_x = nint(abs((x1 - x0)/dx))+1
  n_y = nint(abs((y1 - y0)/dy))+1
  if(de /= 0.)n_z = nint(abs((e1 - e0)/de))+1
  if(de == 0.)n_z=1

  allocate(tune(0:n_x*n_y*n_z, 3, 2))

  dx = (x1-x0)/(n_x-1)
  dy = (y1-y0)/(n_y-1)
  if(n_z /= 1)de = (e1-e0)/(n_z-1)

  ix_dot = index(file_name,'.')
  prefix_name = file_name(1:ix_dot-1)
  call string_trim(prefix_name, prefix_name,ix_pn)
  call file_suffixer (file_name, in_file, '.out', .true.)

  print '(/,1x,3(a6,e12.4,4x))',' x0 = ',x0, ' y0 = ',y0, ' e0 = ',e0
  print '(1x,3(a6,e12.4,4x))',' dx = ',dx, ' dy = ',dy, ' de = ',de
  print '(1x,3(a6,i4,6x,4x))','n_x = ',n_x,'n_y = ',n_y,'n_z = ',n_z
  print '(1x,(a9,i),/)','n_turn = ',n_turn

  ALLOCATE(fft1(1:n_turn))
  ALLOCATE(fft2(1:n_turn))
  ALLOCATE(fft3(1:n_turn))
  ALLOCATE(fft4(1:n_turn))
  ALLOCATE(fft5(1:n_turn))
  ALLOCATE(fft6(1:n_turn))

  fft_turns = n_turn / 2

  ALLOCATE(n1(1,1:fft_turns))
  ALLOCATE(n2(1,1:fft_turns))
  ALLOCATE(n3(1,1:fft_turns))
  ALLOCATE(m1(1,1:fft_turns))
  ALLOCATE(m2(1,1:fft_turns))
  ALLOCATE(m3(1,1:fft_turns))

  ALLOCATE(rftxa(1:fft_turns))
  ALLOCATE(rftxb(1:fft_turns))
  ALLOCATE(rftya(1:fft_turns))
  ALLOCATE(rftyb(1:fft_turns))
  ALLOCATE(rftea(1:fft_turns))
  ALLOCATE(rfteb(1:fft_turns))

  ! Multiplying coordinated by Tinv will convert them to normal
  ! coordinates.
  CALL transfer_matrix_calc(ring, .true., one_turn_mat)
  CALL mat_eigen (one_turn_mat, eval_r, eval_i, evec_r, evec_i, error)
  do i=1,6,2
    T(i,:)   =  evec_r(i,:)
    T(i+1,:) = -evec_i(i,:)
  enddo
  CALL mat_inverse(T,Tinv)

  g= 0
  do i_z = 0,n_z-1
    start_coord%vec(1:6) = 0.
    start_coord%vec(6) = e0 + de * i_z
    WRITE(out_file,'(2A,I3.3,A)') trim(out_file_prefix),".e",int(start_coord%vec(6)*1000),".fm"
    open(unit=13, file=out_file)
    do i_y = 0,n_y-1
      start_coord%vec(3) = y0 + dy * i_y
      do i_x = 0,n_x-1
        start_coord%vec(1) = x0 + dx * i_x
 
        co(0)%vec = orb(0)%vec + start_coord%vec

        do j=1,n_turn
          call track_all(ring, co)
          co(0)%vec = co(ring%n_ele_track)%vec
          if(ring%param%lost)then
            WRITE(*,*) "Particle lost in turn ", j
            exit
          else
            psv(1) = co(0)%vec(1) - orb(0)%vec(1)
            psv(2) = co(0)%vec(2) - orb(0)%vec(2)
            psv(3) = co(0)%vec(3) - orb(0)%vec(3)
            psv(4) = co(0)%vec(4) - orb(0)%vec(4)
            psv(5) = co(0)%vec(5) - orb(0)%vec(5)
            psv(6) = co(0)%vec(6) - orb(0)%vec(6)

            psv_norm = MATMUL(psv, Tinv)

            fft1(j) = psv_norm(1)
            fft2(j) = psv_norm(2)
            fft3(j) = psv_norm(3)
            fft4(j) = psv_norm(4)
            fft5(j) = psv_norm(5)
            fft6(j) = psv_norm(6)
          endif
        end do

        if(.not. ring%param%lost)then
          !interpolated FFT with Hanning window
          do j=1, fft_turns
            hanning = 2*(sin(pi*j/fft_turns)**2)
            n1(1, j)= cmplx(fft1(j),  fft2(j)) * hanning
            n2(1, j)= cmplx(fft3(j),  fft4(j)) * hanning
            n3(1, j)= cmplx(fft5(j), -fft6(j)) * hanning
          enddo

          do j=1, fft_turns
            hanning = 2*(sin(pi*j/fft_turns)**2)
            m1(1, j)= cmplx(fft1(j+fft_turns),  fft2(j+fft_turns)) * hanning
            m2(1, j)= cmplx(fft3(j+fft_turns),  fft4(j+fft_turns)) * hanning
            m3(1, j)= cmplx(fft5(j+fft_turns), -fft6(j+fft_turns)) * hanning
          enddo

          isign= 1

          call fourrow(n1(:,1:fft_turns), isign)    ! NR FFT
          forall(i=1:fft_turns)rftxa(i)=sqrt(n1(1,i)*conjg(n1(1,i)))

          call fourrow(m1(:,1:fft_turns), isign)
          forall(i=1:fft_turns)rftxb(i)=sqrt(m1(1,i)*conjg(m1(1,i)))

          call fourrow(n2(:,1:fft_turns), isign)
          forall(i=1:fft_turns)rftya(i)=sqrt(n2(1,i)*conjg(n2(1,i)))

          call fourrow(m2(:,1:fft_turns), isign)
          forall(i=1:fft_turns)rftyb(i)=sqrt(m2(1,i)*conjg(m2(1,i)))

          call fourrow(n3(:,1:fft_turns), isign)
          forall(i=1:fft_turns)rftea(i)=sqrt(n3(1,i)*conjg(n3(1,i)))

          call fourrow(m3(:,1:fft_turns), isign)
          forall(i=1:fft_turns)rfteb(i)=sqrt(m3(1,i)*conjg(m3(1,i)))

          !tune of x for 1st half etc.
          Q= 1
          do k= 2, fft_turns
            if (rftxa(Q)<rftxa(k)) Q= k
          end do
          d= rftxa(Q)
          b= rftxa(Q+1)
          c= cos(twopi/fft_turns)
          A= (-(d+b*c)*(d-b)+b*sqrt(c*c*(d+b)*(d+b)-2*d*b*(2*c*c-c-1)))/(d*d+b*b+2*d*b*c)
          tune(g, 1, 1)= Q/fft_turns+(1/twopi)*asin(A*sin(twopi/fft_turns))
 
          !tune of x for 2nd half etc.
          Q= 1
          do k= 2, fft_turns
            if (rftxb(Q)<rftxb(k)) Q= k
          end do
          d= rftxb(Q)
          b= rftxb(Q+1)
          c= cos(twopi/fft_turns)
          A= (-(d+b*c)*(d-b)+b*sqrt(c*c*(d+b)*(d+b)-2*d*b*(2*c*c-c-1)))/(d*d+b*b+2*d*b*c)
          tune(g, 1, 2)= Q/fft_turns+(1/twopi)*asin(A*sin(twopi/fft_turns))
   
          !find tune for y first half etc.
          Q= 1
          do k= 2, fft_turns
            if (rftya(Q)<rftya(k)) Q= k
          end do
          d= rftya(Q)
          b= rftya(Q+1)
          c= cos(twopi/fft_turns)
          A= (-(d+b*c)*(d-b)+b*sqrt(c*c*(d+b)*(d+b)-2*d*b*(2*c*c-c-1)))/(d*d+b*b+2*d*b*c)
          tune(g, 2, 1)= Q/fft_turns+(1/twopi)*asin(A*sin(twopi/fft_turns))
   
          !find tune for y second half etc.
          Q= 1
          do k= 2, fft_turns
            if (rftyb(Q)<rftyb(k)) Q= k
          end do
          d= rftyb(Q)
          b= rftyb(Q+1)
          c= cos(twopi/fft_turns)
          A= (-(d+b*c)*(d-b)+b*sqrt(c*c*(d+b)*(d+b)-2*d*b*(2*c*c-c-1)))/(d*d+b*b+2*d*b*c)
          tune(g, 2, 2)= Q/fft_turns+(1/twopi)*asin(A*sin(twopi/fft_turns))
 
          !find tune for z first half etc.
          Q= 1
          do k= 2, fft_turns
            if (rftea(Q)<rftea(k)) Q= k
          end do
          d= rftea(Q)
          b= rftea(Q+1)
          c= cos(twopi/fft_turns)
          A= (-(d+b*c)*(d-b)+b*sqrt(c*c*(d+b)*(d+b)-2*d*b*(2*c*c-c-1)))/(d*d+b*b+2*d*b*c)
          tune(g, 3, 1)= Q/fft_turns+(1/twopi)*asin(A*sin(twopi/fft_turns))
 
          !find tune for z second half etc.
          Q= 1
          do k= 2, fft_turns
            if (rfteb(Q)<rfteb(k)) Q= k
          end do
          d= rfteb(Q)
          b= rfteb(Q+1)
          c= cos(twopi/fft_turns)
          A= (-(d+b*c)*(d-b)+b*sqrt(c*c*(d+b)*(d+b)-2*d*b*(2*c*c-c-1)))/(d*d+b*b+2*d*b*c)
          tune(g, 3, 2)= Q/fft_turns+(1/twopi)*asin(A*sin(twopi/fft_turns))
 
          write(13, '(12e14.5)') start_coord%vec(1), start_coord%vec(3), start_coord%vec(6), &
                                tune(g, :, 1), tune(g, :, 2), tune(g,:,2)-tune(g,:,1)
        endif !.not.ring%param%lost
        ring%param%lost = .false.
      end do !i_x
    end do !i_y
    close(unit=13)
  end do !i_z

  DEALLOCATE(fft1)
  DEALLOCATE(fft2)
  DEALLOCATE(fft3)
  DEALLOCATE(fft4)
  DEALLOCATE(fft5)
  DEALLOCATE(fft6)
  DEALLOCATE(n1)
  DEALLOCATE(n2)
  DEALLOCATE(n3)
  DEALLOCATE(m1)
  DEALLOCATE(m2)
  DEALLOCATE(m3)
  DEALLOCATE(rftxa)
  DEALLOCATE(rftxb)
  DEALLOCATE(rftya)
  DEALLOCATE(rftyb)
  DEALLOCATE(rftea)
  DEALLOCATE(rfteb)
end













