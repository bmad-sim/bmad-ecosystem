program frequency_map

  use bmad
  use bsim_interface
  use nr, only: fourrow
  use z_tune_mod
  
  implicit none
  
type (lat_struct), target :: ring
type (coord_struct), allocatable :: orb(:), co(:), data(:)
type (coord_struct) start_coord
type (ele_struct), pointer :: ele
type (twiss_struct)  twiss1, twiss2

integer i, j, k, s, n, Npts, asum
integer ix
integer n_turn
integer ios, track_state
integer i_z, iq
integer i_a, i_b, n_a, n_b
integer n_z
integer ix_pn, ix_dot
integer version
integer :: fft_turns = 0
integer arg_num, iargc, ind
integer Nthmin, Nrad
integer stat
integer integer_qx, integer_qy

real(rp) rad, th, Ar, Br, Rmax
real(rp) phy_x_set, phy_y_set
real(rp) tune(3,2)
real(rp) x0, y0, e0, x1, y1, e1 
real(rp) :: dx=0., dy=0., de=0.
real(rp) :: rf_frequency
real(rp), ALLOCATABLE :: fft1(:), fft2(:), fft3(:)
real(rp), ALLOCATABLE :: fft4(:), fft5(:), fft6(:)
complex(rp), ALLOCATABLE ::  n1(:,:), m1(:,:),  n2(:,:), m2(:,:), n3(:,:), m3(:,:), n4(:,:)
real(rp), ALLOCATABLE ::  rftxa(:), rftxb(:), rftya(:), rftyb(:), rftea(:), rfteb(:)
real(rp), allocatable :: grid_pts(:,:)
real(rp) hanning
real(rp) U(4,4), V(4,4), Ubar(4,4), Vbar(4,4), G(4,4), tgamma
real(rp) Vinv(4,4)
real(rp) Qx, Qy, Qz, target_tunes(3)

real(rp) one_turn_mat(1:6,1:6)
real(rp) T(1:6,1:6), Tinv(1:6,1:6), P(1:6,1:6), Pinv(1:6,1:6)
real(rp) psv(1:6), psv_norm(1:6)

complex(rp) eigen_val(6), eigen_vec(6,6)

character(40) lattice, out_file_prefix
character(40) :: qtune_mask = ''
character(180) lat_file, out_file
character(80) line, last_line, file_name, prefix_name
character(1) answer
character(60) in_file
character(2) grid_type ! 'xy' or 'rt' for rectangular or polar

logical keep_trying/.true./, error
logical write_orbit/.false./
logical ok, aperture_limits, err
logical :: z_cut = .true.

namelist / parameters /lat_file, out_file_prefix, grid_type, &
     x0, y0, e0, x1, y1, e1, dx, dy, de, qx, qy, qz, qtune_mask, &
     n_turn, aperture_limits, fft_turns, z_cut, Rmax, Nrad, Br, Nthmin

arg_num=iargc()
if(arg_num==0) then
   file_name='freq_map.in'
else
   call getarg(1, file_name)
end if
call string_trim (file_name, file_name, ix)
open (unit= 1, file = file_name, status = 'old', iostat = ios)

if(ios.ne.0)then
   print *
   print *, 'ERROR: CANNOT OPEN FILE: ', trim(file_name)
endif

grid_type = 'xy'
qx = 0
qy = 0
Nrad = 17
Nthmin = 11
Br = 8.0
read(1, nml = parameters)

print *, ' lat_file = ', lat_file

call fullfilename(lat_file, lat_file)
call bmad_parser(lat_file, ring)



if (z_cut) then
  do ix = 1, ring%n_ele_track
     if (ring%ele(ix)%key == key_name_to_key_index('RFCAVITY',.true.)) then
        ! note: picks up first cavity; will not work if first cavity is a 3rd-harmonic, for example!
        rf_frequency = value_of_attribute(ring%ele(ix),'RF_FREQUENCY', err)
        exit
     endif
  enddo
else
   rf_frequency = 1.e-6 ! something small; prevents longitudinal cuts until much larger
endif

call reallocate_coord(co, ring%n_ele_max)
call reallocate_coord(orb, ring%n_ele_max)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
bmad_com%aperture_limit_on = aperture_limits
CALL reallocate_coord(co,ring%n_ele_max)
CALL set_on_off(rfcavity$, ring, on$)
CALL twiss_and_track(ring,co)
CALL calc_z_tune(ring)

IF(qx /= 0. .and. qy /= 0.) THEN

   integer_qx = floor(ring%ele(ring%n_ele_track)%a%phi / (2.*pi))
   integer_qy = floor(ring%ele(ring%n_ele_track)%b%phi / (2.*pi))
   target_tunes(1) = (integer_qx + qx)
   target_tunes(2) = (integer_qy + qy)
   
   if(abs(qz) > 0.) then
      target_tunes(3) = qz
   else
      target_tunes(3) = 0.
   endif
   
   call set_tune3(ring, target_tunes, ok, .false., regex_in=qtune_mask)  !takes tunes that have not not been multiplied by 2pi.

   IF(.not. ok) WRITE(*,*) "Qtune failed"
   CALL twiss_and_track(ring,orb)
   WRITE(*,*) "After QTune: Qx = ", ring%a%tune/twopi, "  Qy = ", ring%b%tune/twopi
ENDIF

n_z = nint(abs((e1 - e0)/de))+1

if( grid_type == 'xy' ) then
  n_a = nint(abs((x1 - x0)/dx))+1
  n_b = nint(abs((y1 - y0)/dy))+1
  Npts = n_z*n_b*n_a
  write(*,'(a)') "xy raster"
  write(*,'(/,1x,3(a6,e12.4,4x))') ' x0 = ',x0, ' y0 = ',y0, ' e0 = ',e0
  write(*,'(1x,3(a6,e12.4,4x))') ' dx = ',dx, ' dy = ',dy, ' de = ',de
  write(*,'(1x,3(a6,i4,6x,4x))') 'n_x = ',n_a,'n_y = ',n_b,'n_z = ',n_z
  write(*,'(1x,(a9,i0),/)') 'n_turn = ',n_turn
elseif( grid_type == 'rt' ) then
  write(*,'(a)') "polar raster"
  n_a = Nthmin
  n_b = Nrad
  asum = 0
  do i=1, Nrad
    asum = asum + i - 1
  enddo
  Npts = n_z*( n_b*(2*Nthmin-1) + 2*asum)
  Ar = Rmax / Nrad * sqrt(1.0d0 + Nrad*Nrad/Br/Br)
  write(*,'(2(a9,i11))') 'Nrad = ', Nrad, 'Nthmin = ', Nthmin
  write(*,'(2(a9,f11.5))') 'Ar = ', Ar, 'Br = ', Br
  write(*,'(3(a9,f11.5))') ' e0 = ',e0, ' de = ',de, 'n_z = ',n_z
  write(*,'(1x,(a9,i0),/)') 'n_turn = ',n_turn
endif
allocate(grid_pts(Npts,6))

ix_dot = index(file_name,'.')
prefix_name = file_name(1:ix_dot-1)
call string_trim(prefix_name, prefix_name,ix_pn)
call file_suffixer (file_name, in_file, '.out', .true.)

ALLOCATE(fft1(1:n_turn))
ALLOCATE(fft2(1:n_turn))
ALLOCATE(fft3(1:n_turn))
ALLOCATE(fft4(1:n_turn))
ALLOCATE(fft5(1:n_turn))
ALLOCATE(fft6(1:n_turn))

if(fft_turns == 0)fft_turns = n_turn / 2

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

! Multiplying coordinated by Tinv will convert them to normal coordinates.

CALL transfer_matrix_calc(ring, one_turn_mat)
if(one_turn_mat(6,5) == 0) then
  call mat_symp_decouple (one_turn_mat, stat, U, V, Ubar, Vbar, G, twiss1, twiss2, tgamma, .true.)
  CALL mat_inverse(V,Vinv)
else
  CALL mat_eigen (one_turn_mat, eigen_val, eigen_vec, error)
  do i=1,6,2
    T(i,:)   =  real(eigen_vec(i,:))
    T(i+1,:) = -aimag(eigen_vec(i,:))
  enddo
  CALL mat_inverse(T,Tinv)
endif

!Write grid to file
open(800, file='grid.dat')
n = 0
do i_z = 0, n_z-1
  start_coord%vec(1:6) = 0.0
  start_coord%vec(6) = e0 + de * i_z
  do i_b = 0, n_b-1
    if( grid_type == 'xy' ) then
      start_coord%vec(3) = y0 + dy * i_b
    elseif( grid_type == 'rt' ) then
      rad = Ar * (i_b+1.0d0) / sqrt(1.0d0 + ((i_b+1.0d0)/Br)**2)
    endif
    if( grid_type == 'xy' ) then
      do i_a = 0,n_a-1
        start_coord%vec(1) = x0 + dx * i_a
        n = n + 1
        grid_pts(n,:) = start_coord%vec(:)
        write(800,'(3es14.5)') start_coord%vec(1), start_coord%vec(3), start_coord%vec(5)
      enddo
    elseif( grid_type == 'rt' ) then
      do i_a = 1,2*(Nthmin+i_b)-1
        th = pi/2.0/(Nthmin+i_b-1.0d0)*(i_a-1.0d0)
        start_coord%vec(1) = rad*cos(th)
        start_coord%vec(3) = rad*sin(th)
        n = n + 1
        grid_pts(n,:) = start_coord%vec(:)
        write(800,'(3es14.5)') start_coord%vec(1), start_coord%vec(3), start_coord%vec(5)
      enddo
    endif
  enddo
enddo
close(800)

do i = 1, Npts
  start_coord%vec = grid_pts(i,:)
  if (start_coord%vec(6) < 0) then
     WRITE(out_file,'(2A,I6.5,A)') trim(out_file_prefix),".e",int(start_coord%vec(6)*1000),".fm"
  else
     WRITE(out_file,'(2A,I6.6,A)') trim(out_file_prefix),".e",int(start_coord%vec(6)*1000),".fm"
  endif
  open(unit=13, file=out_file)
  write(*,'(a,6es10.3)') "Coordinates: ", start_coord%vec(:)
  orb(0)%vec = co(0)%vec + start_coord%vec

  do j=1,n_turn
    call track_all(ring, orb, 0, track_state)
    orb(0)%vec = orb(ring%n_ele_track)%vec
    ! for whether particle z is past RF bucket
    if (track_state /= moving_forward$) then
      n = track_state
      print '(a, i8, 2a)', "Particle lost in turn ", j, ',  At element: ', trim(ring%ele(n)%name)
      print '(a, 6f12.6)', 'Orbit at entrance to element particle lost at: ', orb(n-1)%vec
      exit
    elseif ((z_cut) .and. (abs(orb(0)%vec(5)) > 0.25*(c_light/rf_frequency))) then 
      print '(a, i8)', "Particle outside of RF bucket in turn ", j
      print '(a, f12.6, a, f12.6)', 'z_position: ', orb(0)%vec(5), ',  RF wavelength: ', c_light/rf_frequency
      exit
    else
      psv = orb(0)%vec - co(0)%vec
      if(one_turn_mat(6,5) == 0) then
        psv_norm(1:4) = MATMUL(Vinv, psv(1:4))
        psv_norm(5) = psv(5)
        psv_norm(6) = psv(6)
      else
        psv_norm = MATMUL(psv, Tinv)
      endif

      fft1(j) = psv_norm(1)
      fft2(j) = psv_norm(2)
      fft3(j) = psv_norm(3)
      fft4(j) = psv_norm(4)
      fft5(j) = psv_norm(5)
      fft6(j) = psv_norm(6)
    endif
  enddo

  !if(track_state == moving_forward$ .and. (abs(orb(0)%vec(5)) < 0.25*(c_light/rf_frequency)))then
  if(track_state == moving_forward$)then
     if ((abs(orb(0)%vec(5)) < 0.25*(c_light/rf_frequency)) .or. (.not. z_cut)) then
        !interpolated FFT with Hanning window
        do j=1, fft_turns
           hanning = 2*(sin(pi*j/fft_turns)**2)
           n1(1, j)= cmplx(fft1(j),  fft2(j)) * hanning
           n2(1, j)= cmplx(fft3(j),  fft4(j)) * hanning
           n3(1, j)= cmplx(fft5(j), -fft6(j)) * hanning
           
           ind = n_turn - fft_turns
           m1(1, j)= cmplx(fft1(j+ind),  fft2(j+ind)) * hanning
           m2(1, j)= cmplx(fft3(j+ind),  fft4(j+ind)) * hanning
           m3(1, j)= cmplx(fft5(j+ind), -fft6(j+ind)) * hanning
        enddo
        
        call fourrow(n1(:,1:fft_turns), 1)    ! NR FFT
        call fourrow(m1(:,1:fft_turns), 1)
        call fourrow(n2(:,1:fft_turns), 1)
        call fourrow(m2(:,1:fft_turns), 1)
        call fourrow(n3(:,1:fft_turns), 1)
        call fourrow(m3(:,1:fft_turns), 1)
        
        forall(i=1:fft_turns)
           rftxa(i) = sqrt(n1(1,i)*conjg(n1(1,i)))
           rftxb(i) = sqrt(m1(1,i)*conjg(m1(1,i)))
           rftya(i) = sqrt(n2(1,i)*conjg(n2(1,i)))
           rftyb(i) = sqrt(m2(1,i)*conjg(m2(1,i)))
           rftea(i) = sqrt(n3(1,i)*conjg(n3(1,i)))
           rfteb(i) = sqrt(m3(1,i)*conjg(m3(1,i)))
        end forall
        
        tune(1, 1) = calc_tune(rftxa) !tune of x for 1st half etc.
        tune(1, 2) = calc_tune(rftxb) !tune of x for 2nd half etc.
        tune(2, 1) = calc_tune(rftya) !find tune for y first half etc.
        tune(2, 2) = calc_tune(rftyb) !find tune for y second half etc.
        tune(3, 1) = calc_tune(rftea) !find tune for z first half etc.
        tune(3, 2) = calc_tune(rfteb) !find tune for z second half etc.
        
        if (ring%a%tune/twopi < 0.5) then
           tune(1, 1) = 1. - tune(1, 1)
           tune(1, 2) = 1. - tune(1, 2)
        endif
        
        if (ring%b%tune/twopi < 0.5) then
           tune(2, 1) = 1. - tune(2, 1)
           tune(2, 2) = 1. - tune(2, 2) 
        endif
        
        write(13, '(12e14.5)') start_coord%vec(1), start_coord%vec(3), start_coord%vec(6), &
             tune(:, 1), tune(:, 2), tune(:,2)-tune(:,1)
     else
        write(13, '(3e14.5,9a14)') start_coord%vec(1), start_coord%vec(3), start_coord%vec(6), &
             "*","*","*","*","*","*","*","*","*"
     endif
  endif
  
  track_state = moving_forward$
enddo

contains
  function calc_tune(rft) result(tune)
    real(rp) tune
    real(rp) rft(:)
    integer k, iQ

    real(rp) d, b, c, A

    iQ = 1
    do k= 2, fft_turns
      if (rft(iQ)<rft(k)) iQ= k
    enddo
    d= rft(iQ)
    b= rft(iQ+1)
    c= cos(twopi/fft_turns)
    A= (-(d+b*c)*(d-b)+b*sqrt(c*c*(d+b)*(d+b)-2*d*b*(2*c*c-c-1)))/(d*d+b*b+2*d*b*c)
    tune = real(iQ, rp)/fft_turns+(1.0d0/twopi)*asin(A*sin(twopi/fft_turns))
  end function

end program













