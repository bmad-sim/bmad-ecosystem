program freq_map

use bmad
use nr
use eigen_mod
use nrutil, only: dpc
use z_tune_mod

implicit none

type (lat_struct), target :: ring
type (coord_struct), allocatable :: co(:), orb(:), data(:)
type (coord_struct) start_coord
type (normal_modes_struct) mode
type (ele_struct), pointer :: ele

integer i, j, k, isign, r, s, g
integer ix
integer n_turn
integer ios, track_state
integer i_x, i_y, i_z, iq
integer n_x, n_y, n_z
integer ix_pn, ix_dot
integer version
integer :: fft_turns = 0
integer arg_num, iargc, ind

real(rp) phy_x_set, phy_y_set
real(rp) delta_e/0.0/, chromx, chromy, qp_x, qp_y
real(rp), allocatable :: dk1(:)
real(rp) eps_rel(4), eps_abs(4)
real(rp), allocatable :: tune(:,:,:)
real(rp) d, b, c, A
real(rp) x0, y0, e0, x1, y1, e1 
real(rp) :: dx=0., dy=0., de=0.
real(rp) :: rf_frequency
real(rp), ALLOCATABLE :: fft1(:), fft2(:), fft3(:)
real(rp), ALLOCATABLE :: fft4(:), fft5(:), fft6(:)
complex(dpc), ALLOCATABLE ::  n1(:,:), m1(:,:),  n2(:,:), m2(:,:), n3(:,:), m3(:,:), n4(:,:)
real(rp), ALLOCATABLE ::  rftxa(:), rftxb(:), rftya(:), rftyb(:), rftea(:), rfteb(:)
real(rp) hanning

real(rp) one_turn_mat(1:6,1:6)
real(rp) T(1:6,1:6), Tinv(1:6,1:6), P(1:6,1:6), Pinv(1:6,1:6)
real(rp) psv(1:6), psv_norm(1:6)

complex(rp) eigen_val(6), eigen_vec(6,6)

character*40 lattice, out_file_prefix
character*180 lat_file, out_file
character*80 line, last_line, file_name, prefix_name
character*4 file_type
character*1 answer
character*60 in_file

logical keep_trying/.true./, error
logical write_orbit/.false./
logical ok, aperture_limits, err
logical :: auto_bookkeeper = .true.
logical :: z_cut = .true.

namelist / parameters /lat_file, file_type, out_file_prefix, &
              x0, y0, e0, x1, y1, e1, dx, dy, de, n_turn, aperture_limits,  &
              fft_turns, auto_bookkeeper, z_cut



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

read(1, nml = parameters)

bmad_com%auto_bookkeeper = auto_bookkeeper

print *, ' lat_file = ', lat_file

call fullfilename(lat_file, lat_file)
call bmad_parser(lat_file, ring)

if (z_cut) then
  do ix = 1, ring%n_ele_track
    if (ring%ele(ix)%key == key_name_to_key_index('RFCAVITY',.true.)) then
      rf_frequency = value_of_attribute(ring%ele(ix),'RF_FREQUENCY', err)
      exit
    endif
  enddo
else
  rf_frequency = 1.e-6 ! something small; prevents longitudinal cuts until much larger
endif

call reallocate_coord(orb, ring%n_ele_max)
call reallocate_coord(co, ring%n_ele_max)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
bmad_com%aperture_limit_on = aperture_limits
CALL reallocate_coord(orb,ring%n_ele_max)
CALL set_on_off(rfcavity$, ring, on$)
CALL twiss_and_track(ring,orb)
CALL calc_z_tune(ring)


n_x = nint(abs((x1 - x0)/dx))+1
n_y = nint(abs((y1 - y0)/dy))+1
if(de /= 0.)n_z = nint(abs((e1 - e0)/de))+1
if(de == 0.)n_z=1

allocate(tune(0:n_x*n_y*n_z, 3, 2))

if (dx == 0.) dx = (x1-x0)/(n_x-1)
if (dy == 0.) dy = (y1-y0)/(n_y-1)
if((n_z .ne. 1) .and. (de == 0)) de = (e1-e0)/(n_z-1)

ix_dot = index(file_name,'.')
prefix_name = file_name(1:ix_dot-1)
call string_trim(prefix_name, prefix_name,ix_pn)
call file_suffixer (file_name, in_file, '.out', .true.)

print '(/,1x,3(a6,e12.4,4x))',' x0 = ',x0, ' y0 = ',y0, ' e0 = ',e0
print '(1x,3(a6,e12.4,4x))',' dx = ',dx, ' dy = ',dy, ' de = ',de
print '(1x,3(a6,i4,6x,4x))','n_x = ',n_x,'n_y = ',n_y,'n_z = ',n_z
print '(1x,(a9,i0),/)','n_turn = ',n_turn

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
CALL mat_eigen (one_turn_mat, eigen_val, eigen_vec, error)
do i=1,6,2
  T(i,:)   =  real(eigen_vec(i,:))
  T(i+1,:) = -aimag(eigen_vec(i,:))
enddo
CALL mat_inverse(T,Tinv)

g= 0
do i_z = 0,n_z-1
  start_coord%vec(1:6) = 0.0
  start_coord%vec(6) = e0 + de * i_z
  if (start_coord%vec(6) < 0) then
     WRITE(out_file,'(2A,I6.5,A)') trim(out_file_prefix),".e",int(start_coord%vec(6)*1000),".fm"
  else
     WRITE(out_file,'(2A,I6.6,A)') trim(out_file_prefix),".e",int(start_coord%vec(6)*1000),".fm"
  endif
  open(unit=13, file=out_file)
  do i_y = 0,n_y-1
    start_coord%vec(3) = y0 + dy * i_y
    do i_x = 0,n_x-1
      start_coord%vec(1) = x0 + dx * i_x
     
      write(*,*) " " 
      write(*,'(a,i0,a,i0)') "ix, iy = ", i_x, "   ", i_y
      write(*,'(a,6es10.3)') "Coordinates: ", start_coord%vec(:)
      
      co(0)%vec = orb(0)%vec + start_coord%vec

      do j=1,n_turn
        call track_all(ring, co, 0, track_state)
        co(0)%vec = co(ring%n_ele_track)%vec
        if((track_state /= moving_forward$) .or. (abs(co(0)%vec(5)) > 0.25*(c_light/rf_frequency))) then ! add check for whether particle z is increasing past RF bucket
          WRITE(*,*) "Particle lost in turn ", j
          exit
        else
          psv = co(0)%vec - orb(0)%vec
          psv_norm = MATMUL(psv, Tinv)

          fft1(j) = psv_norm(1)
          fft2(j) = psv_norm(2)
          fft3(j) = psv_norm(3)
          fft4(j) = psv_norm(4)
          fft5(j) = psv_norm(5)
          fft6(j) = psv_norm(6)
        endif
      end do

      if(track_state == moving_forward$ .and. (abs(co(0)%vec(5)) < 0.25*(c_light/rf_frequency)))then
        !interpolated FFT with Hanning window
        do j=1, fft_turns
          hanning = 2*(sin(pi*j/fft_turns)**2)
          n1(1, j)= cmplx(fft1(j),  fft2(j)) * hanning
          n2(1, j)= cmplx(fft3(j),  fft4(j)) * hanning
          n3(1, j)= cmplx(fft5(j), -fft6(j)) * hanning
        enddo

        do j=1, fft_turns
          hanning = 2*(sin(pi*j/fft_turns)**2)
          ind = n_turn - fft_turns
          m1(1, j)= cmplx(fft1(j+ind),  fft2(j+ind)) * hanning
          m2(1, j)= cmplx(fft3(j+ind),  fft4(j+ind)) * hanning
          m3(1, j)= cmplx(fft5(j+ind), -fft6(j+ind)) * hanning
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
        iQ = 1
        do k= 2, fft_turns
          if (rftxa(iQ)<rftxa(k)) iQ= k
        end do
        d= rftxa(iQ)
        b= rftxa(iQ+1)
        c= cos(twopi/fft_turns)
        A= (-(d+b*c)*(d-b)+b*sqrt(c*c*(d+b)*(d+b)-2*d*b*(2*c*c-c-1)))/(d*d+b*b+2*d*b*c)
        tune(g, 1, 1)= real(iQ, rp)/fft_turns+(1/twopi)*asin(A*sin(twopi/fft_turns))
        if (ring%a%tune/twopi < 0.5) tune(g, 1, 1) = 1. - tune(g, 1, 1)
 
        !tune of x for 2nd half etc.
        iQ = 1
        do k= 2, fft_turns
          if (rftxb(iQ)<rftxb(k)) iQ= k
        end do
        d= rftxb(iQ)
        b= rftxb(iQ+1)
        c= cos(twopi/fft_turns)
        A= (-(d+b*c)*(d-b)+b*sqrt(c*c*(d+b)*(d+b)-2*d*b*(2*c*c-c-1)))/(d*d+b*b+2*d*b*c)
        tune(g, 1, 2)= real(iQ, rp)/fft_turns+(1/twopi)*asin(A*sin(twopi/fft_turns))
        if (ring%a%tune/twopi < 0.5) tune(g, 1, 2) = 1. - tune(g, 1, 2)

        !find tune for y first half etc.
        iQ = 1
        do k= 2, fft_turns
          if (rftya(iQ)<rftya(k)) iQ= k
        end do
        d= rftya(iQ)
        b= rftya(iQ+1)
        c= cos(twopi/fft_turns)
        A= (-(d+b*c)*(d-b)+b*sqrt(c*c*(d+b)*(d+b)-2*d*b*(2*c*c-c-1)))/(d*d+b*b+2*d*b*c)
        tune(g, 2, 1)= real(iQ, rp)/fft_turns+(1/twopi)*asin(A*sin(twopi/fft_turns))
        if (ring%b%tune/twopi < 0.5) tune(g, 2, 1) = 1. - tune(g, 2, 1)
 
        !find tune for y second half etc.
        iQ = 1
        do k= 2, fft_turns
          if (rftyb(iQ)<rftyb(k)) iQ= k
        end do
        d= rftyb(iQ)
        b= rftyb(iQ+1)
        c= cos(twopi/fft_turns)
        A= (-(d+b*c)*(d-b)+b*sqrt(c*c*(d+b)*(d+b)-2*d*b*(2*c*c-c-1)))/(d*d+b*b+2*d*b*c)
        tune(g, 2, 2)= real(iQ, rp)/fft_turns+(1/twopi)*asin(A*sin(twopi/fft_turns))
        if (ring%b%tune/twopi < 0.5) tune(g, 2, 2) = 1. - tune(g, 2, 2) 

        !find tune for z first half etc.
        iQ = 1
        do k= 2, fft_turns
          if (rftea(iQ)<rftea(k)) iQ= k
        end do
        d= rftea(iQ)
        b= rftea(iQ+1)
        c= cos(twopi/fft_turns)
        A= (-(d+b*c)*(d-b)+b*sqrt(c*c*(d+b)*(d+b)-2*d*b*(2*c*c-c-1)))/(d*d+b*b+2*d*b*c)
        tune(g, 3, 1)= real(iQ, rp)/fft_turns+(1/twopi)*asin(A*sin(twopi/fft_turns))

        !find tune for z second half etc.
        iQ = 1
        do k= 2, fft_turns
          if (rfteb(iQ)<rfteb(k)) iQ= k
        end do
        d= rfteb(iQ)
        b= rfteb(iQ+1)
        c= cos(twopi/fft_turns)
        A= (-(d+b*c)*(d-b)+b*sqrt(c*c*(d+b)*(d+b)-2*d*b*(2*c*c-c-1)))/(d*d+b*b+2*d*b*c)
        tune(g, 3, 2)= real(iQ, rp)/fft_turns+(1/twopi)*asin(A*sin(twopi/fft_turns))
 
        write(13, '(12e14.5)') start_coord%vec(1), start_coord%vec(3), start_coord%vec(6), &
                              tune(g, :, 1), tune(g, :, 2), tune(g,:,2)-tune(g,:,1)
      endif ! moving_forward
      track_state = moving_forward$
    end do !i_x
  end do !i_y
  close(unit=13)
end do !i_z

end program













