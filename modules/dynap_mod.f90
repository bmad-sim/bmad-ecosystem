module dynap_mod

use bmad, only: rp
use pisa_mod, only: pop_struct, pool_struct

implicit none

type, extends(pop_struct) :: smart_pop_struct
  real(rp), allocatable :: ApC(:)
  real(rp), allocatable :: Q1(:,:)
end type

type mag_struct
  character(1) type  ! 'c' for chromatic, 'h' for harmonic
  character(18) name
  character(5) property  ! element property to set
  real(rp) lb    !lower bound for constraint calculation
  real(rp) ub    !upper bound for constraint calculation
  real(rp) lir   !lower bound for initial population range
  real(rp) uir   !upper bound for initial population range
  real(rp) mutate_delta
end type

contains

subroutine set_magnet_strengths(mags,ring,strengths)
  use bmad

  implicit none

  type(mag_struct) mags(:)
  type(lat_struct) ring
  real(rp) strengths(:)

  integer n_mags, n_loc
  type (ele_pointer_struct), allocatable :: eles(:)
  logical err
  character*18 var_str
  character*30 set_str
  integer i,j

  n_mags = size(mags)

  do i=1, n_mags
    if(mags(i)%name == '') exit
    call lat_ele_locator(mags(i)%name, ring, eles, n_loc, err)
    write(var_str,'(f18.8)') strengths(i)
    set_str = trim(adjustl(mags(i)%property))//'='//trim(adjustl(var_str))

    do j=1, n_loc
      call set_ele_attribute (eles(j)%ele, set_str, ring, err)
      if(err) then
        write(*,*) "Set ele attribute error.  Terminating.", set_str
        !error stop
        call err_exit
      endif
    enddo
  enddo

  if(allocated(eles)) deallocate(eles)
end subroutine

function is_linear_ele(ele)
  use bmad
  implicit none
  logical is_linear_ele
  type(ele_struct) ele

  if( (ele%key == sbend$) .or. &
      (ele%key == sextupole$) .or. &
      (ele%key == multipole$) ) then
    is_linear_ele = .false.
  else
    is_linear_ele = .true.
  endif
end function

function frac(x)
  use bmad, only: rp
  real(rp) frac, x
  frac = x - int(x,kind(x))
end function

subroutine read_initial_population(pool, n_pop, n_chrom, n_harmo, filename, ring, set_chrom_x, set_chrom_y, c_mags)
  use bmad

  implicit none

  type(pool_struct) pool(:)
  integer n_pop, n_chrom, n_harmo
  character(*) filename
  type(lat_struct) ring
  real(rp) set_chrom_x, set_chrom_y
  type(mag_struct) c_mags(:)

  character(1000) line
  integer n_omega
  integer i, j, k
  integer cr_dim
  integer dummy_int, iostat
  integer n_mags
  logical err_flag
  real(rp) feasible
  real(rp) c_mag_str(n_chrom)
  real(rp) h_mag_str(n_harmo)
  real(rp) throw_away(4)
  real(rp) mag_str(n_chrom+n_harmo)
  real(rp), allocatable :: omega(:)
  real(rp) ApC(n_chrom)
  real(rp) Q1(n_chrom,n_chrom-2)

  n_mags = n_chrom + n_harmo
  cr_dim = 2
  n_omega = n_chrom - cr_dim

  allocate(omega(n_omega))

  !- Make chromaticity matrices
  call lattice_bookkeeper(ring)
  call build_chrom_mat(ring, set_chrom_x, set_chrom_y, c_mags, ApC, Q1, err_flag)
  if(err_flag) then
    write(*,'(a)') "Could not build chromaticity matrices when reading initial population.  Aborting."
    stop
  endif

  open(100,file=filename)
  i=1
  do while (i .le. n_pop)
    read(100,'(a)',iostat=iostat) line
    if(iostat .ne. 0) then
      write(*,*) "Seed population smaller than n_pop.  Terminating program."
      error stop
    endif
    line = adjustl(line)
    if( (line(1:1)=="#") .or. (trim(line)=="") ) then
      cycle
    endif
    read(line,*) dummy_int, mag_str, throw_away
    c_mag_str = mag_str(1:n_chrom)
    h_mag_str = mag_str(1+n_chrom:n_harmo+n_chrom)

    call K2_to_omega(c_mag_str,Q1,ApC,omega)

    pool(i)%name = i
    pool(i)%x(1:n_omega) = omega
    pool(i)%x(1+n_omega:n_harmo+n_omega) = h_mag_str
    i = i + 1
  enddo

  close(100)
  if(allocated(omega)) deallocate(omega)
end subroutine

subroutine build_chrom_mat(ring, set_chrom_x, set_chrom_y, c_mags, ApC, Q1, err_flag)
  use bmad
  use f95_lapack
  use calc_ring_mod
  use make_pseudoinverse_mod

  implicit none

  type(lat_struct) ring
  real(rp) set_chrom_x, set_chrom_y
  type(mag_struct) c_mags(:)
  real(rp) ApC(:)
  real(rp) Q1(:,:)
  logical err_flag

  type(coord_struct), allocatable :: co(:)
  type (ele_pointer_struct), allocatable :: eles(:)
  real(rp), allocatable :: A(:,:)
  real(rp), allocatable :: Ap(:,:)
  real(rp), allocatable :: Ident(:,:)
  real(rp), allocatable :: B(:,:), tau(:)

  real(rp), allocatable :: sbend_k2_state(:)
  real(rp), allocatable :: sextupole_state(:)
  real(rp), allocatable :: multipole_state(:)

  integer, parameter :: la_work_size=1000
  real(rp) la_work(la_work_size)

  integer info
  character*16 var_str
  character*25 set_str
  real(rp), parameter :: delta = 100.0d0 !Step size in evaluating numerical derivative 

  integer i, j, k, n_loc
  integer n_chrom, cr_dim, n_omega
  real(rp) nat_chrom_x, nat_chrom_y
  real(rp) chrom_x, chrom_y, chrom_vec(2)

  err_flag = .true.

  n_chrom = size(c_mags)
  cr_dim = 2
  n_omega = n_chrom-2

  allocate(A(cr_dim,n_chrom))
  allocate(Ap(n_chrom,cr_dim))
  allocate(Ident(n_chrom,n_chrom))
  allocate(B(n_chrom,n_chrom))
  allocate(tau(n_chrom))

  call set_on_off(sextupole$, ring, off_and_save$, saved_values=sextupole_state)
  call set_on_off(multipole$, ring, off_and_save$, saved_values=multipole_state)
  call set_on_off(sbend$,     ring, off_and_save$, saved_values=sbend_k2_state, ix_attrib=k2$)

  call calc_ring(ring,4,co,err_flag)
  if(err_flag) then
    write(*,*) "Could not calculate ring without sextupoles & multipoles in build_chrom_mat."
    call early_exit()
    return
  endif

  call chrom_calc(ring, 1.0d-5, nat_chrom_x, nat_chrom_y, err_flag)
  if(err_flag) then
    write(*,*) "Could not calculate natural chromaticity in build_chrom_mat."
    call early_exit()
    return
  endif

  call set_on_off(sextupole$, ring, restore_state$, saved_values=sextupole_state)
  call set_on_off(multipole$, ring, restore_state$, saved_values=multipole_state)
  call set_on_off(sbend$,     ring, restore_state$, saved_values=sbend_k2_state, ix_attrib=k2$)

  deallocate(sextupole_state)
  deallocate(multipole_state)
  deallocate(sbend_k2_state)

  do k=1,n_chrom  !calculate response for family k
    do i=1,n_chrom  !loop over all families, turn family k on, an all others off
      if ( i .eq. k ) then
        write(var_str,'(F14.6)') delta
      else
        write(var_str,'(F14.6)') 0.0d0
      endif
      set_str = trim(adjustl(c_mags(i)%property))//'='//trim(adjustl(var_str))

      call lat_ele_locator(c_mags(i)%name, ring, eles, n_loc, err_flag)
      do j=1, n_loc
        call set_ele_attribute (eles(j)%ele, set_str, ring, err_flag)
        if(err_flag) then
          write(*,*) "Set ele attribute error.  Terminating."
          error stop
        endif
      enddo
      deallocate(eles)
    enddo
    call lattice_bookkeeper(ring)
    call calc_ring(ring,4,co,err_flag)
    if(err_flag) then
      write(*,*) "Could not calculate ring in response loop in build_chrom_mat."
      call early_exit()
      return
    endif

    call chrom_calc(ring, 1.0d-5, chrom_x, chrom_y, err_flag)
    if(err_flag) then
      write(*,*) "Could not calculate chromaticity in response loop in build_chrom_mat."
      call early_exit()
      return
    endif
    A(1,k) = (chrom_x-nat_chrom_x)/delta
    A(2,k) = (chrom_y-nat_chrom_y)/delta
  enddo

  call make_pseudoinverse(A,Ap)

  Ident = 0.0d0
  do i=1,n_chrom
    Ident(i,i) = 1.0d0
  enddo

  B = Ident - matmul(Ap,A)
  call dgeqrf(n_chrom,n_chrom,B,n_chrom,tau,la_work,la_work_size,info)
  if(info .ne. 0) then
    write(*,*) "Call to dgeqrf failed in build_chrom_mat."
    call early_exit()
    return
  endif
  call dorgqr(n_chrom,n_chrom,n_chrom,B,n_chrom,tau,la_work,la_work_size,info)
  if(info .ne. 0) then
    write(*,*) "Call to dorgqr failed in build_chrom_mat."
    call early_exit()
    return
  endif
  Q1=B(:,1:n_omega)
  do i=1,n_chrom-(cr_dim+1)
    Q1(i,i+1:n_omega) = 0.0d0
  enddo

  chrom_vec(1) = nat_chrom_x - set_chrom_x
  chrom_vec(2) = nat_chrom_y - set_chrom_y
  ApC = matmul(Ap,-1*chrom_vec)

  err_flag = .false.

  deallocate(A)
  deallocate(Ap)
  deallocate(Ident)
  deallocate(B)
  deallocate(tau)

  contains
    subroutine early_exit()
      ApC = 0.0d0
      Q1 = 0.0d0
      deallocate(A)
      deallocate(Ap)
      deallocate(Ident)
      deallocate(B)
      deallocate(tau)
    end subroutine
end subroutine

subroutine K2_to_omega(K2,Q1,ApC,omega)
  implicit none

  real(rp) omega(:), K2(:)
  real(rp) ApC(:), Q1(:,:)

  omega = matmul(transpose(Q1),(K2-ApC))
end subroutine

subroutine omega_to_K2(omega,ApC,Q1,K2)
  implicit none

  real(rp) omega(:), K2(:)
  real(rp) ApC(:), Q1(:,:)

  K2 = ApC + matmul(Q1,omega)
end subroutine

subroutine count_feasible_in_pop(pop,n_feasible)
  use bmad
  implicit none

  integer n_feasible
  type(smart_pop_struct) pop(:)
  integer i
  
  n_feasible = 0
  do i=1,size(pop)
    if(pop(i)%name .gt. 0) then
      if(all(pop(i)%c(:) .ge. 0)) then
        n_feasible = n_feasible + 1
      endif
    endif
  enddo
end subroutine

subroutine write_population(pop, generate_feasible_seeds_only, n_chrom, gen_num, filename, prec)
  use bmad
  implicit none

  type(smart_pop_struct) pop(:)
  integer gen_num
  integer generate_feasible_seeds_only
  integer n_chrom
  character(*) filename
  integer, optional :: prec

  character(3) prec_str
  character(3) prec_str2
  character(20) format_str
  integer i, cr_dim, n_omega
  real(rp) feasible
  real(rp) K2(n_chrom)

  cr_dim = 2
  n_omega = n_chrom - cr_dim

  if(present(prec)) then
    write(prec_str,'(i3)') prec
    write(prec_str2,'(i3)') prec+8
    format_str = '(i8,50es'//trim(adjustl(prec_str2))//'.'//trim(adjustl(prec_str))//')'
  else
    format_str = '(i8,50es19.11)'
  endif

  open(22,file=filename,access='append')
  write(22,'(a,i6)') "# Generation ", gen_num
  do i=1,size(pop)
    if(pop(i)%name .gt. 0) then
      if( all(pop(i)%c(:) .ge. 0.0) .or. (generate_feasible_seeds_only .le. 0) ) then
        call omega_to_K2(pop(i)%x(1:n_omega),pop(i)%ApC,pop(i)%Q1,K2)
        if(any(pop(i)%c(:) .lt. 0)) then
          feasible = 0.0
        else
          feasible = 1.0
        endif
        write(22,format_str) pop(i)%name, K2(:), pop(i)%x(1+n_omega:), pop(i)%o(:), feasible
      endif
    endif
  enddo
  write(22,*)
  write(22,*)
  close(22)
end subroutine write_population

end module
