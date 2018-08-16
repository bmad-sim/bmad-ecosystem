module dynap_mod

use bmad, only: rp
use pisa_mod, only: pop_struct, pool_struct
use moga_struct_mod

implicit none

contains

subroutine get_magnet_strengths(mags,ring,strengths)
  use bmad

  implicit none

  type(mag_struct) mags(:)
  type(lat_struct) ring
  real(rp) strengths(:)

  integer n_mags, n_loc
  type (ele_pointer_struct), allocatable :: eles(:)
  logical err
  integer i

  n_mags = size(mags)

  do i=1, n_mags
    if(mags(i)%name == '') exit
    call lat_ele_locator(mags(i)%name, ring, eles, n_loc, err)
    if(err .or. (n_loc .lt. 1)) then
      write(*,*) "Get ele attribute error.  Terminating."
      call err_exit
    endif
    strengths(i) = eles(1)%ele%value(k1$)
  enddo

  if(allocated(eles)) deallocate(eles)
end subroutine

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
    if(n_loc == 0) then
      write(*,*) "Magnet ", mags(i)%name, " not found!"
      call err_exit
    endif
    write(var_str,'(f18.8)') strengths(i)
    set_str = trim(adjustl(mags(i)%property))//'='//trim(adjustl(var_str))

    do j=1, n_loc
      call set_ele_attribute (eles(j)%ele, set_str, err)
      if(err) then
        write(*,*) "Set ele attribute error.  Terminating.", set_str
        !error stop
        call err_exit
      endif
    enddo
  enddo

  if(allocated(eles)) deallocate(eles)
  call lattice_bookkeeper(ring)
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

subroutine read_initial_population(pool, n_pop, n_linear, n_chrom, n_harmo, filename, ring, crm, err_flag)
  use bmad
  use crm_mod

  implicit none

  type(pool_struct) pool(:)
  integer n_pop, n_linear, n_chrom, n_harmo
  character(*) filename
  type(lat_struct) ring
  type(crm_struct) crm
  logical err_flag

  type(lat_struct) ring_working
  character(1000) line
  integer n_omega
  integer i, j, k
  integer dummy_int, iostat
  real(rp) feasible
  real(rp) l_mag_str(n_linear)
  real(rp) c_mag_str(n_chrom)
  real(rp) h_mag_str(n_harmo)
  real(rp) throw_away(4)
  real(rp) mag_str(n_linear+n_chrom+n_harmo)
  real(rp), allocatable :: omega(:)

  n_omega = n_chrom - 2

  allocate(omega(n_omega))

  open(100,file=filename)
  i=1
  do while (i .le. n_pop)
    ring_working = ring
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
    l_mag_str = mag_str(1:n_linear)
    c_mag_str = mag_str(1+n_linear:n_linear+n_chrom)
    h_mag_str = mag_str(1+n_linear+n_chrom:n_linear+n_harmo+n_chrom)

    if(n_linear .gt. 0) then
      call set_magnet_strengths(crm%l_mags,ring_working,l_mag_str)
      call crm_build(ring_working, crm, err_flag)
      if(err_flag) exit
    endif

    call K2_to_omega(c_mag_str,omega,crm)

    pool(i)%name = i
    pool(i)%x(1:n_linear) = l_mag_str
    pool(i)%x(1+n_linear:n_linear+n_omega) = omega
    pool(i)%x(1+n_omega+n_linear:n_linear+n_harmo+n_omega) = h_mag_str
    i = i + 1
  enddo
  close(100)

  call deallocate_lat_pointers(ring_working)
  deallocate(omega)
end subroutine

subroutine K2_to_omega(K2,omega,crm)
  use crm_mod
  implicit none

  real(rp) omega(:), K2(:)
  type(crm_struct) crm

  omega = matmul(crm%Q1t,(K2-crm%ApC))
end subroutine

subroutine omega_to_K2(omega,crm,K2)
  use crm_mod
  implicit none

  real(rp) omega(:), K2(:)
  type(crm_struct) crm

  K2 = crm%ApC + matmul(crm%Q1,omega)
end subroutine

subroutine count_feasible_in_pop(pop,n_feasible)
  use bmad
  implicit none

  integer n_feasible
  type(pop_struct) pop(:)
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

subroutine write_population(pop, generate_feasible_seeds_only, n_linear, n_chrom, gen_num, filename, prec)
  use bmad
  implicit none

  type(pop_struct) pop(:)
  integer gen_num
  integer generate_feasible_seeds_only
  integer n_linear, n_chrom
  character(*) filename
  integer, optional :: prec

  character(3) prec_str
  character(3) prec_str2
  character(20) format_str
  integer i, n_omega
  real(rp) feasible
  real(rp) K2(n_chrom)

  n_omega = n_chrom - 2

  if(present(prec)) then
    write(prec_str,'(i3)') prec
    write(prec_str2,'(i3)') prec+8
    format_str = '(i8,200es'//trim(adjustl(prec_str2))//'.'//trim(adjustl(prec_str))//')'
  else
    format_str = '(i8,200es19.11)'
  endif

  open(22,file=filename,access='append')
  write(22,'(a,i6)') "# Generation ", gen_num
  do i=1,size(pop)
    if(pop(i)%name .gt. 0) then
      if( all(pop(i)%c(:) .ge. 0.0) .or. (generate_feasible_seeds_only .le. 0) ) then
        if(any(pop(i)%c(:) .lt. 0)) then
          feasible = 0.0
        else
          feasible = 1.0
        endif
        write(22,format_str) pop(i)%name, pop(i)%x_phys(:), pop(i)%o(:), feasible
      endif
    endif
  enddo
  write(22,*)
  write(22,*)
  close(22)
end subroutine write_population

end module
