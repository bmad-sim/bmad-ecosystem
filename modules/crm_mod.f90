module crm_mod

use precision_def
use moga_struct_mod

implicit none

real(rp), parameter :: numerical_delta = 1.0d0 !Step size in evaluating numerical derivative 

type crm_struct
  logical stale
  type(mag_struct), pointer :: l_mags(:)
  type(mag_struct), pointer :: c_mags(:)
  real(rp) set_chrom_x
  real(rp) set_chrom_y
  real(rp), allocatable :: ApC(:)
  real(rp), allocatable :: Q1(:,:)
  real(rp), allocatable :: Q1t(:,:)
end type

contains

  subroutine crm_alloc(crm,n_chrom)
    type(crm_struct) crm
    integer n_chrom

    allocate(crm%ApC(n_chrom)) 
    allocate(crm%Q1(n_chrom,n_chrom-2)) 
    allocate(crm%Q1t(n_chrom-2,n_chrom)) 
  end subroutine

  !-

  subroutine crm_build(ring, crm, err_flag)
    use bmad
    use bmad_routine_interface
    use f95_lapack

    implicit none

    type(lat_struct) ring
    type(lat_struct) ring_working
    type(crm_struct) crm
    integer status
    logical err_flag

    type(coord_struct), allocatable :: co(:)
    type (ele_pointer_struct), allocatable :: eles(:)
    type (ele_pointer_struct), allocatable :: slaves(:)
    real(rp), allocatable :: A(:,:)
    real(rp), allocatable :: Ap(:,:)
    real(rp), allocatable :: Ident(:,:)
    real(rp), allocatable :: B(:,:), tau(:)

    integer, parameter :: la_work_size=1000
    real(rp) la_work(la_work_size)

    integer info
    character*16 var_str
    character*25 set_str

    integer i, j, k, n_loc
    integer n_chrom, n_omega
    integer n_slave
    real(rp) init_chrom_x, init_chrom_y
    real(rp) chrom_x, chrom_y, chrom_vec(2)

    ring_working = ring

    n_chrom = size(crm%c_mags)
    n_omega = n_chrom-2

    allocate(A(2,n_chrom))
    allocate(Ap(n_chrom,2))
    allocate(Ident(n_chrom,n_chrom))
    allocate(B(n_chrom,n_chrom))
    allocate(tau(n_chrom))

    do i=1,n_chrom
      if(crm%c_mags(i)%type == 'c') then
        call set_magnet_k2(crm%c_mags(i)%name, ring_working, 0.0d0)
      endif
    enddo

    call twiss_and_track(ring_working,co,status)
    if(status /= ok$) then
      write(*,*) "Could not calculate ring without sextupoles & multipoles in build_chrom_mat."
      call early_exit()
      return
    endif

    call chrom_calc(ring_working, 1.0d-6, init_chrom_x, init_chrom_y, err_flag)
    !write(*,*) "Chromaticity with indicated sextupoles off: ", init_chrom_x, init_chrom_y
    if(err_flag) then
      write(*,*) "Could not calculate initial chromaticity in build_chrom_mat."
      call early_exit()
      return
    endif

    !build chromaticity response matrix
    do k=1,n_chrom  !calculate response for family k
      do i=1,n_chrom  !loop over all families, turn family k on, an all others off
        if ( i .eq. k ) then
          call set_magnet_k2(crm%c_mags(i)%name, ring_working, numerical_delta)
        else
          call set_magnet_k2(crm%c_mags(i)%name, ring_working, 0.0d0)
        endif
      enddo

      call twiss_and_track(ring_working,co,status)
      if(status /= ok$) then
        write(*,*) "Could not calculate ring in response loop in build_chrom_mat."
        call early_exit()
        return
      endif

      call chrom_calc(ring_working, 1.0d-6, chrom_x, chrom_y, err_flag)
      !write(*,*) "second chrom calc: ", chrom_x, chrom_y
      if(err_flag) then
        write(*,*) "Could not calculate chromaticity in response loop in build_chrom_mat."
        call early_exit()
        return
      endif

      A(1,k) = (chrom_x-init_chrom_x)/numerical_delta
      A(2,k) = (chrom_y-init_chrom_y)/numerical_delta
    enddo

    call mat_pseudoinverse(A,Ap)

    chrom_vec(1) = init_chrom_x - crm%set_chrom_x
    chrom_vec(2) = init_chrom_y - crm%set_chrom_y
    crm%ApC = matmul(Ap,-1*chrom_vec)

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
    crm%Q1=B(:,1:n_omega)
    do i=1,n_chrom-(2+1)
      crm%Q1(i,i+1:n_omega) = 0.0d0
    enddo

    crm%Q1t = transpose(crm%Q1)

    err_flag = .false.

    deallocate(A)
    deallocate(Ap)
    deallocate(Ident)
    deallocate(B)
    deallocate(tau)
    deallocate(co)

    call deallocate_lat_pointers(ring_working)

    contains
      subroutine set_magnet_k2(ele_name, ring, str)
        implicit none

        character(*) ele_name
        type(lat_struct) ring
        real(rp) str
        type (ele_pointer_struct), allocatable :: eles(:), slaves(:)
        integer j, n_loc, n_slave
        logical err_flag
        character(2) :: prop = 'k2'
        character(14) var_str
        character(25) set_str

        call lat_ele_locator(ele_name,ring_working,eles,n_loc)
        !write(*,'(a,i3,a,a)') "FOO Found ", n_loc, " magnets named ", ele_name
        write(var_str,'(F14.6)') str
        set_str = trim(adjustl(prop))//'='//trim(adjustl(var_str))
        !if( eles(1)%ele%lord_status .ne. not_a_lord$ ) then
        if(.false. ) then
          call get_slave_list(eles(1)%ele, slaves, n_slave)
          do j=1,n_loc
            call set_ele_attribute (slaves(j)%ele, set_str, err_flag)
          enddo
          deallocate(slaves)
        else
          do j=1,n_loc
            call set_ele_attribute( eles(j)%ele, set_str, err_flag)
          enddo
        endif
        deallocate(eles)
        call lattice_bookkeeper(ring_working)
      end subroutine

      subroutine early_exit()
        crm%ApC = 0.0d0
        crm%Q1 = 0.0d0
        crm%Q1t = 0.0d0
        deallocate(A)
        deallocate(Ap)
        deallocate(Ident)
        deallocate(B)
        deallocate(tau)
      end subroutine
  end subroutine

end module
