module crm_mod

use precision_def
use moga_struct_mod

implicit none

real(rp), parameter :: numerical_delta = 100.0d0 !Step size in evaluating numerical derivative 

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
    use f95_lapack
    use make_pseudoinverse_mod

    implicit none

    type(lat_struct) ring
    type(lat_struct) ring_working
    type(crm_struct) crm
    integer status
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

    integer i, j, k, n_loc
    integer n_chrom, n_omega
    real(rp) nat_chrom_x, nat_chrom_y
    real(rp) chrom_x, chrom_y, chrom_vec(2)

    ring_working = ring

    n_chrom = size(crm%c_mags)
    n_omega = n_chrom-2

    allocate(A(2,n_chrom))
    allocate(Ap(n_chrom,2))
    allocate(Ident(n_chrom,n_chrom))
    allocate(B(n_chrom,n_chrom))
    allocate(tau(n_chrom))

    call set_on_off(sextupole$, ring_working, off_and_save$, saved_values=sextupole_state)
    call set_on_off(multipole$, ring_working, off_and_save$, saved_values=multipole_state)
    call set_on_off(sbend$,     ring_working, off_and_save$, saved_values=sbend_k2_state, ix_attrib=k2$)

    call twiss_and_track(ring_working,co,status)
    if(status /= ok$) then
      write(*,*) "Could not calculate ring without sextupoles & multipoles in build_chrom_mat."
      call early_exit()
      return
    endif

    call chrom_calc(ring_working, 1.0d-5, nat_chrom_x, nat_chrom_y, err_flag)
    if(err_flag) then
      write(*,*) "Could not calculate natural chromaticity in build_chrom_mat."
      call early_exit()
      return
    endif

    call set_on_off(sextupole$, ring_working, restore_state$, saved_values=sextupole_state)
    call set_on_off(multipole$, ring_working, restore_state$, saved_values=multipole_state)
    call set_on_off(sbend$,     ring_working, restore_state$, saved_values=sbend_k2_state, ix_attrib=k2$)

    deallocate(sextupole_state)
    deallocate(multipole_state)
    deallocate(sbend_k2_state)

    do k=1,n_chrom  !calculate response for family k
      do i=1,n_chrom  !loop over all families, turn family k on, an all others off
        if ( i .eq. k ) then
          write(var_str,'(F14.6)') numerical_delta
        else
          write(var_str,'(F14.6)') 0.0d0
        endif
        set_str = trim(adjustl(crm%c_mags(i)%property))//'='//trim(adjustl(var_str))

        call lat_ele_locator(crm%c_mags(i)%name, ring_working, eles, n_loc, err_flag)
        do j=1, n_loc
          call set_ele_attribute (eles(j)%ele, set_str, ring_working, err_flag)
          if(err_flag) then
            write(*,*) "Set ele attribute error.  Terminating."
            error stop
          endif
        enddo
        deallocate(eles)
      enddo
      call lattice_bookkeeper(ring_working)
      call twiss_and_track(ring_working,co,status)
      if(status /= ok$) then
        write(*,*) "Could not calculate ring in response loop in build_chrom_mat."
        call early_exit()
        return
      endif

      call chrom_calc(ring_working, 1.0d-5, chrom_x, chrom_y, err_flag)
      if(err_flag) then
        write(*,*) "Could not calculate chromaticity in response loop in build_chrom_mat."
        call early_exit()
        return
      endif
      A(1,k) = (chrom_x-nat_chrom_x)/numerical_delta
      A(2,k) = (chrom_y-nat_chrom_y)/numerical_delta
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
    crm%Q1=B(:,1:n_omega)
    do i=1,n_chrom-(2+1)
      crm%Q1(i,i+1:n_omega) = 0.0d0
    enddo

    crm%Q1t = transpose(crm%Q1)

    chrom_vec(1) = nat_chrom_x - crm%set_chrom_x
    chrom_vec(2) = nat_chrom_y - crm%set_chrom_y
    crm%ApC = matmul(Ap,-1*chrom_vec)

    err_flag = .false.

    deallocate(A)
    deallocate(Ap)
    deallocate(Ident)
    deallocate(B)
    deallocate(tau)
    deallocate(co)

    call deallocate_lat_pointers(ring_working)

    contains
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
