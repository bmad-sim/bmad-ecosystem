!+
! module naff_mod
!
! This module implements the NAFF algorithm for calculating the spectra
! of periodic data.
!
! Decomposes data of the form D(:) = D1(:) + i D2(:).
!
! freqs contains the frequencies found in the data.
! amps contains the complex amplitudes of these frequencies.
! freqs and amps are sorted with strongest abs(amps) first.
!
! If opt_dump_spectra=<some integer> the FFT spectra will be dumped to a fort.<some integer> file.
! This will also cause other debug information to be printed to stdout.
!
! The steps of NAFF are:
! 1) Estimate peak omega_1 in frequency spectrum using interpolated FFT.
! 2) Refine estimate by using optimizer to maximize <data|e^{i omega_1}>, and also return amplitude of this component.
! 3) Remove e_1 = amp*e^{i omega_i} component from the data.
! 4) Repeat step 1 to estimate the new frequency component.
! 5) Repeat step 2 to refine the new frequency component.
! 6) Use Gram-Schmidt to orthogonalize the new frequency component w.r.t. other components found thus far.
! 7) Delete the orthogonalized basis function from the data.
! 8) Repeat at step 4 until new frequency components are no longer significant.
!-
module naff_mod

use physical_constants

implicit none

contains
  !+
  ! subroutine naff(data1,data2,freqs,amps,n_freqs,opt_dump_spectra)
  !
  ! This subroutine implements the NAFF algorithm for calculating the spectra
  ! of periodic data.
  !
  ! Decomposes data of the form D(:) = D1(:) + i D2(:).
  !
  ! freqs contains the frequencies found in the data.
  ! amps contains the complex amplitudes of these frequencies.
  ! freqs and amps are sorted with strongest abs(amps) first.
  !
  ! If opt_dump_spectra=<some integer> the FFT spectra will be dumped to a fort.<some integer> file.
  ! This will also cause other debug information to be printed to stdout.
  !
  ! The steps of NAFF are:
  ! 1) Estimate peak omega_1 in frequency spectrum using interpolated FFT.
  ! 2) Refine estimate by using optimizer to maximize <data|e^{i omega_1}>, and also return amplitude of this component.
  ! 3) Remove e_1 = amp*e^{i omega_i} component from the data.
  ! 4) Repeat step 1 to estimate the new frequency component.
  ! 5) Repeat step 2 to refine the new frequency component.
  ! 6) Use Gram-Schmidt to orthogonalize the new frequency component w.r.t. other components found thus far.
  ! 7) Delete the orthogonalized basis function from the data.
  ! 8) Repeat at step 4 until new frequency components are no longer significant.
  !
  ! Input:
  !   data1(:)         - real(rp), real part of data.
  !   data2(:)         - real(rp), imaginary part of data.
  !   opt_dump_spectra - integer, optional, spectra written to fort file, also debug info printed.
  ! Output:
  !   freqs(:)         - real(rp), frequency components found in units of 0 to 1.
  !   amps(:)          - complex(rp), amplitudes of frequency components.
  !   n_freqs          - integer, number of components found.
  !-
  subroutine naff(data1,data2,freqs,amps,n_freqs,opt_dump_spectra)
    implicit none

    real(rp) data1(:)
    real(rp) data2(:)
    real(rp) freqs(:)
    complex(rp) amps(:)
    integer, optional :: opt_dump_spectra

    real(rp) r_amp, i_amp
    complex(rp) ak
    real(rp) amp_threshold
    real(rp) f_threshold 
    real(rp) freq
    complex(rp) basis(size(freqs),size(freqs))
    real(rp) temp_array(size(freqs))
    integer sorted(size(freqs))
    integer ix_strongest_freq

    integer i, j, max_freqs, size_data
    integer n_freqs
    logical new_freq

    size_data = size(data1)
    max_freqs = size(freqs)

    amp_threshold = (size_data-10)**0.5
    f_threshold = 2.0d0/(size_data-1)
    
    n_freqs = 0
    do while( n_freqs .lt. max_freqs )
      freq = interpolated_fft(data1(:), data2(:), opt_dump_spectrum=opt_dump_spectra, opt_dump_index=0)
      freq = maximize_projection(freq, data1(:), data2(:))

      new_freq = .true.
      do j=1,n_freqs
        if (abs(freq-freqs(j)) .lt. f_threshold) then
          new_freq = .false.
          i=j
          exit
        endif
      enddo

      if (new_freq) then
        n_freqs = n_freqs + 1
        i = n_freqs
        freqs(i) = freq
        call generate_basis_coeff(basis, freqs, size_data, i)
        call calc_projection(data1, data2, freqs, basis, i, ak)
        amps(i) = ak/sqrt((size_data-1)*1.0d0)
        if ( i .ne. 1 ) then
          if ( abs(amps(i)) .gt. abs(amps(ix_strongest_freq)) ) then
            ix_strongest_freq = i
          endif
        else
          ix_strongest_freq = 1
        endif
      else
        call calc_projection(data1, data2, freqs, basis, i, ak)
      endif

      if( present(opt_dump_spectra) ) then
        write(*,'(a,i3,a,f15.6,a,f15.6)') "naff debug info ( pre-removal projection ", i, ", amplitude and angle): ", abs(ak), " ", atan2(aimag(ak),real(ak))
      endif

      if( abs(ak)/sqrt((size_data-1)*1.0d0) .lt. abs(amps(ix_strongest_freq))/amp_threshold ) then
        exit
      endif

      call remove_component(data1, data2, freqs, basis, i, ak)

      if( present(opt_dump_spectra) ) then
        do j=1,i
          call calc_projection(data1(:), data2(:), freqs(:), basis, j, ak)
          write(*,'(a,i3,a,f15.6,a,f15.6)') "naff debug info (post-removal projection ", j, ", amplitude and angle): ", abs(ak), " ", atan2(aimag(ak),real(ak))
        enddo
      endif
    enddo

    !this inefficient sort algorithm is sufficient, because n_freqs should be at most 10 or so.
    do i=1,n_freqs
      temp_array(i) = abs(amps(i))
    enddo
    do i=1,n_freqs
      sorted(i) = 1
      do j=2,n_freqs
        if( temp_array(j) .gt. temp_array(sorted(i)) ) then
          sorted(i) = j
        endif
      enddo
      temp_array(sorted(i)) = -1.0
    enddo
    freqs(1:n_freqs) = freqs(sorted(1:n_freqs))
    amps(1:n_freqs) = amps(sorted(1:n_freqs))
  end subroutine

  !+
  ! recursive function evaluate_u(basis, tunes, ix, t) result(x)
  !
  ! Evaluates orthogonalized basis function for frequency tunes(ix).
  !-
  recursive function evaluate_u(basis,tunes,ix,t) result(x)
    implicit none
    complex(rp) basis(:,:)
    real(rp) tunes(:)
    integer ix
    real(rp) t

    complex(rp) x
    integer i
    
    x = evaluate_v(tunes(ix),t)
    do i=ix-1,1,-1
      x = x - basis(ix,i) * evaluate_u(basis,tunes,i,t)
    enddo
  end function

  !+
  ! function evaluate_v
  !
  ! Evaluates Euler's Formula: exp(2pi i tune t)
  !-
  function evaluate_v(tune,t) result(a)
    implicit none
    real(rp) tune
    real(rp) t
    complex(rp) a
    a = exp( twopi * (0.0d0,1.0d0) * tune * t)
  end function

  !+
  ! function u_onto_v
  !
  ! Evaluates projection of orthogonalized basis function onto Euler's Formula
  !-
  function u_onto_v(basis,tunes,N,ia,ib) result(x)
    implicit none
    complex(rp) basis(:,:)
    real(rp) tunes(:)
    integer N
    integer ia, ib
    complex(rp) x

    integer i

    ! easyInt takes complex conjugate of first argument.
    x = easyInt(tunes,N,ia,ib) / (N-1)
    do i=ia-1,1,-1
      x = x + conjg(calc_coeff(basis,ia,i)) * easyInt(tunes,N,i,ib) / (N-1)
    enddo
  end function

  !+
  ! function u_onto_u
  !
  ! Evaluates projection of one orthogonalized basis function onto another.
  !-
  function u_onto_u(basis,tunes,N,ia,ib) result(x)
    implicit none
    complex(rp) basis(:,:)
    real(rp) tunes(:)
    integer N
    integer ia, ib

    complex(rp) x

    integer j
    
    x = u_onto_v(basis,tunes,N,ia,ib)
    do j=ib-1,1,-1 
      x = x + calc_coeff(basis,ib,j) * u_onto_v(basis,tunes,N,ia,j) 
    enddo
  end function

  !+
  ! subroutine generate_basis_coeff
  !
  ! Generates basis coefficients for a new frequency component, which orthogonalize the new frequency component
  ! with respect to the other frequency components found thus far.
  !-
  subroutine generate_basis_coeff(basis,tunes,N,ix)
    implicit none
    complex(rp) basis(:,:)
    real(rp) tunes(:)
    integer N,ix
    real(rp) numerator

    integer i, j

    do j=ix-1,1,-1  
      basis(ix,j) = u_onto_v(basis,tunes,N,j,ix)
    enddo
  end subroutine

  !+
  ! function calc_coeff
  !
  ! Accumulates the basis coefficients needed when projecting component ixa onto ixb.
  !-
  function calc_coeff(basis,ixa,ixb) result(coeff)
    implicit none

    complex(rp) basis(:,:)
    integer ixa, ixb
    complex(rp) coeff

    integer comask(size(basis(:,1)),size(basis(:,1)))
    integer ixx, num_terms
    complex(rp) coeff_term
    integer i,j,N, k

    N=size(basis(:,1))

    if( ixa-ixb .lt. 1 ) then
      write(*,*) "Error in calc_coeff indexing: ", ixa, ixb
      stop
    endif
    num_terms = 2**(ixa-ixb-1)
    coeff = 0.0d0
    do ixx=0, num_terms-1
      call get_comask(ixa,ixb,ixx,comask)
      coeff_term = 1.0d0
      do i=1,N
        do j=1,N
          if(comask(i,j) == 1) then
            coeff_term = -coeff_term * basis(i,j)
          endif
        enddo
      enddo
      coeff = coeff + coeff_term
    enddo
  end function

  !+
  ! subroutine get_comask
  !
  ! This private subroutine implements a combinotorical method for calculating the coefficients
  ! necessary to orthogonalize the basis functions.
  !-
  subroutine get_comask(ia,ib,ix,comask)
    implicit none
    integer ia,ib,ix
    integer comask(:,:)

    integer sub_mask_size
    integer submask(ia-ib,ia-ib)
    integer binvec(ia-ib+1)

    sub_mask_size = ia - ib 

    comask = 0
    if( sub_mask_size == 1 ) then
      comask(ia,ib) = 1
    else
      binvec(1) = 1
      binvec(sub_mask_size+1) = 1
      call int2binvec(ix, binvec(2:sub_mask_size))
      call binvec2comask(binvec,submask)
      comask(ib+1:ia,ib:ia-1) = submask
    endif

  end subroutine

  !+
  ! subroutine int2binvec
  !
  ! Convert an integer to its binary representation.  Used only by get_comask.
  !-
  subroutine int2binvec(dec,binvec)
    implicit none
    integer dec, x, i, binvec(:)

    x=dec

    binvec(:) = 0
    do i=size(binvec),1,-1
      if( x .ge. 2**(i-1) ) then
        binvec(size(binvec)-i+1) = 1
        x = x - 2**(i-1) 
      endif
    enddo
  end subroutine

  !+
  ! subroutine binvec2comask
  !
  ! Generates a mask for the basis matrix.  Used for computing coefficients for
  ! orthogonalizing the basis functions.
  !-
  subroutine binvec2comask(binvec,submask)
    implicit none
    integer binvec(:)
    integer submask(:,:)
    integer N
    integer i,j
    integer i_index, j_index
    N = size(submask(1,:))
    submask(:,:) = 0

    i = 1  ! i gives row
    j = 1  ! j gives column
    do while(.true.)
      j = j + 1
      if ( binvec(j) .eq. 1 ) then
        i_index = N-i+1
        j_index = N-j+2
        submask(i_index,j_index) = 1
        if(j .eq. size(binvec)) then
          exit
        endif
        i=j
      endif
    enddo
  end subroutine

  !+
  ! subroutine easyInt
  !
  ! Analytic solution to the integral: int_0^(N-1) conjg(exp(-2pi i nu_i t))*exp(-2pi i nu_j t) dt
  !-
  function easyInt(tunes,N,j,k) result(a)
    implicit none
    real(rp) tunes(:)
    integer N,j,k
    complex(rp) a
    if( j == k ) then
      a = (N-1.0d0)*cmplx(1.0d0,0.0d0)
    else
      a = (0.0d0,-1.0d0)/twopi/(-tunes(j)+tunes(k)) * ( (1.0d0,0.0d0) - exp( (0.0d0,-1.0d0)*twopi*(-tunes(j)+tunes(k))*(N-1.0d0)) )
    endif
  end function

  !+
  ! subroutine calc_projection
  !
  ! Calculates the projection of the data onto an orthogonalized basis function: <data1(:) + i data2(:) | u_i/sqrt(N-1)>
  !-
  subroutine calc_projection(data1, data2, tunes, basis, ix, ak)
    implicit none

    real(rp) data1(:), data2(:)
    real(rp) tunes(:)
    complex(rp) basis(:,:)
    integer ix
    complex(rp) ak

    real(rp) t
    integer i, j, N

    N = size(data1)

    !Calculate component
    ak = 0.0d0
    do i=1, N
      t = i-1.0d0
      ak = ak + conjg(cmplx(data1(i),data2(i))) * evaluate_u(basis,tunes,ix,-t) / sqrt(N-1.0d0)
    enddo
  end subroutine

  !+
  ! subroutine remove_component
  !
  ! Deletes a component from the data: data1 + i data2 = data1 + i data2 - <ak | u_i/sqrt(N-1)>
  !-
  subroutine remove_component(data1, data2, tunes, basis, ix, ak)
    implicit none

    real(rp) data1(:), data2(:)
    real(rp) tunes(:)
    complex(rp) basis(:,:)
    integer ix
    complex(rp) ak
    real(rp) r_amp, i_amp
    complex(rp) component(size(data1))
    
    real(rp) t
    integer i, N

    N = size(data1)

    do i=1, N
      t = i-1.0d0
      component(i) = conjg(ak) * evaluate_u(basis,tunes,ix,-t) / sqrt(N-1.0d0)
    enddo

    do i=1, N
      data1(i) = data1(i) - real(component(i))
      data2(i) = data2(i) - aimag(component(i))
    enddo
  end subroutine

  !+
  ! function maximize_projection
  !
  ! Optimizer that uses Numerical Recipes brent to find a local maximum, which is the frequency that maximizes the projection.
  !
  !-
  function maximize_projection(seed, data1, data2)
    use nr

    implicit none

    real(rp) data1(:), data2(:)
    real(rp) maximize_projection  !fractional frequency. range [0:1]
    real(rp) seed
    
    integer i, N
    real(rp) bracket
    real(rp) :: tol = 1d-8
    real(rp) fmin !throw away

    real(rp) r, window, hsamp

    complex(rp), allocatable :: data(:)

    !

    N = size(data1)
    bracket = 5.0d0/N

    allocate(data(N))
    data(:) = cmplx(data1(:),data2(:))

    !apply window
    hsamp = (N-1)/2.0d0
    do i=1, N
      window = 0.5d0*(1.0d0 + cos(twopi*(i-hsamp-1)/(N-1)))  !hanning
      !!r = 8.0d0
      !!window = exp( -0.5d0*(r*(1.0d0*i-hsamp-1)/(N-1))**2)  !gaussian
      data(i)= data(i) * window
    enddo

    fmin = brent(seed-bracket, seed, seed+bracket, special_projection, tol, maximize_projection)
    deallocate(data)

    contains
    !+
    ! function special_projection
    !
    ! Calculates <data1 + i data2 | exp(i theta)>
    ! Used only by maximize projection.  Uses data global to the module to accomodate stock NR routine.
    !-
    function special_projection(f)
      implicit none

      real(rp), intent(in) :: f !fractional frequency.  range [0:1]
      real(rp) special_projection

      complex(rp) ak
      integer i
      real(rp) theta

      ak = 0
      do i=1, size(data(:))
        theta = twopi*f*(i-1)
        ak = ak + data(i) * exp( (0.0d0,1.0d0) * theta)
      enddo

      special_projection = -abs(ak)
    end function
  end function

  !+
  !  function interpolated_fft
  !
  !  Windows the complex data and used Numerical Recipes four1 to find the peak in the spectrum.
  !  The result is interpolated to improve the accuracy.  Hanning and Gaussian windowing are
  !  available.
  !-
  function interpolated_fft(data1, data2, opt_dump_spectrum, ak, opt_dump_index)
    use nr

    real(rp) data1(:), data2(:)
    integer, optional :: opt_dump_spectrum
    integer, optional :: opt_dump_index
    complex, optional :: ak

    integer dump_spectrum
    integer dump_index

    complex(rp) cplx_data(size(data1))
    real(rp) fft_amp(size(data1))
    real(rp) interpolated_fft
    real(rp) window, r, hsamp
    real(rp) lk, lkm, lkp, A

    integer n_samples
    integer max_ix
    integer i
    integer isign

    dump_spectrum = 0
    if( present(opt_dump_spectrum) ) then
      dump_spectrum = opt_dump_spectrum
    endif

    dump_index = 0
    if( present(opt_dump_index) ) then
      dump_index = opt_dump_index
    endif

    n_samples = size(data1)
    hsamp = (n_samples-1)/2.0d0

    !apply window
    do i=1, n_samples
      !window = 0.5d0*(1.0d0 + cos(twopi*(i-hsamp-1)/(n_samples-1)))  !hanning
      r = 8.0
      window = exp( -0.5d0*(r*(1.0d0*i-hsamp-1.0d0)/(n_samples-1.0d0))**2)  !gaussian
      ! window = 1
      cplx_data(i)= cmplx( data1(i),  data2(i) ) * window
    enddo

    isign = 1
    call four1(cplx_data(:), isign)
    fft_amp(:)=sqrt(cplx_data(:)*conjg(cplx_data(:)))

    if( dump_spectrum > 10 ) then
      do i=1,n_samples
        write(dump_spectrum,*) dump_index, (i-1.0d0)/n_samples, fft_amp(i)
      enddo
      write(dump_spectrum,*)
      write(dump_spectrum,*)
    endif

    max_ix = 2
    do i= 3, n_samples-1
      if ( fft_amp(max_ix) < fft_amp(i) ) then
        max_ix = i
      endif
    end do

    !Gaussian Interpolation (use with gaussian window)
    lk = log(fft_amp(max_ix))
    lkm = log(fft_amp(max_ix-1))
    lkp = log(fft_amp(max_ix+1))
    A = (lkp-lkm) / 2.0d0 / (2.0d0*lk - lkp - lkm)
    interpolated_fft = ( 1.0d0*(max_ix-1) + A ) / n_samples

    ! Parabolic Interpolation (use with hanning window)
    ! lk = fft_amp(max_ix)
    ! lkm = fft_amp(max_ix-1)
    ! lkp = fft_amp(max_ix+1)
    ! A = (lkp-lkm) / 2.0 / (2.0*lk - lkp - lkm)
    ! interpolated_fft = 1.0d0*max_ix/n_samples + A/n_samples
  end function

end module naff_mod







