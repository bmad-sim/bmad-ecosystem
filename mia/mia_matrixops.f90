module mia_matrixops

  !Contains matrix operations (SVD, FFT, etc) for MIA
  use mia_types

contains

  subroutine svd_fft(data)
    !
    !Calls SVD and FFT funtions
    !
    type(data_set) :: data 
    call svd(data)
    call fft(data)
  end subroutine svd_fft

  subroutine svd(data)

    !		      
    !This routine executes an SVD of the position history matrix using
    !LAPACK subroutine dgesdd.
    !

    implicit none
    type(data_set) data   !Data set
    integer :: i, j, info, lwork
    real(rp), allocatable :: work(:), A(:,:), temp(:,:)
    integer, allocatable :: iwork(:)
    real(rp) :: q(1,1)

    allocate(data%tau_mat(NUM_TURNS, 2*NUM_BPMS))    
    allocate(data%lambda(2*NUM_BPMS))
    allocate(data%pi_mat(2*NUM_BPMS, 2*NUM_BPMS)) 
    allocate(iwork(16*NUM_BPMS))
    allocate(A(NUM_TURNS, 2*NUM_BPMS))    
    allocate(temp(2*NUM_BPMS, 2*NUM_BPMS))
    !Allocates the pi, tau and lambda matrices based on turns and active
    !processors

    info = 0
    do i=1, 2*NUM_BPMS
       data%lambda(i) = 0.0_rp
       do j=1, 2*NUM_BPMS
          data%pi_mat(i,j) = 0.0_rp
          data%tau_mat(i,j) = 0.0
       enddo
    enddo

    A = data%poshis(1:NUM_TURNS,1:2*NUM_BPMS)
!I don't know what this is:
!    call out(data%poshis, "poshis")
    lwork = -1
    allocate(work(7883))
    !The first call to dgesdddetermines the ideal length of temporary array 
    !'work'; the SVD is not done until the second call.
    call dgesdd('S', NUM_TURNS, 2*NUM_BPMS, A, NUM_TURNS, data%lambda, &
         data%tau_mat, NUM_TURNS, data%pi_mat, 2*NUM_BPMS, work, lwork, &
         iwork, info)
    !Reallocate work to the ideal size
    lwork = work(1)
    deallocate(work)
    allocate(work(lwork))
    !Now do the SVD...
    call dgesdd('S', NUM_TURNS, 2*NUM_BPMS, A, NUM_TURNS, data%lambda, &
         data%tau_mat, NUM_TURNS, data%pi_mat, 2*NUM_BPMS, work, lwork, &
         iwork, info)

    if (.not.(info==0)) then
       Print *, "Error in column", info
    endif

    call transpose(data%pi_mat,2*NUM_BPMS, 2*NUM_BPMS)

    deallocate(work)
    deallocate(iwork)
  end subroutine svd

  subroutine transpose(matrix, row, col)

    real(rp), allocatable :: temp(:,:), matrix(:,:)
    integer :: row, col, i, j

    allocate(temp(col, row))
    do i=1, row
       do j=1, col
          temp(j,i) = matrix(i,j)
       enddo
    enddo
    deallocate (matrix)
    allocate(matrix(col,row))
    matrix = temp
  end subroutine transpose

  subroutine gullotine(arr, length)
    !
    !Truncates an array to a certain length.
    !
    real(rp), allocatable :: arr(:), temp(:)
    integer :: i, j, length

    allocate(temp(length))
    do i=1, length
       temp(i) = arr(i)
    enddo

    deallocate(arr)
    allocate(arr(length))
    arr = temp

  end subroutine gullotine

  subroutine powerof2(n)
    !
    !Determines the highest power of two which is less than or equal to n
    !Used to truncate arrays with NUM_TURNS not equal to a power of two
    !for FFT.
    !
    integer :: i, power, length, n
    logical :: goOn
    length =  NUM_TURNS
    goOn = .true.
    power = 0
    do while (goOn)
       length = length / 2
       power = power + 1
       if (length < 2) then
          goOn = .false.
       endif
    enddo
    if (length > 0) then
       n = 2**power
    else 
       n = NUM_TURNS
    endif
  end subroutine powerof2

  subroutine fft(data)
    !
    !calls fft routine
    !
    implicit none
    type(data_set) data         !Data set
    integer :: n, &             !Number of turns
         i, &                   !Counter
         fr_peak                !Frequency peak 
    real(rp), allocatable :: a(:), & !Col of tau_mat to be analyzed
         p(:)                   !Column of phi_spec   
    integer :: n_eigen          !Number of eigenvectors to do fft on
                                !It saves time to do the calculation only for
                                !the top 5-6 eigenmodes.

    n_eigen = 6
    call powerof2(n)            !Make n the next smaller power of 2
                                !if n is not a power of 2 already.
    allocate(data%spectrum(n,2*data%bpmproc))
    allocate(data%phi_spec(n,2*data%bpmproc))
    allocate(p(n))
    allocate(data%fr_peak(2*data%bpmproc))


    do i = 1, n_eigen

       allocate(a(n))
       a(:) = data%tau_mat(1:n,i)

       call polar_fft_f77(a,p,n,fr_peak)
       data%phi_spec(:,i) = p(:)
       data%spectrum(:,i) = a(:)
       data%fr_peak(i) = fr_peak
       deallocate(a)
    end do

    deallocate(p)

  end subroutine fft

  subroutine polar_fft_f77(amp,p,n,fr_peak)

    !  Receives a and n, returns amplitude a, phase (radians) and fr_peak
    implicit none

    integer ::  n,&        !?
         fr_peak, &    !Peak frequency
         i, isgn    !

    integer,save :: m = 0
    real(rp), save :: mf, nf
    real(rp) :: amp(n), &   !Amplitude
         
         p(n)          !?
    integer,save :: np = 0     !?
    real(rp) :: ar(n),ai(n), &  !ar = amp and ai = 0??
         max_amp       !Maximum amplitude

    ! Initialize
    if(n /= np) then
       np = n
       nf = n;
       mf = log(nf)/log(2.0_rp) + 0.001
       m = mf
    endif

    do i=1,n
       ar(i) = amp(i)
       ai(i) = 0
    enddo

    call rect_fft_f77(ar,ai,m)

    max_amp = -1.0_rp
    fr_peak = -1
    do i=1,n
       amp(i) = 2.0_rp*sqrt(ar(i)*ar(i)+ai(i)*ai(i))/n
       if(i >= 0 .AND. amp(i) > max_amp .and. i < n/2) then  ! find peak AC
          max_amp = amp(i)
          fr_peak = i-1   ! because FORTRAN arrays start with 1
       endif
       !	  if(ar(i).ne.0.0) then
       p(i) = atan2(ai(i),ar(i))
       !	    if(ar(i).lt.0.0) p(i) = isgn(ai(i))*pi+p(i)
       !	  else
       !	    p(i) = isgn(ai(i))*pi/2.0
       !	  endif
    enddo

    return

  end subroutine polar_fft_f77


  subroutine rect_fft_f77(ar,ai,m)
    !
    ! Radix-2, in-place FFT, after Cooley, Lewis and Welch.
    ! Complex ar,ai are input and output.
    ! Very fast if m doesn't change between calls.
    !
    implicit none

    real(rp) :: ar(*),ai(*)
    integer m

    integer, save :: n
    integer, save :: mp = 0
    real(rp), save ::  wr(32768),wi(32768)
    integer, save :: nv2,nm1
    real(rp) ::  ur,ui,tr,ti
    integer L,Le,Le1
    integer i,j,k,ip

    ! Initialize
    if(mp /= m) then
       mp = m
       n = 2**m
       do L=1,m
          Le = 2**L
          Le1 = Le/2
          !****  w(Le1) = exp(-j*pi/Le)
          wr(Le1) = cos(pi/Le1)
          wi(Le1) = sin(-pi/Le1)
       enddo
       nv2 = n/2
       nm1 = n-1
    endif

    ! Bit Reverse
    j = 1
    do i=1,nm1
       if(i < j) then
          tr = ar(j)
          ti = ai(j)
          ar(j) = ar(i)
          ai(j) = ai(i)
          ar(i) = tr
          ai(i) = ti
       endif
       k = nv2
       do while(k < j)
          j = j-k
          k = k/2
       enddo
       j = j+k
    enddo

    ! 	butterfly stages
    do L=1,m
       Le = 2**L
       Le1 = Le/2
       ur = 1.0_rp
       ui = 0.0_rp
       do j=1,Le1
          do i=j,n,Le
             ip = i+Le1
             !****  t = a(ip)*u   
             tr = ar(ip)*ur-ai(ip)*ui
             ti = ar(ip)*ui+ai(ip)*ur
             !****  a(ip) = a(ip)-t   
             ar(ip) = ar(i)-tr
             ai(ip) = ai(i)-ti
             !****  a(i) = a(i)+t   
             ar(i) = ar(i)+tr
             ai(i) = ai(i)+ti
          enddo
          !****  u = u*w  
          tr = ur*wr(Le1)-ui*wi(Le1)
          ti = ur*wi(Le1)+ui*wr(Le1)
          ur = tr
          ui = ti
       enddo
    enddo

    return

  end subroutine rect_fft_f77

  subroutine match_tau_column(nset,data)
    !
    !finds the columns that match in tau matrix, and finds if they are the 
    !horizontal pair, or the vertical pair.
    !
    integer :: nset               !Number of files
    type (data_set) :: data(*)    !Data set
    integer :: i, q, c, cset, &   !Counters; cset = current set
         col_counta, col_countb,& !Counters for columns found in A and B modes
         col_a_p, &               !Horizontal positive
         col_a_n, &               !Horizontal negative
         col_b_n, col_b_p, &      !Vertical positive and negative
         place, &                 !Counter--from 1 to nset
         colu1(5), colu2(5), &    !Column numbers of first and second match
         col_a_(2), col_b_(2), &  !Contains columns with A or B mode matches
         colu(4) , &              !Contains columns with Eigen modes 
         af(2), bf(2), &          !Which set A and B modes were found in
         lambdas, &               !Number of lambas assumed not to be noise
         pair_cap                 !Number of mode pairs to find before
    !moving on to the next file.
    real(rp) :: ave(2), &         !2*average of filtered noise (threshold)
         noise_sum(2), &          !Sum of filtered noise (last 5 lambda values)
         sum_a(4), &              !Sum of odd rows of pi matrix
         sum_b(4)                 !Sum of even rows of pi matrix
    logical :: mode(2), &         !True if a file has not been found for that
         more                     ! mode yet. 1 is A, 2 is B

    mode(1) = .true.
    mode(2) = .true.
    noise_sum(1) = 0.0_rp
    noise_sum(2) = 0.0_rp
    more = .false.

    !Number of lambdas used only affects if data is rejected for
    !having too much noise.
    lambdas = floor(2 * 0.6 * data(1)%bpmproc)   !Average 40% of lambda values

    !Last 5 lambda values are assumed to be noise. Nothing 2x the
    !average of these values is used to determine the eigenmode pairs.
    !(ave(cset) is the threshold)
    do cset = 1, nset
       call findNoise(data(cset)%lambda, ave(cset))
       ave(cset) = 2.0*ave(cset)
!       print *, ave(cset)
    end do
!    do cset = 1,nset
!       do i = lambdas, 2*NUM_BPMS
!          noise_sum(cset) = noise_sum(cset) + data(cset)%lambda(i)
!       end do
!       ave(cset) = 2.0*(noise_sum(cset)/(2*data(1)%bpmproc - lambdas))
!    end do

    !Compares all BPMs to each other to find eigen mode matches.
    ! Stops at the first match.
    !*Change to find more matches?
    do cset = 1, nset
       pair_cap = 0
       do i = 1, 2*data(1)%bpmproc-1
          do q = i+1, 2*data(1)%bpmproc-1
            !Finds a match if the difference between frequency peaks is greater
            !than 2 and both lambdas are above threshold.
             if   ( iabs(data(cset)%fr_peak(i)-data(cset)%fr_peak(q)) < 2 &
                  .and. data(cset)%lambda(i) >= ave(cset) &
                  .and. data(cset)%lambda(q) >= ave(cset)) then
                write (*, '(1x,2a,i2,a,i2,a,i2, a/)', advance = "no") &
                     "Potential", " Eigen Mode match for file ", cset,  &
                     " found in columns:", i, "  and ", q, "."

                ! Pair cap is the number of pairs to find before moving on
                ! to the next file.
                ! pair_cap = pair_cap + 1
!                call logic_get('N', 'Y', 'Use this match?', more)
                !more = .not. more
!                if (.not. more) then
                   colu1(cset) = i
                   colu2(cset) = q
                   pair_cap = pair_cap + 1
                   goto 101
!                endif
                if (nset == 2 .and. pair_cap == 1) then
                   go to 101
                end if
             end if
          end do
       end do
       if (pair_cap == 0) then
         Print *, "No eigen mode match found. Try again with a different file."
      endif
101    continue  
    end do

    colu(1) = colu1(1)
    colu(3) = colu1(2)
    colu(2) = colu2(1)
    colu(4) = colu2(2)

    col_counta = 0
    col_countb = 0

    do c = 1,4 

       sum_a(c) = 0.0_rp
       sum_b(c) = 0.0_rp

       !Sum_a(c) is the sum of odd columns of the pi matrix (x values),
       !sum_b(c) is the sum of even columns (y values).
       if (nset==1) then
          call splitsum(data(1)%pi_mat(:, colu(c)), sum_a(c), sum_b(c))
       else
          place = ceiling(1.0*c/nset)
          call splitsum(data(place)%pi_mat(:, colu(c)), sum_a(c), sum_b(c))
       endif

      !Assign horizontal or vertical based on which direction has a greater sum
       if (sum_a(c) > sum_b(c)) then
          col_counta = col_counta + 1
          col_a_(col_counta) = colu(c)
          af(col_counta) = ceiling(1.0*c/nset)
       else if (sum_a(c) < sum_b(c)) then 
          col_countb = col_countb + 1
          col_b_(col_countb) = colu(c)
          bf(col_countb) = ceiling(1.0*c/nset)
       end if

       !Mode() accounts for col_countx remaining 2 when the second file
       !is analyzed (MIA would otherwise think the second file was
       !in both modes).
       if (col_counta==2 .and. mode(1)) then
          Print *, "Excitation for file ", ceiling(1.0*c/nset), &
               "is in the A (x) mode."
          mode(1) = .false.
          xfile = data(ceiling(1.0*c/nset))%filename
       else if (col_countb==2 .and. mode(2)) then
          Print *, "Excitation for file ", ceiling(1.0*c/nset), &
               "is in the B (y) mode."
          yfile = data(ceiling(1.0*c/nset))%filename
          mode(2) = .false.
       endif

    end do

    !Check for files that only have excitation in a single mode
    if (col_counta > 2 .or. col_countb > 2) then
!       Print *, "counta: ", col_counta
!       Print *, "countb: ", col_countb
       Print *, "Excitation for both files appears to be in the same plane."
       Print *, "Try again with different files."
       STOP
    endif

    if (nset == 1) then 
       PRINT *, "nset==1"
       af(2) = 1
       bf(1) = 1
       bf(2) = 1
    end if

    !    write (*,'(1x,a,i2,a,i2)')"Horizontal Match: ", &
    !    col_a_(1),",",col_a_(2), "Vertical Match: ", col_b_(1),",", col_b_(2)

    if (abs(data(af(1))%phi_spec(data(af(1))%fr_peak(col_a_(1)),col_a_(1)) &
         - data( af(2))%phi_spec(data(af(2))%fr_peak(col_a_(2)),col_a_(2))) &
         < pi) then
       if (data(af(1))%phi_spec(data(af(1))%fr_peak(col_a_(1)),col_a_(1)) < &
            data(af(2))%phi_spec(data(af(2))%fr_peak(col_a_(2)),col_a_(2))) &
            then
          data_struc%col_a_p = col_a_(1)
          data_struc%col_a_n = col_a_(2)
       else 
          data_struc%col_a_p = col_a_(2)
          data_struc%col_a_n = col_a_(1)
       end if

    else 

       if (data(af(1))%phi_spec(data(af(1))%fr_peak(col_a_(1)),col_a_(1)) < &
            data(af(2))%phi_spec(data(af(2))%fr_peak(col_a_(2)),col_a_(2))) then
          data_struc%col_a_p = col_a_(2)
          data_struc%col_a_n = col_a_(1)
       else 
          data_struc%col_a_p = col_a_(1)
          data_struc%col_a_n = col_a_(2)
       end if
    end if

    if (af(1)==af(2)) then
       data_struc%set_num_a = af(1)
    else 
       print*, "af(1) is not = af(2)!"
    end if

    if (abs(data(bf(1))%phi_spec(data(bf(1))%fr_peak(col_b_(1)),col_b_(1)) &
         -  data(bf(2))%phi_spec(data(bf(2))%fr_peak(col_b_(2)),col_b_(2)))&
         < pi) then
       if (data(bf(1))%phi_spec(data(bf(1))%fr_peak(col_b_(1)),col_b_(1)) &
            < data(bf(2))%phi_spec(data(bf(2))%fr_peak(col_b_(2)),col_b_(2)))&
            then
          data_struc%col_b_p = col_b_(1)
          data_struc%col_b_n = col_b_(2)
       else 
          data_struc%col_b_p = col_b_(2)
          data_struc%col_b_n = col_b_(1)
       end if

    else 

       if (data(bf(1))%phi_spec(data(bf(1))%fr_peak(col_b_(1)),col_b_(1)) < &
            data(bf(2))%phi_spec(data(bf(2))%fr_peak(col_b_(2)),col_b_(2))) then
          data_struc%col_b_p = col_b_(2)
          data_struc%col_b_n = col_b_(1)
       else
          data_struc%col_b_p = col_b_(1)
          data_struc%col_b_n = col_b_(2)
       end if

    end if

    if (bf(1)==bf(2)) then
       data_struc%set_num_b = bf(1)
    else 
       print*, "bf(1) is not = bf(2)!"
    end if

    !    write (*,'(1x,a,i2)') &
    !         "Horizontal Positive: ", data_struc%col_a_p, &
    !         "Horizontal Negative: ", data_struc%col_a_n, &
    !         "Vertical Positive: ", data_struc%col_b_p,&
    !         "Vertical Negative: ", data_struc%col_b_n

  end subroutine match_tau_column

  subroutine findNoise(lambda, ave)

    real(rp) :: lambda(:), current, sum, ave
    integer :: i, j, k
    logical :: found

    sum = 0
    j = 0
    found = .false.
    do i = (2*NUM_BPMS-1), 6, -1
       if (lambda(i) > lambda(i+1)*10 .and. .not. found) then
          k = i
          current = lambda(i)
          found = .true.
          j = 1
          sum = lambda(i)
       endif

       if (found .and. lambda(i) < 10*current) then
          j = j+1
          sum = sum + lambda(i)
       end if
    end do

    if (.not. found) then
       j = 1
       sum = lambda(12) !Arbitrary, but will be near average...
    endif
    ave = sum / j
!    Print *, j
!    Print *, sum
!    Print *, ave

  end subroutine findNoise

  subroutine splitsum(pi_mat, sum_a, sum_b)

    !
    !Sums odd and even columns of the pi matrix separately.
    !

    integer :: i                   !Counter
    real(rp) :: sum_a, sum_b,  &   !Horizontal (a) and Vertical(b) sums
         pi_mat(:)                 !Pi matrix

    !Sums x and y values from the pi matrix (horizontal are odd colums 
    !and vertical are even columns).

    sum_a = 0
    sum_b = 0

    do i = 1, NUM_BPMS
       sum_a = sum_a + abs(pi_mat(2*i-1))             !Odd values (x)
       sum_b = sum_b + abs(pi_mat(2*i))               !Even values (y)
    end do

  end subroutine splitsum


  subroutine tune(data)
    !
    !Finds the tune of the machine and phi(t)
    !
    type(data_set) data(*)
    integer :: n, i, j, x(2), y(2), k, f, bin(2,3), n_set, turns
    real(rp) :: detM, M(3,3), C(3), M_inv(3,3), smallDet, MC(3), temp, &
         smallM(2,2), x_max, temp_tune(2), temp_phi_t(2)
    REAL :: tu, switch 

    call powerof2(turns)
    do n_set = 1, 2
       bin(1,2) = data(n_set)%fr_peak(data_struc%col_a_p)
       bin(2,2) = data(n_set)%fr_peak(data_struc%col_b_p)
       do n = 1, 2            !Only test first two modes.
          bin(n,1) = bin(n,2)-1
          bin(n,3) = bin(n,2)+1
          ! Print *, "fr_peak: ", bin(n,2)
          if (bin(n,2)==0) cycle  

          do i = 1, 3
             do j = 1, 3
                M(i,j) = bin(n,j)**(3-i)
             enddo
             C(i) = data(n_set)%spectrum(bin(n,i),n)
          enddo

          detM = 0.0_rp
          call det3x3(M, detM)

          !Find the inverse of M:
          do i = 1, 3  
             if (i==1) then 
                y(1) = 2
             else
                y(1) = 1
             endif
             if (i==3) then
                y(2) = 2
             else 
                y(2) = 3
             endif
             do j = 1, 3
                if (j==2) then
                   x(1) = 3
                   x(2) = 1
                else if(j==1) then
                   x(1) = 2
                   x(2) = 3
                else
                   x(1) = 1
                   x(2) = 2
                endif
                if (i==2) then
                   temp = x(2)
                   x(2) = x(1)
                   x(1) = temp
                endif
                call small_M(M, x, y, smallM)
                call det2x2(smallM, smallDet)
                M_inv(i,j) = smallDet / detM
             enddo
          enddo

          !Solve A = M_inverse * C:
          do k = 1, 3
             MC(k) = 0.
             do f = 1, 3
                MC(k) = MC(k) + M_inv(f,k)*C(f)
             enddo

          enddo

          x_max = -MC(2) / (2*MC(1))
          temp_tune(n) = ((x_max / turns) * FREQ)
          if (temp_tune(n) < FREQ) then
             temp_tune(n) = FREQ - temp_tune(n)
          endif

          temp_phi_t(n) = 2*pi*temp_tune(n) / FREQ

          data_struc%phi_t(n_set) = data_struc%phi_t(n_set) + &
               temp_phi_t(n)/2
          data_struc%tune(n_set) = data_struc%tune(n_set) + temp_tune(n)/2

       enddo
    enddo

    if (data_struc%set_num_a /= 1) then
       switch = data_struc%tune(1)
       data_struc%tune(1) = data_struc%tune(2)
       data_struc%tune(2) = switch
       switch = data_struc%phi_t(1)
       data_struc%phi_t(1) = data_struc%phi_t(2)
       data_struc%phi_t(2) = switch
    endif

    data_struc%tune(1) = mod(data_struc%tune(1), FREQ)
    data_struc%tune(2) = mod(data_struc%tune(2), FREQ)

  end subroutine tune


  subroutine small_M(M, x, y, smallM)
    !
    !Makes a 2x2 matrix from a 3x3 matrix, given a set of two x 
    !and two y coordinates to take values from.
    !
    integer :: i, j, x(2), y(2)
    real(rp) :: M(3,3), smallM(2,2)
    do i = 1, 2
       do j = 1, 2
          smallM(i,j) = M(x(i), y(j))
       enddo
    enddo
  end subroutine small_M

  subroutine det2x2(mat, detM)
    !
    !Finds the determinant of a 2x2 matrix
    !
    real(rp) :: mat(2,2), detM

    detM = mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1)

  end subroutine det2x2

  subroutine det3x3(mat, detM)
    !
    !Finds determinant of a 3x3 matrix
    !
    real(rp) :: mat(3,3), detM

    detM = mat(1,1)*mat(2,2)*mat(3,3) - &
         mat(1,1)*mat(3,2)*mat(2,3) - &
         mat(2,1)*mat(1,2)*mat(3,3) + &
         mat(2,1)*mat(3,2)*mat(1,3) + &
         mat(3,1)*mat(1,2)*mat(2,3) - &
         mat(3,1)*mat(2,2)*mat(1,3)

  end subroutine det3x3

  subroutine regen(data)
    !
    !Should regenerate data file read in...doesn't quit work.
    !
    type(data_set) data
    real(rp), allocatable :: temp(:,:), tlambda(:,:)
    integer :: i, j

    allocate(temp(NUM_TURNS, 2*NUM_BPMS))
    allocate(tlambda(NUM_TURNS, 2*NUM_BPMS))

    do i=1, NUM_TURNS
       do j=1, 2*NUM_BPMS
          if(i==j) then
             tlambda(i,j) = sqrt(data%lambda(i))
          else
             tlambda(i,j) = 0.
          endif
          temp(i,j) = 0.
       enddo
    enddo

    temp = matmul(tlambda, data%pi_mat)
    temp = matmul(data%tau_mat, temp)

    call out(temp, "regen ")
    deallocate(temp)
    deallocate (tlambda)

  end subroutine regen

  subroutine wls()
    !
    !Matrix svd and regeneration experiment by Willard Louis Sigma
    !
    INTEGER i,j
    INTEGER, PARAMETER :: M=2,N=3

    REAL(rp) :: a(M,N) = RESHAPE((/1.0,4.0,2.0,5.0,3.0,6.0/),(/M,N/))
    REAL(rp) :: b(N,M) = RESHAPE((/1.0,2.0,3.0,4.0,5.0,6.0/),(/N,M/))
    REAL(rp) :: c(M,M)
    REAL(rp) :: e(M,M)    
    real(rp) :: lam(M,M)
    REAL(rp) :: d(M,M) = RESHAPE((/1.0,5.0,2.0,7.0/),(/M,M/))
    REAL(rp) :: tau(M,M), lambda(M), pi(M,M)

    write(*,*) 'Matrix [a]'
    do i=1,M
       write(*,1000) (a(i,j),j=1,N)
    enddo
    write(*,*)

    write(*,*) 'Matrix [b]'
    do i=1,N
       write(*,1000) (b(i,j),j=1,M)
    enddo
    write(*,*)

    write(*,*) 'Matrix [d]'
    do i=1,M
       write(*,1000) (d(i,j),j=1,M)
    enddo
    write(*,*) 

    c = matmul(a, b)
    write(*,*) 'Matrix [c] = [a] x [b]'
    do i = 1,M
       write(*,1000) (c(i,j),j=1,M)
    enddo

    e = matmul(c, d)
    write(*,*) 'Matrix [e] = [c] x [d]'
    do i = 1,M
       write(*,1000) (e(i,j),j=1,M)
    enddo

    tau = e
    call svdcmp(tau, lambda, pi)
    write(*,*) 'tau'
    do i = 1,M
       write(*,1000) (tau(i,j),j=1,M)
    enddo
    write(*,*) 'Lambda'
    do i = 1,M
       write(*,1000) (lambda(i))
    enddo
    write(*,*) 'pi'
    do i = 1,M
       write(*,1000) (pi(i,j),j=1,M)
    enddo

    lam = RESHAPE((/0.0,0.0,0.0,0.0/),(/M,M/))
    lam(1,1) = lambda(1)
    lam(2,2) = lambda(2)
    e = matmul(lam, pi)
    e = matmul(tau, e)
    write(*,*) 'Matrix [e] = [c] x [d]'
    do i = 1,M
       write(*,1000) (e(i,j),j=1,M)
    enddo

1000 FORMAT(1x,1P10E14.6)
  end subroutine wls

  subroutine out(temp, name)
    !
    !Output function for a position history matrix.
    !
    implicit none
    real(rp) :: temp(:,:)
    integer :: i,j, openstatus
    integer :: bpm
    character(30) :: filename
    character(6) ::  name

    filename = "./data/" // trim(name) // ".out"

!    call wls

    open (unit = 27, file = trim(filename), &
         action = "write", position = "rewind",&
         iostat = openstatus)
    if (openstatus > 0) print *, "*** Cannot open output file ***",&
         openstatus
    do i=1, NUM_BPMS
       write(27,*) "BPM# ", i
       do j=1, NUM_TURNS
          write (27,*) temp(j,i), temp(j, 2*i)
       enddo
    enddo
  end subroutine out

end module mia_matrixops
