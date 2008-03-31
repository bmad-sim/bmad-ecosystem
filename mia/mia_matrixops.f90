module mia_matrixops

  !Contains matrix operations (SVD and FFT) for MIA
  use mia_types

!  integer :: NUM_BPMS, PI
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
    !this routine executes the SVD and uses the nr and precision_def modules
    !

    implicit none
    type(data_set) data   !Data set
    integer :: i, j

    allocate(data%tau_mat(1:NUM_TURNS, 1:2*NUM_BPMS))    
    allocate(data%lambda(1:2*NUM_BPMS))		     
    allocate(data%pi_mat(1:2*NUM_BPMS, 1:2*NUM_BPMS)) 
    !allocates the pi, tau and lambda matrices based on turns and active
    !processors

    do i=1, 2*NUM_BPMS
       data%lambda(i) = 0.0_rp
       do j=1, 2*NUM_BPMS
          data%pi_mat(i,j) = 0.0_rp
       enddo
    enddo

    data%tau_mat = data%poshis  
    !sets matrix to be used in SVD to be the position history matrix

    call svdcmp(data%tau_mat, data%lambda, data%pi_mat)

  end subroutine svd

  subroutine fft(data)

    !
    !calls fft routine
    !

    implicit none
    type(data_set) data         !Data set
    integer :: n, &             !Number of turns
         i, &                   !Counter
         fr_peak                !Frequency peak 
    real(rp), allocatable :: a(:), & !Takes col of tau_mat, becomes col of spectrum
         p(:)                   !Column of phi_spec   

    n = data%numturns

    !  allocate(data%tau_vec(n,2*data%bpmproc))
    allocate(data%spectrum(n,2*data%bpmproc))
    allocate(data%phi_spec(n,2*data%bpmproc))
    allocate(p(n))
    allocate(data%fr_peak(2*data%bpmproc))

    do i = 1, 2*data%bpmproc

       allocate(a(n))
       a(:) = data%tau_mat(:,i)

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
    real(rp) :: fr_peak(28), &    !Peak frequency
         ave(2), &                !2*average of filtered noise (threshold)
         noise_sum(2), &          !Sum of filtered noise (last 5 lambda values)
         sum_a(4), &              !Sum of odd rows of pi matrix
         sum_b(4)                 !Sum of even rows of pi matrix


    noise_sum(1) = 0.0_rp
    noise_sum(2) = 0.0_rp
    lambdas = floor(2 * 0.6 * data(1)%bpmproc)   !Average 40% of lambda values
    Print *, "Using ", lambdas, "lambda values"

   !Filters noise--last 5 lambda values are assumed to be noise. Nothing 2x the
   !average of these values is used. (ave(cset) is the threshold)

    !*Make this more general--don't assume 14 BPMs
    do cset = 1,nset
       do i = lambdas, 2*data(1)%bpmproc
          noise_sum(cset) = noise_sum(cset) + data(cset)%lambda(i)
       end do
       ave(cset) = 2.0*(noise_sum(cset)/(2*data(1)%bpmproc - lambdas))
    end do

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
                pair_cap = pair_cap + 1
                !Pair cap is the number of pairs to find before moving on
                !to the next file.
             !   PRINT *, "q:", q
                colu1(cset) = i        
                colu2(cset) = q

                if (nset == 2 .and. pair_cap == 1) then
                   go to 101
                end if
             end if
          end do
       end do
101    continue  
    end do

    colu(1) = colu1(1)
    colu(3) = colu1(2)
    colu(2) = colu2(1)
    colu(4) = colu2(2)

    col_counta = 1
    col_countb = 1

    do c = 1,4 

       sum_a(c) = 0.0_rp
       sum_b(c) = 0.0_rp

       !Sum_a(c) is the sum of odd columns of the pi matrix (x values),
       !sum_b(c) is the sum of even columns (y values).
       if (nset==1) then
          call splitsum(data(1)%pi_mat(:, colu(c)), sum_a(c), sum_b(c))
       elseif (nset==2) then
          place = ceiling(1.0*c/nset)
          call splitsum(data(place)%pi_mat(:, colu(c)), &
               sum_a(c), sum_b(c))
       endif

      !Assign horizontal or vertical based on which direction has a greater sum
       if (sum_a(c) > sum_b(c)) then
          col_a_(col_counta) = colu(c)
          af(col_counta) = ceiling(1.0*c/nset)
          col_counta = col_counta + 1
       else if (sum_a(c) < sum_b(c)) then 
          col_b_(col_countb) = colu(c)
          bf(col_countb) = ceiling(1.0*c/nset)
          col_countb = col_countb + 1
       end if

       if (col_counta==2) then
          Print *, "Excitation for file ", ceiling(c/2.0), &
               "is in the A (x) mode."
      !    data_struc%A_mode = ceiling(c/2.0)
       else if (col_countb==2) then
          Print *, "Excitation for file ", ceiling(c/2.0), &
               "is in the B (y) mode."
      !    data_struc%B_mode = ceiling(c/2.0)
       endif

    end do

    !Check for files that only have excitation in a single mode
    if (col_counta > 3 .or. col_countb > 3) then
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

    write (*,'(1x,a,i2,a,i2)')"Horizontal Match: ", &
         col_a_(1),",",col_a_(2), "Vertical Match: ", col_b_(1),",", col_b_(2)

    if ( abs(data( af(1))%phi_spec( data( af(1))%fr_peak(col_a_(1)), &
         col_a_(1)) - &
         data( af(2))%phi_spec( data(af(2))%fr_peak(col_a_(2)), &
         col_a_(2))) < pi) then
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
         -  data(bf(2))%phi_spec(data(bf(2))%fr_peak(col_b_(2)),col_b_(2)))< pi)then
       if (data(bf(1))%phi_spec(data(bf(1))%fr_peak(col_b_(1)),col_b_(1)) &
            < data(bf(2))%phi_spec(data(bf(2))%fr_peak(col_b_(2)),col_b_(2))) then
          data_struc%col_b_p = col_b_(1)
          data_struc%col_b_n = col_b_(2)
       else 
          data_struc%col_b_p = col_b_(2)
          data_struc%col_b_n = col_b_(1)
       end if

    else 

       if (data(bf(1))%phi_spec(data(bf(1))%fr_peak(col_b_(1)),col_b_(1))< &
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

    write (*,'(1x,a,i2)') &
         "Horizontal Positive: ", data_struc%col_a_p, &
         "Horizontal Negative: ", data_struc%col_a_n, &
         "Vertical Positive: ", data_struc%col_b_p,&
         "Vertical Negative: ", data_struc%col_b_n

  end subroutine match_tau_column

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
    integer :: n, i, j, x(2), y(2), k, f, bin(2,3), n_set
    real(rp) :: detM, M(3,3), C(3), M_inv(3,3), smallDet, MC(3), temp, &
         smallM(2,2), x_max, temp_tune(2), temp_phi_t(2)
    REAL :: tu, switch 



    do n_set = 1, 2
       bin(1,2) = data(n_set)%fr_peak(data_struc%col_a_p)
       bin(2,2) = data(n_set)%fr_peak(data_struc%col_b_p)
       do n = 1, 2            !Only test first two modes.
                               !Increments by two because values are paired
          
 
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
          temp_tune(n) = ((x_max / 1024) * FREQ)
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

end module mia_matrixops
