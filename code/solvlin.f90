!+
! logical Function Solvlin(mat, v, z, dim, size)
!
! Function to solve for 'z' in the linear equation, 'Mat*z = v', where 'Mat'
! and 'v' are given.
!
! Input:
!     mat(size,*) -- Real: Matrix to be inverted
!     v(*) -- Real: known input vector
!     dim -- Integer: dimension of vectors and matrix
!     size -- Integer: array size of matrix to be inverted
!
! Output:
!     z(*) -- Real: vector solution to the linear equation
!-

  logical function solvlin(mat,v,z,dim,size)
  implicit none
  integer dim, size
  integer i, n, m, p, max1, np1
  integer nsub
  real zero
  parameter(ZERO=0.0)
  real mat(size,*), z(*), v(*)
  real y(100), a(100,100), mult, sum
  real amax, ahold, yhold
  integer order(100), imax, tempord
  logical shiftFlag

  solvlin = .false.   !linear equations not yet solved

! **** INITIALIZE VARIABLES ****
  do n = 1, dim
     z(n) = ZERO	!set solution vector to 0.0
     y(n) = v(n)  !copy the vector v
     order(n) = n	!initial order of rows
     do m = 1, dim
  	a(m,n) = mat(m,n)   !copy the matrix mat
   	   enddo !m
  enddo !n

! **** PERFORM GAUSSIAN ELIMINATION ****
  max1 = dim - 1
  do n = 1, max1
      np1 = n + 1
      amax = abs( a(n,n) )
      imax = n
      shiftFlag = .FALSE.
      do i = n+1, dim
  	if( abs( a(i,n) ) > amax ) then
  	    amax = a(i,n)
  	    imax = i
  	    shiftFlag = .TRUE.
  	endif
      enddo !i
      if(shiftFlag) then
  	tempord = order(imax)
  	order(imax) = order(n)
  	order(n) = tempord
  	do nsub = n, dim
  	    ahold = a(n,nsub)
  	    a(n,nsub) = a(imax, nsub)
  	    a(imax,nsub) = ahold
  	enddo !nsub
  	yhold = y(n)
  	y(n) = y(imax)
  	y(imax) = yhold
      endif
      do m = np1, dim
  	if(a(n,n)==0.0) return
  	mult = - a(m,n) / a(n,n)
  	do p = np1, dim
  	   a(m,p) = mult * a(n,p) + a(m,p)
  	enddo !p
  	y(m) = mult * y(n) + y(m)
      enddo !m
  enddo !n

! **** BACK SUBSTITUTE FOR THE SOLUTION ****
  if(a(dim,dim)==0.0) return
  z(dim) = y(dim) / a(dim,dim)
  do n = max1, 1, -1
      sum = y(n)
      do m = n+1, dim
  	sum = sum - a(n,m)*z(m)
      enddo !m
      z(n) = sum / a(n,n)
  enddo !n

  solvlin = .true.   !linear equations solved

  return
  end
