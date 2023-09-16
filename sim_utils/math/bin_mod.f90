module bin_mod

use utilities_mod
use precision_def
use output_mod

implicit none

type bin_struct                       ! 1D bins
  real(rp), allocatable :: count(:)   ! Counts (or weight) in each bin
  real(rp) :: min                     ! Bounds for the bins
  real(rp) :: max                     !
  real(rp) :: delta                   ! Size of a bin
  integer  :: n                       ! Number of bins
end type

type general_bin_struct               ! Similar to bin_struct, but for multiple dimensions
  real(rp), allocatable :: count(:)   ! Counts (or weight) in each bin
  real(rp) :: min(3)                  ! Bounds for the bins
  real(rp) :: max(3)                  !
  real(rp) :: delta(3)                ! Size of a bin
  integer  :: dim = 3                 ! Number of dimensions
  integer  :: n(3) = 1                ! number of bins in each dimension
end type

contains

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function bin_data_density (bin_data, x, order) result (r)
! 
! Calculate the density of binned data at an arbitrary location x.
! Zero is returned if x is out of bounds of the binned data.
!
! Input: 
!   bin_data    -- bin_struct: binned data
!   x           -- real(rp): position to query
!   order       -- integer, optional: interpolation order: 0 or 1 only. Default: 1
!
! Output:
!   r           -- real(rp): density at x
!-

function bin_data_density (bin_data, x, order) result (r)

type(bin_struct) :: bin_data
real(rp) :: x, x1, r, r1, r2, rel_x
integer, optional :: order
integer :: ix1, ix2, ord
character(30), parameter :: r_name = 'density'

!

ord = integer_option(1, order)

! Does the same thing: ord = merge(order, 1, present(order))
! Get nearest index

ix1 = bin_index(x, bin_data%min, bin_data%delta)

r1 = count_at_index(bin_data, ix1)/bin_data%delta

! Zeroth order interpolation 
if (ord == 0) then
  r = r1
  return
endif

if (ord /= 1) then
  call out_io (s_error$, r_name, 'Order must be 0 or 1')
  call err_exit
endif

! Get next nearest bin
x1 = bin_x_center(ix1, bin_data%min, bin_data%delta)
rel_x = abs(x-x1)/bin_data%delta
ix2 = merge(ix1-1, ix1+1, x<x1)
r2 = count_at_index(bin_data, ix2)/bin_data%delta

! Linear interpolation
r = (1-rel_x)*r1 + rel_x*r2

end function

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function bin_data (data, weight, min, max, n_bins) result (binned_data)
!  
! Bin centers are at [ x_min + 1/2*delta, x_min + 3/2*delta, ..., x_max - 1/2*delta ] 
! with delta  = (max-min)/n_bins
!
! Input: 
!   data(:)    -- real(rp): 1D data to bin
!   weight(:)  -- real(rp), optional: 1D weights for each data. Default: 1
!   min        -- real(rp), optional: minimum considered. Default: minval(data)
!   max        -- real(rp), optional: maximum considered. Default: maxval(data)
!   n_bins     -- integer,  optional: number of bins. Default: 2*size(data)^(1/3) (Rice's rule)
!
! Output:
!   binned_data   -- bin_struct: binned data.
!-

function bin_data (data, weight, min, max, n_bins) result (binned_data)

implicit none
type(bin_struct) :: binned_data
real(rp) :: data(:)
real(rp), optional :: weight(:), min, max
integer, optional :: n_bins
integer :: i, ix
character(30), parameter :: r_name = 'bin'

!

binned_data%min = real_option(minval(data), min)
binned_data%max = real_option(maxval(data), max)
binned_data%n = integer_option(n_bins_automatic(size(data)), n_bins)

! Allocate and initialize 
if (allocated(binned_data%count)) deallocate(binned_data%count)
allocate(binned_data%count(binned_data%n))
binned_data%count = 0

! Step size
binned_data%delta = (binned_data%max - binned_data%min)/binned_data%n

! Populate
do i=lbound(data, 1), ubound(data, 1)
  ix = bin_index(data(i), binned_data%min, binned_data%delta) 
  if (ix < 1 .or. ix > binned_data%n) cycle ! Out of range
  if (present(weight))then
    binned_data%count(ix) = binned_data%count(ix) + weight(i)
  else
    binned_data%count(ix) = binned_data%count(ix) + 1.0_rp
  endif
enddo
  
end function bin_data

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function bin_data_density_2d (bin_data, x, y, order) result(r0)
! 
! Calculate the density of 2D binned data at an arbitrary location (x, y)
! Zero is returned if x or y is out of bounds of the binned data
!
! Input: 
!   bin_data    -- general_bin_struct: binned data. Must have %dim == 2
!   x           -- real(rp): x position to query
!   y           -- real(rp): y position to query
!   order       -- integer, optional: interpolation order: 0 or 1 only
!                                     Default: 1
! Output:
!   r0          -- real(rp): density at (x,y)
!-
function bin_data_density_2d (bin_data, x, y, order) result(r0)

type(general_bin_struct) :: bin_data
real(rp) :: x, y, x1, y1, r(4), r0, rel_x, rel_y
integer, optional :: order
integer :: i, ix(4), iy(4), ord, shift_x, shift_y
character(30), parameter :: r_name = 'density_2d'

!

ord = integer_option(1, order)

if (bin_data%dim /=2) then
  call out_io (s_error$, r_name, '%dim must be 2')
  call err_exit
endif

! Get nearest bin 
ix(1) = bin_index(x, bin_data%min(1), bin_data%delta(1))
iy(1) = bin_index(y, bin_data%min(2), bin_data%delta(2))

! Zeroth order interpolation 
if (ord == 0) then
  r0 = general_bin_count(bin_data, ix(1), iy(1))/(bin_data%delta(1)*bin_data%delta(2))
  return
endif

if (ord /= 1) then
  call out_io (s_error$, r_name, 'Order must be 0 or 1')
  call err_exit
endif

! Get nearest bin center
x1  = bin_x_center(ix(1), bin_data%min(1), bin_data%delta(1))
y1  = bin_x_center(iy(1), bin_data%min(2), bin_data%delta(2))

! Get indices of surrounding bins
rel_x = (x - x1)/bin_data%delta(1)
rel_y = (y - y1)/bin_data%delta(2)
shift_x = merge(1, -1, rel_x>0)
shift_y = merge(1, -1, rel_y>0)
rel_x  = abs(rel_x)
rel_y  = abs(rel_y)

ix(2) = ix(1) + shift_x
iy(2) = iy(1) 

ix(3) = ix(1) 
iy(3) = iy(1) + shift_y

ix(4) = ix(2)
iy(4) = iy(3)

! Get counts
do i=1, 4
  r(i) = general_bin_count(bin_data, ix(i), iy(i))
enddo

! Bi-Linear interpolation

r0 = (1-rel_x)*(1-rel_y)*r(1) + (rel_x)*(1-rel_y)*r(2) + (1-rel_x)*(rel_y)*r(3) + (rel_x)*(rel_y)*r(4)
! Scale for density
r0 = r0/(bin_data%delta(1)*bin_data%delta(2))

end function bin_data_density_2d

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function bin_2d(data1, data2, weight, min1, max1, min2, max2, n_bins1, n_bins2)
!          result (bin_data)
!
! Similiar to bin(...), but for two dimensions.
!
! Input: 
!   data1(:)   -- real(rp): data to bin in dimension 1
!   data2(:)   -- real(rp): data to bin in dimension 2
!   weight(:)  -- real(rp), optional: 1D weights for each data. 
!                                     Default: 1
!   min1       -- real(rp), optional: minimum considered for data1. Default: minval(data1)
!   min2       -- real(rp), optional: minimum considered for data2. Default: minval(data2)
!   max1       -- real(rp), optional: maximum considered for data1. Default: maxval(data1)
!   max2       -- real(rp), optional: maximum considered for data2. Default: maxva2(data1)
!   n_bins1    -- integer,  optional: number of bins for dimension 1.
!                                        Default: 2*size(data1)^(1/6) 
!   n_bins2    -- integer,  optional: number of bins for dimension 2.
!                                        Default: 2*size(data2)^(1/6) 
! Output:
!   bin_data   -- general_bin_struct: binned data, with %dim==2
!- 

function bin_2d(data1, data2, weight, min1, max1, min2, max2, n_bins1, n_bins2) result (bin_data)
implicit none
type(general_bin_struct) :: bin_data
real(rp) :: data1(:), data2(:)
real(rp), optional :: weight(:), min1, max1, min2, max2
integer, optional :: n_bins1, n_bins2
integer :: i, ix(2), index
character(30), parameter :: r_name = 'bin_2d'

!

bin_data%dim = 2

bin_data%min(1) = real_option(minval(data1), min1)
bin_data%min(2) = real_option(minval(data2), min2)

bin_data%max(1) = real_option(maxval(data1), max1)
bin_data%max(2) = real_option(maxval(data2), max2)

! Rice rule 
bin_data%n(1) = integer_option(ceiling((2.0_rp*size(data1))**(1.0/6.0)), n_bins1)
bin_data%n(2) = integer_option(ceiling((2.0_rp*size(data2))**(1.0/6.0)), n_bins2)

! Allocate and initialize 
!if (allocated(bin_data%count)) deallocate(bin_data%count)
if (allocated(bin_data%count)) deallocate(bin_data%count)
allocate(bin_data%count(bin_data%n(1)*bin_data%n(2) ))
bin_data%count = 0

! Step size
do i=1, bin_data%dim
  bin_data%delta(i) = (bin_data%max(i) - bin_data%min(i))/bin_data%n(i)
enddo

if ( (lbound(data1,1) /= lbound(data2,1)) .or. (ubound(data1,1) /= ubound(data2,1)) ) then
  call out_io (s_error$, r_name, 'Array labels must be the same for binning' )
  call err_exit
endif

! Populate
do i=lbound(data1, 1), ubound(data1, 1)
  ix(1) = bin_index(data1(i), bin_data%min(1), bin_data%delta(1)) 
  ix(2) = bin_index(data2(i), bin_data%min(2), bin_data%delta(2)) 
  
  ! check for out of range
  if (.not. general_bin_index_in_bounds(bin_data, ix(1), ix(2)) ) cycle

  index = general_bin_index(bin_data, ix(1), ix(2))
  
  if (present(weight))then
    bin_data%count(index) = bin_data%count(index) + weight(i)
  else
    bin_data%count(index) = bin_data%count(index) + 1.0_rp
  endif
enddo

end function

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function bin_index(x, bin1_x_min, bin_delta) result (ix_bin)
!
! Helper function to locate the appropriate histogram bin index.
!
! Input:
!   x           -- real(rp): Input value to bin.
!   bin1_x_min  -- real(rp): Minimum value of bin with index 1.
!   bin_delta   -- real(rp): Bin width.
!
! Output:
!   ix_bin      -- integer: Index of bin x is in.
!-

function bin_index (x, bin1_x_min, bin_delta) result(ix_bin)

real(rp) :: x, bin1_x_min, bin_delta
integer :: ix_bin

ix_bin = ceiling((x-bin1_x_min)  / bin_delta)

end function bin_index

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
! Function bin_x_center (ix_bin, bin1_x_min, bin_delta) result(x_center)
!
! Helper function to locate the center of a histogram bin.
!
! Input:
!   ix_bin      -- integer: Index of bin under question.
!   bin1_x_min  -- real(rp): Minimum value of bin with index 1.
!   bin_delta   -- real(rp): Bin width.
!
! Output:
!   ix_bin      -- int

function bin_x_center (ix_bin, bin1_x_min, bin_delta) result(x_center)

real(rp) :: x_center, bin1_x_min, bin_delta
integer :: ix_bin

!

x_center = (ix_bin-0.5_rp) * bin_delta + bin1_x_min

end function bin_x_center

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
! Wrapper to return zero if out of bounds

function count_at_index(bin_data, index ) result(c)
type(bin_struct) bin_data
real(rp) :: c
integer :: index
!
if (index < 1 .or. index > size(bin_data%count)) then
  c = 0
else 
  c = bin_data%count(index)
endif
end function

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
! Function to automatically select the number of bins

function n_bins_automatic(n_data) result(n)
implicit none
integer :: n_data, n
! Rice rule
n = ceiling((2.0_rp*n_data)**(1.0/3.0))
! Sturges' rule: n = ceiling(1.0_rp + log(1.0_rp*n_data)/0.693147_rp)
end function

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
! Function for looking up an index in the 1D count array

function general_bin_index(bin_data, ix1, ix2, ix3) result (index)
implicit none
type(general_bin_struct) :: bin_data
integer ix1, index
integer, optional :: ix2, ix3
index = ix1
if (present(ix2)) index = index + bin_data%n(1)*(ix2-1)
if (present(ix3)) index = index + bin_data%n(1)*bin_data%n(2)*(ix3-1)
end function

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
! Function for getting the count at a general index. Count will be 0 if out of bounds

function general_bin_count(bin_data, ix1, ix2, ix3) result (count)
implicit none
type(general_bin_struct) :: bin_data
integer ix1
integer, optional :: ix2, ix3
real(rp) :: count
count = 0
if (present(ix3)) then
  if (general_bin_index_in_bounds(bin_data, ix1, ix2, ix3) ) &
    count = bin_data%count(general_bin_index(bin_data, ix1, ix2, ix3))
else if (present(ix2) ) then
  if (general_bin_index_in_bounds(bin_data, ix1, ix2) ) &
    count = bin_data%count(general_bin_index(bin_data, ix1, ix2))
else
  if (general_bin_index_in_bounds(bin_data, ix1) ) &
    count = bin_data%count(general_bin_index(bin_data, ix1))
endif
end function

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
! Function for checking bounds

function general_bin_index_in_bounds(bin_data, ix1, ix2, ix3) result(in_bounds)
implicit none
type(general_bin_struct) :: bin_data
integer ix1
integer, optional :: ix2, ix3
logical :: in_bounds
in_bounds = .true.
if (present(ix3)) then
  if (ix3 < 1 .or. ix3 > bin_data%n(3) ) in_bounds = .false.
endif
if (present(ix2)) then
  if (ix2 < 1 .or. ix2 > bin_data%n(2) ) in_bounds = .false.
endif
if (ix1 < 1 .or. ix1 > bin_data%n(1) ) in_bounds = .false.
end function

end module
