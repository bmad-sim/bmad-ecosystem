program real_8_example
use pointer_lattice   ! Read in structure definitions, etc.
implicit none

type (real_8) r8      ! Define a real_8 variable named r8
real(dp) x            ! Define a double precision number

!

nice_taylor_print = .true.    ! Shorten "call print" output
call init (only_2d0, 3, 0)    ! Initialize: #Vars = 2, Order = 3

x = 0.1d0
call alloc(r8)          ! Initialize memory for r8
r8 = x                  ! This will make r8 act as a real

print '(/,a)', 'r8 is now acting as a real:'
call print (r8)

r8 = 0.7d0 + dz_8(1) + 2*dz_8(2)**3 ! Init r8 as a Taylor series
print '(/,a)', 'r8 is now acting as a Taylor series:'
call print(r8)

r8 = r8**4  ! Raise the Taylor series to the 4th power
print '(/,a)', 'This is r8^4:'
call print (r8)

call kill(r8)
end program
