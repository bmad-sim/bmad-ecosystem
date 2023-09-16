!+
! Function virtual_memory_usage() result (usage)
!
! Routine to return the amount of virtual memory currently used by the process calling this routine.
! This routine is useful for memory debugging.
! This routine works with Linux and will definitely not work with Windows nor OSX.
!
! Output:
!   usage   -- integer: Virtual memory usage in KB. Will be set to -1 or -2 if there is a read error.
!-

function virtual_memory_usage() result (usage)

use sim_utils_interface, dummy => virtual_memory_usage

implicit none

integer usage, iu, ios
character(64) line

!

usage = -1

iu = lunget()
open (iu, file = '/proc/self/status', action = 'READ')
do
  read(iu, '(a)', iostat = ios) line
  if (ios /= 0) exit
  if (line(1:7) /= 'VmSize:') cycle
  read (line(8:), *, iostat = ios) usage
  if (ios /= 0) usage = -2
  exit
enddo

close(iu)

end function

