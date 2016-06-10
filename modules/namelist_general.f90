module namelist_general

! provides:
! lat_file, use_hybrid, use_line, periodicity

use bmad

implicit none

character*100 lat_file
logical use_hybrid

namelist / general /    lat_file, &
                        use_hybrid

end module
