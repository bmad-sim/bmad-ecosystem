subroutine break_wall_into_segments (wall, seg_len_max)

use synrad_struct
use synrad_interface

implicit none

type (wall_struct) wall

integer i_seg, ip, n_seg, i, n

real(rp) seg_len_max, wall_len, rr

! If the distance between two wall points is too large then break the interval
! into more than 1 segment.
! x_seg and s_seg are the (x,s) coords of the center of the segment

! First count the number of segments needed and allocate

n = 0
do ip = 1, wall%n_pt_tot
  wall_len = sqrt((wall%pt(ip)%x - wall%pt(ip-1)%x)**2 + &
                             (wall%pt(ip)%s - wall%pt(ip-1)%s)**2)
  n_seg = 1 + wall_len / seg_len_max
  n = n + n_seg
enddo

if (allocated(wall%seg)) deallocate(wall%seg)
allocate (wall%seg(n))

! Now fill in the information.

i_seg = 0
do ip = 1, wall%n_pt_tot
  wall_len = sqrt((wall%pt(ip)%x - wall%pt(ip-1)%x)**2 + &
                             (wall%pt(ip)%s - wall%pt(ip-1)%s)**2)
  n_seg = 1 + wall_len / seg_len_max
  wall%pt(ip)%n_seg = n_seg
  wall%pt(ip)%ix_seg = i_seg
  do i = i_seg+1, i_seg+n_seg
    wall%seg(i)%ix_pt = ip
    rr = (i - i_seg - 0.5) / n_seg
    wall%seg(i)%x_mid = wall%pt(ip-1)%x * (1 - rr) + wall%pt(ip)%x * rr
    wall%seg(i)%s_mid = wall%pt(ip-1)%s * (1 - rr) + wall%pt(ip)%s * rr
    rr = float(i - i_seg) / n_seg
    wall%seg(i)%x = wall%pt(ip-1)%x * (1 - rr) + wall%pt(ip)%x * rr
    wall%seg(i)%s = wall%pt(ip-1)%s * (1 - rr) + wall%pt(ip)%s * rr
    if (i /= 1) wall%seg(i)%len = &
                      sqrt((wall%seg(i)%s - wall%seg(i-1)%s)**2 + &
                      (wall%seg(i)%x - wall%seg(i-1)%x)**2)
  end do
  i_seg = i_seg + n_seg
enddo

wall%n_seg_tot = i_seg

end subroutine
