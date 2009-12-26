!+
! Program to translate a bmad lattice to a format for
! Andy Wolski to use in his ORM analysis.
!
! Note: The original program is from [cesr.orm.code]bmad_lat_to_wolski.f90
!-

program lat_to_wolski_orm_format

use bmad
use cesr_basic_mod

implicit none

type (lat_struct), target :: ring
type (ele_struct), pointer :: ele

real(rp) l, angle, e1, e2, k1, k2, ks, tilt, hkick, vkick

integer i, ix, j, j1, j2

character(40) db_lat, lat_name
character(100) lat_file, line, file

!

db_lat = ''
call choose_cesr_lattice (lat_name, lat_file, db_lat, ring) 

! lattice file

file = trim(lat_name) // '.bmad'
print *, 'File: ', trim(file)

open (1, file = file, status = 'new', carriagecontrol = 'list', recl = 500)

write (1, '(2a17, 3x, 32a14)') 'Name', 'Type',  'S',  'L', &
      'Angle', 'E1', 'E2', 'K1', 'K2', 'Ks', 'Tilt', 'HKick', 'VKick', &
      'M11', 'M12 ...'


do i = 1, ring%n_ele_track

  ele => ring%ele(i)

  l = ele%value(l$)
  hkick = ele%value(hkick$)
  vkick = ele%value(vkick$)

  angle = 0
  e1 = 0
  e2 = 0
  k1 = 0
  k2 = 0
  ks = 0
  tilt = 0

  select case (ele%key)
  case (solenoid$)
    ks = ele%value(ks$)
    tilt = ele%value(tilt$)
  case (sol_quad$)
    ks = ele%value(ks$)
    k1 = ele%value(k1$)
    tilt = ele%value(tilt$)
  case (quadrupole$)
    k1 = ele%value(k1$)
    tilt = ele%value(tilt$)
  case (sbend$)
    angle = ele%value(angle$)
    e1 = ele%value(e1$)
    e2 = ele%value(e2$)
  case (sextupole$)
    k2 = ele%value(k2$)
    tilt = ele%value(tilt$)
  end select

  write (1, '(a, 1x, a, 32es14.5)') ele%name, key_name(ele%key), ele%s, l, &
            angle, e1, e2, k1, k2, ks, tilt, hkick, vkick, &
            ((ele%mat6(j1,j2), j2 = 1, 4), j1 = 1, 4) 

enddo

close (1)

! steering file

file = 'steerings.list'
open (1, file = file, status = 'new', carriagecontrol = 'list')

do i = ring%n_ele_track+1, ring%n_ele_max

  ele => ring%ele(i)

  if (ele%lord_status /= overlay_lord$) cycle
  if (ele%component_name /= 'HKICK' .and. ele%component_name /= 'VKICK') cycle
  if (ele%name(1:1) /= 'H' .and. ele%name(1:1) /= 'V') cycle

  do j = ele%ix1_slave, ele%ix2_slave
    ix = ring%control(j)%ix_slave
    j2 = 17 * (j - ele%ix1_slave) + 1
    write (line(j2:), '(a)') ring%ele(ix)%name
  enddo

  write (1, '(3a)') ele%name, ' ', trim(line)

enddo

close(1)

end program
