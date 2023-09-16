!+
! Subroutine write_beam_floor_positions (file_name, beam, ele, new_file)
!
! Routine to write a file of beam positions in global floor coordinates.
!
! Input:
!   file_name     -- character(*): Name of file.
!   beam          -- beam_struct: Beam to write
!   ele           -- ele_struct: Element that the beam is at.
!   new_file      -- logical, optional: New file or append? Default = True.
!-

subroutine write_beam_floor_positions (file_name, beam, ele, new_file)

use bmad_routine_interface, dummy => write_beam_floor_positions

implicit none

type (beam_struct), target :: beam
type (bunch_struct), pointer :: bunch
type (coord_struct), pointer :: p
type (ele_struct) :: ele

real(rp) floor_vec(6)
integer iu, ib, ip

logical, optional :: new_file
logical error, append

character(*) file_name
character(200) full_name
character(*), parameter :: r_name = 'write_beam_floor_positions'

!

iu = lunget()
call fullfilename (file_name, full_name)

if (logic_option(.true., new_file)) then
  open (iu, file = full_name)
else
  open (iu, file = full_name, access = 'append')
endif

!

do ib = 1, size(beam%bunch)
  bunch => beam%bunch(ib)

  write (iu, '(a, i8)')     '# ix_ele       =', ele%ix_ele
  write (iu, '(a, i8)')     '# ix_bunch     =', ib
  write (iu, '(a, i8)')     '# n_particle   =', size(bunch%particle)
  write (iu, '(a, a)')      '# species      = ', trim(species_name(bunch%particle(1)%species))
  write (iu, '(a, es16.9)') '# bunch_charge =', bunch%charge_tot
  write (iu, '(a, es16.9)') '# z_center     =', bunch%z_center
  write (iu, '(a, es16.9)') '# t_center     =', bunch%t_center
  write (iu, '((a, a7), (a12, 5a18, 6x), 2x, (a10, 2a15, 5x), 2x, (a13, 2a15, 2x), 2x, a)') '#', 'Ix', 'X', 'Y', 'Z', &
          'Px (eV)', 'Py (eV)', 'Pz (eV)', 'Time', 'S-pos', 'Charge', 'Spin_x', 'Spin_y', 'Spin_z', 'State' 
  do ip = 1, size(bunch%particle)
    p => bunch%particle(ip)
    floor_vec = orbit_to_floor_phase_space(p, ele)
    write (iu, '(i8, 6es18.10, 2x, 3es15.7, 2x, 3es15.7, 2x, a)') ip, floor_vec(1:3), &
          floor_vec(4:6)*p%p0c*(1+p%vec(6)),  p%t, p%s, p%charge, p%spin, trim(coord_state_name(p%state))
  enddo
enddo

end subroutine
