!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine write_power_header (iu, file, gen_params)

  use sr_struct
  use sr_interface

  implicit none

  type (general_param_struct) gen_params
  character*(*) file
  integer iu

!

  open (unit = iu, file = file, carriagecontrol = 'list')
  write (iu, *) 'Lattice: ', gen_params%lattice
  write (iu, *) 'I_beam    =', gen_params%i_beam,    ' ! Amps/beam'
  write (iu, *) 'Epsilon_y =', gen_params%epsilon_y, ' ! Vertical emittance'

  write (iu, '(3(/,2x,a))') &
'          Segment                                  |        Source', &
'Ix Name         S_seg X_seg  P/len   P/Area  P_tot | N E/P  Name           S',&
'                 (m)   (m)   (W/m)  (kW/m^2)   (W) |                      (m)'

end subroutine
