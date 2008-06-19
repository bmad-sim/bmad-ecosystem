!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine write_power_header (iu, file, gen_params)

  use sr_struct
  use sr_interface

  implicit none

  type (synrad_param_struct) gen_params
  character(*) file
  integer iu

!

  open (unit = iu, file = file, carriagecontrol = 'list')
  write (iu, *) 'Lattice: ', gen_params%lat_file
  write (iu, *) 'I_beam    =', gen_params%i_beam,    ' ! Amps/beam'
  write (iu, *) 'Epsilon_y =', gen_params%epsilon_y, ' ! Vertical emittance'

  write (iu, '(3(/,2x,a))') &
'          Segment                                  ', &
'  Ix  Name          S_seg      X_seg     P/len      P/Area     P_tot     Phot/sec      A Beta    B Beta    Ele Type       Relevant               Ele',&
'                     (m)        (m)      (W/m)     (W/mm^2)      (W)      (1/s)         (m)        (m)     at s_mid       Attribute              Name'

end subroutine
