subroutine check_aperture ( u, hit )

  use sr_mod
  use cesrv_struct

  implicit none

  type (universe_struct) u
  logical hit

  integer i
  real(rp) epsx, epsy, sige, n_sigma

  !

  hit  = .false.
  epsx = u%ring%a%emit
  epsy = u%ring%b%emit
  sige = u%energy_data%d1%d(1)%model

  type '(a, $)', ' Number of sigma to use?'
  read *, n_sigma
  type *, ' epsx is : ',epsx
  type *, ' epsy is : ',epsy
  type *, ' dE/E is : ',sige

  call bmad_parser2( 'u:[cesr.bmad.layout]aperture.bmad',u%ring )

  do i=1, u%ring%n_ele_max
    if ( u%ring%ele(i)%value(x1_limit$) == 0 )  cycle
  enddo

  type *, 'apertures loaded'

end subroutine
