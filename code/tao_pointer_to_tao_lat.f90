!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Function tao_pointer_to_tao_lat (u, lat_type) result (tao_lat)
!
! Routine to set a pointer to a tao_lat.
!
! Also see:
!   tao_pointer_to_universe
!
! Input:
!   u         -- Tao_universe_struct: Universe to work with
!   lat_type  -- Integer, optional: model$ (default), design$, or base$.
!
! Output:
!   tao_lat   -- tao_lattice_struct, pointer: Tao_lat pointer. 
!                   Points to u%model, u%design, or u%base
!-

function tao_pointer_to_tao_lat (u, lat_type) result (tao_lat)

use tao_struct

implicit none

type (tao_universe_struct), target :: u
type (tao_lattice_struct), pointer :: tao_lat
integer, optional :: lat_type
character(28) :: r_name = 'tao_pointer_to_tao_lat'

!

if (.not. present(lat_type)) then
  tao_lat => u%model
  return
endif

select case (lat_type)
case (design$)
  tao_lat => u%design
case (model$)
  tao_lat => u%model
case (base$)
  tao_lat => u%base
case default
  call err_exit
end select

end function tao_pointer_to_tao_lat

