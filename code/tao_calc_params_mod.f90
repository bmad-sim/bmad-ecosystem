module tao_calc_params_mod

use tao_mod
use macro_utils_mod
use beam_mod


contains

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine tao_calc_params (u, ix_ele)
!
! Finds all the lattice and particle parameters at element ix_ele
!
! Input:
!  u               -- tao_universe_struct
!  ix_ele          -- Integer: Find at this ele
!
! Output:
!  u%model%ele%                      -- ele_struct
!      x,y,z                 -- twiss parameters
!      vec0
!      mat6
!      c_mat
!  u%maco_beam%macro_params  -- macro_bunch_params_struct
!  u%beam%params             -- bunch_struct_params
!
!-

subroutine tao_calc_params (u, ix_ele)

implicit none

type (tao_universe_struct), target :: u
type (ring_struct), pointer :: lat

integer ix_ele

lat => u%model

if (s%global%track_type .eq. "single" .or. &
      u%model%param%lattice_type .eq. circular_lattice$) then
  call make_mat6 (lat%ele_(ix_ele), lat%param, u%model_orb(ix_ele-1), &
                    u%model_orb(ix_ele), .true.)
  call twiss_propagate1 (lat%ele_(ix_ele-1), lat%ele_(ix_ele))
elseif (s%global%track_type .eq. "beam") then
  call make_mat6 (lat%ele_(ix_ele), lat%param, u%model_orb(ix_ele-1), &
                    u%model_orb(ix_ele), .true.)
  call twiss_propagate1 (lat%ele_(ix_ele-1), lat%ele_(ix_ele))
  call calc_bunch_params (u%beam%beam%bunch(s%global%bunch_to_plot), &
                                lat%ele_(ix_ele), u%beam%params)
elseif (s%global%track_type .eq. "macro") then
  call make_mat6 (lat%ele_(ix_ele), lat%param, u%model_orb(ix_ele-1), &
                    u%model_orb(ix_ele), .true.)
  call twiss_propagate1 (lat%ele_(ix_ele-1), lat%ele_(ix_ele))
  call calc_macro_bunch_params (u%macro_beam%beam%bunch(s%global%bunch_to_plot), &
                                lat%ele_(ix_ele), u%macro_beam%params)
endif

end subroutine tao_calc_params
    

end module tao_calc_params_mod
