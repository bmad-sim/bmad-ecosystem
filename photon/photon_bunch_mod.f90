module photon_bunch_mod

use photon_target_mod
use track1_photon_mod

implicit none

contains

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine track_photon_bunch (bunch, branch, ix_start, ix_end)
!
! Routine to track a photon_bunch through a single lattice element.
!
! Input:
!   bunch       -- bunch_struct: Starting bunch position.
!   branch      -- branch_struct: lattice branch to track through.
!   ix_start    -- integer, optional: Starting element index. Default is zero.
!                    Tracking starts at the downstream end of the element.
!   ix_end      -- integer, optional: Ending element index. Default is zero.
!                    Tracking ends at the downstream end of the element.
!
! Output:
!   bunch       -- Bunch_struct: Ending bunch position.
!-

subroutine track_photon_bunch (bunch, branch, ix_start, ix_end)

type (bunch_struct) bunch
type (branch_struct) branch
type (coord_struct) orb
type (ele_struct), pointer :: ele
type (ele_struct), pointer :: targ_ele
type (surface_grid_struct), pointer :: gr
type (photon_target_struct), pointer :: target

integer, optional :: ix_start, ix_end
integer i, ie, ie2, npart0, i_start, i_end, lb(2), ub(2), sz(2)
integer ix, iy, nx, ny

logical setup_grid_target

character(*), parameter :: r_name = 'track_photon_bunch'

! Deterministic targeting?

i_start = integer_option(0, ix_start)
i_end   = integer_option(branch%n_ele_track, ix_end)
npart0 = size(bunch%particle)

ie = i_start
do 
  ie = ie + 1
  if (ie > i_end) exit

  ele => branch%ele(ie)

  setup_grid_target = grid_target_found(ele)
  if (setup_grid_target) then
    target => ele%photon%target
    target%deterministic_grid = .true.
    targ_ele => pointer_to_ele(ele%branch%lat, target%ele_loc)
    gr => targ_ele%photon%surface%grid
    lb = lbound(gr%pt); ub = ubound(gr%pt)
    sz = ub + 1 - lb 
    gr => targ_ele%photon%surface%grid
    gr%pt = surface_grid_pt_struct()
  endif

  ! Track particles

  if (setup_grid_target) then

    ! Track to target and gather statistics.

    do i = 1, size(bunch%particle)
      do ix = lb(1), ub(1)
      do iy = lb(2), ub(2)

        target%ix_grid = ix
        target%iy_grid = iy
        orb = bunch%particle(i)

        do ie2 = ie, targ_ele%ix_ele - 1
          call track1 (orb, branch%ele(ie2), branch%param, orb)
          if (orb%state /= alive$) exit
        enddo

        if (orb%state == alive$) then
          call photon_add_to_detector_statistics (bunch%particle(i), orb, targ_ele, nx, ny)
          if (nx /= ix .or. ny /= iy) then
            call out_io (s_fatal$, r_name, 'BAD TARGETING')
            if (global_com%exit_on_error) call err_exit
          endif
        endif

        bunch%particle(i) = orb

      enddo
      enddo
    enddo

    ! Generate photons from target grid.
    ! Don't need to do this if there is a grid target

    ie = targ_ele%ix_ele - 1
    if (grid_target_found(targ_ele)) then      
      call reallocate_bunch (bunch, sz(1)*sz(2))
    else
      call reallocate_bunch (bunch, npart0)
    endif

    do i = 1, size(bunch%particle)
      ix = modulo(i-1, sz(1)) + lb(1)
      iy = modulo((i-1)/sz(1), sz(2)) + lb(2)
      bunch%particle(i)%vec(1) = ix * gr%dr(1) + gr%r0(1)
      bunch%particle(i)%vec(3) = iy * gr%dr(2) + gr%r0(2)
      bunch%particle(i)%field(1) = abs(gr%pt(ix,iy)%E_x)
      bunch%particle(i)%field(2) = abs(gr%pt(ix,iy)%E_y)
      bunch%particle(i)%phase(1) = atan2(aimag(gr%pt(ix,iy)%E_x), real(gr%pt(ix,iy)%E_x))
      bunch%particle(i)%phase(2) = atan2(aimag(gr%pt(ix,iy)%E_y), real(gr%pt(ix,iy)%E_y))
    enddo

  ! Simple case without a grid target.
  ! Just track through a single element

  else
    do i = 1, size(bunch%particle)
      call track1 (bunch%particle(i), ele, branch%param, bunch%particle(i))
    enddo
  endif

enddo

!-------------------------------------------------------------------------------------
contains

function grid_target_found (ele) result (target_found)

type (ele_struct) ele
logical target_found

!

target_found = .false.
if (.not. associated(ele%photon)) return
if (ele%photon%target%type == grided$) target_found = .true.

end function grid_target_found

end subroutine track_photon_bunch

end module
