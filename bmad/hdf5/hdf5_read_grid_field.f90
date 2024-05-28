!+
! Subroutine hdf5_read_grid_field (file_name, ele, g_field, err_flag, pmd_header, combine)
!
! Routine to read an hdf5 file that encodes an array of grid_field structures.
!
! Input:
!   file_name     -- character(*): File to create.
!   ele           -- ele_struct: Element associated with the map.
!   g_field(:)    -- grid_field_struct, pointer: Array with possible existing data.
!   combine       -- logical, optional: If False (the default), discard existing data in input g_field.
!                     If True, output g_field(:) array combines old and new data.
!
! Ouput:
!   g_field(:)    -- grid_field_struct, pointer: Grid field array.
!   err_flag      -- logical: Set True if there is an error. False otherwise.
!   pmd_header    -- pmd_header_struct, optional: Extra info like file creation date.
!-

subroutine hdf5_read_grid_field (file_name, ele, g_field, err_flag, pmd_header, combine)

use hdf5_openpmd_mod
use bmad_interface, dummy => hdf5_read_grid_field

implicit none

type (grid_field_struct), pointer :: g_field(:)
type (grid_field_struct), pointer :: gf
type (ele_struct) ele
type (pmd_header_struct), optional :: pmd_header
type (pmd_header_struct) :: pmd_head
type(H5O_info_t) :: infobuf 

integer(hid_t) f_id, f2_id, z_id, complex_t
integer i, j, k, it, n0(3), n1(3), iver, h5_err, n_grid, storage_type, max_corder
integer n_links, h5_stat
integer(hsize_t) idx
integer(hsize_t) this_size

logical, optional :: combine
logical err_flag, err

character(*) file_name
character(40) master_name
character(100) c_name, name
character(*), parameter :: r_name = 'hdf5_read_grid_field'

!

err_flag = .true.
call hdf5_open_file (file_name, 'READ', f_id, err);  if (err) return
call pmd_init_compound_complex(complex_t)

call hdf5_read_attribute_alloc_string (f_id, 'externalFieldPath', pmd_head%basePath, err, .true.);         if (err) return
call hdf5_read_attribute_alloc_string (f_id, 'openPMD', pmd_head%openPMD, err, .true.);                    if (err) return
call hdf5_read_attribute_alloc_string (f_id, 'openPMDextension', pmd_head%openPMDextension, err, .true.);  if (err) return

call hdf5_read_attribute_alloc_string (f_id, 'software', pmd_head%software, err, .false.)
call hdf5_read_attribute_alloc_string (f_id, 'softwareVersion', pmd_head%softwareVersion, err, .false.)
call hdf5_read_attribute_alloc_string (f_id, 'date', pmd_head%date, err, .false.)

if (present(pmd_header)) pmd_header = pmd_head

! Deallocate if needed

if (associated(g_field) .and. .not. logic_option(.false., combine)) then
  call unlink_fieldmap(grid_field = g_field)
endif

! No "%T" case.

it = index(pmd_head%basePath, '/%T/')
if (it == 0) then
  call allocate_this_field (g_field, n_grid)
  z_id = hdf5_open_group(f_id, pmd_head%basePath, err, .true.); if (err) return
  call read_this_field (z_id, ele, g_field(n_grid), err)
  call H5Gclose_f(z_id, h5_err)
  return
endif

! Count grids

f2_id = hdf5_open_group(f_id, pmd_head%basePath(1:it), err, .true.)
call H5Gget_info_f (f2_id, storage_type, n_links, max_corder, h5_err)
do idx = 0, n_links-1
  call H5Lget_name_by_idx_f (f2_id, '.', H5_INDEX_NAME_F, H5_ITER_INC_F, idx, c_name, h5_err, this_size)
  call to_f_str(c_name, name)
  call H5Oget_info_by_name_f(f2_id, name, infobuf, h5_stat)
  if (infobuf%type /= H5O_TYPE_GROUP_F) cycle    ! Ignore non-group elements.
  if (.not. is_integer(name)) then               ! This assumes basepath uses "/%T/" to designate different grids.
    call out_io (s_warn$, r_name, 'NAME OF DIRECTORY IN PATH IS NOT AN INTEGER: ' // quote(name))
    cycle
  endif
  call allocate_this_field (g_field, n_grid)
  z_id = hdf5_open_group(f2_id, name, err, .true.);         if (err) return
  call read_this_field(z_id, ele, g_field(n_grid), err);     if (err) return
  call h5gclose_f(z_id, h5_err)
enddo
call H5Gclose_f(f2_id, h5_err)

!

call pmd_kill_compound_complex(complex_t)
call h5fclose_f(f_id, h5_err)

err_flag = .false.

!--------------------------------------------------------------------------------
contains

subroutine allocate_this_field (g_field, n_grid)

type (grid_field_struct), pointer :: g_field(:), g_temp(:)
integer n_grid

!

if (.not. associated(g_field)) then
  allocate (g_field(1))
  n_grid = 1
  return
endif

n_grid = size(g_field) + 1
g_temp => g_field
allocate (g_field(n_grid))
g_field(1:n_grid-1) = g_temp

end subroutine allocate_this_field

!--------------------------------------------------------------------------------
! contains

subroutine read_this_field(root_id, ele, gf, err_flag)

type (ele_struct) ele
type (grid_field_struct), target :: gf
type (grid_field_pt1_struct), allocatable, target :: gpt(:,:,:)
type (grid_field_pt1_struct), pointer :: gptr(:,:,:)

real(rp) field_scale, total_scale, rho

integer indx(3), lb(3), g_size(3), ub(3)
integer(hid_t) root_id, z_id

character(40) name
character(8) component_name(3), B_name(3), E_name(3)
logical err_flag, error, b_field_here, e_field_here

!

err_flag = .true.

call hdf5_read_attribute_string (root_id, 'masterParameter', name, error, .false.)
if (name == '') then
  gf%master_parameter = 0
  call hdf5_read_attribute_real(root_id, 'fieldScale', gf%field_scale, error, .false., dflt_value = 1.0_rp)
else
  gf%master_parameter = attribute_index(ele, upcase(name))
  call hdf5_read_attribute_real(root_id, 'componentFieldScale', gf%field_scale, error, .false.)
endif

call hdf5_read_attribute_string (root_id, 'gridGeometry', name, error, .false.)
select case (name)
case ('rectangular')
  gf%geometry = xyz$
  indx = [1, 2, 3]
case ('cylindrical')
  gf%geometry = rotationally_symmetric_rz$
  indx = [1, 3, 2]
case default
  call out_io (s_error$, r_name, 'gridGeometry name not recognized: ' // name)
  return
end select

call hdf5_read_attribute_string (root_id, 'eleAnchorPt', name, error, .false.); if (error) return
call match_word (name, anchor_pt_name(1:), gf%ele_anchor_pt, error, .false.);   if (error) return

call hdf5_read_attribute_int (root_id, 'harmonic', gf%harmonic, error, .false.)
call hdf5_read_attribute_real (root_id, 'RFphase', gf%phi0_fieldmap, error, .false.)
if (gf%harmonic /= 0) then
  call hdf5_read_attribute_real (root_id, 'fundamentalFrequency', ele%value(rf_frequency$), error, .false.)
  if (ele%key == lcavity$) then
    gf%phi0_fieldmap = gf%phi0_fieldmap / gf%harmonic
  else
    gf%phi0_fieldmap = 0.25_rp - gf%phi0_fieldmap / gf%harmonic
  endif
endif

call hdf5_read_attribute_int (root_id, 'interpolationOrder', gf%interpolation_order, error, .false., dflt_value = 1)
call hdf5_read_attribute_int (root_id, 'gridLowerBound', lb, error, .true.);          if (error) return
call hdf5_read_attribute_int (root_id, 'gridSize', g_size, error, .true.);            if (error) return
call hdf5_read_attribute_real (root_id, 'gridOriginOffset', gf%r0, error, .true.);    if (error) return
call hdf5_read_attribute_real (root_id, 'gridSpacing', gf%dr, error, .true.);         if (error) return
call hdf5_read_attribute_real (root_id, 'gridCurvatureRadius', rho, error, .false.)
if (error) call hdf5_read_attribute_real (root_id, 'curvedRefFrame', rho, error, .false.)  ! Old style

ub = lb + g_size - 1
gf%curved_ref_frame = (rho /= 0)
gf%dr = gf%dr(indx)

!

allocate (gf%ptr)
gf%ptr%file = file_name

select case (gf%geometry)
case (xyz$)
  component_name = xyz_axislabels
  B_name = [character(6):: 'Bx', 'By', 'Bz']     ! Old style
  E_name = [character(6):: 'Ex', 'Ey', 'Ez']     ! Old style
  allocate (gf%ptr%pt(lb(1):ub(1), lb(2):ub(2), lb(3):ub(3)))
  gptr => gf%ptr%pt
case (rotationally_symmetric_rz$)
  component_name = rthetaz_axislabels
  B_name = [character(6):: 'Br', 'Btheta', 'Bz']       ! Old style
  E_name = [character(6):: 'Er', 'Etheta', 'Ez']       ! Old style
  allocate(gf%ptr%pt(lb(1):ub(1), lb(3):ub(3), lb(2):ub(2)))
  allocate (gpt(lb(1):ub(1), lb(2):ub(2), lb(3):ub(3)))
  gptr => gpt
end select

b_field_here = .false.
e_field_here = .false.

if (hdf5_exists(root_id, 'magneticField', error, .false.)) then
  b_field_here = .true.
  z_id = hdf5_open_group(root_id, 'magneticField', err, .true.);  if (err) return
  do i = 1, 3
    if (hdf5_exists(z_id, component_name(i), err, .false.)) then
      call pmd_read_complex_dataset(z_id, trim(component_name(i)), complex_t, 1.0_rp, component_name, gptr%B(i), error)
    else
      gptr%B(i) = 0
    endif
    if (gf%geometry == rotationally_symmetric_rz$) gf%ptr%pt(:,:,lb(2))%B(i) = gptr(:,lb(2),:)%B(i)
  enddo
  call H5Gclose_f(z_id, h5_err)
endif

if (hdf5_exists(root_id, 'electricField', error, .false.)) then
  e_field_here = .true.
  z_id = hdf5_open_group(root_id, 'electricField', err, .true.);  if (err) return
  do i = 1, 3
    if (hdf5_exists(z_id, component_name(i), err, .false.)) then
      call pmd_read_complex_dataset(z_id, trim(component_name(i)), complex_t, 1.0_rp, component_name, gptr%E(i), error)
    else
      gptr%E(i) = 0
    endif
    if (gf%geometry == rotationally_symmetric_rz$) gf%ptr%pt(:,:,lb(2))%E(i) = gptr(:,lb(2),:)%E(i)
  enddo
  call H5Gclose_f(z_id, h5_err)
endif

! Old style

do i = 1, 3
  if (hdf5_exists(root_id, B_name(i), error, .false.)) then
    call pmd_read_complex_dataset(root_id, trim(B_name(i)), complex_t, 1.0_rp, component_name, gptr%B(i), error)
    if (gf%geometry == rotationally_symmetric_rz$) gf%ptr%pt(:,:,lb(2))%B(i) = gptr(:,lb(2),:)%B(i)
    b_field_here = .true.
  endif

  if (hdf5_exists(root_id, E_name(i), error, .false.)) then
    call pmd_read_complex_dataset(root_id, trim(E_name(i)), complex_t, 1.0_rp, component_name, gptr%E(i), error)
    if (gf%geometry == rotationally_symmetric_rz$) gf%ptr%pt(:,:,lb(2))%E(i) = gptr(:,lb(2),:)%E(i)
    e_field_here = .true.
  endif
enddo

!

if (e_field_here .and. b_field_here) then
  gf%field_type = mixed$
elseif (e_field_here) then
  gf%field_type = electric$
elseif (b_field_here) then
  gf%field_type = magnetic$
else
  gf%field_type = mixed$
endif

err_flag = .false.

end subroutine read_this_field

end subroutine hdf5_read_grid_field
