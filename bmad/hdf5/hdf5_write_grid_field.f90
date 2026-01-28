!+
! Subroutine hdf5_write_grid_field (file_name, ele, g_field, err_flag)
!
! Routine to create an hdf5 file encoding an array of grid_field structures.
! Note: Conventionally, the file name should have an ".h5" suffix.
!
! Input:
!   file_name     -- character(*): File to create.
!   ele           -- ele_struct: Element associated with the map.
!   g_field(:)    -- grid_field_struct: Grid field.
!
! Output:
!   err_flag      -- logical: Set True if there is an error. False otherwise.
!-

subroutine hdf5_write_grid_field (file_name, ele, g_field, err_flag)

use hdf5_openpmd_mod
use bmad_interface, dummy => hdf5_write_grid_field

implicit none

type (grid_field_struct), target :: g_field(:)
type (grid_field_struct), pointer :: gf
type (grid_field_pt1_struct), allocatable, target :: gpt(:,:,:)
type (grid_field_pt1_struct), pointer :: gptr(:,:,:)
type (ele_struct) ele

integer i, j, k, n, ix, im, ig, igf, h5_err, indx(3), n3
integer(hid_t) f_id, r_id, b_id, b2_id, b3_id, z_id, complex_t
logical err_flag, err

character(*) file_name
character(*), parameter :: r_name = 'dhf5_write_grid_field'
character(40) this_path
character(28) date_time, root_path, sub_path
character(8) :: component_name(3)

!

err_flag = .true.
call hdf5_open_file (file_name, 'WRITE', f_id, err);  if (err) return
call pmd_init_compound_complex(complex_t)

call date_and_time_stamp (date_time, .true., .true.)
root_path = '/ExternalFieldMesh/'
sub_path = '%T/'

call hdf5_write_attribute_string(f_id, 'dataType',          'Bmad:grid_field', err)
call hdf5_write_attribute_string(f_id, 'openPMD',           '2.0.0', err)
call hdf5_write_attribute_string(f_id, 'openPMDextension',  'BeamPhysics;SpeciesType', err)
call hdf5_write_attribute_string(f_id, 'externalFieldPath',  trim(root_path) // trim(sub_path), err)
call hdf5_write_attribute_string(f_id, 'software',          'Bmad', err)
call hdf5_write_attribute_string(f_id, 'softwareVersion',   '1.0', err)
call hdf5_write_attribute_string(f_id, 'date',              date_time, err)

call H5Gcreate_f(f_id, trim(root_path), b_id, h5_err)

do igf = 1, size(g_field)
  gf => g_field(igf)

  write (sub_path, '(i0)') igf
  call H5Gcreate_f(b_id, trim(sub_path), b2_id, h5_err)

  gptr => gf%ptr%pt
  select case (gf%geometry)
  case (xyz$)
    component_name = xyz_axislabels
    indx = [1, 2, 3]
    call hdf5_write_attribute_string(b2_id,  'gridGeometry', 'rectangular', err)

  case (rotationally_symmetric_rz$)
    component_name = rthetaz_axislabels
    indx = [1, 3, 2]
    allocate(gpt(lbound(gptr, 1):ubound(gptr, 1), lbound(gptr, 3):ubound(gptr, 3), lbound(gptr, 2):ubound(gptr, 2)))
    gptr => gpt
    n3 = lbound(gf%ptr%pt, 3)
    do i = 1, 3
      gptr(:,n3,:)%B(i) = gf%ptr%pt(:,:,n3)%B(i)
      gptr(:,n3,:)%E(i) = gf%ptr%pt(:,:,n3)%E(i)
    enddo
    call hdf5_write_attribute_string(b2_id, 'gridGeometry', 'cylindrical', err)
  end select

  im = gf%master_parameter 
  if (im == 0) then
    call hdf5_write_attribute_real(b2_id,    'fieldScale',           gf%field_scale, err)
  else
    call hdf5_write_attribute_string(b2_id,  'masterParameter',      attribute_name(ele, im), err)
    call hdf5_write_attribute_real(b2_id,    'fieldScale',           gf%field_scale * ele%value(im), err)
  endif
  call hdf5_write_attribute_real(b2_id,    'componentFieldScale',  gf%field_scale, err)

  
  if (gf%harmonic /= 0) then
    call hdf5_write_attribute_real(b2_id, 'fundamentalFrequency', ele%value(rf_frequency$), err)
    if (ele%key == lcavity$) then
      call hdf5_write_attribute_real(b2_id, 'RFphase', gf%harmonic * gf%phi0_fieldmap, err)
    else
      call hdf5_write_attribute_real(b2_id, 'RFphase', gf%harmonic * (0.25_rp - gf%phi0_fieldmap), err)      
    endif
  endif

  call hdf5_write_attribute_string(b2_id,  'axisLabels', component_name(3:1:-1), err) ! Write labels in reverse since using Fortran
  call hdf5_write_attribute_string(b2_id,  'eleAnchorPt',          downcase(anchor_pt_name(gf%ele_anchor_pt)), err)
  call hdf5_write_attribute_real(b2_id,    'gridOriginOffset',     gf%r0, err)
  call hdf5_write_attribute_real(b2_id,    'gridSpacing',          gf%dr(indx), err)
  call hdf5_write_attribute_int(b2_id,     'harmonic',             gf%harmonic, err)
  call hdf5_write_attribute_int(b2_id,     'interpolationOrder',   gf%interpolation_order, err)
  call hdf5_write_attribute_int(b2_id,     'gridLowerBound',       lbound(gptr), err)
  call hdf5_write_attribute_int(b2_id,     'gridSize',             shape(gptr), err)

  if ((ele%key == sbend$ .or. ele%key == rf_bend$) .and. gf%curved_ref_frame) then
    call hdf5_write_attribute_real(b2_id,     'gridCurvatureRadius',       ele%value(rho$), err)
  else
    call hdf5_write_attribute_real(b2_id,     'gridCurvatureRadius',       0.0_rp, err)
  endif

  !

  if (gf%field_type == magnetic$ .or. gf%field_type == mixed$) then
    call H5Gcreate_f(b2_id, 'magneticField/', b3_id, h5_err)
    do i = 1, 3
      call pmd_write_complex_to_dataset (b3_id, component_name(i), complex_t, component_name(i), unit_tesla, gptr%B(i), err)
    enddo
    call H5Gclose_f(b3_id, h5_err)
  endif

  if (gf%field_type == electric$ .or. gf%field_type == mixed$) then
    call H5Gcreate_f(b2_id, 'electricField/', b3_id, h5_err)
    do i = 1, 3
      call pmd_write_complex_to_dataset (b3_id, component_name(i), complex_t, component_name(i), unit_V_per_m, gptr%E(i), err)
    enddo
    call H5Gclose_f(b3_id, h5_err)
  endif

  call H5Gclose_f(b2_id, h5_err)
enddo

call H5Gclose_f(b_id, h5_err)

call pmd_kill_compound_complex(complex_t)
call H5Fclose_f(f_id, h5_err)

err_flag = .false.

end subroutine hdf5_write_grid_field

