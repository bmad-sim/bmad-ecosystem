!+
! Module wavefront_hdf5_mod
!
! Genesis 1.3 Version 4 field file input and output for wavefront_struct.
!
! The format is the one written by Wavefront.write_genesis4 and read by
! Wavefront.from_genesis4 in openPMD-beamphysics, which is in turn the format written by
! Genesis's WriteFieldHDF5 (src/IO/writeFieldHDF5.cpp). At the file root:
!
!   gridpoints    -- integer, one element: Number of transverse grid points, nx = ny.
!   gridsize      -- double, one element:  Transverse grid spacing [m], dx = dy.
!   refposition   -- double, one element:  Pulse reference position along the lattice [m].
!   wavelength    -- double, one element:  Central wavelength [m].
!   slicecount    -- integer, one element: Number of slices, nz.
!   slicespacing  -- double, one element:  Slice spacing [m].
!
! and then one group per slice, named slice000001 through slice<nz>, each holding two
! datasets field-real and field-imag of nx*ny doubles.
!
! Two conversions matter.
!
! Units. Genesis stores the field as a dimensionless amplitude dfl in sqrt(W), related to
! the electric field in V/m by E = dfl * sqrt(2 * Z0) / dx, where Z0 = mu_0 * c is the
! impedance of free space. Note this scales with dx, so it is not simply a change of units.
!
! Index order. The per-slice datasets are flat arrays of nx*ny values with the x index
! running fastest, that is, (x1,y1), (x2,y1), ... (xn,y1), (x1,y2), ... The Python
! implementation has to transpose to write that, since NumPy is row major. Fortran does
! not: a column major (nx, ny) array flattens with the first index fastest already, so
! reshape is all that is needed and it is a copy, not a permutation.
!
! Genesis's format holds one polarisation component only, so a write must select either Ex
! or Ey, and a read fills Ex.
!-

module wavefront_hdf5_mod

use wavefront_mod
use hdf5_interface

implicit none

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine wavefront_write_genesis4 (wf, file_name, err_flag, polarization)
!
! Routine to write a wavefront to a Genesis 1.3 Version 4 field file.
!
! Genesis requires a square transverse grid with equal spacings, so nx must equal ny and
! dx must equal dy. Genesis also carries a single polarisation component, so if both Ex and
! Ey are present the polarization argument must say which one to write.
!
! Input:
!   wf            -- wavefront_struct: Wavefront to write.
!   file_name     -- character(*): File to create. Conventionally ends in ".fld.h5".
!   polarization  -- character(*), optional: 'x' or 'y'. Default is whichever single
!                      component is present; required if both are present.
!
! Output:
!   err_flag      -- logical: Set True on error, False otherwise.
!-

subroutine wavefront_write_genesis4 (wf, file_name, err_flag, polarization)

type (wavefront_struct), target :: wf
complex(wf_rp), pointer :: field(:,:,:)
integer(hid_t) f_id, g_id
integer n_grid(3), nx, ny, nz, iz, h5_err
real(rp) dfl_scale
real(rp), allocatable :: work(:)
logical err_flag, err
character(*) file_name
character(*), optional :: polarization
character(*), parameter :: r_name = 'wavefront_write_genesis4'
character(20) group_name

!

err_flag = .true.

call wavefront_check (wf, err);  if (err) return

! Select the component to write.

if (present(polarization)) then
  select case (polarization)
  case ('x')
    if (.not. allocated(wf%Ex)) then
      call out_io (s_error$, r_name, 'WAVEFRONT HAS NO Ex COMPONENT TO WRITE.')
      return
    endif
    field => wf%Ex

  case ('y')
    if (.not. allocated(wf%Ey)) then
      call out_io (s_error$, r_name, 'WAVEFRONT HAS NO Ey COMPONENT TO WRITE.')
      return
    endif
    field => wf%Ey

  case default
    call out_io (s_error$, r_name, 'POLARIZATION MUST BE "x" OR "y". GOT: ' // polarization)
    return
  end select

elseif (allocated(wf%Ex) .and. allocated(wf%Ey)) then
  call out_io (s_error$, r_name, 'WAVEFRONT HAS BOTH Ex AND Ey BUT THE GENESIS4 FORMAT HOLDS ONE', &
                                 'COMPONENT ONLY. SPECIFY polarization = "x" OR "y".')
  return

elseif (allocated(wf%Ex)) then
  field => wf%Ex

else
  field => wf%Ey
endif

! Genesis restrictions on the grid.

n_grid = shape(field)
nx = n_grid(1); ny = n_grid(2); nz = n_grid(3)

if (nx /= ny) then
  call out_io (s_error$, r_name, 'THE GENESIS4 FORMAT REQUIRES nx = ny. GOT: \2i6\ ', i_array = [nx, ny])
  return
endif

if (wf%dx /= wf%dy) then
  call out_io (s_error$, r_name, 'THE GENESIS4 FORMAT REQUIRES dx = dy. GOT: \2es14.6\ ', r_array = [wf%dx, wf%dy])
  return
endif

!

call hdf5_open_file (file_name, 'WRITE', f_id, err);  if (err) return

call hdf5_write_dataset_int  (f_id, 'gridpoints',   [nx],              err);  if (err) return
call hdf5_write_dataset_real (f_id, 'gridsize',     [wf%dx],           err);  if (err) return
call hdf5_write_dataset_real (f_id, 'refposition',  [wf%ref_position], err);  if (err) return
call hdf5_write_dataset_real (f_id, 'wavelength',   [wf%wavelength],   err);  if (err) return
call hdf5_write_dataset_int  (f_id, 'slicecount',   [nz],              err);  if (err) return
call hdf5_write_dataset_real (f_id, 'slicespacing', [wf%dz],           err);  if (err) return

! E [V/m] -> dfl [sqrt(W)].

dfl_scale = wf%dx / sqrt(2 * (mu_0_vac * c_light))
allocate (work(nx*ny))

do iz = 1, nz
  write (group_name, '(a, i0.6)') 'slice', iz
  call H5Gcreate_f (f_id, trim(group_name), g_id, h5_err)
  if (h5_err < 0) then
    call out_io (s_error$, r_name, 'CANNOT CREATE GROUP: ' // trim(group_name))
    return
  endif

  work = dfl_scale * reshape(real(field(:,:,iz), rp), [nx*ny])
  call hdf5_write_dataset_real (g_id, 'field-real', work, err);  if (err) return

  work = dfl_scale * reshape(aimag(field(:,:,iz)), [nx*ny])
  call hdf5_write_dataset_real (g_id, 'field-imag', work, err);  if (err) return

  call H5Gclose_f (g_id, h5_err)
enddo

call H5Fclose_f (f_id, h5_err)

err_flag = .false.

end subroutine wavefront_write_genesis4

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine wavefront_read_genesis4 (wf, file_name, err_flag)
!
! Routine to read a wavefront from a Genesis 1.3 Version 4 field file.
!
! The field is read into Ex; Ey is left unallocated. This matches
! Wavefront.from_genesis4, and reflects that the format carries no polarisation state.
!
! Input:
!   file_name   -- character(*): File to read.
!
! Output:
!   wf          -- wavefront_struct: Wavefront read from the file.
!   err_flag    -- logical: Set True on error, False otherwise.
!-

subroutine wavefront_read_genesis4 (wf, file_name, err_flag)

type (wavefront_struct) wf
integer(hid_t) f_id, g_id
integer nx, nz, iz, ivec(1), h5_err
real(rp) rvec(1), dx, dz, wavelength, ref_position, e_scale
real(rp), allocatable :: re_work(:), im_work(:)
logical err_flag, err
character(*) file_name
character(*), parameter :: r_name = 'wavefront_read_genesis4'
character(20) group_name

!

err_flag = .true.

call hdf5_open_file (file_name, 'READ', f_id, err);  if (err) return

call hdf5_read_dataset_int (f_id, 'gridpoints', ivec, err, 'gridpoints');  if (err) return
nx = ivec(1)

call hdf5_read_dataset_int (f_id, 'slicecount', ivec, err, 'slicecount');  if (err) return
nz = ivec(1)

call hdf5_read_dataset_real (f_id, 'gridsize', rvec, err, 'gridsize');  if (err) return
dx = rvec(1)

call hdf5_read_dataset_real (f_id, 'slicespacing', rvec, err, 'slicespacing');  if (err) return
dz = rvec(1)

call hdf5_read_dataset_real (f_id, 'wavelength', rvec, err, 'wavelength');  if (err) return
wavelength = rvec(1)

call hdf5_read_dataset_real (f_id, 'refposition', rvec, err, 'refposition');  if (err) return
ref_position = rvec(1)

if (nx < 1 .or. nz < 1) then
  call out_io (s_error$, r_name, 'FILE HAS A NON-POSITIVE GRID SIZE. gridpoints, slicecount: \2i6\ ', &
                                                                        i_array = [nx, nz])
  return
endif

call wavefront_init (wf, nx, nx, nz, dx, dx, dz, wavelength, 'x', ref_position)

! dfl [sqrt(W)] -> E [V/m].

e_scale = sqrt(2 * (mu_0_vac * c_light)) / dx
allocate (re_work(nx*nx), im_work(nx*nx))

do iz = 1, nz
  write (group_name, '(a, i0.6)') 'slice', iz
  g_id = hdf5_open_group (f_id, trim(group_name), err, .true.);  if (err) return

  ! The dataset extent has to be checked before reading rather than trusted, because
  ! hdf5_read_dataset_real goes through H5LTread_dataset, which takes no buffer size and so
  ! writes however much the file says the dataset holds. A file whose gridpoints disagrees
  ! with the actual slice datasets would otherwise run off the end of re_work, and on a
  ! forgiving allocator that reads as working.

  if (.not. slice_size_ok(g_id, 'field-real', nx*nx, group_name)) return
  if (.not. slice_size_ok(g_id, 'field-imag', nx*nx, group_name)) return

  call hdf5_read_dataset_real (g_id, 'field-real', re_work, err, trim(group_name) // '/field-real')
  if (err) return
  call hdf5_read_dataset_real (g_id, 'field-imag', im_work, err, trim(group_name) // '/field-imag')
  if (err) return

  wf%Ex(:,:,iz) = e_scale * reshape(cmplx(re_work, im_work, wf_rp), [nx, nx])

  call H5Gclose_f (g_id, h5_err)
enddo

call H5Fclose_f (f_id, h5_err)

call wavefront_check (wf, err);  if (err) return

err_flag = .false.

end subroutine wavefront_read_genesis4

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function slice_size_ok (g_id, dataset_name, n_expect, group_name) result (is_ok)
!
! Routine to check that a slice dataset holds exactly the expected number of values, so that
! reading it cannot overrun the buffer. See the call site for why this is not optional.
!
! Input:
!   g_id          -- integer(hid_t): ID of the slice group.
!   dataset_name  -- character(*): 'field-real' or 'field-imag'.
!   n_expect      -- integer: Number of values the dataset must hold.
!   group_name    -- character(*): Group name, for the error message.
!
! Output:
!   is_ok         -- logical: True if the dataset is present and the right size.
!-

function slice_size_ok (g_id, dataset_name, n_expect, group_name) result (is_ok)

integer(hid_t) g_id
integer n_expect
integer(hsize_t) n_have
logical is_ok, err
character(*) dataset_name, group_name
character(*), parameter :: r_name = 'wavefront_read_genesis4'
type (hdf5_info_struct) info

!

is_ok = .false.

info = hdf5_object_info (g_id, dataset_name, err, .true.);  if (err) return

! data_dim is always three long, with the unused trailing dimensions zero, so a rank 1
! dataset of n values reports [n, 0, 0]. Take the product over the dimensions in use.

n_have = product(info%data_dim, mask = (info%data_dim > 0))

if (n_have /= n_expect) then
  call out_io (s_error$, r_name, &
        'DATASET ' // trim(group_name) // '/' // trim(dataset_name) // ' HOLDS \i0\ VALUES', &
        'BUT gridpoints IMPLIES \i0\. THE FILE IS INCONSISTENT.', &
        i_array = [int(n_have), n_expect])
  return
endif

is_ok = .true.

end function slice_size_ok

end module wavefront_hdf5_mod
