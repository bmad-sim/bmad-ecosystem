!+
! Subroutine bmad_and_xsif_parser (lat_file, lat, make_mats6, digested_read_ok, use_line, err_flag)
!
! Subroutine to parse either a Bmad, XSIF (eXtended Standard Input Format), or digested lattice file.
!
! This routine will assume that the lattice file uses Bmad syntax except when:
!   * lat_file is prefixed by the string 'xsif::' 
!   * lat_file has a '.xsif' suffix. Example: lat_file = 'xsif::/nfs/user/dcs/this_lat'
!   * lat_file has '.digested' in the name.
!
! Input:
!   lat_file   -- Character(*): Name of the Bmad or xsif file.
!   make_mats6 -- Logical, optional: Compute the 6x6 transport matrices for the
!                   Elements? Default is True.
!                   This argument ignored if lat_file corresponds to a digested file.
!   use_line   -- Character(*), optional: If present and not blank, override the use 
!                   statement in the lattice file and use use_line instead.
!                   This argument ignored if lat_file corresponds to a digested file.
!
! Output:
!   lat         -- lat_struct: Structure holding the lattice information.
!   digested_read_ok -- Logical, optional: Set True if the digested file was
!                        successfully read. False otherwise.
!   err_flag    -- Logical, optional: Set True if there is an error. False otherwise
!-

subroutine bmad_and_xsif_parser (lat_file, lat, make_mats6, digested_read_ok, use_line, err_flag)

use bmad_struct
use bmad_interface, except_dummy => bmad_and_xsif_parser

implicit none

type (lat_struct), target :: lat

character(*) :: lat_file
character(*), optional :: use_line

logical, optional :: make_mats6, digested_read_ok, err_flag

integer ix, inc_version

! Digested file?

if (index(lat_file, '.digested') /= 0) then
  call read_digested_bmad_file (lat_file, lat, inc_version, err_flag)
  if (present(digested_read_ok)) digested_read_ok = (inc_version == bmad_inc_version$)
  return
endif

! look for 'xsif::' prefix. Allow lat_file to have leading blanks.

ix = index(lat_file, 'xsif::')
if (ix /= 0) then
  if (lat_file(1:ix-1) == '') then
    call xsif_parser (lat_file(ix+6:), lat, make_mats6, digested_read_ok, use_line, err_flag)
    return
  endif
endif

!

ix = index(lat_file, '.xsif')
if (ix /= 0 .and. len_trim(lat_file) == ix+4) then
  call xsif_parser (lat_file, lat, make_mats6, digested_read_ok, use_line, err_flag)
else 
  call bmad_parser (lat_file, lat, make_mats6, digested_read_ok, use_line, err_flag)
endif

end subroutine
