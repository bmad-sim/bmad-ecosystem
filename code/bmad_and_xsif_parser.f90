!+
! Subroutine bmad_and_xsif_parser (lat_file, lat, make_mats6, digested_read_ok, use_line)
!
! Subroutine to parse either a Bmad or XSIF (extended standard input format) lattice file.
!
! This routine will assume that the lattice file uses Bmad syntax except when 
! lat_file is prefixed by the string 'xsif::'. Example:
!   lat_file = 'xsif::/nfs/user/dcs/this_lat'
!
! Note: The presence of an LCavity element in the input file (even if it is not
! used in the lattice) will make lat%lattice_type = linear_lattice$.
!
! Modules needed:
!   use bmad
!
! Input:
!   file       -- Character(*): Name of the Bmad or xsif file.
!   make_mats6 -- Logical, optional: Compute the 6x6 transport matrices for the
!                   Elements? Default is True.
!   use_line   -- Character(*), optional: If present and not blank, override the use 
!                   statement in the lattice file and use use_line instead.
!
! Output:
!   lat         -- lat_struct: Structure holding the lattice information.
!   bmad_status  -- Bmad status common block.
!     %ok            -- Set True if parsing is successful. False otherwise.
!-

subroutine bmad_and_xsif_parser (lat_file, lat, make_mats6, digested_read_ok, use_line)

use bmad_struct
use bmad_interface, except_dummy => bmad_and_xsif_parser

implicit none

type (lat_struct), target :: lat

character(*) :: lat_file
character(*), optional :: use_line

logical, optional :: make_mats6, digested_read_ok

integer ix

! look for 'xsif::' prefix. Allow lat_file to have leading blanks.

ix = index(lat_file, 'xsif::')
if (ix /= 0) then
  if (lat_file(1:ix-1) /= '') ix = 0
endif

if (ix == 0) then
  call bmad_parser (lat_file, lat, make_mats6, digested_read_ok, use_line)
else
  call xsif_parser (lat_file(ix+6:), lat, make_mats6, digested_read_ok, use_line)
endif

end subroutine
