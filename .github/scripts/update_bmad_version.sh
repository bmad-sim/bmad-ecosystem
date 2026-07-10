#!/bin/bash
#
# Usage: update_bmad_version.sh <version>
#
# Stamps bmad/modules/bmad_version_mod.f90 with the given version string.
# Run by the release workflow so the reported version matches the release tag.

set -euo pipefail

VERSION="${1:?Usage: update_bmad_version.sh <version>}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
VERSION_FILE="$REPO_ROOT/bmad/modules/bmad_version_mod.f90"

cat >"$VERSION_FILE" <<EOF
!+
! Module bmad_version_mod
!
! This file is generated automatically via a GitHub action. Do not modify by hand.
!-

module bmad_version_mod
character(*), parameter :: bmad_version_date = "$VERSION"
end module
EOF

echo "Stamped $VERSION_FILE with version $VERSION"
