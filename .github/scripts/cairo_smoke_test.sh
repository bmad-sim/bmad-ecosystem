#!/bin/bash
# Smoke test for the quick_plot Cairo backend (ACC_PLOT_PACKAGE="cairo").
#
# Runs Tao on a lattice with an on-screen plot page and writes PDF and PNG hardcopy.
# The plot window program (qp_cairo_viewer) is run with SDL's "dummy" video driver so no
# display is needed, and it saves what its window shows (QP_CAIRO_VIEWER_SNAPSHOT).
# Checks that all three outputs are produced and look like the right file types.

echo -e "## Cairo plotting smoke test\n"

source ./util/dist_source_me > /dev/null

set -e

out="${RUNNER_TEMP:-/tmp}/cairo_smoke_test"
mkdir -p "$out"

viewer="$ACC_EXE/qp_cairo_viewer"
if [ ! -x "$viewer" ]; then
  echo "FAIL: plot window program not built: $viewer"
  exit 1
fi

export SDL_VIDEODRIVER=dummy
export QP_CAIRO_VIEWER_SNAPSHOT="$out/window.png"
export QP_CAIRO_DEBUG=1

cd bmad-doc/tao_examples/cesr
echo '```'
# Give Tao time to start, draw, and write the files, then quit.
(sleep 15; echo quit) | tao -command "write pdf $out/tao.pdf; write gif $out/tao.png"
echo '```'

fail=0
check () {   # check <file> <text expected in the first 8 bytes>
  local file=$1 magic=$2
  if [ ! -s "$file" ]; then
    echo "FAIL: missing or empty: $file"
    fail=1
  elif ! head -c 8 "$file" | grep -q "$magic"; then
    echo "FAIL: $file does not start with $magic"
    fail=1
  else
    echo "OK: $file ($(wc -c < "$file") bytes)"
  fi
}

check "$out/tao.pdf" "%PDF"
check "$out/tao.png" "PNG"
check "$out/window.png" "PNG"

if pgrep -f qp_cairo_viewer > /dev/null; then
  echo "FAIL: qp_cairo_viewer still running after Tao exited"
  fail=1
fi

exit $fail
