#!/bin/bash
#
# Check of the quick_plot Cairo backend (ACC_PLOT_PACKAGE="cairo") with Tao:
#   1. With no display: Tao must carry on (the plot window draws off-screen) and write PDF, PS
#      and PNG files.
#   2. With a virtual X display (Xvfb): the plot must be drawn into the GTK window.
# Run after Bmad has been built with ACC_PLOT_PACKAGE="cairo".

set -e

echo -e "## Cairo plotting check\n"
echo '```'

source ./util/dist_source_me > /dev/null

OUT=$(mktemp -d)

# No library or program may link X11, PGPLOT or PLplot directly. (GTK itself may use X11, so
# only each file's own NEEDED entries are checked, not everything ldd would list.)
echo "* Linked libraries"
bad=0
for f in production/lib/*.so production/bin/*; do
  [ -f "$f" ] || continue
  needed=$(readelf -d "$f" 2>/dev/null | grep NEEDED | grep -iE 'libX11|pgplot|plplot|csirocsa|qsastime' || true)
  if [ -n "$needed" ]; then echo "  FAIL  $f links: $needed"; bad=1; fi
done
[ $bad == 0 ] && echo "  OK    no X11, PGPLOT or PLplot libraries are linked" || exit 1

cd bmad-doc/tao_examples/cesr

# A plot with negative tick labels, written in every hardcopy format.
cmds="place r22 orbit; write pdf $OUT/tao.pdf; write ps $OUT/tao.ps; write gif $OUT/tao.png; quit -force"

echo "* Tao with no display"
env -u DISPLAY -u WAYLAND_DISPLAY tao -command "$cmds" < /dev/null > "$OUT/no_display.log" 2>&1 || {
  cat "$OUT/no_display.log"; echo "Tao failed"; exit 1;
}
grep -q "No display available" "$OUT/no_display.log" || {
  cat "$OUT/no_display.log"; echo "Expected the no-display message"; exit 1;
}

echo "* Tao with a virtual X display"
mkdir "$OUT/x"
QP_GTK_SNAPSHOT="$OUT/window.png" xvfb-run -a -s "-screen 0 1280x1024x24" \
  tao -command "place r22 orbit; write pdf $OUT/x/tao.pdf; quit -force" < /dev/null > "$OUT/xvfb.log" 2>&1 || {
  cat "$OUT/xvfb.log"; echo "Tao failed"; exit 1;
}

python3 - "$OUT" <<'EOF'
import re, struct, subprocess, sys, zlib

out = sys.argv[1]
ok = True

def check(cond, what):
    global ok
    print(('  OK    ' if cond else '  FAIL  ') + what)
    ok = ok and cond

def png_info(file):
    """Return (width, height, number of non-white pixels) of an 8-bit RGB(A) PNG."""
    data = open(file, 'rb').read()
    assert data[:8] == b'\x89PNG\r\n\x1a\n', 'not a PNG file'
    pos, idat = 8, b''
    while pos < len(data):
        n, typ = struct.unpack('>I4s', data[pos:pos+8])
        chunk = data[pos+8:pos+8+n]
        if typ == b'IHDR': w, h, depth, ctype = struct.unpack('>IIBB', chunk[:10])
        if typ == b'IDAT': idat += chunk
        pos += 12 + n
    assert depth == 8 and ctype in (2, 6), 'unexpected PNG format'
    bpp = 3 if ctype == 2 else 4
    raw, stride = zlib.decompress(idat), w * bpp
    prev, rows, i = bytearray(stride), [], 0
    for _ in range(h):
        f, line = raw[i], bytearray(raw[i+1:i+1+stride])
        i += 1 + stride
        for x in range(stride):
            a = line[x-bpp] if x >= bpp else 0
            b, c = prev[x], prev[x-bpp] if x >= bpp else 0
            if f == 1: line[x] = (line[x] + a) & 255
            elif f == 2: line[x] = (line[x] + b) & 255
            elif f == 3: line[x] = (line[x] + (a + b) // 2) & 255
            elif f == 4:
                p = a + b - c
                pa, pb, pc = abs(p - a), abs(p - b), abs(p - c)
                line[x] = (line[x] + (a if pa <= pb and pa <= pc else b if pb <= pc else c)) & 255
        rows.append(line)
        prev = line
    ink = sum(1 for r in rows for x in range(0, stride, bpp) if min(r[x:x+3]) < 200 and (bpp == 3 or r[x+3] > 128))
    return w, h, ink

def pdf_text(file):
    return subprocess.run(['pdftotext', file, '-'], capture_output=True, text=True).stdout

print('* No display')
check(open(out + '/tao.pdf', 'rb').read(5) == b'%PDF-', 'tao.pdf is a PDF file')
check(open(out + '/tao.ps', 'rb').read(4) == b'%!PS', 'tao.ps is a PostScript file')
text = pdf_text(out + '/tao.pdf')
check('Horizontal Orbit' in text, 'PDF has the graph title')
check(re.search('^−[0-9.]+$', text, re.M) is not None, 'PDF tick labels use the Unicode minus sign')
check(re.search('^-[0-9.]+$', text, re.M) is None, 'PDF tick labels have no hyphen minus')
w, h, ink = png_info(out + '/tao.png')
check(w > 100 and h > 100 and ink > 1000, 'tao.png has a plot (%dx%d, %d dark pixels)' % (w, h, ink))

print('* Virtual X display')
log = open(out + '/xvfb.log').read()
check('No display available' not in log, 'the plot window opened')
w, h, ink = png_info(out + '/window.png')
check(w > 100 and h > 100 and ink > 1000, 'the window shows a plot (%dx%d, %d dark pixels)' % (w, h, ink))
check('Horizontal Orbit' in pdf_text(out + '/x/tao.pdf'), 'PDF written with a display has the graph title')

sys.exit(0 if ok else 1)
EOF

echo '```'
echo -e "\n## Cairo plotting check passed\n"
