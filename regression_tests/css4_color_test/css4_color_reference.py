"""
Verify quick_plot CSS4 color module data against the reference source.

Uses generate_css4_module.py as the single source of truth for color definitions.
Verifies that the compiled Fortran module's RGB arrays match expected values by
running a small Fortran program that dumps the palette.

Usage:
    python3 css4_color_reference.py [--swatch output.png]
"""

import os
import subprocess
import sys

from PIL import Image

# Import canonical color data from the generator (single source of truth)
sys.path.insert(0, os.path.dirname(__file__))
from generate_css4_module import ORIGINAL_COLORS, CSS4_COLORS, farthest_point_first


def get_ordered_palette():
    """Get the full color palette in index order.

    Returns
    -------
    list of tuple
        List of (name, r, g, b) for indices 0 through 154.
    """
    # Original 17 colors (indices 0-16)
    palette = [(c[0], c[1], c[2], c[3]) for c in
               [ORIGINAL_COLORS[i] for i in range(17)]]

    # CSS4 colors in farthest-point-first order (indices 17-154)
    seed_rgb = [(c[1], c[2], c[3]) for c in palette]
    candidates = [(name, r, g, b) for name, (r, g, b) in CSS4_COLORS.items()]
    ordered = farthest_point_first(seed_rgb, candidates)
    palette.extend(ordered)

    return palette


def verify_against_fortran(bin_dir):
    """Run the Fortran color dump program and compare output to reference.

    Parameters
    ----------
    bin_dir : str
        Path to directory containing compiled test executables.

    Returns
    -------
    bool
        True if all colors match.
    """
    exe = os.path.join(bin_dir, "css4_color_test")
    if not os.path.isfile(exe):
        print(f"  Fortran test executable not found: {exe}")
        print("  Skipping Fortran verification (run build first).")
        return True

    proc = subprocess.run([exe], capture_output=True, text=True)
    if proc.returncode != 0:
        print(f"  css4_color_test failed with exit code {proc.returncode}")
        return False

    palette = get_ordered_palette()
    errors = 0

    for line in proc.stdout.splitlines():
        parts = line.split()
        if len(parts) < 5:
            continue
        idx = int(parts[0])
        r, g, b = int(parts[2]), int(parts[3]), int(parts[4])
        if idx >= len(palette):
            continue
        name, er, eg, eb = palette[idx]
        if r != er or g != eg or b != eb:
            print(f"  MISMATCH index {idx} ({name}): "
                  f"got ({r},{g},{b}) expected ({er},{eg},{eb})")
            errors += 1

    if errors == 0:
        print(f"  All {len(palette)} colors verified OK")
    else:
        print(f"  {errors} mismatches found")
    return errors == 0


def generate_swatch(output_path, swatch_w=60, swatch_h=30, cols=10):
    """Generate a reference color swatch PNG using PIL.

    Parameters
    ----------
    output_path : str
        Output PNG file path.
    swatch_w : int
        Width of each color swatch in pixels.
    swatch_h : int
        Height of each color swatch in pixels.
    cols : int
        Number of columns in the grid.
    """
    palette = get_ordered_palette()
    n_colors = len(palette)
    rows = (n_colors + cols - 1) // cols

    width = cols * swatch_w
    height = rows * swatch_h

    img = Image.new("RGB", (width, height), (255, 255, 255))

    for idx, (name, r, g, b) in enumerate(palette):
        col = idx % cols
        row = idx // cols
        x0 = col * swatch_w
        y0 = row * swatch_h
        for y in range(y0, min(y0 + swatch_h, height)):
            for x in range(x0, min(x0 + swatch_w, width)):
                img.putpixel((x, y), (r, g, b))

    img.save(output_path)
    print(f"Reference swatch written to: {output_path}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Verify CSS4 color module data.")
    parser.add_argument("--swatch", help="Generate reference swatch PNG at this path.")
    parser.add_argument("--bin-dir", help="Path to compiled test binaries for Fortran verification.")
    args = parser.parse_args()

    palette = get_ordered_palette()
    print(f"Reference palette: {len(palette)} colors (17 original + {len(palette)-17} CSS4)")

    if args.swatch:
        generate_swatch(args.swatch)

    if args.bin_dir:
        print("\nVerifying Fortran module against reference...")
        ok = verify_against_fortran(args.bin_dir)
        sys.exit(0 if ok else 1)
    else:
        print("\nTo verify against compiled Fortran, pass --bin-dir <path>")
        print("Printing first 20 CSS4 colors (by distinctiveness):")
        for i, (name, r, g, b) in enumerate(palette[17:37]):
            print(f"  {17+i:3d}: {name:24s} ({r:3d}, {g:3d}, {b:3d})")
