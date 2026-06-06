"""
Verify quick_plot CSS4 color module data against the reference source.

Contains the canonical color definitions and ordering algorithm used to
generate qp_css4_colors_mod.f90. Also verifies that the compiled Fortran
module's RGB arrays match expected values.

Usage:
    python3 css4_color_reference.py [--swatch output.png]
"""

import os
import subprocess
import sys

from PIL import Image
from PIL.ImageColor import colormap as _pil_colormap

# Original quick_plot colors (indices 0-16)
# Indices 3, 8, 12 updated to match CSS4 spec values.
# Old values: green=(0,255,0) [now "lime"], orange=(255,127,0), purple=(127,0,255)
ORIGINAL_COLORS = {
    0: ("white", 255, 255, 255),
    1: ("black", 0, 0, 0),
    2: ("red", 255, 0, 0),
    3: ("green", 0, 128, 0),
    4: ("blue", 0, 0, 255),
    5: ("cyan", 0, 255, 255),
    6: ("magenta", 255, 0, 255),
    7: ("yellow", 255, 255, 0),
    8: ("orange", 255, 165, 0),
    9: ("yellow_green", 127, 255, 0),
    10: ("light_green", 0, 255, 127),
    11: ("navy_blue", 0, 127, 255),
    12: ("purple", 128, 0, 128),
    13: ("reddish_purple", 255, 0, 127),
    14: ("dark_grey", 85, 85, 85),
    15: ("light_grey", 170, 170, 170),
    16: ("transparent", 0, 0, 0),
}

# CSS4 named colors from PIL (the authoritative source).
CSS4_COLORS = {
    name: (int(hexval[1:3], 16), int(hexval[3:5], 16), int(hexval[5:7], 16))
    for name, hexval in _pil_colormap.items()
}


def color_distance_sq(c1, c2):
    """Compute squared Euclidean distance in RGB space.

    Parameters
    ----------
    c1 : tuple
        RGB tuple (r, g, b).
    c2 : tuple
        RGB tuple (r, g, b).

    Returns
    -------
    int
        Squared distance.
    """
    return (c1[0] - c2[0])**2 + (c1[1] - c2[1])**2 + (c1[2] - c2[2])**2


def farthest_point_first(seed_colors, candidate_colors):
    """Order candidate colors by maximum distinctiveness from already-selected set.

    Uses the farthest-point-first (greedy) algorithm: at each step, select the
    candidate that maximizes its minimum distance to all previously selected colors.

    Parameters
    ----------
    seed_colors : list of tuple
        Initial set of (r, g, b) tuples already "selected" (the original 17 colors).
    candidate_colors : list of tuple
        List of (name, r, g, b) tuples to order.

    Returns
    -------
    list of tuple
        Ordered list of (name, r, g, b) tuples.
    """
    selected_rgb = list(seed_colors)
    remaining = list(candidate_colors)
    ordered = []

    while remaining:
        best_idx = -1
        best_min_dist = -1

        for i, (name, r, g, b) in enumerate(remaining):
            min_dist = min(color_distance_sq((r, g, b), s) for s in selected_rgb)
            if min_dist > best_min_dist:
                best_min_dist = min_dist
                best_idx = i

        chosen = remaining.pop(best_idx)
        ordered.append(chosen)
        selected_rgb.append((chosen[1], chosen[2], chosen[3]))

    return ordered


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
