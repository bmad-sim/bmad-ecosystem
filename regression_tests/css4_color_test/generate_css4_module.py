#!/usr/bin/env python3
"""Generate qp_css4_colors_mod.f90 with colors ordered by distinctiveness.

Uses "farthest point first" algorithm to order CSS4 colors so that the most
visually distinct colors get lower indices. This ensures that devices with
limited color support (e.g., PGPLOT with 99 indices) get the most useful
colors natively, while near-whites and subtle variants end up at higher indices.

The first 17 colors (indices 0-16) are the original quick_plot colors and are
unchanged. CSS4 colors start at index 17.
"""

import numpy as np
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
            # Find minimum distance to any already-selected color
            min_dist = min(color_distance_sq((r, g, b), s) for s in selected_rgb)
            if min_dist > best_min_dist:
                best_min_dist = min_dist
                best_idx = i

        chosen = remaining.pop(best_idx)
        ordered.append(chosen)
        selected_rgb.append((chosen[1], chosen[2], chosen[3]))

    return ordered


def generate_fortran_module(ordered_css4):
    """Generate the Fortran module source code.

    Parameters
    ----------
    ordered_css4 : list of tuple
        Ordered list of (name, r, g, b) tuples for indices 17-154.

    Returns
    -------
    str
        Complete Fortran module source.
    """
    n_css4 = len(ordered_css4)
    max_index = 16 + n_css4
    continuous_start = max_index + 1

    lines = []
    lines.append("module qp_css4_colors_mod")
    lines.append("")
    lines.append("implicit none")
    lines.append("")
    lines.append("! CSS4 (CSS Color Level 4) named colors for quick_plot.")
    lines.append("!")
    lines.append("! Index scheme:")
    lines.append("!   0-15:    Original quick_plot named colors")
    lines.append("!   16:      Transparent (special)")
    lines.append(f"!   17-{max_index}:  CSS4 named colors ({n_css4} colors)")
    lines.append(f"!   {continuous_start}+:    Continuous color mapping")
    lines.append("!")
    lines.append("! Colors are ordered by visual distinctiveness (farthest-point-first algorithm)")
    lines.append("! so that devices with limited color indices get the most useful colors first.")
    lines.append("! Colors that duplicate original quick_plot names (indices 0-15) sort to the")
    lines.append("! end since their RGB distance to the seeds is zero.")
    lines.append("")
    lines.append(f"integer, parameter :: qp_n_css4_colors = {n_css4}")
    lines.append(f"integer, parameter :: qp_max_color_index$ = {max_index}")
    lines.append(f"integer, parameter :: qp_continuous_color_start$ = {continuous_start}")
    lines.append("integer, parameter :: qp_custom_color$ = -2   ! Sentinel: use qp_custom_rgb instead of palette lookup")
    lines.append("")
    lines.append("! Storage for arbitrary RGB color specified via hex (#RRGGBB) or RGB(r,g,b) string.")
    lines.append("! Set by qp_string_to_enum when it returns qp_custom_color$.")
    lines.append("integer, save :: qp_custom_rgb(3) = [0, 0, 0]  ! R, G, B values (0-255)")
    lines.append("")

    # Build full RGB arrays (indices 0 to max_index)
    all_colors = []  # (name, r, g, b) for index 0..max_index
    for i in range(17):
        orig = ORIGINAL_COLORS[i]
        all_colors.append((orig[0], orig[1], orig[2], orig[3]))
    for name, r, g, b in ordered_css4:
        all_colors.append((name, r, g, b))

    # Generate RED array
    lines.append("! RGB values for all discrete colors (indices 0:{})".format(max_index))
    lines.append("! Ordered: original 16 + transparent, then CSS4 by distinctiveness.")
    lines.append(f"integer, parameter :: qp_color_red(0:{max_index}) = [ &")
    for i in range(0, len(all_colors), 10):
        chunk = all_colors[i:i+10]
        vals = ", ".join(f"{c[1]:3d}" for c in chunk)
        if i + 10 >= len(all_colors):
            lines.append(f"  {vals} ]")
        else:
            comment = f"  ! {i}-{i+len(chunk)-1}"
            lines.append(f"  {vals}, &{comment}")
    lines.append("")

    # Generate GREEN array
    lines.append(f"integer, parameter :: qp_color_green(0:{max_index}) = [ &")
    for i in range(0, len(all_colors), 10):
        chunk = all_colors[i:i+10]
        vals = ", ".join(f"{c[2]:3d}" for c in chunk)
        if i + 10 >= len(all_colors):
            lines.append(f"  {vals} ]")
        else:
            comment = f"  ! {i}-{i+len(chunk)-1}"
            lines.append(f"  {vals}, &{comment}")
    lines.append("")

    # Generate BLUE array
    lines.append(f"integer, parameter :: qp_color_blue(0:{max_index}) = [ &")
    for i in range(0, len(all_colors), 10):
        chunk = all_colors[i:i+10]
        vals = ", ".join(f"{c[3]:3d}" for c in chunk)
        if i + 10 >= len(all_colors):
            lines.append(f"  {vals} ]")
        else:
            comment = f"  ! {i}-{i+len(chunk)-1}"
            lines.append(f"  {vals}, &{comment}")
    lines.append("")

    # Generate name array
    lines.append(f"character(24), parameter :: qp_css4_color_name(17:{max_index}) = [ character(24) :: &")
    for i, (name, r, g, b) in enumerate(ordered_css4):
        idx = 17 + i
        if i == len(ordered_css4) - 1:
            lines.append(f"  '{name}' ]")
        elif (i + 1) % 5 == 0:
            lines.append(f"  '{name}', &  ! index {idx}")
        else:
            lines.append(f"  '{name}', &")
    lines.append("")

    # Generate integer parameter constants for each CSS4 color
    lines.append("! Integer constants for CSS4 color indices.")
    lines.append('! Note: The color "tan" uses the constant css_tan$ (not tan$) to avoid conflict')
    lines.append("! with the tangent operator (tan$) defined in bmad_struct.f90.")
    lines.append("")

    # Group constants 3 per line for readability
    const_entries = []
    for i, (name, r, g, b) in enumerate(ordered_css4):
        idx = 17 + i
        const_name = "css_tan" if name == "tan" else name
        const_entries.append(f"{const_name}$ = {idx}")

    for i in range(0, len(const_entries), 3):
        chunk = const_entries[i:i+3]
        lines.append(f"integer, parameter :: {', '.join(chunk)}")

    lines.append("")
    lines.append("end module qp_css4_colors_mod")
    lines.append("")

    return "\n".join(lines)


def main():
    """Generate the reordered CSS4 color module."""
    # Seed colors: the original 17 quick_plot colors
    seed_rgb = [(c[1], c[2], c[3]) for c in ORIGINAL_COLORS.values()]

    # Prepare CSS4 candidates as (name, r, g, b)
    candidates = [(name, r, g, b) for name, (r, g, b) in CSS4_COLORS.items()]

    print(f"Ordering {len(candidates)} CSS4 colors by distinctiveness...")
    ordered = farthest_point_first(seed_rgb, candidates)

    print("\nFirst 20 colors (most distinctive, get lowest indices):")
    for i, (name, r, g, b) in enumerate(ordered[:20]):
        print(f"  {17+i:3d}: {name:24s} ({r:3d}, {g:3d}, {b:3d})")

    print("\nLast 20 colors (least distinctive, highest indices):")
    for i, (name, r, g, b) in enumerate(ordered[-20:]):
        idx = 17 + len(ordered) - 20 + i
        print(f"  {idx:3d}: {name:24s} ({r:3d}, {g:3d}, {b:3d})")

    # Generate and write the Fortran module
    source = generate_fortran_module(ordered)
    outpath = "../../sim_utils/plot/qp_css4_colors_mod.f90"
    with open(outpath, "w") as f:
        f.write(source)
    print(f"\nWrote: {outpath}")
    print(f"Total colors: {len(ordered)} CSS4 + 17 original = {len(ordered)+17}")
    print(f"Max color index: {16 + len(ordered)}")
    print(f"Continuous color start: {17 + len(ordered)}")

    # Print the "boundary" color (index 99) for device with 99-color limit
    if len(ordered) > 82:
        print(f"\nWith a 99-color device (indices 0-99):")
        print(f"  Native colors: 0-99 ({100} total, including {83} CSS4)")
        print(f"  Color at index 99: {ordered[82][0]}")
        print(f"  Colors needing fallback: {len(ordered) - 83} (indices 100-{16+len(ordered)})")


if __name__ == "__main__":
    main()
