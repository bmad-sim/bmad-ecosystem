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

# Full CSS4 named colors (excluding the 10 that share names with original quick_plot colors:
# black, blue, cyan, green, magenta, orange, purple, red, white, yellow
# All excluded names now have identical RGB values to their CSS4 definitions.)
CSS4_COLORS = {
    "aliceblue": (240, 248, 255),
    "antiquewhite": (250, 235, 215),
    "aqua": (0, 255, 255),
    "aquamarine": (127, 255, 212),
    "azure": (240, 255, 255),
    "beige": (245, 245, 220),
    "bisque": (255, 228, 196),
    "blanchedalmond": (255, 235, 205),
    "blueviolet": (138, 43, 226),
    "brown": (165, 42, 42),
    "burlywood": (222, 184, 135),
    "cadetblue": (95, 158, 160),
    "chartreuse": (127, 255, 0),
    "chocolate": (210, 105, 30),
    "coral": (255, 127, 80),
    "cornflowerblue": (100, 149, 237),
    "cornsilk": (255, 248, 220),
    "crimson": (220, 20, 60),
    "darkblue": (0, 0, 139),
    "darkcyan": (0, 139, 139),
    "darkgoldenrod": (184, 134, 11),
    "darkgray": (169, 169, 169),
    "darkgreen": (0, 100, 0),
    "darkgrey": (169, 169, 169),
    "darkkhaki": (189, 183, 107),
    "darkmagenta": (139, 0, 139),
    "darkolivegreen": (85, 107, 47),
    "darkorange": (255, 140, 0),
    "darkorchid": (153, 50, 204),
    "darkred": (139, 0, 0),
    "darksalmon": (233, 150, 122),
    "darkseagreen": (143, 188, 143),
    "darkslateblue": (72, 61, 139),
    "darkslategray": (47, 79, 79),
    "darkslategrey": (47, 79, 79),
    "darkturquoise": (0, 206, 209),
    "darkviolet": (148, 0, 211),
    "deeppink": (255, 20, 147),
    "deepskyblue": (0, 191, 255),
    "dimgray": (105, 105, 105),
    "dimgrey": (105, 105, 105),
    "dodgerblue": (30, 144, 255),
    "firebrick": (178, 34, 34),
    "floralwhite": (255, 250, 240),
    "forestgreen": (34, 139, 34),
    "fuchsia": (255, 0, 255),
    "gainsboro": (220, 220, 220),
    "ghostwhite": (248, 248, 255),
    "gold": (255, 215, 0),
    "goldenrod": (218, 165, 32),
    "gray": (128, 128, 128),
    "greenyellow": (173, 255, 47),
    "grey": (128, 128, 128),
    "honeydew": (240, 255, 240),
    "hotpink": (255, 105, 180),
    "indianred": (205, 92, 92),
    "indigo": (75, 0, 130),
    "ivory": (255, 255, 240),
    "khaki": (240, 230, 140),
    "lavender": (230, 230, 250),
    "lavenderblush": (255, 240, 245),
    "lawngreen": (124, 252, 0),
    "lemonchiffon": (255, 250, 205),
    "lightblue": (173, 216, 230),
    "lightcoral": (240, 128, 128),
    "lightcyan": (224, 255, 255),
    "lightgoldenrodyellow": (250, 250, 210),
    "lightgray": (211, 211, 211),
    "lightgreen": (144, 238, 144),
    "lightgrey": (211, 211, 211),
    "lightpink": (255, 182, 193),
    "lightsalmon": (255, 160, 122),
    "lightseagreen": (32, 178, 170),
    "lightskyblue": (135, 206, 250),
    "lightslategray": (119, 136, 153),
    "lightslategrey": (119, 136, 153),
    "lightsteelblue": (176, 196, 222),
    "lightyellow": (255, 255, 224),
    "lime": (0, 255, 0),
    "limegreen": (50, 205, 50),
    "linen": (250, 240, 230),
    "maroon": (128, 0, 0),
    "mediumaquamarine": (102, 205, 170),
    "mediumblue": (0, 0, 205),
    "mediumorchid": (186, 85, 211),
    "mediumpurple": (147, 112, 219),
    "mediumseagreen": (60, 179, 113),
    "mediumslateblue": (123, 104, 238),
    "mediumspringgreen": (0, 250, 154),
    "mediumturquoise": (72, 209, 204),
    "mediumvioletred": (199, 21, 133),
    "midnightblue": (25, 25, 112),
    "mintcream": (245, 255, 250),
    "mistyrose": (255, 228, 225),
    "moccasin": (255, 228, 181),
    "navajowhite": (255, 222, 173),
    "navy": (0, 0, 128),
    "oldlace": (253, 245, 230),
    "olive": (128, 128, 0),
    "olivedrab": (107, 142, 35),
    "orangered": (255, 69, 0),
    "orchid": (218, 112, 214),
    "palegoldenrod": (238, 232, 170),
    "palegreen": (152, 251, 152),
    "paleturquoise": (175, 238, 238),
    "palevioletred": (219, 112, 147),
    "papayawhip": (255, 239, 213),
    "peachpuff": (255, 218, 185),
    "peru": (205, 133, 63),
    "pink": (255, 192, 203),
    "plum": (221, 160, 221),
    "powderblue": (176, 224, 230),
    "rebeccapurple": (102, 51, 153),
    "rosybrown": (188, 143, 143),
    "royalblue": (65, 105, 225),
    "saddlebrown": (139, 69, 19),
    "salmon": (250, 128, 114),
    "sandybrown": (244, 164, 96),
    "seagreen": (46, 139, 87),
    "seashell": (255, 245, 238),
    "sienna": (160, 82, 45),
    "silver": (192, 192, 192),
    "skyblue": (135, 206, 235),
    "slateblue": (106, 90, 205),
    "slategray": (112, 128, 144),
    "slategrey": (112, 128, 144),
    "snow": (255, 250, 250),
    "springgreen": (0, 255, 127),
    "steelblue": (70, 130, 180),
    "tan": (210, 180, 140),
    "teal": (0, 128, 128),
    "thistle": (216, 191, 216),
    "tomato": (255, 99, 71),
    "turquoise": (64, 224, 208),
    "violet": (238, 130, 238),
    "wheat": (245, 222, 179),
    "whitesmoke": (245, 245, 245),
    "yellowgreen": (154, 205, 50),
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
    lines.append(f"!   17-{max_index}:  CSS4 named colors ({n_css4} colors, excluding 10 that duplicate existing names)")
    lines.append(f"!   {continuous_start}+:    Continuous color mapping")
    lines.append("!")
    lines.append("! Colors are ordered by visual distinctiveness (farthest-point-first algorithm)")
    lines.append("! so that devices with limited color indices get the most useful colors first.")
    lines.append("!")
    lines.append("! Excluded CSS4 names (conflict with existing quick_plot names):")
    lines.append("!   black, blue, cyan, green, magenta, orange, purple, red, white, yellow")
    lines.append("")
    lines.append(f"integer, parameter :: qp_n_css4_colors = {n_css4}")
    lines.append(f"integer, parameter :: qp_max_color_index$ = {max_index}")
    lines.append(f"integer, parameter :: qp_continuous_color_start$ = {continuous_start}")
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
