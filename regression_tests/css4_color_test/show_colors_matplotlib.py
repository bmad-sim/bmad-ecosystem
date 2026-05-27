#!/usr/bin/env python3
"""Display what PGPLOT would render vs the expected CSS4 colors.

Shows three columns per row:
  - Expected CSS4 color (from the spec / matplotlib)
  - What PGPLOT would actually render (exact if within device limit,
    nearest-color fallback if above)
  - Color name and index

Usage
-----
    python3 show_colors_matplotlib.py                  # default device limit 99
    python3 show_colors_matplotlib.py --device-limit 127
    python3 show_colors_matplotlib.py --original       # show original 16 colors
"""

import os
import sys

import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

# Add parent path for module import
sys.path.insert(0, os.path.dirname(__file__))
from css4_color_reference import ORIGINAL_COLORS, CSS4_COLORS, get_ordered_palette

# Derive name lists from the canonical source
_palette = get_ordered_palette()
QP_ORIGINAL_NAMES = [name for (name, r, g, b) in _palette[:16]]
QP_CSS4_NAMES = [name for (name, r, g, b) in _palette[17:]]


def get_qp_module_rgb():
    """Get RGB values from the reference palette.

    Returns
    -------
    list of tuple
        RGB values for indices 0-154, each as (r, g, b) with values 0-255.
    """
    palette = get_ordered_palette()
    return [(r, g, b) for (name, r, g, b) in palette]


def nearest_color_index(ix_color, all_rgb, device_limit):
    """Find the nearest color within device range (Euclidean RGB distance).

    This replicates the PGPLOT qp_set_color_basic fallback logic.

    Parameters
    ----------
    ix_color : int
        The requested color index.
    all_rgb : list of tuple
        All RGB values indexed 0-154.
    device_limit : int
        Maximum color index supported by the device.

    Returns
    -------
    int
        Index of the nearest color within device range.
    """
    r, g, b = all_rgb[ix_color]
    best_dist = float('inf')
    best_idx = 0
    for ix in range(device_limit + 1):
        ri, gi, bi = all_rgb[ix]
        dist = (r - ri)**2 + (g - gi)**2 + (b - bi)**2
        if dist < best_dist:
            best_dist = dist
            best_idx = ix
    return best_idx


def pgplot_rendered_color(ix_color, all_rgb, device_limit):
    """Get the RGB that PGPLOT would actually render for a given color index.

    Parameters
    ----------
    ix_color : int
        The requested color index (0-154).
    all_rgb : list of tuple
        All RGB values indexed 0-154.
    device_limit : int
        Maximum color index supported by the PGPLOT device.

    Returns
    -------
    tuple
        (r, g, b) that PGPLOT would actually display (0-255).
    """
    if ix_color <= device_limit:
        return all_rgb[ix_color]
    else:
        nearest_idx = nearest_color_index(ix_color, all_rgb, device_limit)
        return all_rgb[nearest_idx]


def plot_pgplot_comparison(device_limit=99):
    """Show expected CSS4 color vs what PGPLOT actually renders.

    Parameters
    ----------
    device_limit : int
        Maximum color index the PGPLOT device supports (from pgqcol).
    """
    all_rgb = get_qp_module_rgb()

    # Build list: index, name, expected_rgb, pgplot_rgb
    entries = []
    for i, name in enumerate(QP_ORIGINAL_NAMES):
        expected = all_rgb[i]
        rendered = pgplot_rendered_color(i, all_rgb, device_limit)
        entries.append((i, name, expected, rendered))
    # Skip transparent (16)
    for i, name in enumerate(QP_CSS4_NAMES):
        idx = i + 17
        expected = all_rgb[idx]
        rendered = pgplot_rendered_color(idx, all_rgb, device_limit)
        entries.append((idx, name, expected, rendered))

    n = len(entries)
    per_page = 46
    n_pages = (n + per_page - 1) // per_page

    for page in range(n_pages):
        start = page * per_page
        end = min(start + per_page, n)
        page_entries = entries[start:end]
        n_rows = len(page_entries)

        fig, ax = plt.subplots(1, 1, figsize=(12, n_rows * 0.4 + 1.5))
        ax.set_xlim(0, 12)
        ax.set_ylim(0, n_rows + 1)
        ax.set_axis_off()
        ax.set_title(
            f"PGPLOT rendering simulation (device limit={device_limit}) — "
            f"page {page+1}/{n_pages}",
            fontsize=12, fontweight='bold', pad=10
        )

        for i, (idx, name, expected, rendered) in enumerate(page_entries):
            y = n_rows - i - 0.5

            exp_norm = (expected[0]/255, expected[1]/255, expected[2]/255)
            ren_norm = (rendered[0]/255, rendered[1]/255, rendered[2]/255)

            # Expected color swatch
            rect_exp = mpatches.FancyBboxPatch(
                (0.3, y - 0.35), 2.2, 0.7,
                boxstyle="round,pad=0.02",
                facecolor=exp_norm, edgecolor='gray', linewidth=0.5
            )
            ax.add_patch(rect_exp)

            # PGPLOT rendered swatch
            rect_pg = mpatches.FancyBboxPatch(
                (2.7, y - 0.35), 2.2, 0.7,
                boxstyle="round,pad=0.02",
                facecolor=ren_norm, edgecolor='gray', linewidth=0.5
            )
            ax.add_patch(rect_pg)

            # Color name and index
            label = f"{idx:3d}: {name}"
            ax.text(5.2, y, label, va='center', fontsize=8,
                    fontfamily='monospace')

            # Show if above device limit
            if idx > device_limit:
                nearest_idx = nearest_color_index(idx, all_rgb, device_limit)
                nearest_name = ""
                if nearest_idx < 17:
                    nearest_name = QP_ORIGINAL_NAMES[nearest_idx]
                else:
                    nearest_name = QP_CSS4_NAMES[nearest_idx - 17]
                ax.text(9.0, y, f"→ {nearest_idx}: {nearest_name}",
                        va='center', fontsize=7, fontfamily='monospace',
                        color='red')

        # Column headers
        ax.text(1.4, n_rows + 0.4, "Expected", ha='center', fontsize=9,
                fontweight='bold')
        ax.text(3.8, n_rows + 0.4, "PGPLOT", ha='center', fontsize=9,
                fontweight='bold')
        ax.text(5.2, n_rows + 0.4, "Name", ha='left', fontsize=9,
                fontweight='bold')
        ax.text(9.0, n_rows + 0.4, "Fallback", ha='left', fontsize=9,
                fontweight='bold', color='red')

        plt.tight_layout()

    plt.show()


def plot_original_colors():
    """Show the original 16 quick_plot colors."""
    all_rgb = get_qp_module_rgb()

    fig, ax = plt.subplots(1, 1, figsize=(8, 4))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 17)
    ax.set_axis_off()
    ax.set_title("Original quick_plot colors (indices 0-15)", fontsize=12,
                 fontweight='bold')

    for i, name in enumerate(QP_ORIGINAL_NAMES):
        y = 16 - i - 0.5
        rgb = all_rgb[i]
        norm = (rgb[0] / 255.0, rgb[1] / 255.0, rgb[2] / 255.0)

        rect = mpatches.FancyBboxPatch(
            (0.5, y - 0.4), 3.0, 0.8,
            boxstyle="round,pad=0.02",
            facecolor=norm, edgecolor='gray', linewidth=0.5
        )
        ax.add_patch(rect)
        ax.text(4.0, y, f"{i:2d}: {name}", va='center', fontsize=9,
                fontfamily='monospace')

    plt.tight_layout()
    plt.show()


def hls_to_rgb(h, l, s):
    """Convert HLS to RGB (0-1 scale).

    Parameters
    ----------
    h : float
        Hue in degrees (0-360).
    l : float
        Lightness (0-1).
    s : float
        Saturation (0-1).

    Returns
    -------
    tuple
        (r, g, b) each 0-1.
    """
    import colorsys
    # colorsys expects h in 0-1
    r, g, b = colorsys.hls_to_rgb(h / 360.0, l, s)
    return (r, g, b)


def plplot_continuous_color(t):
    """Compute the PLPLOT continuous colormap color at position t.

    PLPLOT uses: plscmap1l(.false., [0, 1], [240, 0], [0.6, 0.6], [0.8, 0.8])
    which is HLS interpolation from hue=240 to hue=0, sat=0.6, light=0.8.
    The .false. means hue interpolates through lower values: 240→180→120→60→0.

    Parameters
    ----------
    t : float
        Position in colormap (0.0 to 1.0).

    Returns
    -------
    tuple
        (r, g, b) each 0-1.
    """
    hue = 240.0 * (1.0 - t)  # 240 at t=0, 0 at t=1
    lightness = 0.6
    saturation = 0.8
    return hls_to_rgb(hue, lightness, saturation)


def pgplot_continuous_color(t):
    """Compute the improved PGPLOT continuous color at position t.

    Now matches PLPLOT: HLS with hue 240→0, saturation=0.8, lightness=0.6.
    Uses pgscr to define the color on-the-fly at a scratch index.

    Parameters
    ----------
    t : float
        Position in colormap (0.0 to 1.0).

    Returns
    -------
    tuple
        (r, g, b) each 0-1.
    """
    # Same HLS formula as PLPLOT
    return plplot_continuous_color(t)


def plot_continuous_colors():
    """Show PLPLOT vs PGPLOT continuous colormaps side-by-side."""
    n_steps = 256
    fig, axes = plt.subplots(2, 1, figsize=(12, 3))

    # PLPLOT continuous colorbar
    plplot_colors = np.zeros((1, n_steps, 3))
    for i in range(n_steps):
        t = i / (n_steps - 1)
        plplot_colors[0, i, :] = plplot_continuous_color(t)

    axes[0].imshow(plplot_colors, aspect='auto', extent=[0, 1, 0, 1])
    axes[0].set_yticks([])
    axes[0].set_title("PLPLOT continuous colormap (HLS: hue 240→0, sat=0.8, light=0.6)",
                      fontsize=11, fontweight='bold')
    axes[0].set_xlabel("real_color (0.0 → 1.0)")

    # PGPLOT continuous colorbar (now matches PLPLOT via pgscr)
    pgplot_colors = np.zeros((1, n_steps, 3))
    for i in range(n_steps):
        t = i / (n_steps - 1)
        pgplot_colors[0, i, :] = pgplot_continuous_color(t)

    axes[1].imshow(pgplot_colors, aspect='auto', extent=[0, 1, 0, 1])
    axes[1].set_yticks([])
    axes[1].set_title("PGPLOT continuous colormap (now matches PLPLOT via pgscr on-the-fly)",
                      fontsize=11, fontweight='bold')
    axes[1].set_xlabel("real_color (0.0 → 1.0)")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Simulate PGPLOT CSS4 color rendering vs expected")
    parser.add_argument('--original', action='store_true',
                        help="Show original 16 quick_plot colors only")
    parser.add_argument('--continuous', action='store_true',
                        help="Show PLPLOT vs PGPLOT continuous colormaps")
    parser.add_argument('--device-limit', type=int, default=99,
                        help="PGPLOT device color index limit (default: 99)")
    args = parser.parse_args()

    if args.original:
        plot_original_colors()
    elif args.continuous:
        plot_continuous_colors()
    else:
        plot_pgplot_comparison(device_limit=args.device_limit)
