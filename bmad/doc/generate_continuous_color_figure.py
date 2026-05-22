#!/usr/bin/env python3
"""Generate the continuous colormap figure for quick-plot.tex.

Produces plot-continuous-color.pdf showing the HLS-based continuous colormap
used by qp_continuous_color (hue 240->0, saturation 0.8, lightness 0.6).
"""

import colorsys

import matplotlib.pyplot as plt
import numpy as np


def hls_to_rgb(hue_deg, lightness, saturation):
    """Convert HLS (degrees, 0-1, 0-1) to RGB (0-1, 0-1, 0-1).

    Parameters
    ----------
    hue_deg : float
        Hue in degrees (0-360).
    lightness : float
        Lightness (0-1).
    saturation : float
        Saturation (0-1).

    Returns
    -------
    tuple
        (r, g, b) each 0-1.
    """
    return colorsys.hls_to_rgb(hue_deg / 360.0, lightness, saturation)


def main():
    """Generate the continuous color figure."""
    n_steps = 512
    colors = np.zeros((1, n_steps, 3))
    for i in range(n_steps):
        t = i / (n_steps - 1)
        hue = 240.0 * (1.0 - t)  # 240 (blue) at t=0, 0 (red) at t=1
        colors[0, i, :] = hls_to_rgb(hue, 0.6, 0.8)

    fig, ax = plt.subplots(figsize=(8, 1.2))
    ax.imshow(colors, aspect='auto', extent=[0, 1, 0, 1])
    ax.set_yticks([])
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xlabel("real_color (0.0 → 1.0)", fontsize=10)
    ax.set_title("qp_continuous_color: HLS colormap (hue 240°→0°, S=0.8, L=0.6)",
                 fontsize=11)
    plt.tight_layout()
    plt.savefig("plot-continuous-color.pdf", bbox_inches='tight')
    print("Wrote plot-continuous-color.pdf")


if __name__ == "__main__":
    main()
