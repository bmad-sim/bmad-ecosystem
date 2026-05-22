"""
Generate CSS4 color reference swatch and verify quick_plot color module data.

This script:
1. Generates a reference PNG showing all CSS4 colors available in quick_plot.
2. Verifies the RGB values in qp_css4_colors_mod.f90 match the W3C CSS4 spec.

Usage:
    python3 css4_color_reference.py
"""

import os
import struct
import zlib

# CSS4 hex definitions (W3C CSS Color Level 4 specification)
CSS4_HEX = {
    'aliceblue': 'F0F8FF', 'antiquewhite': 'FAEBD7', 'aqua': '00FFFF',
    'aquamarine': '7FFFD4', 'azure': 'F0FFFF', 'beige': 'F5F5DC',
    'bisque': 'FFE4C4', 'blanchedalmond': 'FFEBCD', 'blueviolet': '8A2BE2',
    'brown': 'A52A2A', 'burlywood': 'DEB887', 'cadetblue': '5F9EA0',
    'chartreuse': '7FFF00', 'chocolate': 'D2691E', 'coral': 'FF7F50',
    'cornflowerblue': '6495ED', 'cornsilk': 'FFF8DC', 'crimson': 'DC143C',
    'darkblue': '00008B', 'darkcyan': '008B8B', 'darkgoldenrod': 'B8860B',
    'darkgray': 'A9A9A9', 'darkgreen': '006400', 'darkgrey': 'A9A9A9',
    'darkkhaki': 'BDB76B', 'darkmagenta': '8B008B', 'darkolivegreen': '556B2F',
    'darkorange': 'FF8C00', 'darkorchid': '9932CC', 'darkred': '8B0000',
    'darksalmon': 'E9967A', 'darkseagreen': '8FBC8F', 'darkslateblue': '483D8B',
    'darkslategray': '2F4F4F', 'darkslategrey': '2F4F4F', 'darkturquoise': '00CED1',
    'darkviolet': '9400D3', 'deeppink': 'FF1493', 'deepskyblue': '00BFFF',
    'dimgray': '696969', 'dimgrey': '696969', 'dodgerblue': '1E90FF',
    'firebrick': 'B22222', 'floralwhite': 'FFFAF0', 'forestgreen': '228B22',
    'fuchsia': 'FF00FF', 'gainsboro': 'DCDCDC', 'ghostwhite': 'F8F8FF',
    'gold': 'FFD700', 'goldenrod': 'DAA520', 'gray': '808080',
    'greenyellow': 'ADFF2F', 'grey': '808080', 'honeydew': 'F0FFF0',
    'hotpink': 'FF69B4', 'indianred': 'CD5C5C', 'indigo': '4B0082',
    'ivory': 'FFFFF0', 'khaki': 'F0E68C', 'lavender': 'E6E6FA',
    'lavenderblush': 'FFF0F5', 'lawngreen': '7CFC00', 'lemonchiffon': 'FFFACD',
    'lightblue': 'ADD8E6', 'lightcoral': 'F08080', 'lightcyan': 'E0FFFF',
    'lightgoldenrodyellow': 'FAFAD2', 'lightgray': 'D3D3D3', 'lightgreen': '90EE90',
    'lightgrey': 'D3D3D3', 'lightpink': 'FFB6C1', 'lightsalmon': 'FFA07A',
    'lightseagreen': '20B2AA', 'lightskyblue': '87CEFA', 'lightslategray': '778899',
    'lightslategrey': '778899', 'lightsteelblue': 'B0C4DE', 'lightyellow': 'FFFFE0',
    'lime': '00FF00', 'limegreen': '32CD32', 'linen': 'FAF0E6',
    'maroon': '800000', 'mediumaquamarine': '66CDAA', 'mediumblue': '0000CD',
    'mediumorchid': 'BA55D3', 'mediumpurple': '9370DB', 'mediumseagreen': '3CB371',
    'mediumslateblue': '7B68EE', 'mediumspringgreen': '00FA9A',
    'mediumturquoise': '48D1CC', 'mediumvioletred': 'C71585', 'midnightblue': '191970',
    'mintcream': 'F5FFFA', 'mistyrose': 'FFE4E1', 'moccasin': 'FFE4B5',
    'navajowhite': 'FFDEAD', 'navy': '000080', 'oldlace': 'FDF5E6',
    'olive': '808000', 'olivedrab': '6B8E23', 'orangered': 'FF4500',
    'orchid': 'DA70D6', 'palegoldenrod': 'EEE8AA', 'palegreen': '98FB98',
    'paleturquoise': 'AFEEEE', 'palevioletred': 'DB7093', 'papayawhip': 'FFEFD5',
    'peachpuff': 'FFDAB9', 'peru': 'CD853F', 'pink': 'FFC0CB',
    'plum': 'DDA0DD', 'powderblue': 'B0E0E6', 'rebeccapurple': '663399',
    'rosybrown': 'BC8F8F', 'royalblue': '4169E1', 'saddlebrown': '8B4513',
    'salmon': 'FA8072', 'sandybrown': 'F4A460', 'seagreen': '2E8B57',
    'seashell': 'FFF5EE', 'sienna': 'A0522D', 'silver': 'C0C0C0',
    'skyblue': '87CEEB', 'slateblue': '6A5ACD', 'slategray': '708090',
    'slategrey': '708090', 'snow': 'FFFAFA', 'springgreen': '00FF7F',
    'steelblue': '4682B4', 'tan': 'D2B48C', 'teal': '008080',
    'thistle': 'D8BFD8', 'tomato': 'FF6347', 'turquoise': '40E0D0',
    'violet': 'EE82EE', 'wheat': 'F5DEB3', 'whitesmoke': 'F5F5F5',
    'yellowgreen': '9ACD32'
}

# Quick_plot CSS4 color names in index order (17:154)
# Read from module to stay in sync with reordering
def _read_css4_names_from_module():
    """Read CSS4 color names from the Fortran module file.

    Returns
    -------
    list of str
        Color names in index order (17 to 154).
    """
    import re
    module_path = os.path.join(os.path.dirname(__file__),
                               '..', '..', 'sim_utils', 'plot',
                               'qp_css4_colors_mod.f90')
    with open(module_path) as f:
        content = f.read()
    # Extract names from the qp_css4_color_name array
    pattern = r"qp_css4_color_name\(17:154\)\s*=\s*\[\s*character\(24\)\s*::\s*&(.*?)\]"
    match = re.search(pattern, content, re.DOTALL)
    if not match:
        raise ValueError("Could not parse qp_css4_color_name from module")
    block = match.group(1)
    names = re.findall(r"'(\w+)'", block)
    return names

QP_CSS4_NAMES = _read_css4_names_from_module()

# Original quick_plot colors (indices 0:15)
QP_ORIGINAL_NAMES = [
    'white', 'black', 'red', 'green', 'blue', 'cyan', 'magenta', 'yellow',
    'orange', 'yellow_green', 'light_green', 'navy_blue', 'purple',
    'reddish_purple', 'dark_grey', 'light_grey'
]

QP_ORIGINAL_RGB = [
    (255, 255, 255), (0, 0, 0), (255, 0, 0), (0, 128, 0),
    (0, 0, 255), (0, 255, 255), (255, 0, 255), (255, 255, 0),
    (255, 165, 0), (127, 255, 0), (0, 255, 127), (0, 127, 255),
    (128, 0, 128), (255, 0, 127), (85, 85, 85), (170, 170, 170)
]


def parse_fortran_module(filepath):
    """Parse the RGB arrays from qp_css4_colors_mod.f90."""
    with open(filepath, 'r') as f:
        content = f.read()

    def extract_array(name):
        """Extract integer values from a Fortran parameter array."""
        import re
        # Find the array declaration
        pattern = rf'{name}\(0:154\)\s*=\s*\[\s*&(.*?)\]'
        match = re.search(pattern, content, re.DOTALL)
        if not match:
            raise ValueError(f"Could not find array {name}")
        array_text = match.group(1)
        # Remove comments, continuation chars, and extract numbers
        values = []
        for line in array_text.split('\n'):
            line = re.sub(r'!.*', '', line)  # Remove comments
            line = line.replace('&', '').strip()
            if line:
                for token in line.split(','):
                    token = token.strip()
                    if token:
                        values.append(int(token))
        return values

    red = extract_array('qp_color_red')
    green = extract_array('qp_color_green')
    blue = extract_array('qp_color_blue')
    return red, green, blue


def verify_module_data(module_path):
    """Verify module RGB data matches CSS4 spec."""
    red, green, blue = parse_fortran_module(module_path)

    print(f"Parsed {len(red)} red, {len(green)} green, {len(blue)} blue values")
    assert len(red) == 155, f"Expected 155 red values, got {len(red)}"
    assert len(green) == 155, f"Expected 155 green values, got {len(green)}"
    assert len(blue) == 155, f"Expected 155 blue values, got {len(blue)}"

    # Verify original 16 colors (indices 0-15)
    print("\nVerifying original 16 colors...")
    orig_errors = 0
    for i, (name, (er, eg, eb)) in enumerate(zip(QP_ORIGINAL_NAMES, QP_ORIGINAL_RGB)):
        if red[i] != er or green[i] != eg or blue[i] != eb:
            print(f"  MISMATCH index {i} ({name}): "
                  f"got ({red[i]},{green[i]},{blue[i]}) expected ({er},{eg},{eb})")
            orig_errors += 1
    if orig_errors == 0:
        print("  All 16 original colors OK")

    # Verify CSS4 colors (indices 17-154)
    print("\nVerifying 138 CSS4 colors...")
    css4_errors = 0
    for i, name in enumerate(QP_CSS4_NAMES):
        idx = i + 17
        h = CSS4_HEX[name]
        er = int(h[0:2], 16)
        eg = int(h[2:4], 16)
        eb = int(h[4:6], 16)
        if red[idx] != er or green[idx] != eg or blue[idx] != eb:
            print(f"  MISMATCH index {idx} ({name}): "
                  f"got ({red[idx]},{green[idx]},{blue[idx]}) expected ({er},{eg},{eb})")
            css4_errors += 1

    if css4_errors == 0:
        print("  All 138 CSS4 colors OK")

    total_errors = orig_errors + css4_errors
    print(f"\nTotal errors: {total_errors}")
    return total_errors == 0


def write_png(filename, width, height, pixels):
    """
    Write a minimal PNG file.

    Parameters
    ----------
    filename : str
        Output file path.
    width : int
        Image width in pixels.
    height : int
        Image height in pixels.
    pixels : list of list of tuple
        pixels[y][x] = (r, g, b), each 0-255.
    """
    def make_chunk(chunk_type, data):
        chunk = chunk_type + data
        return (struct.pack('>I', len(data)) + chunk +
                struct.pack('>I', zlib.crc32(chunk) & 0xffffffff))

    # PNG signature
    signature = b'\x89PNG\r\n\x1a\n'

    # IHDR
    ihdr_data = struct.pack('>IIBBBBB', width, height, 8, 2, 0, 0, 0)
    ihdr = make_chunk(b'IHDR', ihdr_data)

    # IDAT
    raw_data = b''
    for row in pixels:
        raw_data += b'\x00'  # filter byte
        for r, g, b in row:
            raw_data += struct.pack('BBB', r, g, b)
    compressed = zlib.compress(raw_data)
    idat = make_chunk(b'IDAT', compressed)

    # IEND
    iend = make_chunk(b'IEND', b'')

    with open(filename, 'wb') as f:
        f.write(signature + ihdr + idat + iend)


def generate_reference_swatch(output_path):
    """
    Generate a reference color swatch PNG showing all quick_plot colors.

    Layout: grid of color swatches, 10 columns.
    """
    swatch_w = 60
    swatch_h = 30
    cols = 10
    n_colors = 16 + 138  # original + CSS4
    rows = (n_colors + cols - 1) // cols

    width = cols * swatch_w
    height = rows * swatch_h

    # Build color list: original 16 + CSS4 138
    all_colors = []
    for name, rgb in zip(QP_ORIGINAL_NAMES, QP_ORIGINAL_RGB):
        all_colors.append((name, rgb))
    for name in QP_CSS4_NAMES:
        h = CSS4_HEX[name]
        rgb = (int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16))
        all_colors.append((name, rgb))

    # Create pixel array
    pixels = [[(255, 255, 255)] * width for _ in range(height)]

    for idx, (name, (r, g, b)) in enumerate(all_colors):
        col = idx % cols
        row = idx // cols
        x0 = col * swatch_w
        y0 = row * swatch_h
        for y in range(y0, min(y0 + swatch_h, height)):
            for x in range(x0, min(x0 + swatch_w, width)):
                pixels[y][x] = (r, g, b)

    write_png(output_path, width, height, pixels)
    print(f"Reference swatch written to: {output_path}")

    # Also write a text file with the expected colors for programmatic comparison
    txt_path = output_path.replace('.png', '_expected.txt')
    with open(txt_path, 'w') as f:
        f.write("# Quick_plot CSS4 color test - expected RGB values\n")
        f.write("# index  name                     R    G    B\n")
        for i, (name, rgb) in enumerate(zip(QP_ORIGINAL_NAMES, QP_ORIGINAL_RGB)):
            f.write(f"{i:5d}  {name:24s} {rgb[0]:3d}  {rgb[1]:3d}  {rgb[2]:3d}\n")
        f.write(f"   16  {'transparent':24s}   0    0    0\n")
        for i, name in enumerate(QP_CSS4_NAMES):
            idx = i + 17
            h = CSS4_HEX[name]
            r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
            f.write(f"{idx:5d}  {name:24s} {r:3d}  {g:3d}  {b:3d}\n")
    print(f"Expected values written to: {txt_path}")


if __name__ == '__main__':
    script_dir = os.path.dirname(os.path.abspath(__file__))
    module_path = os.path.join(script_dir, '../../sim_utils/plot/qp_css4_colors_mod.f90')

    # Verify module data
    print("=" * 60)
    print("Verifying qp_css4_colors_mod.f90 RGB data against CSS4 spec")
    print("=" * 60)
    ok = verify_module_data(module_path)

    # Generate reference swatch
    print("\n" + "=" * 60)
    print("Generating reference color swatch")
    print("=" * 60)
    swatch_path = os.path.join(script_dir, 'css4_reference_swatch.png')
    generate_reference_swatch(swatch_path)

    if ok:
        print("\n*** ALL CHECKS PASSED ***")
    else:
        print("\n*** ERRORS FOUND - see above ***")
        exit(1)
