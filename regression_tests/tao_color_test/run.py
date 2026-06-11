"""Regression test for CSS4 named colors, hex, and RGB color formats in Tao.

Verifies that Tao correctly accepts and stores color values in all supported
formats: CSS4 named colors, hex (#RRGGBB, #RGB), and RGB(r,g,b).
Also tests normalization (case, whitespace) and label_color.
"""

import os
import subprocess
import sys

bin_dir = sys.argv[1]

# Use the small lattice from tao_test (relative to this test directory)
lat_file = os.path.join(os.path.dirname(__file__), "..", "tao_test", "small.bmad")
curve = "r13.g.a"

# Colors to test: (format_label, color_string, expected_readback)
test_cases = [
    # CSS4 named colors (normalization: all stored lowercase)
    ("css4_name", "teal", "teal"),
    ("css4_name_mixed_case", "Teal", "teal"),
    ("css4_name_upper", "TEAL", "teal"),
    ("css4_name_salmon", "salmon", "salmon"),
    ("css4_name_gold", "gold", "gold"),
    # Original color names (also normalized to lowercase)
    ("original_name", "red", "red"),
    ("original_name_blue", "blue", "blue"),
    ("original_name_upper", "RED", "red"),
    # Hex colors (normalization: digits uppercased)
    ("hex_6digit", "#FF5733", "#FF5733"),
    ("hex_6digit_lower", "#ff5733", "#FF5733"),
    ("hex_6digit_mixed", "#fF5733", "#FF5733"),
    ("hex_3digit", "#F53", "#F53"),
    ("hex_3digit_lower", "#abc", "#ABC"),
    # RGB format (normalization: strip spaces, uppercase prefix)
    ("rgb_format", "RGB(0,128,255)", "RGB(0,128,255)"),
    ("rgb_format_spaces", "RGB(255, 0, 128)", "RGB(255,0,128)"),
    ("rgb_lower_prefix", "rgb(10,20,30)", "RGB(10,20,30)"),
    ("rgb_mixed_prefix", "Rgb(1,2,3)", "RGB(1,2,3)"),
    # Invalid RGB (should be rejected)
    ("rgb_four_values", "RGB(1,2,3,4)", None),
]

results = []

for label, color_str, expected in test_cases:
    # Build the tao command: set the color, then show the curve to read it back
    cmd = (
        f"set curve {curve} line%color = {color_str};"
        f"show curve {curve};"
        f"quit"
    )
    exe = f'{bin_dir}tao -noplot -lat {lat_file} -command "{cmd}"'
    proc = subprocess.run(exe, shell=True, capture_output=True, text=True)

    # Parse the line%color from show output
    actual = ""
    for line in proc.stdout.splitlines():
        if "line%color" in line and "z_color" not in line and "=" in line:
            # Format: "line%color           = teal"
            actual = line.split("=", 1)[1].strip()
            break

    # Check for errors
    has_error = "[ERROR" in proc.stdout or proc.returncode != 0

    if expected is None:
        # Expect rejection (error message)
        results.append(f'"{label}" STR "{has_error}"')
    elif has_error:
        results.append(f'"{label}" STR "ERROR"')
    else:
        results.append(f'"{label}" STR "{actual}"')

# Test that an invalid color produces an error (doesn't crash)
cmd = f"set curve {curve} line%color = not_a_color;quit"
exe = f'{bin_dir}tao -noplot -lat {lat_file} -command "{cmd}"'
proc = subprocess.run(exe, shell=True, capture_output=True, text=True)
# Should get an error message but not crash (exit 0)
has_error_msg = "ERROR" in proc.stdout
results.append(f'"invalid_color_rejected" STR "{has_error_msg}"')

# Test label_color via python plot_graph command
label_color_tests = [
    ("label_color_css4", "teal", "teal"),
    ("label_color_upper", "GOLD", "gold"),
    ("label_color_hex", "#FF0000", "#FF0000"),
]

for label, color_str, expected in label_color_tests:
    cmd = (
        f"set plot r13 x%label_color = {color_str};"
        f"python plot_graph r13.g;"
        f"quit"
    )
    exe = f'{bin_dir}tao -noplot -lat {lat_file} -command "{cmd}"'
    proc = subprocess.run(exe, shell=True, capture_output=True, text=True)

    # Parse x.label_color from python output
    actual = ""
    for line in proc.stdout.splitlines():
        if "x.label_color" in line:
            # Format: "x.label_color;ENUM;T;teal"
            parts = line.split(";")
            if len(parts) >= 4:
                actual = parts[3].strip()
            break

    has_error = "[ERROR" in proc.stdout or proc.returncode != 0
    if has_error:
        results.append(f'"{label}" STR "ERROR"')
    else:
        results.append(f'"{label}" STR "{actual}"')

# Write output.now
with open("output.now", "w") as f:
    for r in results:
        f.write(r + "\n")

print("tao_color_test: All tests completed.")
for r in results:
    print(" ", r)
