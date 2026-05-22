"""Regression test for CSS4 named colors, hex, and RGB color formats in Tao.

Verifies that Tao correctly accepts and stores color values in all supported
formats: CSS4 named colors, hex (#RRGGBB, #RGB), and RGB(r,g,b).
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
    ("css4_name", "teal", "teal"),
    ("css4_name_mixed_case", "Teal", "teal"),
    ("css4_name_salmon", "salmon", "salmon"),
    ("css4_name_gold", "gold", "gold"),
    ("original_name", "red", "red"),
    ("original_name_blue", "blue", "blue"),
    ("hex_6digit", "#FF5733", "#FF5733"),
    ("hex_6digit_lower", "#ff5733", "#FF5733"),
    ("hex_3digit", "#F53", "#FF5533"),
    ("rgb_format", "RGB(0,128,255)", "RGB(0,128,255)"),
    ("rgb_format_spaces", "RGB(255, 0, 128)", "RGB(255,0,128)"),
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
        if "line%color" in line and "z_color" not in line:
            # Format: "line%color           = teal"
            actual = line.split("=", 1)[1].strip()
            break

    # Check for errors
    has_error = "[ERROR" in proc.stdout or proc.returncode != 0

    if has_error:
        results.append(f'"{label}" STR ERROR')
    else:
        results.append(f'"{label}" STR {actual}')

# Test that an invalid color produces an error (doesn't crash)
cmd = f"set curve {curve} line%color = not_a_color;quit"
exe = f'{bin_dir}tao -noplot -lat {lat_file} -command "{cmd}"'
proc = subprocess.run(exe, shell=True, capture_output=True, text=True)
# Should get an error message but not crash (exit 0)
has_error_msg = "ERROR" in proc.stdout
results.append(f'"invalid_color_rejected" STR {has_error_msg}')

# Write output.now
with open("output.now", "w") as f:
    for r in results:
        f.write(r + "\n")

print("tao_color_test: All tests completed.")
for r in results:
    print(" ", r)
