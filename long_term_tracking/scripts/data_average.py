# Python script to do averaging, over some number of rows of data, in files created by the long_term_tracking program.
# The output is a table of averaged values.
# This script is meant to be modified to suit the current needs.

import sys
import numpy as np


datf = open(sys.argv[1], 'r')

# Read header parameters (potentially useful).

T = True      # Since "T" and "F" are used in the data file parameter header
F = False

for line in datf:
  if line[0] == '#':
    if '=' not in line: continue    # Not a parameter
    line = line.replace('%', '_')   # EG: ltt%ramping_on -> ltt_ramping_on
    exec(line[1:].strip())

  else:
    break

# Now do a centered average.  
# With n_ave = 5, the output will not contain rows corresponding to the first two rows and the last two
# rows of the input data file.

n_ave = 5     # Averaging window. The algorithm assumes this is an odd number
n_here = -1

rows = [None] * n_ave
out_name = sys.argv[1] + '.ave'
outf = open(out_name, 'w')
datf.seek(0)

for line in datf:
  if 'Turn' in line: outf.write(line)
  if line.strip() == '': continue
  if line[0] == '#': continue

  cols = np.array(list(map(float, line.split())))

  if n_here == n_ave-1:
    rows = rows[1:] + [cols]
  else:
    n_here += 1
    rows[n_here] = cols

  if n_here == n_ave-1:
    sum_row = sum(rows) / n_ave
    outf.write(' '.join(list(map(str, sum_row))) + '\n')

print ('Created: ' + out_name)
