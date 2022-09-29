#!/usr/bin/env python3





def ix(n):
  return '{0:0>2}'.format(n)



def fixline(line, ixcell, prefix):
  newline = line
  n = ixcell*2
  newline = newline.replace('FX', prefix)
  newline = newline.replace('X0', str(ix(n-2)))
  newline = newline.replace('X1', str(ix(n-1)))
  newline = newline.replace('X2', str(ix(n)))
  newline = newline.replace('ZZ', str(ix(ixcell)))
  return newline


def parse_template(file):
  template = []
  with open(file, 'r') as f:
    for line in f:
      template.append(line)
  return template

      
  
def write_section(template, outfile, n1, n2, prefix):
    for n in range(n1,n2+1):
      print('writing cell', n, 'with prefix', prefix)
      for line in template:
        x = fixline(line, n, prefix)
        outfile.write(x)


def write_cell_lines(outfile, name, n1, n2, prefix):
  print('writing line', name)
  outfile.write('!--- '+name+'\n')
  outfile.write(name+': line = (\n')
  for n in range(n1,n2):
    outfile.write('  '+prefix+'.CELL'+str(ix(n))+',\n')
  # final cell
  outfile.write  ('  '+prefix+'.CELL'+str(ix(n2))+')\n')
  


def write_cell_geometry(outfile, n, factor, prefix):
  print('writing geometry for cell', n, 'with prefix', prefix)
  ix1 = str(ix(2*n-1))
  ix2 = str(ix(2*n))
  q1=prefix+'.Qua'+ix1
  q2=prefix+'.Qua'+ix2
  d1=prefix+'.Pip'+ix1
  d2=prefix+'.Pip'+ix2
  patch1=prefix+'.patch'+ix1
  patch2=prefix+'.patch'+ix2
  #patch3=prefix+'.patch'+ix2+'a'
  #patch4=prefix+'.patch'+ix2+'b'
  outfile.write('!--- Geometry for '+prefix+' cell '+str(n)+'\n')
  outfile.write(q1+'[x_offset] = FF.Qua01[x_offset]*'+str(factor)+'\n')
  outfile.write(q2+'[x_offset] = FF.Qua02[x_offset]*'+str(factor)+'\n')
  outfile.write(patch1+'[x_pitch] = FF.patch1[x_pitch]*'+str(factor)+'\n')
  outfile.write(patch2+'[x_pitch] = FF.patch2[x_pitch]*'+str(factor)+'\n')
  #outfile.write(patch3+'[x_pitch] = FF.patch3[x_pitch]*'+str(factor)+'\n')
  #outfile.write(patch4+'[x_pitch] = FF.patch4[x_pitch]*'+str(factor)+'\n')
  
  #drift_stretch 
  outfile.write(d1+'[L] = FF.Pip01[L] + FF.stretch*'+ str(1-factor)+'\n')
  outfile.write(d2+'[L] = FF.Pip02[L] + FF.stretch*'+ str(1-factor)+'\n')
  
  # Move quads by stretch
  outfile.write(q1+'[offset] = -FF.quad_padding - FF.stretch*'+str(1-factor)+'\n')
  outfile.write(q2+'[offset] = +FF.quad_padding + FF.stretch*'+str(1-factor)+'\n')
  outfile.write('\n')


# Simple function f(0) = 0, f(1) = 1, f'(0) = f'(1) = 0
def f1(x):
  return 3*x**2 - 2*x**3

# Higher order merge
def f2(x):
  a5 = -1.45
  a7 = -7.32
  return x**2 * (3-2*x) + a5*x**2 * (1-x)**2 * (1-2*x) + a7*x**3*(1-x)**3 * (1-2*x)

def f3(x):
  a0 = 1
  a1 = 2
  a2 = 3.816
  a3 = 7.820
  a4 = 19.810
  return 0.5 + (x-0.5)*a0 + (x-0.5)*a1 * (x*(1-x)) + (x-0.5)*a2 * (x*(1-x))**2 + (x-0.5)*a3 * (x*(1-x))**3 + (x-0.5)*a4 * (x*(1-x))**4

def f4(x):
  # Scott Berg March 30, 2017
  return   x - (0.5 - x)*(1 - x)*x*(1.788 + 3.954*(1 - x)*x + 6.58*(1 - x)**2*x**2)


def write_geometry(outfile, n1, n2, prefix, reverse=False):
  for n in range(n1, n2+1):
    factor = 1 - f4(n/(n2-n1 +2))
    if reverse:
      factor = 1 - factor
    write_cell_geometry(outfile, n, factor, prefix)


# BPM for cell n
def bpm(n, prefix=''):
  ix1 = ix(n)
  ix2 = ix(2*n-1)
  ele = '.BLK'+ix1
  s = prefix+'.BPM'+ix1+': instrument, superimpose, ref = '+prefix+ele +'\n'
  print(s)
  return s

def write_bpms(filename, locell, ncells, prefix):
  f=open(filename, 'w')
  for n in range(locell, ncells+1):
    f.write( bpm(n, prefix))
  f.close()



N_ARC = 16
N_ARC_EXTRA = 0
N_TRANSITION  = 12*2
N_FA = N_ARC 
N_FB = N_ARC
N_TA = N_TRANSITION
N_ZA = 13
N_ZB = 14

TEMPLATE = 'cell/template_cell.bmad'
template = parse_template(TEMPLATE)

# FA
fout = open('fa.cells.bmad', 'w')
write_section(template, fout, 1, N_FA+N_ARC_EXTRA, 'FA')
write_cell_lines(fout, 'FA.arc_cells', 1+N_ARC_EXTRA, N_ARC+N_ARC_EXTRA, 'FA')
fout.close()

# TA
fout = open('../tx/ta.cells.bmad', 'w')
write_section(template, fout, 1, N_TA, 'TA')
write_cell_lines(fout, 'TA.transition_cells', 1,  N_TA, 'TA') 
fout.close()

fout = open('../tx/ta.geometry.bmad', 'w')
write_geometry(fout, 1, N_TA, 'TA')
fout.close()


# FB
fout = open('fb.cells.bmad', 'w')
write_section(template, fout, 1, N_FB+N_ARC_EXTRA+1, 'FB')
write_cell_lines(fout, 'FB.arc_cells', 1, N_FB, 'FB')
fout.close()

# TB
fout = open('../tx/tb.cells.bmad', 'w')
write_section(template, fout, 1, N_TA, 'TB')
write_cell_lines(fout, 'TB.transition_cells', 1,  N_TA, 'TB') 
fout = open('../tx/tb.geometry.bmad', 'w')
write_geometry(fout, 1, N_TA, 'TB', reverse = True)
fout.close()

# ZA
template = parse_template('cell/template_straight_cell.bmad')
fout = open('../zx/za.cells.bmad', 'w')
write_section(template, fout, 1, N_ZA, 'ZA')
write_cell_lines(fout, 'ZA.straight_cells', 1, N_ZA, 'ZA')
fout.close()

# ZB
fout = open('../zx/zb.cells.bmad', 'w')
write_section(template, fout, 1, N_ZB, 'ZB')
write_cell_lines(fout, 'ZB.straight_cells', 1, N_ZB, 'ZB')
fout.close()


# BPMs
write_bpms('fa.bpms.bmad', 0, N_FA, 'FA')
write_bpms('fb.bpms.bmad', 1, N_FB+1, 'FB')
write_bpms('../tx/ta.bpms.bmad', 1, N_TA, 'TA')
write_bpms('../tx/tb.bpms.bmad', 1, N_TA, 'TB')
write_bpms('../zx/za.bpms.bmad', 1, N_ZA, 'ZA')
write_bpms('../zx/zb.bpms.bmad', 1, N_ZB, 'ZB')
