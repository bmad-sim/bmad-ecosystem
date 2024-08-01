#!/usr/bin/env python3
# coding: utf-8

# # Generate documentation for Tao python commands.
# This will auto-generate documentation for the tao python commands. 
# It extracts a formatted comment above each case statement.

# Note: To run tests, do the following:
#   1) python generate_tex_for_python_commands.py   # Run this file
#   2) cd to pytao directory
#   3) python generate_interface_commands.py
#   4) python run_tests.py

import json
import os
import re
from string import Formatter

#F90_FILE = os.path.join(os.environ['ACC_ROOT_DIR'], 'tao/code/tao_python_cmd.f90')
#F90_FILE = 'tao_python_cmd.f90'
F90_FILE = '../code/tao_python_cmd.f90'

TEXFILE = 'python-interface-commands.tex'
JSONFILE =  'python-interface-commands.json'

assert os.path.exists(F90_FILE)
LINES = open(F90_FILE).readlines()


# In[9]:

def get_params(command_str):
    return [fn for _, fn, _, _ in Formatter().parse(command_str) if fn is not None]

# Collect commands with the header like:
#      !%% beam ----------------
command = None
inside = False
jsonlines = {}
texdoc = []
where = ''

for i, line in enumerate(LINES):
    if line.startswith('!%%'):   
        command = line.split()[1]
        inside = True
        jsonlines[command] = {'description':''}
        texdoc.append({'command':command, 'overview':'', 'notes':''})
        thisdoc = texdoc[-1]
        where = 'overview'

    elif inside and line.startswith('!'):
        line = line[1:]
        if line[0]==' ':
            line = line[1:]
        # Extract command syntax from the next line
        if 'Command syntax' in line:
            where = 'syntax'
            cmd_str = LINES[i+1].strip('!').strip()
            thisdoc['command_str'] = LINES[i+1][1:]
            j = 2
            # Allow continued line
            while cmd_str.endswith('^^'):
                cmd_str+= LINES[i+j].strip('!').strip()
                thisdoc['command_str'] += LINES[i+j][1:]
                j+=1 
            jsonlines[command]['command_str'] = cmd_str
            jsonlines[command]['parameters'] = get_params(cmd_str)
        #
        if where == 'syntax' and j == 0: where = 'post-syntax'
        if where == 'syntax': j -= 1
        if line.strip() == 'Notes': where = 'notes'
        if line.strip() == 'Parameters': where = 'parameters'
        if line.strip() == 'Returns': where = 'other'
        #
        jsonlines[command]['description']+=line
        if where == 'overview': thisdoc['overview'] += line
        if where == 'post-syntax': thisdoc['notes'] += line

    else:
        inside = False

# Write to LaTeX file

re_pe = re.compile(r'((^| )\S+%\S+( |\n|$))')

with open(TEXFILE, 'w') as f:
    f.write('% WARNING: this is automatically generated. DO NOT EDIT.\n')
    for doc in texdoc:
        cmd = doc["command"]
        cmd2 = cmd.replace('_', '.')
        f.write(f'''
%% python {cmd} ------------------------------------
\\subsection{{python {cmd}}}
\\index{{python!{cmd}}}
\\label{{p:{cmd2}}}

''')
        f.write(doc['overview'].replace('&', '\\&'))
        f.write('\\begin{example}\n')
        f.write(doc['command_str'].replace('{','\\{').replace('}','\}'))
        f.write('\\end{example}\n')
        f.write('\\begin{verbatim}\n')
        noff = 0
        #notes = doc['notes'].replace('{','\\{').replace('}','\}').replace('^', '\\^{}') \
        #                                        .replace('&', '\\&').replace('#', '\\#').replace('$', '\\$')
        notes = doc['notes']
        #for match in re.finditer(re_pe, notes):
        #  m = match.group()
        #  ix0 = 0
        #  ix1 = len(m)
        #  if m[0] == ' ': ix0 = 1
        #  if m[-1] == ' ' or m[-1] == '\n': ix1 -= 1
        #  iz0 = match.start() + noff
        #  iz1 = match.end() + noff
        #  notes = notes[:iz0] + m[:ix0] + '\\vn{' + m[ix0:ix1] + '}' + m[ix1:] + notes[iz1:]
        #  noff += 5  # Number of characters added
        #
        f.write(notes.strip() + '\n')
        f.write('\\end{verbatim}\n')

print('Written:', TEXFILE)

# Write to json file

json.dump(jsonlines, open(JSONFILE, 'w'))
print('Written:', JSONFILE)

