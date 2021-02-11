#!/usr/bin/env python
# coding: utf-8

# # Generate documentation for Tao python commands
# 
# This will auto-generate documentation for the tao python commands. 
# 
# It extracts a formatted comment above each case statement

# In[8]:


import json
import os
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


# In[10]:


# Collect commands with the header like:
#      !%% beam ----------------
command = None
inside = False
doclines = {}
for i, line in enumerate(LINES):
    if line.startswith('!%%'):   
        command = line.split()[1]
        inside = True
        doclines[command] = {'description':''}
    elif inside and line.startswith('!'):
        line = line[1:]
        if line[0]==' ':
            line = line[1:]
        # Extract command syntax from the next line
        if 'Command syntax' in line:
            cmd_str = LINES[i+1].strip('!').strip()
            j = 2
            # Allow continued line
            while cmd_str.endswith('^^'):
                cmd_str+= LINES[i+j].strip('!').strip()
                j+=1 
            doclines[command]['command_str'] = cmd_str
            doclines[command]['parameters'] = get_params(cmd_str)
        doclines[command]['description']+=line
    else:
        inside = False
doclines['datum_create']['command_str']


# In[11]:


def tex_from_lines(lines, label):
    text = '\\subsection{python '+label+'}\n'
    text += '\\begin{example}\n'
    text +=  lines
    text += '\n\\end{example}\n'
    return text
#print(tex_from_lines(doclines['beam'], 'X'))


# In[12]:


# Individual command 
def tex_from_lines(lines, label):
    text = '\\item['+label+'] \\Newline'
    text += '\\begin{example}\n'
    text += lines
    text += '\n\\end{example}\n'
    return text
#print(tex_from_lines(doclines['beam'], 'X'))


# In[13]:


# Write to file



with open(TEXFILE, 'w') as f:
    f.write('% WARNING: this is automatically generated. DO NOT EDIT.\n')
    f.write('\\begin{description}\n')
    for k, lines in doclines.items():
        tex = tex_from_lines(lines['description'], k)
        f.write(tex)
    f.write('\\end{description}\n')  

print('Written:', TEXFILE)


# In[14]:


json.dump(doclines, open(JSONFILE, 'w'))
print('Written:', JSONFILE)

