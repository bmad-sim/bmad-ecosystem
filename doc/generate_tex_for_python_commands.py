#!/usr/bin/env python
# coding: utf-8

# # Generate documentation for Tao python commands
# 
# This will auto-generate documentation for the tao python commands. 
# 
# It extracts a formatted comment above each case statement

# In[19]:


import json
import os

#F90_FILE = os.path.join(os.environ['ACC_ROOT_DIR'], 'tao/code/tao_python_cmd.f90')
F90_FILE = '../code/tao_python_cmd.f90'

TEXFILE = 'python-interface-commands.tex'
JSONFILE =  'python-interface-commands.json'

assert os.path.exists(F90_FILE)
LINES = open(F90_FILE).readlines()


# In[20]:


# Collect commands with the header like:
#      !%% beam ----------------
command = None
inside = False
doclines = {}
for line in LINES:
    if line.startswith('!%%'):   
        command = line.split()[1]
        inside = True
        doclines[command] = ''
    elif inside and line.startswith('!'):
        line = line[1:]
        if line[0]==' ':
            line = line[1:]
        doclines[command]+=line
    else:
        inside = False
doclines['beam']


# In[21]:


def tex_from_lines(lines, label):
    text = '\\subsection{python '+label+'}\n'
    text += '\\begin{example}\n'
    text +=  lines
    text += '\n\\end{example}\n'
    return text
#print(tex_from_lines(doclines['beam'], 'X'))


# In[22]:


# Individual command 
def tex_from_lines(lines, label):
    text = '\\item['+label+'] \\Newline'
    text += '\\begin{example}\n'
    text += lines
    text += '\n\\end{example}\n'
    return text
#print(tex_from_lines(doclines['beam'], 'X'))


# In[23]:


# Write to file



with open(TEXFILE, 'w') as f:
    f.write('% WARNING: this is automatically generated. DO NOT EDIT.\n')
    f.write('\\begin{description}\n')
    for k, lines in doclines.items():
        tex = tex_from_lines(lines, k)
        f.write(tex)
    f.write('\\end{description}\n')  

print('Written:', TEXFILE)


# In[24]:


json.dump(doclines, open(JSONFILE, 'w'))
print('Written:', JSONFILE)

