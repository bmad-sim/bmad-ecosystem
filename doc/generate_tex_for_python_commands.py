#!/usr/bin/env python
# coding: utf-8

# # Generate documentation for Tao python commands
# 
# This will auto-generate documentation for the tao python commands. It extracts the headers above each case statement

# In[1]:


import os

#F90_FILE = os.path.join(os.environ['ACC_ROOT_DIR'], 'tao/code/tao_python_cmd.f90')
F90_FILE = '../code/tao_python_cmd.f90'

OUTFILE = 'python-interface-commands.tex'

assert os.path.exists(F90_FILE)
LINES = open(F90_FILE).readlines()


# In[2]:


from itertools import groupby

# Collect by command
CLINE = '!----------------------------------------------------------------------\n'

HEADERS = [list(group) for k, group in groupby(LINES, lambda x: x == CLINE) if not k][1:-3]
HEADERS[-1]


# In[3]:



def convert_header(lines):
    """
    Convert lines that 
    
    """
    
    comment = []
    for line in lines:
        if line.startswith('!'):
            comment.append(line[1:-1])
        if line.strip().startswith('case'):
            #print(line)
            sline = line.split("'")
            if len(sline)==1:
                return None
            command = sline[1]
            return {command:comment}
    return None

convert_header(HEADERS[-1])


# In[4]:


doclines = {}
for res in map(convert_header, HEADERS):
    if res:
        doclines.update(res)
    else:
        break
doclines.keys()


# In[5]:


def tex_from_lines(lines, label):
    text = '\\subsection{python '+label+'}\n'
    text += '\\begin{example}\n'
    text += ('\n'.join(lines))
    text += '\n\\end{example}\n'
    return text
#print(tex_from_lines(doclines['beam'], 'X'))


# In[6]:


# Individual command 
def tex_from_lines(lines, label):
    text = '\\item['+label+'] \\Newline'
    text += '\\begin{example}\n'
    text += ('\n'.join(lines))
    text += '\n\\end{example}\n'
    return text
#print(tex_from_lines(doclines['beam'], 'X'))


# In[7]:


# Write to file



with open(OUTFILE, 'w') as f:
    f.write('% WARNING: this is automatically generated. DO NOT EDIT.\n')
    f.write('\\begin{description}\n')
    for k, lines in doclines.items():
        tex = tex_from_lines(lines, k)
        f.write(tex)
    f.write('\\end{description}\n')  

print('Written:', OUTFILE)


# In[ ]:




