#!/usr/bin/env python

import os, sys
import shutil
import re, random
from uuid import uuid4


#def main( formula_file, input_file=0 ):
def main(argv):
  ''' Usage:
  python3 convert_regex.py mode formula_file input_files
 
  mode is: -check or -convert

  The input files cannot be binary
  
  The target strings in the input_file will be updated based on the formula_file
  with format:
  old_string1   new_string1
  old_string2   new_string2
  ...           ...
  Example: python3 convert.py -convert formula *txt

  '''
  print(argv)
  
  if len(argv) < 3:
    print('At least 3 arguments required!')  
  else:
    mode = argv[0]
    formula_file = argv[1]
    file_list = argv[2:]
    
    if mode=='-check':
      target_files = check_files(formula_file, file_list)
      print('Target files:', target_files)
    elif mode == '-convert':
      target_files = check_files(formula_file, file_list)
      for i in range(len(target_files)):
        convert_one_file(formula_file, target_files[i])
    else:
      print('The 2nd argument has to be "check" or "convert"!!')

def check_files(formula_file, file_list):

    list1 = parse2(formula_file)
    print('Conversion scheme:',list1)
    
    working_dir = os.getcwd()
    #print('Checking files in...',working_dir)
    print('Checking files in the file list:',file_list)
    target_files = []
    #local_files = os.listdir()
    local_files = file_list
    for i in range(len(local_files)):
      file_now = local_files[i]
      if (file_now == 'convert_regex.py' ):   # skip this python file
        continue
      if (file_now == formula_file):     # skip the formula file
        continue
      if ( not os.path.isfile(file_now) ): # skip directories
        continue
      matches = []
      with open(file_now,'r') as f:
        try:
          for line in f:
            for item1, item2 in list1:
              match = re.findall(item1+'(?!\w)',line,flags = re.IGNORECASE)
              if match: 
                matches.append(match)
      
        except:
          print('error reading', file_now, ', skipping')
      if matches:
        print(file_now, ':', len(matches))      
        target_files.append(file_now)
    

    return target_files

def parse2(formula_file):

  list1=[]
  with open(formula_file,'r') as f:
    for line in f:
      print(line)
      list1.append(line.split()[0:2]) 
  return list1

def convert_one_file(formula_file, input_file):
  list1 = parse2(formula_file)
  newfile = input_file+'_converted'
  with open(input_file,'r') as f2:
    with open(newfile,'w') as f3: 
      for line in f2:
        garbage_str_list = []
        count = 0
        for item1, item2 in list1: 
          # In each line, replace list1[0] 
          #(Case insensitive, and not followed by any a-z, A-Z, 0-9) 
          # with some unique garbage_str
          # THEN replace the unique garbage with list1[1]
          garbage_str = str(uuid4())
          garbage_str_list.append(garbage_str) # save the unique garbage string
          #print(garbage_str)
          #line = re.sub(item1+'(?!\w)', garbage_str, line, flags = re.IGNORECASE)
          line = re.sub('(?<!\w)'+item1+'(?!\w)', garbage_str, line, flags = re.IGNORECASE)
          count = count + 1
        count = 0
        for item1, item2 in list1: 
          #print(garbage_str_list[count])
          line = re.sub(garbage_str_list[count], item2, line, flags = re.IGNORECASE)
          count = count + 1
        f3.write(line)
  
  # overwrite input_file with the temp.txt
  # comment this out to prevent overwrite
  shutil.move(newfile, input_file)

if __name__ == "__main__":
  args = sys.argv[1:]
  if len(args) < 3:
    print(main.__doc__)
  else:
    main(sys.argv[1:])
