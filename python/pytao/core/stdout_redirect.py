import os, sys, select

# the pipe would fail for some reason if I didn't write to stdout at some point
# so I write a space, then backspace (will show as empty in a normal terminal)
#sys.stdout.write(' \b')

pipe_out, pipe_in = os.pipe()

# save a copy of stdout
stdout = os.dup(1)

# save a copy of stderr
stderr = os.dup(2)

def set():
  # replace stdout with our write pipe
  os.dup2(pipe_in, 1)
#  os.dup2(pipe_in, 2)
  
def unset():
  # Restore to terminal
  os.dup2(stdout, 1)
#  os.dup2(stderr, 2)

# check if we have more to read from the pipe
def more_data():
  r, _, _ = select.select([pipe_out], [], [], 0)
  return bool(r)

# read the whole pipe
def read():
  out = ''
  while more_data():
    out += os.read(pipe_out, 1024).decode('utf-8')
  return out

