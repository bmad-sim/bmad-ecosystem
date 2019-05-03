import tkinter as tk
from tkinter import messagebox
import sys

#---------------------------------------------------------------
# List window 

class tao_list_window(tk.Toplevel):

  def __init__(self, root, title, list, *args, **kwargs):
    tk.Toplevel.__init__(self, root, *args, **kwargs)
    self.title(title)


#---------------------------------------------------------------
# Root window 

class tao_root_window(tk.Tk):

  def __init__(self, *args, **kwargs):
    tk.Tk.__init__(self, *args, **kwargs)

    self.title("Tao")
    self.geometry('350x200')
    self.protocol("WM_DELETE_WINDOW", self.quit_cmd)

    # Menu bar

    menubar = tk.Menu(self)

    file_menu = tk.Menu(menubar)
    file_menu.add_command(label = 'Read...', command = self.read_cmd)
    file_menu.add_command(label = 'Write...', command = self.write_cmd)
    file_menu.add_command(label = 'Reinit...', command = self.reinit_cmd)
    file_menu.add_separator()
    file_menu.add_command(label = 'Quit', command = self.quit_cmd, accelerator = 'Ctrl+Q')
    menubar.add_cascade(label = 'File', menu = file_menu)

    window_menu = tk.Menu(menubar)
    window_menu.add_command(label = 'Optimizer...', command = self.optimizer_cmd)
    window_menu.add_command(label = 'Plotting...', command = self.plotting_cmd)
    window_menu.add_command(label = 'Wave...', command = self.wave_cmd)
    menubar.add_cascade(label = 'Window', menu = window_menu)

    self.config(menu=menubar)

    # Init

    if len(sys.argv) > 1:
      pass

    else:
      pass


    init_frame = tk.Frame(self, width = 20, height = 10)
    init_frame.pack()
    beam_str = tk.Entry(init_frame)
    beam_str.pack()
    b = tk.Button(init_frame, text = "Tao Init", command = self.tao_init)
    b.pack()
      
    # Bindings

    self.bind_all("<Control-q>", self.quit_cmd)
    self.bind_all('<Return>', self.return_event)

  #-------------------------

  def tao_init(self):
    window = tao_list_window(self, 'Tao Init', [])


  #-------------------------
  # Menu bar callbacks

  def read_cmd(self):
    print ('Read called')

  def write_cmd(self):
    print ('Write called')

  def reinit_cmd(self):
    print ('Reinit called')

  def quit_cmd(self, event = ''):
    result = messagebox.askquestion("Quit", "Are You Sure?", icon='warning')
    if result == 'yes':
      sys.exit(0)
    else:
      return

  def optimizer_cmd(self):
    print ('Optimizer called')
    win = tk.Toplevel(self)
    win.title('Optimizer')

  def plotting_cmd(self):
    print ('Plotting called')

  def wave_cmd(self):
    print ('Wave called')

  # Other callbacks

  def return_event(self, event):
    print("You hit return.")



#---------------------------------------------------------------

if __name__ == "__main__":
  root = tao_root_window()
  root.mainloop()
