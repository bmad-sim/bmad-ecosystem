from tkinter import *

root = Tk()

def key(event):
    print ("pressed", repr(event.char))

def callback(event):
    frame.focus_set()
    print ("clicked at" + str(event.x) + ' ' + str(event.y))

frame = Frame(root, width=100, height=100)
frame.bind("<Key>", key)
frame.bind("<Button-1>", callback)
frame.pack()

root.mainloop()
