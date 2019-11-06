from .main import tao_root_window
if __name__ == "__main__":
    root = tao_root_window()
    if root.do_mainloop: 
        root.mainloop()
