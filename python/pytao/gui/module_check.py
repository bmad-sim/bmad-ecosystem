import pkgutil
import sys
def module_check():
    '''
    Checks for the required modules (listed in the required list).
    If any modules are missing, information on how to install the missing modules
    is printed and execution is halted.
    '''
    missing = [] # Missing modules
    upgrade = False
    required = ["tkinter",
            "pexpect",
            "matplotlib"]
    for mod in required:
        if pkgutil.find_loader(mod) == None:
            missing.append(mod)
    # Check matplotlib version
    if "matplotlib" not in missing:
        import matplotlib
        v = matplotlib.__version__
        if v[0] != '3':
            print("Warning: matplotlib version 3 or higher is required"
                    + " (your version:" + v + ")")
            upgrade = True

    if missing != []:
        print("Warning: some required python modules were not found.")
        print("Missing modules:")
        for mod in missing:
            print(mod)
        print("Please install these modules using your system package manager,"
                + "pip, macports, or whatever you use to manager your python modules")
        print("Linux system package managers:")
        print("      sudo apt-get install python-" + mod)
        print("      sudo pacman -S python-" + mod)
        print("Using pip:")
        print("      pip install " + missing[0])
        print("Using macports (Mac users):")
        print("      sudo port install py36-tkinter")
        print("Note: it is generally a bad idea to use sudo pip install,"
                + " as this may break parts of your system python installation"
                + " that the os depends on to function properly.    "
                + "It is a better idea to set up a virtual environment "
                + "and install your modules there.")
        print("Note: it is a good idea to search for the exact name of the "
                + "modules you need before trying to install them.  "
                + "Use the following commands to do so.")
        print("Debian based distros (incl. Ubuntu):")
        print("      apt-search " + mod)
        print("Arch based distros (incl. Manjaro):")
        print("      pacman -Ss " + mod)
        print("Using pip:")
        print("      pip search " + missing[0])
        print("Using macports (Mac users):")
        print("      port search tkinter")
        sys.exit()
    if upgrade:
        sys.exit()
