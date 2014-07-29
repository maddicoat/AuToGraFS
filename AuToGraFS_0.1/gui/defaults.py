""" This is a module to handle generic ASE (gui) defaults from a ~/.ase/gui.py configuration file, if it exists.
It is imported when opening ag and can then be modified at runtime, if necessary.
syntax for each entry:

gui_default_settings['key'] = value
"""

gui_default_settings = {
    'gui_graphs_string' : 'i, e - E[-1]',   # default for the graph command in the gui
    'covalent_radii' : None
    }

def read_defaults():
    import os.path
    name = os.path.expanduser('~/.ase/gui.py')
    config = gui_default_settings
    if os.path.exists(name):
        execfile(name)
    return config
