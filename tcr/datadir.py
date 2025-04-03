"""
This module defines the default data directory for the library.
It is constructed by joining the directory of the current file with
a subdirectory named `data`. This folder is intended to store any
data files required by the library.
"""

import os
DATA_DIR = os.path.join(os.path.split(__file__)[0], 'data')
