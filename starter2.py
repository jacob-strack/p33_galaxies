
import os
import sys
if not os.path.exists('dtools/starter1.py'):
    print("Need to get the submodule")
    print("git submodule init")
    print("git submodule update")
    sys.exit(-1)
from dtools.starter1 import *

