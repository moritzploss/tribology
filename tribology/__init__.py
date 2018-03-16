"""

Python init file for tribology package

"""

__all__ = ['tribology',
           'boundary_element',
           'constants',
           'data_import',
           'dowson_hamrock',
           'hertz',
           'lubrication',
           'roller_bearings']

from tribology.tribology import *
from tribology.boundary_element import *
from tribology.constants import *
from tribology.data_import import *
from tribology.dowson_hamrock import *
from tribology.hertz import *
from tribology.lubrication import *
from tribology.roller_bearings import *

import sys
packpath = os.sep.join(os.path.realpath(__file__).split(os.sep)[:-1])
sys.path.insert(0, packpath)

subs_dirs = os.walk(packpath)
for sub_dir in subs_dirs:
    if not sub_dir[0].endswith('__pycache__'):
        sys.path.insert(0, sub_dir[0])

pass