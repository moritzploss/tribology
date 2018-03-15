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
           'roller_bearings',
           'p3can']

from tribology.tribology import *
from tribology.boundary_element import *
from tribology.constants import *
from tribology.data_import import *
from tribology.dowson_hamrock import *
from tribology.hertz import *
from tribology.lubrication import *
from tribology.roller_bearings import *

import sys
sys.path.insert(0, "tribology")
sys.path.insert(0, "tribology/p3can")
from tribology.p3can.p3can import *
