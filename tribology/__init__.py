"""
Methods related to general tribology
"""

__all__ = ['tribology',
           'tribology_boundary_element',
           'tribology_constants',
           'tribology_dowson_hamrock',
           'tribology_hertz',
           'tribology_lubrication']

import tribology.tribology
from . import tribology_boundary_element
from . import tribology_constants
from . import tribology_dowson_hamrock
from . import tribology_hertz
from . import tribology_lubrication
