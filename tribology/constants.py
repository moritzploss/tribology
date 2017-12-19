# -*- coding: utf-8 -*-

"""

Relevant constants for tribology research. Use, for example, as follows:

.. code-block:: python

   # get young's modulus of steel in units of MPa
   import tribology as tr
   youngs_modulus_steel = tr.YoungsMod.STEEL.value

"""

from enum import Enum


class RadBall(Enum):
    """

    Scalar constants for ball radii in units of mm. Available values are:

    :TQInch: 3/4 inch
    :HInch: 1/2 inch
    :Inch: 1 inch

    """
    TQInch = 25.4 * 3 / 8
    HInch = 25.4 / 4
    Inch = 25.4 / 2


class YoungsMod(Enum):
    """

    Scalar constants for Young's modulus in units of MPa. Available values are:

    :STEEL:
    :GLASS:
    :SiN:

    """
    STEEL = 210000
    GLASS = 70000
    SiN = 315000


class PoissonRatio(Enum):
    """

    Scalar constants for unitless Poisson ratios. Available values are:

    :STEEL:
    :GLASS:
    :SiN:

    """
    STEEL = 0.3
    GLASS = 0.22
    SiN = 0.26


class MatDens(Enum):
    """

    Scalar constants for material density in :math:`\\text{kg m}^{-3}`.
    Available values are:

    :STEEL:
    :GLASS:

    """
    STEEL = 7800
    SiN = 3200


class LubeDens(Enum):
    """

    Scalar constants for lubricant density in :math:`\\text{g ml}^{-1}`.
    Available values are:

    """
    MINERAL_OIL_GENERIC = 0.862


class LubeVisc(Enum):
    """

    Tuple containing kinematic viscosities (in unit of cSt) of lubricants
    at 40 and 100 :math:`^{\\circ}\\text{C}`. The first value in the tuple
    corresponds to the viscosity at 40 :math:`^{\\circ}\\text{C}`, the second
    to 100 :math:`^{\\circ}\\text{C}`. The **TEMPS** value contains the
    reference temperatures (40, 100). Available values are:

    :TEMPS: Tuple containing temperature values of 40 and 100

    """
    TEMPS = (40, 100)


class PressVisc(Enum):
    """

    Pressure-viscosity coefficients in :math:`\\text{Pa}^{-1}` at room
    temperature. Available values are:

    :ESTER_OIL_GENERIC:
    :MINERAL_OIL_GENERIC:

    """
    ESTER_OIL_GENERIC = 15 * 10 ** (-9)
    MINERAL_OIL_GENERIC = 30 * 10 ** (-9)
